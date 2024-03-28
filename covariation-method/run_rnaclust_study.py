# Copyright (C) <2024>  <The Ohio State University>       

# This program is free software: you can redistribute it and/or modify                              
# it under the terms of the GNU General Public License as published by 
# the Free Software Foundation, either version 3 of the License, or    
# (at your option) any later version.                                                                                       
# This program is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the implied warranty of           
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       
# GNU General Public License for more details.                                                                             

# You should have received a copy of the GNU General Public License 
# along with this program.  If not, see <https://www.gnu.org/licenses/>.


import os
from utils import *


def run_rnaclust_one_rRNA_type(acc_barrnap_loc, rrna_to_analyze, \
                                genome_seq_dict, annot_locations, \
                                max_dist_upstream, max_dist_downstream, \
                                num_nts_into_rrna_3p, num_nts_into_rrna_5p, \
                                leader_trailer_regions, \
                                acc, species):
    """
    """
    # find locations of interest
    annots_of_interest = find_rrna_annots_of_interest(acc_barrnap_loc, "Name=" + rrna_to_analyze + "_rRNA;product=" + rrna_to_analyze + " ribosomal RNA")

    for annot_row in annots_of_interest:
        counter, product, leader_region, trailer_region = extract_leader_trailer_region(annot_row, \
                                                                genome_seq_dict, \
                                                                max_dist_upstream, max_dist_downstream, \
                                                                num_nts_into_rrna_3p, num_nts_into_rrna_5p)
        
        guid = acc + "|" + str(counter)


        if len(leader_region) != max_dist_upstream + num_nts_into_rrna_5p or \
            len(trailer_region) != max_dist_downstream + num_nts_into_rrna_3p:
            print("Leader/trailer region too short! Target Leader Length: {} Actual: {} Target Trailer Length: {} Actual: {}".format( \
                        max_dist_upstream + num_nts_into_rrna_5p, len(leader_region), \
                        max_dist_downstream + num_nts_into_rrna_3p, len(trailer_region)))
            continue

        leader_trailer_regions.append([acc, counter, species, leader_region, trailer_region])

    return leader_trailer_regions


def make_fasta(seq_loc, leader_trailer_regions):
    """
    """
    seq_counter = 0
    with open(seq_loc, "w") as f:
        for acc, counter, species, leader, trailer in leader_trailer_regions:
            guid = species.replace(" ","_") + "-" + str(counter)
            sequence = ""

            for s in leader:
                sequence += s

            for i in range(10):
                sequence += "N"

            for s in trailer:
                sequence += s

            seq_counter += 1

            f.write(">" + str(seq_counter) + "_" + guid + "\n")
            f.write(sequence + "\n")
    f.close()

    return


def run_single_mmseqs_rnaclust(rnaclust_fldr, genus, leader_trailer_regions, \
                                rnaclust_exe, rnaclust_local_bin, \
                                max_dist_upstream, num_nts_into_rrna_5p):
    """
    """
    root_fldr = rnaclust_fldr + genus + "/"
    check_mk_fldrs([root_fldr])

    seq_loc_1 = root_fldr + "step1.fa"
    make_fasta(seq_loc_1, leader_trailer_regions)

    mmseqs_prior = root_fldr + "mmseqs"
    command = "mmseqs" + " easy-cluster" + \
                " " + seq_loc_1 + \
                " " + mmseqs_prior + \
                " tmp" + \
                " --min-seq-id " + str(30.0 / 100.0) + \
                " -c " + str(50.0 / 100.0) + \
                " --cov-mode 1" + \
                " > /dev/null" 
    os.system(command)

    mmseqs_fasta_loc = mmseqs_prior + "_rep_seq.fasta"
    mmseqs_cluster_location = mmseqs_prior + "_cluster.tsv"

    cluster_seqid_dict = {}
    with open(mmseqs_cluster_location, "r") as f:
        for row in f:
            cluster, seqid = row.replace("\n","").split("\t")
            if cluster not in cluster_seqid_dict:
                cluster_seqid_dict[cluster] = []
            cluster_seqid_dict[cluster].append(seqid)
    f.close()

    seqid_dict = {}
    with open(seq_loc_1, "r") as f:
        subj, seq = "", ""
        for row in f:
            if row.startswith(">"):
                if len(subj) != 0:
                    seqid_dict[subj] = seq
                subj = row[1:].replace("\n","")
                seq = ""
            else:
                seq = seq + row.replace("\n","")
        if len(subj) != 0:
            seqid_dict[subj] = seq
    f.close()

    num_Clusters_over_10 = 0
    for cluster in cluster_seqid_dict:
        seqids = cluster_seqid_dict[cluster]
        if len(seqids) >= 10:
            num_Clusters_over_10 += 1

    print("The # of clusters are:", len(cluster_seqid_dict))
    print("# of clusters with at least 10 seqs:", num_Clusters_over_10)

    for cluster in cluster_seqid_dict:
        cluster_dir = root_fldr + cluster + "/"
        os.makedirs(cluster_dir, exist_ok=True)

        seqids = cluster_seqid_dict[cluster]
        with open(cluster_dir + cluster + "__" + "cluster-seqs.out", "w") as f:
            for seqid in seqids:
                f.write(seqid + "\n")
        f.close()

        fasta_loc = cluster_dir + "seqs.fa"
        with open(fasta_loc, "w") as f:
            for seqid in seqids:
                seq = seqid_dict[seqid]

                f.write(">" + seqid + "\n")
                f.write(seq + "\n")
        f.close()

        # run mmseqs to cluster 100% seqs together to save compute time
        mmseqs_prior = cluster_dir + "mmseqs"
        command = "mmseqs" + " easy-cluster" + \
                    " " + fasta_loc + \
                    " " + mmseqs_prior + \
                    " tmp" + \
                    " --min-seq-id " + str(100.0 / 100.0) + \
                    " -c " + str(100.0 / 100.0) + \
                    " --cov-mode 1" + \
                    " > /dev/null" 
        os.system(command)

        mmseqs_fasta_loc = mmseqs_prior + "_rep_seq.fasta"
        mmseqs_cluster_location = mmseqs_prior + "_cluster.tsv"        

        rnaclust_loc = cluster_dir + "rnaclust/"
        if os.path.isdir(rnaclust_loc) == True:
            os.system("rm -rf " + rnaclust_loc)
        
        # create anchor file on-the-fly
        seqids_in_file = []
        with open(mmseqs_fasta_loc, "r") as f:
            for row in f:
                if row.startswith(">"):
                    subj = row[1:].split()[0].replace("\n","")
                    seqids_in_file.append(subj)
        f.close()

        anchor_loc = cluster_dir + "anchor_file.txt"
        with open(anchor_loc, "w") as f:
            for seqid in seqids_in_file:
                init_pos = max_dist_upstream + num_nts_into_rrna_5p
                final_pos = init_pos + 10
                f.write(seqid.replace("-","_") + "\t" + str(init_pos) + "\t" + str(final_pos) + "\t" + "mature_rrna" + "\n")
        f.close()

        command = rnaclust_exe + \
                    " --fasta " + mmseqs_fasta_loc + \
                    " --cpu " + str(12) + \
                    " --localbin " + rnaclust_local_bin + \
                    " --dir " + rnaclust_loc + \
                    " --rnasoup" + \
                    " --mlocarna-opts '--anchor-constraints " + anchor_loc + "'"
        os.system(command)



    """
    rnaclust_loc = root_fldr + "rnaclust/"
    if os.path.isdir(rnaclust_loc) == True:
        os.system("rm -rf " + rnaclust_loc)
    
    command = rnaclust_exe + \
                " --fasta " + mmseqs_fasta_loc + \
                " --cpu " + str(12) + \
                " --localbin " + rnaclust_local_bin + \
                " --dir " + rnaclust_loc + \
                " --rnasoup"
    os.system(command)
    """

    return


def run_rnaclust(genus_species_dict, species_genome_dict, \
                min_number_species_per_genus, \
                min_number_lt_seqs_to_run, \
                max_dist_upstream, max_dist_downstream, \
                num_nts_into_rrna_3p_16S, num_nts_into_rrna_5p_16S, \
                num_nts_into_rrna_3p_23S, num_nts_into_rrna_5p_23S, \
                barrnap_annotation_loc, assembly_dl_root, \
                rnaclust_16s, rnaclust_23s, \
                rnaclust_local_bin, rnaclust_exe):
    """
    """
    for genus in genus_species_dict:
        print("Running clade:", genus)

        genus_species = genus_species_dict[genus]
        if len(genus_species) < min_number_species_per_genus:
            continue

        leader_trailer_regions_16S = []
        leader_trailer_regions_23S = []
        for species in genus_species:
            acc = species_genome_dict[species]

            acc_barrnap_loc = barrnap_annotation_loc + acc + "-barrnap.out"
            fasta_loc = assembly_dl_root + acc + "_genomic.fna"
            gff_loc = assembly_dl_root + acc + "_genomic.gff"

            # load genome dict
            genome_seq_dict = load_fa_to_dict(fasta_loc)

            # extract all annotation locations
            annot_locations = parse_gff_file(gff_loc)

            leader_trailer_regions_16S = run_rnaclust_one_rRNA_type(acc_barrnap_loc, "16S", \
                                            genome_seq_dict, annot_locations, \
                                            max_dist_upstream, max_dist_downstream, \
                                            num_nts_into_rrna_3p_16S, num_nts_into_rrna_5p_16S, \
                                            leader_trailer_regions_16S, \
                                            acc, species)
            leader_trailer_regions_23S = run_rnaclust_one_rRNA_type(acc_barrnap_loc, "23S", \
                                            genome_seq_dict, annot_locations, \
                                            max_dist_upstream, max_dist_downstream, \
                                            num_nts_into_rrna_3p_23S, num_nts_into_rrna_5p_23S, \
                                            leader_trailer_regions_23S, \
                                            acc, species)


        if len(leader_trailer_regions_16S) >= min_number_lt_seqs_to_run:
            run_single_mmseqs_rnaclust(rnaclust_16s, genus, leader_trailer_regions_16S, \
                                rnaclust_exe, rnaclust_local_bin, \
                                max_dist_upstream, num_nts_into_rrna_5p_16S)
        if len(leader_trailer_regions_23S) >= min_number_lt_seqs_to_run:
            run_single_mmseqs_rnaclust(rnaclust_23s, genus, leader_trailer_regions_23S, \
                                rnaclust_exe, rnaclust_local_bin,  \
                                max_dist_upstream, num_nts_into_rrna_5p_23S)
        


    return
