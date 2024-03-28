# Copyright (C) <2022>  <The Ohio State University>       

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
import time


def check_in_taxid_lineage(gtdb_species, clade_of_interest, tree):
    """
    checks if gtdb_species has clade_of_interest in its parent lineage
    """
    if clade_of_interest == "satio":
        satio_species = ["s__archaeoglobus fulgidus", \
                        "s__aquifex aeolicus", \
                        "s__bacillus subtilis", \
                        "s__chlamydia trachomatis", \
                        "s__deinococcus radiodurans", \
                        "s__escherichia coli", \
                        "s__haemophilus influenzae", \
                        "s__helicobacter pylori", \
                        "s__methanocaldococcus jannaschii", \
                        "s__methanothermobacter thermautotrophicus", \
                        "s__mycobacterium tuberculosis", \
                        "s__thermotoga maritima", \
                        "s__treponema pallidum"]
        if gtdb_species.lower() in satio_species:
            return True
        else:
            return False

    parent_lineage = tree.ascend(gtdb_species.lower())
    for node in parent_lineage:
        if node.taxid == clade_of_interest.lower().split("-")[1]:
            return True

    return False


def species_selection(clade_of_interest, sp_clusters, bact_metadata, \
                        rep_genome_only, tree):
    """
    """
    

    # get list of representative accessions per species, for if rep_genome_only is True
    rep_genome_accs = set()
    with open(sp_clusters, "r") as f:
        next(f)
        for row in f:
            row = row.replace("\n","").split("\t")
            gtdb_assembly_acc, species_name = row[:2]

            rep_genome_accs.add(gtdb_assembly_acc)
    f.close()

    # parse over bacterial metadata, identify all genomes that are tagged with "ncbi complete"
    # if rep_genome_only is true, then only grab the rep genome for applicable species
    # Only grab species if it falls underneath clade_of_interest in GTDB taxonomy tree
    species_genome_dict = {}
    accs_added_counter = 0
    with open(bact_metadata, "r") as f:
        next(f)
        for row in f:
            row = row.replace("\n","").split("\t")
            gtdb_assembly_acc = row[0]
            ncbi_genome_annotation = row[45]

            gtdb_taxonomy_lineage = row[16]
            gtdb_species = gtdb_taxonomy_lineage.split(";")[-1]

            if ncbi_genome_annotation == "Complete Genome":
                if (rep_genome_only == True and gtdb_assembly_acc in rep_genome_accs) or rep_genome_only == False:
                    species_in_clade = check_in_taxid_lineage(gtdb_species, clade_of_interest, tree)

                    if species_in_clade == True:
                        # add accession to lookup dictionary!
                        if gtdb_species not in species_genome_dict:
                            species_genome_dict[gtdb_species] = set()

                        gtdb_assembly_acc = gtdb_assembly_acc.replace("RS_", "").replace("GB_", "")
                        species_genome_dict[gtdb_species].add(gtdb_assembly_acc)
                        accs_added_counter += 1
    f.close()

    print("# of species with >=1 acc:", len(species_genome_dict))
    print("total # of accs in consideration:", accs_added_counter)

    return species_genome_dict


def get_assembly_dict(assembly_report_loc):
    """
    """
    assembly_ftp_dict = {}
    with open(assembly_report_loc, "r") as f:
        for row in f:
            if row[0] != "#":
                row = row.replace("\n","").split("\t")
                assembly_acc = row[0]
                asm_acc = row[17]
                ftp_path = row[19]

                # add both assembly_acc and ftp_path to dictionary lookup
                assembly_ftp_dict[assembly_acc] = ftp_path
                assembly_ftp_dict[asm_acc] = ftp_path
    f.close()

    return assembly_ftp_dict



def download_from_assembly_ftp(assembly_bio_data_loc, assembly_acc, ftp_address, ending):
    """
    """
    time.sleep(2) # pause 2 seconds - try not to piss off NCBI...
    #t_stop = time.time() + 60*2 # allow to run for 2 minutes then kill

    web_assembly_report_loc = ftp_address + "/" + ftp_address.split("/")[-1] + ending
    
    assembly_report_dl = assembly_bio_data_loc + assembly_acc + ending
    if ".gbff" in ending:
        assembly_report_dl = assembly_report_dl.replace(".gbff", ".gff")

    command = "timeout 120 wget " + web_assembly_report_loc + " -c -o /dev/null"
    dl_success = os.system(command)
    if dl_success != 0:
        return False

    os.system("mv " + ftp_address.split("/")[-1] + ending + " " + assembly_report_dl)
    
    if ".gz" in ending:
        os.system("gunzip -f " + assembly_report_dl)		   

    return True


def genome_selection(clade_of_interest, \
                        sp_clusters, bact_metadata, \
                        assembly_dl_root, assembly_data_summary, \
                        rep_genome_only, tree, \
                        refseq_assembly_summary_loc):
    """
    """
    refseq_assembly_ftp_dict = get_assembly_dict(refseq_assembly_summary_loc)

    # get dictionary of species : set(genome_accs) meeting conditions
    species_genome_dict = species_selection(clade_of_interest, \
                                                sp_clusters, bact_metadata, \
                                                rep_genome_only, tree)
    
    successful_dl_accs = set()
    for gtdb_species in species_genome_dict:
        counter = 0
        assembly_accs = species_genome_dict[gtdb_species]

        for assembly_acc in assembly_accs:
            if assembly_acc in refseq_assembly_ftp_dict:
                ftp_address = refseq_assembly_ftp_dict[assembly_acc]

                fasta_download_status = download_from_assembly_ftp(assembly_dl_root, assembly_acc, ftp_address, "_genomic.fna.gz")
                gff_download_status = download_from_assembly_ftp(assembly_dl_root, assembly_acc, ftp_address, "_genomic.gff.gz")


                if fasta_download_status == True and gff_download_status == True:
                    successful_dl_accs.add(assembly_acc)
                    counter += 1
                    if counter == 1:
                        break

    # write out successfully downloaded accessions for downstream activities
    with open(assembly_data_summary, "w") as f:
        for gtdb_species in species_genome_dict:
            assembly_accs = species_genome_dict[gtdb_species]

            remaining_accs = assembly_accs.intersection(successful_dl_accs)
            
            if len(remaining_accs) == 0:
                print("NOTICE - no accs could be downloaded for species:", gtdb_species)
            else:
                acc_out_str = ""
                for acc in remaining_accs:
                    acc_out_str = acc_out_str + acc + ";"
                acc_out_str = acc_out_str[:-1]
                f.write(gtdb_species + "\t" + acc_out_str + "\n")
    f.close()

    return