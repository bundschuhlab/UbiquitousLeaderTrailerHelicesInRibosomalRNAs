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
from Bio.Seq import Seq
from ete3 import Tree
import csv
import numpy as np

def check_mk_fldrs(fldrs_in):
    """
    """
    for fldr in fldrs_in:
        os.makedirs(fldr, exist_ok = True)            

    return


def load_in_species_accs(assembly_data_summary):
    """
    """
    species_genome_dict = {}
    genome_species_dict = {}
    with open(assembly_data_summary, "r") as f:
        for row in f:
            row = row.replace("\n","").split("\t")
            species, accs = row

            accs = accs.split(";")
            species_genome_dict[species] = accs[0]

            for acc in accs:
                genome_species_dict[acc] = species
    f.close()

    return species_genome_dict, genome_species_dict


def get_genus_species_dict(min_number_species_per_genus, tree, \
                        species_genome_dict, genome_species_dict, \
                        taxonomy_level_of_interest):
    """
    """
    genus_species_dict = {}
    for species in species_genome_dict:
        genus_taxid = ""
        parent_lineage = tree.ascend(species.lower())
        for node in parent_lineage:
            node_taxid = node.taxid

            if node_taxid.startswith(taxonomy_level_of_interest):
                genus_taxid = node_taxid
            elif taxonomy_level_of_interest == "a__":
                genus_taxid = "all"

        if genus_taxid != "":
            if genus_taxid not in genus_species_dict:
                genus_species_dict[genus_taxid] = []
            genus_species_dict[genus_taxid].append(species)

    genus_run_counter = 0
    for genus in genus_species_dict:
        if len(genus_species_dict[genus]) >= min_number_species_per_genus:
            genus_run_counter += 1
    
    print("# of species clusters to run:", genus_run_counter)

    return genus_species_dict


def load_fa_to_dict(fa_in):
    """
    """
    seqid_dict = {}
    with open(fa_in, "r") as f:
        subj, seq = "", ""
        for row in f:
            if row.startswith(">"):
                if len(subj) != 0:
                    seqid_dict[subj] = seq
                subj = row[1:].replace("\n","").split()[0]
                seq = ""
            else:
                seq = seq + row.replace("\n","").upper()
        
        if len(subj) != 0:
            seqid_dict[subj] = seq
    f.close()

    return seqid_dict


def parse_gff_file(gff_file):
    """
    """
    annot_locations = []
    with open(gff_file, "r") as f:
        for row in f:
            if row.startswith("##FASTA"):
                break
            if row.startswith("#"):
                continue
            
            row = row.replace("\n","").split("\t")

            acc = row[0]
            annot_type = row[2]
            annot_start = int(row[3])
            annot_end = int(row[4])
            strand_direction = row[6]
            definition = row[8]

            annot_locations.append([acc, annot_type, \
                                    annot_start, annot_end, \
                                    strand_direction, \
                                    definition])
    f.close()

    return annot_locations


def find_rrna_annots_of_interest(annot_locations, desired_product):
    """
    """
    annots_of_interest = []
    counter = 0

    with open(annot_locations, "r") as f:
        for row in f:
            if row.startswith("#"):
                continue
            
            row = row.replace("\n","").split("\t")
            genome_seqid = row[0]
            position_start = int(row[3])
            position_end = int(row[4])
            strand_direction = row[6]
            feature_string = row[8]

            if feature_string == desired_product:
                counter += 1

                product = desired_product.split("product=")[1]
                annots_of_interest.append([counter, product, genome_seqid, \
                                                strand_direction, position_start, position_end])
    f.close()

    return annots_of_interest


def extract_leader_trailer_region(annot_row, genome_seq_dict, \
                                    max_dist_upstream, max_dist_downstream, \
                                    num_nts_into_rrna_3p, num_nts_into_rrna_5p):
    """
    """
    # find the closest gene upstream and downstream of the annot
    counter, product, acc, strand_direction, annot_start, annot_end = annot_row
    
    # initialize counters to maxes from config
    dist_from_start = max_dist_upstream
    dist_from_end = max_dist_downstream

    # just used fixed thersholds for now!
    # control for ends of genome
    if annot_start < dist_from_start:
        #print("Adjusted distance from start!")
        dist_from_start = annot_start
    if annot_end + dist_from_end > len(genome_seq_dict[acc]):
        #print("Adjusted distance from end!")
        dist_from_end = len(genome_seq_dict[acc]) - annot_end

    # Just use fixed threshold for now!!
    leader_region, trailer_region = extract_genome_leader_trailer_seqs(genome_seq_dict, acc, \
                                        strand_direction, annot_start, annot_end, \
                                        dist_from_start, dist_from_end, \
                                        num_nts_into_rrna_3p, num_nts_into_rrna_5p)

    return counter, product, leader_region, trailer_region


def extract_genome_leader_trailer_seqs(genome_seq_dict, acc, \
                                        strand_direction, annot_start, annot_end, \
                                        dist_from_start, dist_from_end, \
                                        num_nts_into_rrna_3p, num_nts_into_rrna_5p):
    """
    """
    genome_seq = genome_seq_dict[acc]

    if strand_direction == "+":
        leader_region = genome_seq[annot_start - dist_from_start : annot_start + num_nts_into_rrna_5p]
        trailer_region = genome_seq[annot_end - num_nts_into_rrna_3p : annot_end + dist_from_end]
    elif strand_direction == "-":
        leader_region = genome_seq[annot_end - num_nts_into_rrna_5p : annot_end + dist_from_start]
        trailer_region = genome_seq[annot_start - dist_from_end : annot_start + num_nts_into_rrna_3p]

        # take reverse complement 
        dna = Seq(leader_region)
        leader_region = str(dna.reverse_complement())
        dna = Seq(trailer_region)
        trailer_region = str(dna.reverse_complement())

    leader_region = leader_region.replace("T","U")
    trailer_region = trailer_region.replace("T", "U")

    return leader_region, trailer_region
