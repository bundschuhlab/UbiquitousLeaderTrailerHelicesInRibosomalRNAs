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


# variables
min_number_species_per_genus = 5
min_number_lt_seqs_to_run = 5

num_nts_into_rrna_3p_16S = -2 
num_nts_into_rrna_5p_16S = 15

num_nts_into_rrna_3p_23S = 1 
num_nts_into_rrna_5p_23S = -1

max_dist_upstream = 500 
max_dist_downstream = 300

bp_cutoff_23s = 1000

taxonomy_level_of_interest = "a__"

# fldr paths
rnaclust_16s = "16S_RNAClust/"
rnaclust_23s = "23S_RNAClust/"

rnaclust_lt_bp_summary = "results/rnaclust_summary.csv"
tree_loc = "results/trees/"

# data dl paths
force_redownload = False
temp_data_root = "../leader_trailer/resources/temp/"
genbank_assembly_summary_loc = temp_data_root + "assembly_summary_genbank.txt"
refseq_assembly_summary_loc = temp_data_root + "assembly_summary_refseq.txt"
taxdmp_loc = temp_data_root + "taxdmp/"
assembly_dl_root = "../leader_trailer/resources/assembly_dls/"
sp_clusters = "../leader_trailer/resources/temp/sp_clusters.tsv"
gtdb_full_tree_loc = "../codon_freqs/resources/bac120.tree"

assembly_data_summary = "resources/tot-taxonomy-genome-summary.txt"
barrnap_annotation_loc = "../leader_trailer_v2/resources/barrnap_annotations/"

# tool paths
rnaclust_local_bin = "/home/fung/Documents/tools_resources/rnaclust_bin"
rnaclust_exe = "perl /home/fung/Documents/tools_resources/rnaclust_bin/RNAclust.pl"
