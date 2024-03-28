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


tasks = [
    "sample_download_data", \
    "annotate_genomes",  \
    "rnafold_leader_trailer", \
    "analyze_folding_results"
]


##################################################
# Key choices each time you run the pipeline
rrna_to_analyze = "16S"
cl_of_interest = "d__bacteria" # d__archaea 
clade_of_interest = rrna_to_analyze + "-" + cl_of_interest
rep_genome_only = True
threads = 10

max_dist_upstream = 500 
max_dist_downstream = 300 

num_rand_replicates = 100
num_structs_to_generate = 100
num_middle_insert_Ns = 10

sw_length = "10_8"
min_leader_length = 40
min_trailer_length = 40

if rrna_to_analyze == "16S":
    num_nts_into_rrna_3p = -2 
    num_nts_into_rrna_5p = 15 

elif rrna_to_analyze == "23S":
    num_nts_into_rrna_3p = 1 
    num_nts_into_rrna_5p = -1 

# results paths
result_fldr = "results_" + clade_of_interest + "/"
assembly_data_summary = result_fldr + "taxonomy-genome-summary.txt"

bio_vs_scramb_structs = result_fldr + clade_of_interest + "-" + "bio_vs_scramb_structs-" + str(sw_length) + "nt.csv"
tree_file_fldr = result_fldr + "tree_files-sw-" + str(sw_length) + "nt/"

logging_file = result_fldr + "logging.txt"
genome_bio_rand_structs = result_fldr + rrna_to_analyze + "_structs/"
genome_bio_rand_signals = result_fldr + rrna_to_analyze + "_signals/"

struct_locs = result_fldr + "RNAfold/"
bio_vs_rand_hists = result_fldr + "bio_vs_scramb_hists/"

# support paths
barrnap_annotation_loc = "resources/barrnap_annotations/"

# data dl paths
force_redownload = False
temp_data_root = "../leader_trailer/resources/temp/"
genbank_assembly_summary_loc = temp_data_root + "assembly_summary_genbank.txt"
refseq_assembly_summary_loc = temp_data_root + "assembly_summary_refseq.txt"
taxdmp_loc = temp_data_root + "taxdmp/"
assembly_dl_root = "../leader_trailer/resources/assembly_dls/"
sp_clusters = "../leader_trailer/resources/temp/sp_clusters.tsv"

gtdb_full_tree_loc = "../codon_freqs/resources/bac120.tree"
bact_metadata = temp_data_root + "bac120_metadata_r207.tsv"

# tool paths
barrnap_exe = "barrnap"
rnafold_exe = "RNAfold"
rnaclust_exe = "perl /home/fung/Documents/tools_resources/rnaclust_bin/RNAclust.pl"
rnaclust_local_bin = "/home/fung/Documents/tools_resources/rnaclust_bin"
