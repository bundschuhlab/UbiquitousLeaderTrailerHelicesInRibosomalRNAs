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


from config import *
from utils import *
from download_ncbi_ftp import load_ncbi_data
from run_rnaclust_study import run_rnaclust
from analyze_rnaclust import analyze_rnaclust_runs
from gen_trees import gen_tree_of_results


# load in gtdb taxonomy tree
tree = load_ncbi_data(temp_data_root, taxdmp_loc, genbank_assembly_summary_loc, \
                    refseq_assembly_summary_loc, force_redownload)

# check if all fldrs exist, if not, make them
check_mk_fldrs([rnaclust_16s, rnaclust_23s, "results/", tree_loc])

# get list of species/accs that have successful data download
species_genome_dict, genome_species_dict = load_in_species_accs(assembly_data_summary)

# get lookup of clusters in the dataset to species
# A cluster is some taxonomic rank (e.g., genus, class)
genus_species_dict = get_genus_species_dict(min_number_species_per_genus, tree, \
                                        species_genome_dict, genome_species_dict, \
                                        taxonomy_level_of_interest)

# parse over each species cluster, run LT/RNAClust analysis
run_rnaclust(genus_species_dict, species_genome_dict, \
                min_number_species_per_genus, \
                min_number_lt_seqs_to_run, \
                max_dist_upstream, max_dist_downstream, \
                num_nts_into_rrna_3p_16S, num_nts_into_rrna_5p_16S, \
                num_nts_into_rrna_3p_23S, num_nts_into_rrna_5p_23S, \
                barrnap_annotation_loc, assembly_dl_root, \
                rnaclust_16s, rnaclust_23s, \
                rnaclust_local_bin, rnaclust_exe)

# parse over rnaclust runs, analyze rnaclust consensus structures
analyze_rnaclust_runs(genus_species_dict, species_genome_dict, \
                rnaclust_16s, rnaclust_23s, \
                rnaclust_lt_bp_summary)

# Create trees
gen_tree_of_results("16S", tree_loc, rnaclust_16s, \
                        sp_clusters, gtdb_full_tree_loc, \
                        species_genome_dict, genome_species_dict)
gen_tree_of_results("23S", tree_loc, rnaclust_23s, \
                        sp_clusters, gtdb_full_tree_loc, \
                        species_genome_dict, genome_species_dict)
