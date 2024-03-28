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
from utils import check_mk_fldrs, load_in_species_accs
from download_ncbi_ftp import load_ncbi_data
from genome_selection_activities import genome_selection
from annotate_genomes import run_barrnap
from rnafold_modeling import run_shuffling_method_core
from analyze_folding_results import process_folding_for_signals

def main():
    """
    Main drives the entire pipeline, call the pipeline by adjusting paths in
    config.py and calling "python3 main.py"
    """
    # This function downloads data from NCBI and GTDB related to taxonomy and genomes
    # This function also loads in a Tree object to navigate the gtdb taxonomy tree
    tree = load_ncbi_data(temp_data_root, taxdmp_loc, genbank_assembly_summary_loc, \
                        refseq_assembly_summary_loc, force_redownload)

    # Check if necessary folders exist, and if not, create them
    check_mk_fldrs([result_fldr, barrnap_annotation_loc, struct_locs, \
                        bio_vs_rand_hists, \
                        tree_file_fldr, genome_bio_rand_structs, \
                        genome_bio_rand_signals])
    
    # Initilize logging
    logging = open(logging_file, "a+")

    # This function downloads RefSeq genomes for every selected species of interest
    if "sample_download_data" in tasks:
        print ("Downloading genome data and annotations")
        genome_selection(clade_of_interest, \
                        temp_data_root + "sp_clusters.tsv", \
                        bact_metadata, \
                        assembly_dl_root, assembly_data_summary, \
                        rep_genome_only, tree, \
                        refseq_assembly_summary_loc)

    # This function obtains 2 Python dictionaries of downloaded genomes, one indexed by
    # RefSeq accession and the other by GTDB taxonomy species name 
    species_genome_dict, genome_species_dict = load_in_species_accs(assembly_data_summary, cl_of_interest, tree)

    # This function runs Barrnap over the downloaded genomes to obtain 16S and 23S annotations
    if "annotate_genomes" in tasks:
        print("Annotate each of the genomes to standardize gene/RNA annotations")
        run_barrnap(species_genome_dict, genome_species_dict, \
                        barrnap_exe, assembly_dl_root, \
                        barrnap_annotation_loc, \
                        threads, logging, \
                        cl_of_interest)

    # This function runs the core shuffling method pipeline over each 16S and 23S 
    # pre-rRNA leader trailer to generate shuffled sequences and 100 structures
    # for each sequence (1 biological and 100 shuffled per rRNA)
    if "rnafold_leader_trailer" in tasks:
        print("Extracting leader/trailer regions, running RNAFold on bio + scrambled")
        run_shuffling_method_core(species_genome_dict, genome_species_dict, \
                            struct_locs, barrnap_annotation_loc, \
                            max_dist_upstream, max_dist_downstream, \
                            num_middle_insert_Ns, num_rand_replicates, \
                            threads, rnafold_exe, bio_vs_rand_hists, \
                            bio_vs_scramb_structs, assembly_dl_root, \
                            num_structs_to_generate, \
                            rrna_to_analyze, genome_bio_rand_structs, \
                            num_nts_into_rrna_3p, num_nts_into_rrna_5p, \
                            result_fldr, genome_bio_rand_signals, \
                            min_leader_length, min_trailer_length)

    # This function analyzes the generated LT fold structures, and counts up
    # a signal for each structure. For each original pre-rRNA, the biological
    # signals are compared to the scrambled sequences to extract a z-score
    if "analyze_folding_results" in tasks:
        print("Analyzing folding results to find bio signal from scrambled noise")
        process_folding_for_signals(species_genome_dict, genome_species_dict, \
                                    genome_bio_rand_structs, bio_vs_scramb_structs, \
                                    genome_bio_rand_signals, threads, \
                                    result_fldr)

    return


main()