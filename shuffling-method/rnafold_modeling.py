import subprocess
import RNA
from utils import load_fa_to_dict, check_mk_fldrs
from Bio.Seq import Seq
import random
import statistics
import matplotlib.pyplot as plt
import csv
import os
from multiprocessing import Process, Queue, Pool, Manager
from fasta_dinucleotide_shuffle import dinucleotide_shuffle

# set random seed for consistency
random_seed = 2023
random.seed(random_seed)


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


def extract_leader_trailer_region(annot_row, genome_seq_dict, \
                                    annot_locations, \
                                    max_dist_upstream, max_dist_downstream, \
                                    num_nts_into_rrna_3p, num_nts_into_rrna_5p):
    """
    """
    # find the closest gene upstream and downstream of the annot
    counter, product, acc, strand_direction, annot_start, annot_end = annot_row
    
    # initialize counters to maxes from config
    dist_from_start = max_dist_upstream
    dist_from_end = max_dist_downstream
    r_annot_start, r_annot_end = annot_start, annot_end

    leader_next_annot, trailer_next_annot = "", ""

    # just used fixed thersholds for now!
    # control for ends of genome
    if annot_start < dist_from_start:
        dist_from_start = annot_start
    if annot_end + dist_from_end > len(genome_seq_dict[acc]):
        dist_from_end = len(genome_seq_dict[acc]) - annot_end

    if strand_direction == "+":
        for g_acc, g_annot_type, g_annot_start, g_annot_end, \
                g_strand_direction, g_definition in annot_locations:
            
            if g_annot_type == "rRNA" or acc != g_acc:
                continue
            
            if g_annot_start < r_annot_start and g_annot_end < r_annot_start:
                if r_annot_start - g_annot_start < dist_from_start:
                    dist_from_start = r_annot_start - g_annot_start
                    leader_next_annot = g_definition
                if r_annot_start - g_annot_end < dist_from_start:
                    dist_from_start = r_annot_start - g_annot_end
                    leader_next_annot = g_definition

            if g_annot_start > r_annot_end and g_annot_end > r_annot_end:
                if g_annot_start - r_annot_end < dist_from_end:
                    dist_from_end = g_annot_start - r_annot_end
                    trailer_next_annot = g_definition
                if g_annot_end - r_annot_end < dist_from_end:
                    dist_from_end = g_annot_end - r_annot_end
                    trailer_next_annot = g_definition
                    

    elif strand_direction == "-":
        for g_acc, g_annot_type, g_annot_start, g_annot_end, \
                g_strand_direction, g_definition in annot_locations:
            
            if g_annot_type == "rRNA" or acc != g_acc:
                continue
            
            # check if either start or end is closest to leader (pos end)
            if g_annot_start > r_annot_end and g_annot_end > r_annot_end:
                if g_annot_start - r_annot_end < dist_from_start:
                    dist_from_start = g_annot_start - r_annot_end
                    leader_next_annot = g_definition
                if g_annot_end - r_annot_end < dist_from_start:
                    dist_from_start = g_annot_end - r_annot_end
                    leader_next_annot = g_definition

            if g_annot_start < r_annot_start and g_annot_end < r_annot_start:
                if r_annot_start - g_annot_start < dist_from_end:
                    dist_from_end = r_annot_start - g_annot_start
                    trailer_next_annot = g_definition
                if r_annot_start - g_annot_end < dist_from_end:
                    dist_from_end = r_annot_start - g_annot_end
                    trailer_next_annot = g_definition

    # extract genome leader/trailer regions
    leader_region, trailer_region = extract_genome_leader_trailer_seqs(genome_seq_dict, acc, \
                                        strand_direction, annot_start, annot_end, \
                                        dist_from_start, dist_from_end, \
                                        num_nts_into_rrna_3p, num_nts_into_rrna_5p)

    return counter, product, leader_region, trailer_region


def run_rnafold_struct_only(leader_region, trailer_region, \
                    num_middle_insert_Ns, \
                    num_structs_to_generate, \
                    header, struct_file):
    """
    """
    nuc_seq = ""
    constraint_seq = ""

    # add leader region
    nuc_seq += leader_region
    for s in leader_region:
        constraint_seq += "."

    # add middle region
    mid_region = ""
    for i in range(num_middle_insert_Ns):
        mid_region += "N"

    nuc_seq += mid_region
    for s in mid_region:
        constraint_seq += "x"

    # add trailer region
    nuc_seq += trailer_region
    for s in trailer_region:
        constraint_seq += "." 

    # create model details
    md = RNA.md()

    # activate unique multibranch loop decomposition
    md.uniq_ML = 1
    
    # create fold_compound data structure (required for all subsequently applied  algorithms)
    fc = RNA.fold_compound(nuc_seq, md)

    # add constraints for non-pairing
    fc.constraints_add(constraint_seq, RNA.CONSTRAINT_DB_DEFAULT) 

    # compute MFE and MFE structure
    (mfe_struct, mfe) = fc.mfe()
    
    # rescale Boltzmann factors for partition function computation
    fc.exp_params_rescale(mfe)

    # compute partition function to fill DP matrices
    fc.pf()

    # predict multiple structures
    struct_counter = 0
    for struct in fc.pbacktrack(num_structs_to_generate):
        # Struct = dot/parantheses notation
        # struct_energy computes dG (kcal/mol)
        struct_energy = fc.eval_structure(struct)

        struct_counter += 1

        # if struct counter gets to be more than 1 - e.g., only write out sequence for first entry
        if struct_counter > 1:
            leader_region = ""
            trailer_region = ""
            mid_region = ""

        struct_file.write(str(leader_region) + "\t" + str(trailer_region) + "\t" + \
                str(mid_region) + "\t" + \
                header + "\t" + str(struct_counter) + "\t" + struct + "\n")
    
    return
        

def run_folding(acc, counter, leader_region, trailer_region, \
                        genome_bio_rand_structs, \
                        num_middle_insert_Ns, num_rand_replicates, \
                        num_structs_to_generate):
    """
    """
    # open a file - blank if existing - then open in append mode
    struct_file = open(genome_bio_rand_structs + acc + "-" + str(counter), "w")
    struct_file.write("Leader \t Trailer \t Joiner \t Header \t Header_Num \t Structure \n")
    struct_file = open(genome_bio_rand_structs + acc + "-" + str(counter), "a")

    # gen bio structures
    run_rnafold_struct_only(leader_region, trailer_region, \
                    num_middle_insert_Ns, \
                    num_structs_to_generate, \
                    "BIO", struct_file)

    # gen scrambled structures
    # scramble the sequence with dinucleotide shuffling
    scrambled_leaders = dinucleotide_shuffle(leader_region, random_seed, num_rand_replicates)
    scrambled_trailers = dinucleotide_shuffle(trailer_region, random_seed, num_rand_replicates)

    for i in range(num_rand_replicates):
        scrambled_leader = scrambled_leaders[i]
        scrambled_trailer = scrambled_trailers[i]

        run_rnafold_struct_only(scrambled_leader, scrambled_trailer, \
                        num_middle_insert_Ns, \
                        num_structs_to_generate, \
                        "SCRAMB_" + str(i), struct_file)

    return


def run_single_genome_rnafold(fasta_loc, \
                                max_dist_upstream, max_dist_downstream, \
                                acc_struct_locs, num_middle_insert_Ns, num_rand_replicates, \
                                threads, rnafold_exe, acc_bio_vs_rand_hists, \
                                species, num_structs_to_generate, acc, gff_loc, acc_barrnap_loc, \
                                num_nts_into_rrna_3p, num_nts_into_rrna_5p, \
                                rrna_to_analyze, genome_bio_rand_structs, \
                                min_leader_length, min_trailer_length):
    """
    """
    print("Now running ", acc, species)
    
    # load genome dict
    genome_seq_dict = load_fa_to_dict(fasta_loc)

    # extract all annotation locations
    annot_locations = parse_gff_file(gff_loc)

    # find locations of interest
    annots_of_interest = find_rrna_annots_of_interest(acc_barrnap_loc, "Name=" + rrna_to_analyze + "_rRNA;product=" + rrna_to_analyze + " ribosomal RNA")

    for annot_row in annots_of_interest:
        counter, _, _, _, _, _ = annot_row
        guid = acc + "|" + str(counter)

        # extract leader/trailer region
        counter, product, leader_region, trailer_region = extract_leader_trailer_region(annot_row, \
                                                                genome_seq_dict, \
                                                                annot_locations, \
                                                                max_dist_upstream, max_dist_downstream, \
                                                                num_nts_into_rrna_3p, num_nts_into_rrna_5p,)
        
        # disclude LTs whose minimum lengths are less than the arbitrary threshold
        if len(leader_region) < min_leader_length or len(trailer_region) < min_trailer_length:
            continue

        # run folding, generate 100 structs. save them
        run_folding(acc, counter, leader_region, trailer_region, \
                        genome_bio_rand_structs, \
                        num_middle_insert_Ns, num_rand_replicates, \
                        num_structs_to_generate)
    
    # Clear the variables to reclaim memory
    genome_seq_dict = {}
    annot_locations = []
    annots_of_interest = []

    return


def determine_already_done_genomes(genome_species_dict, barrnap_annotation_loc, \
                                    genome_bio_rand_compar, genome_bio_rand_structs):
    """
    """
    # Determine if an acc has been logged or not
    already_done_accs = set()
    for filename in os.listdir(genome_bio_rand_structs):
        acc = filename.split("-")[0]
        already_done_accs.add(acc)
    return already_done_accs


def run_shuffling_method_core(species_genome_dict, genome_species_dict, \
                            struct_locs, barrnap_annotation_loc, \
                            max_dist_upstream, max_dist_downstream, \
                            num_middle_insert_Ns, num_rand_replicates, \
                            threads, rnafold_exe, bio_vs_rand_hists, \
                            bio_vs_scramb_structs, assembly_dl_root, \
                            num_structs_to_generate, \
                            rrna_to_analyze, genome_bio_rand_structs, \
                            num_nts_into_rrna_3p, num_nts_into_rrna_5p, \
                            result_fldr, genome_bio_rand_signals, \
                            min_leader_length, min_trailer_length):
    """
    """
    # Set up the parallel tasks for multithreading
    m = Manager()
    q = m.Queue()
    p = Pool(threads)

    # Generate structures
    # load a list of already completed accs (from previous runs)
    already_done_accs = determine_already_done_genomes(genome_species_dict, barrnap_annotation_loc, \
                                                        genome_bio_rand_structs, genome_bio_rand_structs)
    
    # Initialize the LTs to-do for downstream parallelization
    tasks = []
    for acc in genome_species_dict:
        species = genome_species_dict[acc]
        print(acc, species)

        if acc in already_done_accs:
            continue

        acc_barrnap_loc = barrnap_annotation_loc + acc + "-barrnap.out"

        fasta_loc = assembly_dl_root + acc + "_genomic.fna"
        gff_loc = assembly_dl_root + acc + "_genomic.gff"

        acc_struct_locs = struct_locs + acc + "/"
        acc_bio_vs_rand_hists = bio_vs_rand_hists + acc + "/"
        check_mk_fldrs([acc_struct_locs, acc_bio_vs_rand_hists])

        tasks.append(p.apply_async(
            run_single_genome_rnafold, (fasta_loc, \
                                    max_dist_upstream, max_dist_downstream, \
                                    acc_struct_locs, num_middle_insert_Ns, num_rand_replicates, \
                                    threads, rnafold_exe, acc_bio_vs_rand_hists, \
                                    species, num_structs_to_generate, \
                                    acc, gff_loc, acc_barrnap_loc, \
                                    num_nts_into_rrna_3p, num_nts_into_rrna_5p, \
                                    rrna_to_analyze, genome_bio_rand_structs, \
                                    min_leader_length, min_trailer_length, \
                                    )))

    print("# of genomes already done:", len(already_done_accs))
    print("# of genomes for tasking:", len(tasks))

    # Kick off the parallel tasks 
    [r.get() for r in tasks]

    return