import os
import csv
from utils import load_fa_to_dict


def process_aln(aln_file):
    """
    """
    msa_seq = ""
    with open(aln_file, "r") as f:
        for row in f:
            row = row.replace("\n","")

            if row.startswith("#A1"):
                row = row.split("#A1")[1].replace(" ","")
                msa_seq += row
    f.close()

    # find location of the mature rRNA
    leader_end, trailer_start = "", ""
    joiner_positions = []
    for i in range(len(msa_seq)):
        s = msa_seq[i]
        if s == "0" and leader_end == "":
            leader_end = i
        if s == "0":
            trailer_start = i
            joiner_positions.append(i)

    return msa_seq, leader_end, trailer_start, joiner_positions


def process_mlocarna(mlocarna):
    """
    """
    sequence = ""
    species_to_grab = ""
    with open(mlocarna, "r") as f:
        for row in f:
            row = row.replace("\n","")

            if species_to_grab == "" and "s__" in row:
                species_to_grab = row.split()[0]

            if species_to_grab != "" and row.startswith(species_to_grab):
                row = row.replace(species_to_grab, "")
                seq = row.replace(" ","")
                sequence += seq
    f.close()

    # find location of the mature rRNA
    leader_end, trailer_start = "", ""
    for i in range(len(sequence)):
        s = sequence[i]
        if s == "N" and leader_end == "":
            leader_end = i
        if s == "N":
            trailer_start = i

    return sequence, leader_end, trailer_start


def correct_alirna_ps(alirna_ps, new_alirna_ps, joiner_positions):
    """
    """
    # make sure joiner length is correct
    if len(joiner_positions) != 10:
        raise Exception("unallowed joiner length - something went wrong?", new_alirna_ps, joiner_positions)

    # find MSA sequence
    msa_seq = ""
    with open(alirna_ps, "r") as f:
        grab_seq = False
        for row in f:
            row = row.replace("\n","")
            if row.startswith("/sequence { ("):
                grab_seq = True
            elif grab_seq == True and row.startswith(") } def"):
                grab_seq = False
            elif grab_seq == True:
                seq = row.replace("\\","").replace(" ","")
                msa_seq += seq
    f.close()

    new_msa_seq = ""
    num_nucs_added = 0
    for i in range(len(msa_seq)):
        s = msa_seq[i]

        # if row length for nucs is more than 255, add a \ and next row
        if num_nucs_added != 0 and num_nucs_added % 255 == 0:
            new_msa_seq += "\\"
            new_msa_seq += "\n"

        if i not in joiner_positions:
            # Not an N, so just add the seq
            new_msa_seq += s
        else:
            # Should be an N, so thus should be a dash
            if s != "-":
                raise Exception("Something went wrong with joiner in PS file!", joiner_positions, i, s)
            # Add a hash tag instead for the mature rRNA
            new_msa_seq += "#"
        num_nucs_added += 1

    # at end of new MSA sequence, add a \ and next row
    new_msa_seq += "\\"
    new_msa_seq += "\n"

    # Now, build the new PS file
    with open(new_alirna_ps, "w") as f:
        grab_seq = False
        added_msa = False
        with open(alirna_ps, "r") as r:
            for row in r:
                if row.startswith("/sequence { ("):
                    f.write(row)
                    grab_seq = True
                elif grab_seq == True and row.startswith(") } def"):
                    grab_seq = False
                elif grab_seq == True:
                    if added_msa == False:
                        f.write(new_msa_seq)
                        added_msa = True
                if grab_seq == False:
                    f.write(row)
        r.close()
    f.close()

    return


def process_alirna_ps(alirna_ps):
    """
    """
    grab_base_pairs = False
    base_pairs = []
    with open(alirna_ps, "r") as f:
        for row in f:
            row = row.replace("\n","")
            if row.startswith("/pairs ["):
                grab_base_pairs = True
            elif row.startswith("] def") and grab_base_pairs == True:
                grab_base_pairs = False
            elif grab_base_pairs == True:
                first_bp = int(row.split()[0].replace("[",""))
                last_bp = int(row.split()[1].replace("]",""))
                base_pairs.append([first_bp, last_bp])
    f.close()

    # determine # of base pairs with covariation spots
    grab_annotations = False
    base_pair_w_covariation = []
    with open(alirna_ps, "r") as f:
        for row in f:
            row = row.replace("\n","")
            if row == "":
                continue

            if row.startswith("% Start Annotations"):
                grab_annotations = True
            elif grab_annotations == True and "End Annotations" in row:
                grab_annotations = False
            elif grab_annotations == True and "colorpair" in row:
                first_bp, last_bp, covariation, disagreements, _ = row.split()
                base_pair_w_covariation.append([int(first_bp), int(last_bp), float(covariation), float(disagreements)])

    # counts up what fraction of base_pair_w_covariation are "Red", e.g., covariation = 0.0
    num_red = 0
    for first_bp, last_bp, covariation, disagreements in base_pair_w_covariation:
        if covariation == 0.0:
            num_red += 1
    
    if len(base_pair_w_covariation) > 0:
        fract_red = num_red / len(base_pair_w_covariation)
    else:
        fract_red = "NA"

    return base_pairs, base_pair_w_covariation, fract_red


def count_num_leader_trailer_bps(base_pairs, leader_end, trailer_start, base_pair_w_covariation):
    """
    """
    number_leader_trailer_bps = 0
    leader_trailer_bp_guids = []
    for bp in base_pairs:
        first, second = bp

        if first <= leader_end and second >= trailer_start:
            number_leader_trailer_bps += 1
            leader_trailer_bp_guids.append(str(first) + "|" + str(second))
    
    # see which leader trailer BPs are NOT red
    num_lt_red = 0
    tot_num = 0
    for first_bp, last_bp, covariation, disagreements in base_pair_w_covariation:
        guid = str(first_bp) + "|" + str(last_bp)
        if guid in leader_trailer_bp_guids:
            tot_num += 1

            if covariation == 0.0:
                num_lt_red += 1
    
    if tot_num > 0:
        fract_lt_red = num_lt_red / tot_num
    else:
        fract_lt_red = "NA"

    return number_leader_trailer_bps, fract_lt_red


def count_sw_lengths(base_pairs, sequence, leader_end, trailer_start, \
                        sliding_window_size, min_num_bps_to_count_window):
    """
    """
    tot_seq_len = len(sequence)

    leader_trailer_base_pairings = {}
    for start, stop in base_pairs:
        # Check to see if bp start/stop are flipped
        if start > stop:
            start, stop = stop, start

        if start <= leader_end and stop >= trailer_start:
            leader_trailer_base_pairings[start] = stop

    num_windows_meeting_threshold = 0
    # Run the leader
    for i in range(leader_end - sliding_window_size):
        window_range = range(i, i + sliding_window_size)
        window_qualifies = False
        # get list of all trailer positions pairing w/ these leader positions
        trailers_pairing_to_this_leader_sw = []
        for pos in window_range:
            if pos in leader_trailer_base_pairings:
                trailer = leader_trailer_base_pairings[pos]
                trailers_pairing_to_this_leader_sw.append(trailer)
        
        # determine the maximum # of trailer bps within the window
        for rooting_trailer in trailers_pairing_to_this_leader_sw:
            # use each trailer as a position to grab following trailers
            max_trailer_pos = rooting_trailer + sliding_window_size

            num_trailers_within_range = sum(1 for trailer in trailers_pairing_to_this_leader_sw \
                                                if trailer >= rooting_trailer and trailer <= max_trailer_pos)

            if num_trailers_within_range >= min_num_bps_to_count_window:
                window_qualifies = True

        if window_qualifies == True:
            num_windows_meeting_threshold += 1

    # Run the trailer
    for i in range(trailer_start, tot_seq_len - sliding_window_size):
        window_range = range(i, i + sliding_window_size)
        window_qualifies = False
        # get list of all leader positions pairing w/ these trailer positions
        leaders_pairing_to_this_trailer_sw = []
        for pos in window_range:
            if pos in leader_trailer_base_pairings:
                leader = leader_trailer_base_pairings[pos]
                leaders_pairing_to_this_trailer_sw.append(leader)
        
        # determine the maximum # of trailer bps within the window
        for rooting_leader in leaders_pairing_to_this_trailer_sw:
            # use each trailer as a position to grab following trailers
            max_leader_pos = rooting_leader + sliding_window_size

            num_leaders_within_range = sum(1 for leader in leaders_pairing_to_this_trailer_sw \
                                                if leader >= rooting_leader and trailer <= max_leader_pos)

            if num_leaders_within_range >= min_num_bps_to_count_window:
                window_qualifies = True

        if window_qualifies == True:
            num_windows_meeting_threshold += 1

    return num_windows_meeting_threshold    


def determine_cluster_dissimilarity(cluster, seqid_dict):
    """
    """
    # Determine the identity of the cluster (note: rnaclust appends IDN_)
    core_seqid = ""
    for seqid in seqid_dict:
        if seqid.endswith(cluster.replace("-","_")):
            core_seqid = seqid
        """
        if seqid.startswith("id"):
            seqid_corrected = seqid.replace(seqid.split("_")[0] + "_","").replace("-","_")
        else:
            seqid_corrected = seqid
        """

    if core_seqid == "":
        print("Unable to find cluster:", cluster)
        return "NA"

    # Make query
    with open("temp_query.fa", "w") as f:
        f.write(">" + core_seqid + "\n")
        f.write(seqid_dict[core_seqid] + "\n")
    f.close()

    # parse over other sequences, make a FASTA, conduct alignment, mine output
    min_pident_pqcov = 1.0
    for seqid in seqid_dict:
        seq = seqid_dict[seqid]

        with open("temp_subj.fa", "w") as f:
            f.write(">" + seqid + "\n")
            f.write(seq + "\n")
        f.close()

        # conduct mapping
        command = "mmseqs easy-search" + \
                    " temp_query.fa temp_subj.fa temp_results.m8" + \
                    " tmp --search-type 3"
        os.system(command)

        with open("temp_results.m8", "r") as f:
            hit_result = False
            for row in f:
                row = row.replace("\n","").split("\t")
                hit_result = True
                query, subj, pident, length, mismatch, gapopen, \
                    qstart, qend, sstart, send, evalue, bitscore = row

                signal = float(pident) * (float(length) / len(seqid_dict[core_seqid]))
                if signal < min_pident_pqcov:
                    min_pident_pqcov = signal
        f.close()

        if hit_result == False:
            print("No alignment for:", seqid, core_seqid)
            min_pident_pqcov = 0.0

    return min_pident_pqcov


def analyze_rrna_rnaclust_structs(rnaclust_fldr, genus, rrna_type):
    """
    """
    with open("results/" + rrna_type + "-summary.csv", "w") as ff:
        out = csv.writer(ff)
        out.writerow(["Tax Clade", "Cluster ID", "Number of Seqs in Cluster", \
                        "Fraction of Bases Paired", \
                        "Fraction of all covariation bps that are RED", \
                        "Fraction of all LT covariation bps that are RED", \
                        "Length of RNAClust Structure", \
                        "Number of LT Base Pairs", \
                        "Number of SWs Meeting Conditions (8 in 10)", \
                        "Number of SWs Meeting Conditions (7 in 9)", \
                        "Number of SWs Meeting Conditions (6 in 8)"])
        
        for genus in os.listdir(rnaclust_fldr):
            root_fldr = rnaclust_fldr + genus + "/"
            mmseqs_prior = root_fldr + "mmseqs"
            mmseqs_cluster_location = mmseqs_prior + "_cluster.tsv"

            cluster_seqid_dict = {}
            with open(mmseqs_cluster_location, "r") as f:
                for row in f:
                    cluster, seqid = row.replace("\n","").split("\t")
                    if cluster not in cluster_seqid_dict:
                        cluster_seqid_dict[cluster] = []
                    cluster_seqid_dict[cluster].append(seqid)
            f.close()

            compiled_rnaclust_struct_loc = "results/" + rrna_type + "-rnaclust_structs/"
            os.makedirs(compiled_rnaclust_struct_loc, exist_ok=True)

            for cluster in cluster_seqid_dict:
                cluster_dir = root_fldr + cluster + "/"
                
                rnaclust_cluster_dir = cluster_dir + "rnaclust/"
                rnaclust_seqs = rnaclust_cluster_dir + "seqs.fasta"
                mloc_result_loc = rnaclust_cluster_dir + "maligs.mloc/"
                mlocarna = rnaclust_cluster_dir + "mlocarna.out"

                # check if mloc results directory exists - if it doesn't means too few seqs 
                if os.path.isdir(mloc_result_loc) == False:
                    continue
                    
                rnaclust_struct = mloc_result_loc + "results/alirna.ps"
                rnaclust_aln = mloc_result_loc + "results/result.aln"

                if os.path.isfile(rnaclust_struct) == False:
                    continue

                os.system("cp " + rnaclust_struct + " " + compiled_rnaclust_struct_loc + cluster + ".ps")

                # compile a signal from the seqs (# of seqs and similarity)
                seqid_dict = load_fa_to_dict(rnaclust_seqs)
                num_seqs = len(seqid_dict)

                # process ALN file to find positions of leader/joiner/trailer
                sequence, leader_end, trailer_start, joiner_positions = process_aln(rnaclust_aln)

                # correct the rnaclust PS structure
                correct_alirna_ps(rnaclust_struct, compiled_rnaclust_struct_loc + cluster + ".ps", joiner_positions)

                #sequence, leader_end, trailer_start = process_mlocarna(mlocarna)
                base_pairs, base_pair_w_covariation, fract_red = process_alirna_ps(rnaclust_struct)

                fract_base_pairs = len(base_pairs) / (len(sequence) / 2)

                number_leader_trailer_bps, fract_lt_red = count_num_leader_trailer_bps(base_pairs, leader_end, trailer_start, base_pair_w_covariation)
                num_sws_meeting_threshold_10_8 = count_sw_lengths(base_pairs, sequence, leader_end, trailer_start, 10, 8)
                num_sws_meeting_threshold_9_7 = count_sw_lengths(base_pairs, sequence, leader_end, trailer_start, 9, 7)
                num_sws_meeting_threshold_8_6 = count_sw_lengths(base_pairs, sequence, leader_end, trailer_start, 8, 6)

                out.writerow([genus, cluster, num_seqs, fract_base_pairs, \
                                fract_red, fract_lt_red, len(sequence), \
                                number_leader_trailer_bps, \
                                num_sws_meeting_threshold_10_8, \
                                num_sws_meeting_threshold_9_7, \
                                num_sws_meeting_threshold_8_6])
    ff.close()

    return


def analyze_rnaclust_runs(genus_species_dict, species_genome_dict, \
                rnaclust_16s, rnaclust_23s, \
                rnaclust_lt_bp_summary):
    """
    """
    # parse over 16S rnaclust results
    analyze_rrna_rnaclust_structs(rnaclust_16s, "bact", "16S")

    # parse over 23S rnaclust results
    analyze_rrna_rnaclust_structs(rnaclust_23s, "bact", "23S")

    return


