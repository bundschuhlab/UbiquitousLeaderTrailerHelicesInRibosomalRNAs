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


import statistics
import numpy as np
import os
import csv
from multiprocessing import Process, Queue, Pool, Manager


def calculate_signal_sliding_window_3(leader, trailer, mid_region, struct):
    """
    sliding window approach - calculate the # of sliding windows with e.g., 6 out of 8
    consecutive bps 
    """
    base_pairing_positions = []
    for i in range(len(struct)):
        s = struct[i]
        if s == "(":
            base_pairing_positions.append([i, ""])
        elif s == ")":
            for j in range(1, len(base_pairing_positions) + 1):
                start, stop = base_pairing_positions[-j]
                if stop == "":
                    base_pairing_positions[-j][1] = i
                    break
    
    # parse over the base pairing positions, identify those that pair between leader and trailer
    last_leader_position = len(leader)
    first_trailer_position = len(leader) + len(mid_region)
    tot_seq_len = len(leader) + len(trailer) + len(mid_region)

    leader_trailer_base_pairings = {}
    for start, stop in base_pairing_positions:
        if start <= last_leader_position and stop >= first_trailer_position:
            leader_trailer_base_pairings[start] = stop

    sliding_window_size = 10
    min_num_bps_to_count_window = 8
    num_windows_meeting_threshold = 0
    # Run the leader
    for i in range(last_leader_position - sliding_window_size):
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
    for i in range(first_trailer_position, tot_seq_len - sliding_window_size):
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


def get_signal_from_file(filename, genome_bio_rand_signals, acc, counter):
    """
    """
    header_signal_dict = {}
    with open(filename, "r") as f:
        for row in f:
            try:
                first_spot = row.split("\t")[0].strip()
                if first_spot == "Leader":
                    continue
                leader_r, trailer_r, mid_r, header, struct_counter, struct = row.replace("\n","").split("\t")
                # Sequence only written out for the first entry in each set
                if int(struct_counter) == 1:
                    leader_region = leader_r
                    trailer_region = trailer_r
                    mid_region = mid_r
                signal = calculate_signal_sliding_window_3(leader_region, trailer_region, mid_region, struct)

                if header not in header_signal_dict:
                    header_signal_dict[header] = []
                header_signal_dict[header].append(signal)
            except:
                print("Error handling filename:", filename, row.split("\t"))
                return

    bio_signal = 0
    rand_signals = []

    for header in header_signal_dict:
        if header == "BIO":
            bio_signal = np.average(header_signal_dict[header])
        else:
            rand_signals.append(np.average(header_signal_dict[header]))

    if len(rand_signals) < 3:
        return "ERROR", bio_signal, "ERROR", "ERROR", "", "", ""
    mean = statistics.mean(rand_signals)
    stdev = statistics.stdev(rand_signals)

    if stdev != 0.0:
        z_score = (bio_signal - mean) / stdev
    else:
        z_score = "ERROR"

    with open(genome_bio_rand_signals + acc + "-" + counter, "w") as f:
        out = csv.writer(f)
        out.writerow([acc, counter, z_score, bio_signal, \
                        mean, stdev, len(leader_region), \
                        len(trailer_region), len(mid_region)])
    f.close()

    return 


def process_folding_for_signals(species_genome_dict, genome_species_dict, \
                                    genome_bio_rand_structs, bio_vs_scramb_structs, \
                                    genome_bio_rand_signals, threads, \
                                    result_fldr):
    """
    """
    # Set up parallelization
    m = Manager()
    q = m.Queue()
    p = Pool(threads)

    # Initalize parallelization for counting up signals
    tasks = []
    for filename in os.listdir(genome_bio_rand_structs):
        acc = filename.split("-")[0]
        if acc not in genome_species_dict:
            continue
        species = genome_species_dict[acc]
        counter = filename.split("-")[1]

        print(acc, species)
        tasks.append(p.apply_async(
            get_signal_from_file, (genome_bio_rand_structs + filename, \
                                    genome_bio_rand_signals, acc, counter, )))
    
    [r.get() for r in tasks]

    # Compile results across all LTs
    genome_rrna_signals = []
    for filename in os.listdir(genome_bio_rand_signals):
        with open(genome_bio_rand_signals + filename, "r") as f:
            reader = csv.reader(f)
            for row in reader:
                acc, counter, z_score, bio_signal, mean, stdev, leader_len, trailer_len, joiner_len = row
                genome_rrna_signals.append([acc, counter, z_score, bio_signal, mean, stdev, leader_len, trailer_len, joiner_len])
        f.close()
    
    # write out results
    all_leaders = []
    all_trailers = []
    all_sum = []
    with open(bio_vs_scramb_structs, "w") as f:
        out = csv.writer(f)
        out.writerow(["Acc", "Counter", "Z Score", "Bio Signal", "Scrambled Avg", "Scrambled StDev", \
                            "Leader Length", "Trailer Length", "Joiner Length", "GC Fract"])

        for acc, counter, z_score, bio_signal, mean, stdev, leader_len, trailer_len, joiner_len in genome_rrna_signals:
            # get the leader + trailer sequences
            with open(genome_bio_rand_structs + acc + "-" + counter, "r") as r:
                for row in r:
                    first_spot = row.split("\t")[0].strip()
                    if first_spot == "Leader":
                        continue
                    leader_r, trailer_r, mid_r, header, struct_counter, struct = row.replace("\n","").split("\t")
                    # Sequence only written out for the first entry in each set
                    if int(struct_counter) == 1:
                        leader_region = leader_r
                        trailer_region = trailer_r
                        mid_region = mid_r
                        break
            r.close()

            # compute GC content
            total_str = leader_region + trailer_region
            num_G = total_str.count("G")
            num_C = total_str.count("C")
            gc_fract = (num_G + num_C) / len(total_str)
            
            out.writerow([acc, counter, z_score, bio_signal, mean, stdev, leader_len, trailer_len, joiner_len, gc_fract])

            all_leaders.append(leader_len)
            all_trailers.append(trailer_len)
            all_sum.append(leader_len + trailer_len)
    f.close()

    return