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


from ete3 import Tree
import os
import csv

def gen_tree_of_results(rrna_type, tree_loc, rnaclust_fldr, \
                        sp_clusters, gtdb_full_tree_loc, \
                        species_genome_dict, genome_species_dict):
    """
    """
    rrna_tree_loc = "results/trees/" + rrna_type + "/"
    os.makedirs(rrna_tree_loc, exist_ok=True)

    adj_species_genome_dict = {}
    for species in species_genome_dict:
        acc = species_genome_dict[species]

        adj_species = species.replace(" ","_").lower()
        adj_species_genome_dict[adj_species] = acc

    all_accs = set()
    acc_signal_dict = {}
    acc_taxonomy_dict = {}
    with open("results/" + rrna_type + "-summary.csv", "r") as f:
        reader = csv.reader(f)
        next(reader, None)

        for row in reader:
            genus, cluster, num_seqs, fract_base_pairs, \
                                fract_red, fract_lt_red, seq_len, \
                                number_leader_trailer_bps, \
                                num_sws_meeting_threshold_10_8, \
                                num_sws_meeting_threshold_9_7, \
                                num_sws_meeting_threshold_8_6 = row

            root_fldr = rnaclust_fldr + genus + "/"
            cluster_dir = root_fldr + cluster + "/"

            covariation_thresholds = 0.75
            sw_threshold = 20

            if float(fract_red) > covariation_thresholds:
                # covariation = poor support of structure
                if float(num_sws_meeting_threshold_8_6) >= sw_threshold:
                    call = "NOCOV,STEM"
                else:
                    call = "NOCOV,NOSTEM"
            else:
                # covariation = good support of structure
                if float(num_sws_meeting_threshold_8_6) >= sw_threshold:
                    call = "COV,STEM"
                else:
                    call = "COV,NOSTEM"

            with open(cluster_dir + cluster + "__" + "cluster-seqs.out", "r") as ff:
                for row in ff:
                    seqid = row.replace("\n","")
                    prior = seqid.split("_")[0]
                    corrected_seqid = seqid.split(prior + "_")[1]
                    
                    counter = seqid.split("-")[-1]
                    taxid = seqid.split(prior + "_")[1][:-(1+len(counter))]

                    if taxid.lower() not in adj_species_genome_dict:
                        print("unable to ID", seqid, taxid)
                        continue
                    acc = adj_species_genome_dict[taxid.lower()]

                    all_accs.add(acc)
                    acc_signal_dict[acc + "-" + counter] = call
                    acc_taxonomy_dict[acc + "-" + counter] = genome_species_dict[acc] + "-" + counter + " [" + str(call) + "]"
            ff.close()
    f.close()

    acc_fullnode_dict = {}
    with open(sp_clusters, "r") as f:
        next(f)
        for row in f:
            data = row.split("\t")[0]
            
            if "GB_" in data:
                acc = data.split("GB_")[1]
            elif "RS_" in data:
                acc = data.split("RS_")[1]
            
            acc_fullnode_dict[acc] = data
    f.close()

    t = Tree(gtdb_full_tree_loc, format=1, quoted_node_names = True)
    bact_tree_folder = rrna_tree_loc + "bact/"
    os.makedirs(bact_tree_folder, exist_ok=True)

    # extract the region of interest (p__Bacteroidota)
    clade_tree_search = t.search_nodes(name = "d__Bacteria")[0] 
    downselected_node = clade_tree_search.detach()

    # now prune the tree to only those nodes of interest
    desired_nodes = set()
    acc_node_map = {}
    for acc in all_accs:
        node_name = acc_fullnode_dict[acc]

        acc_node_map[acc] = node_name
        desired_nodes.add(node_name)

    internal_labels = set()
    clean_external_nodes = set()
    for node in downselected_node.search_nodes():
        if node.is_leaf() == False:
            if "__" in node.name:
                # check if leaves in set
                node_leaves = node.get_leaf_names()
                for n in node_leaves:
                    if n in desired_nodes:
                        internal_labels.add(node.name)
        else:
            if node.name in desired_nodes:
                clean_external_nodes.add(node.name)

    downselected_node.prune(list(clean_external_nodes.union(internal_labels)), preserve_branch_length = True)

    # add lower nodes for multiple rRNA copies per organism
    for guid in acc_signal_dict:
        acc, counter = guid.split("-")
        node = acc_node_map[acc]

        tree_search = downselected_node.search_nodes(name = node)[0]
        tree_search.add_child(name = guid)

    downselected_node.write(format = 1, outfile = bact_tree_folder + "itol.tree")

    # create a color gradient from min to max observed z-score
    from colour import Color
    light_red = Color("#e3cecb") 
    dark_red = Color("#F91607")

    with open(bact_tree_folder + "tree-labels.txt", "w") as f:
        f.write("TREE_COLORS" + "\n")
        f.write("SEPARATOR TAB" + "\n")
        f.write("DATA" + "\n")

        for acc in acc_signal_dict:
            #node_name = acc_node_map[acc]
            node_name = acc
            acc_signal = acc_signal_dict[acc]

            if acc_signal == "NOCOV,STEM":
                color = "#FF7276" # light red
            elif acc_signal == "NOCOV,NOSTEM":
                color = "#A8D8FF" # light blue
            elif acc_signal == "COV,STEM":
                color = "#80393C" # dark red
            elif acc_signal == "COV,NOSTEM":
                color = "#06038D" # dark blue

            f.write(node_name + "\tlabel\t" + str(color) + "\n")
            f.write(node_name + "\tbranch\t" + str(color) + "\tnormal\t1\n")
    f.close()

    with open(bact_tree_folder + "tree-names.txt", "w") as f:
        f.write("LABELS" + "\n")
        f.write("SEPARATOR TAB" + "\n")
        f.write("DATA" + "\n")

        for acc in acc_signal_dict:
            taxonomy = acc_taxonomy_dict[acc]
            node_name = acc
            f.write(node_name + "\t" + taxonomy + "\n")
    f.close()    


    return
