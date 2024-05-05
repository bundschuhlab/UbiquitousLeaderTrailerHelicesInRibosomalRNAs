import pandas
import os
import numpy as np
import math
import csv
from plot_utils import gen_histogram, gen_barplot_taxonomy, \
                    gen_histogram_cut_axes, \
                    gen_barplot, gen_colored_scatterplot

rrna_types = ["16S", "23S"]
for rrna_type in rrna_types:
    print("rRNA type:", rrna_type)


    plot_loc = "plots_" + rrna_type + "/"
    os.makedirs(plot_loc, exist_ok = True)

    data_loc = "combined-" + rrna_type + "-concat.xlsx"

    data = pandas.read_excel(data_loc, sheet_name='main')

    # First, gen the single-shuffled z-score histograms
    z_scores = []
    for z_score in data['SINGLE: Z Score']:
        if math.isnan(z_score) == False:
            z_scores.append(float(z_score))

    x_min, x_max = -2.5, 15.0
    y_max = 0.15

    bins = np.linspace(x_min, x_max, num=int(((x_max - x_min) / 0.5)+1))
    gen_histogram(plt_loc = plot_loc + rrna_type + "-singleseq-zscores.png", \
                    title = "Sequence-Shuffled " + rrna_type + " Z-Scores", \
                    xname = "Z-Scores", \
                    yname = "Frequency", \
                    data = z_scores, \
                    bins = bins, \
                    x_lims = [x_min, x_max], \
                    y_lims = [0, y_max], \
                    vert_lines = [0, 1, 2])

    # Second, generate RNAclust SW signal plots
    rnaclust_numsws = []
    for score in data['RNACLUST: Num SWs 8_10']:
        if math.isnan(score) == False:
            rnaclust_numsws.append(float(score))

    x_min, x_max = 0, 115
    y_max = 0.31
    bins = np.linspace(x_min, x_max, num=int(((x_max - x_min) / 5)+1))
    gen_histogram(plt_loc = plot_loc + rrna_type + "-rnaclust-sws.png", \
                    title = "Sequence-Covariation " + rrna_type + " Number SWs", \
                    xname = "Number of SSWs", \
                    yname = "Frequency", \
                    data = rnaclust_numsws, \
                    bins = bins, \
                    x_lims = [x_min, x_max], \
                    y_lims = [0, y_max], \
                    vert_lines = [15])

    # Third, generate RNAclust covariation plots
    rnaclust_covariation = []
    for score in data['RNACLUST: Covariation']:
        if score != "NA" and math.isnan(score) == False:
            rnaclust_covariation.append(float(score))

    x_min, x_max = 0, 1.0
    y_max = 0.26
    bins = np.linspace(x_min, x_max, num=int(((x_max - x_min) / 0.05)+1))
    gen_histogram(plt_loc = plot_loc + rrna_type + "-rnaclust-covariation.png", \
                    title = "Sequence-Covariation " + rrna_type + " Covariation", \
                    xname = "Fraction LT bps with Covariation", \
                    yname = "Frequency", \
                    data = rnaclust_covariation, \
                    bins = bins, \
                    x_lims = [x_min, x_max], \
                    y_lims = [0, y_max], \
                    vert_lines = [0.25])

    # Find fraction of each taxonomy 
    singleseq_threshold = 1
    rnaclust_sws_threshold = 15
    rnaclust_covariation_threshold = 0.25

    # Load a lookup of species : lineage array
    species_lineage_dict = {}
    with open("combined_lineages.csv", "r") as f:
        reader = csv.reader(f)
        next(reader, None)
        for row in reader:
            s__,g__,f__,o__,c__,p__,d__ = row
            species_lineage_dict[s__] = [g__,f__,o__,c__,p__,d__]
    f.close()

    tax_signals_dict = {}
    taxa_domain_dict = {}
    for index, row in data.iterrows():
        species = row["Species"]
        g__,f__,o__,c__,p__,d__ = species_lineage_dict[species.lower()]
        taxa = p__ # taxid to split barplot on
        taxa_domain_dict[p__] = d__

        z_score = row["SINGLE: Z Score"]
        rnaclust_sws = row["RNACLUST: Num SWs 8_10"]
        rnaclust_covariation = row["RNACLUST: Covariation"]

        singleseq_call = None
        if math.isnan(z_score) == False:
            if float(z_score) >= singleseq_threshold:
                singleseq_call = True
            else:
                singleseq_call = False
        
        rnaclust_call = None
        if math.isnan(rnaclust_sws) == False:
            if float(rnaclust_sws) >= rnaclust_sws_threshold:
                rnaclust_call = True
            else:
                rnaclust_call = False

        if singleseq_call == True or rnaclust_call == True:
            composite_call = True
        elif singleseq_call == None or rnaclust_call == None:
            composite_call = None
        else:
            composite_call = False
        
        if taxa not in tax_signals_dict:
            tax_signals_dict[taxa] = []
        tax_signals_dict[taxa].append(composite_call)
    
    y_max = 1.0

    gen_barplot_taxonomy(plt_loc = plot_loc + rrna_type + "-taxonomy-fractstems.png", \
                    title = "Phyla Fraction of " + rrna_type + " LTs with Stem", \
                    xname = "Phyla", \
                    yname = "Fraction of LTs with Stem", \
                    tax_signals_dict = tax_signals_dict, \
                    taxa_domain_dict = taxa_domain_dict)
    
    # get per-species fractions
    accs, stem_predicts = [], []
    for acc in data["Acc"]:
        accs.append(acc)
    for lt_stem_call in data["Both False?"]:
        lt_stem_call = bool(lt_stem_call)
        if lt_stem_call == False:
            stem_predicts.append(True)
        elif lt_stem_call == True:
            stem_predicts.append(False)
        else:
            raise Exception("unknown lt stem call:", lt_stem_call)

    # get acc - stem/no stem counts
    acc_stem_fracts = {}
    for i in range(len(accs)):
        acc = accs[i]
        predict = stem_predicts[i]

        if acc not in acc_stem_fracts:
            acc_stem_fracts[acc] = [0, 0]

        if predict == True:
            acc_stem_fracts[acc][0] += 1
        else:
            acc_stem_fracts[acc][1] += 1

    # get data to plot fracts of each acc that contain stem
    acc_fract_stems = []
    num_0_fract, species_w_0_fract = 0, []
    num_not0_not1_fract = 0
    for acc in acc_stem_fracts:
        num_with, num_without = acc_stem_fracts[acc]
        fract = num_with / (num_with + num_without)

        acc_fract_stems.append(fract)
        if fract == 0.0:
            num_0_fract += 1
            species_w_0_fract.append(acc)
        elif fract != 0.0 and fract != 1.0:
            num_not0_not1_fract += 1

    print("# of species with 0 fract LT:", num_0_fract)
    print("# of species with NON 0 and NON 1 fract LT:", num_not0_not1_fract)

    x_min, x_max = 0, 1.0
    y_max = 5000
    bins = np.linspace(x_min, x_max, num=6)

    gen_histogram_cut_axes(plt_loc = plot_loc + rrna_type + "-species-lt-fract.png", \
                    title = rrna_type, \
                    xname = "Fraction of Species LTs with Stem", \
                    yname = "Number of Species", \
                    data = acc_fract_stems, \
                    bins = bins, \
                    x_lims = [x_min, x_max], \
                    y_lims = [0, y_max])

    # Look at # of LTs per species, computing enrichment ratio of those with LT stem vs not
    acc_num_lts_per_species = {}
    for i in range(len(accs)):
        acc = accs[i]

        if acc not in acc_num_lts_per_species:
            acc_num_lts_per_species[acc] = 0
        acc_num_lts_per_species[acc] += 1

    # for ALL species, compute the normalized fraction of LTs per species
    tot_num_species = len(acc_num_lts_per_species)
    number_lt_per_species_ALL = {}
    for i in range(11):
        number_lt_per_species_ALL[i] = 0.0

    for acc in acc_num_lts_per_species:
        num_lts = acc_num_lts_per_species[acc]

        if num_lts not in number_lt_per_species_ALL:
            continue

        number_lt_per_species_ALL[num_lts] += 1 / tot_num_species

    # Now, just for species whose fract LT = 0, compute normalized fract
    num_species_0fract = len(species_w_0_fract)
    number_lt_per_species_FRACT0 = {}
    for i in range(11):
        number_lt_per_species_FRACT0[i] = 0.0

    num_lt_from_0fract_species = 0
    for acc in species_w_0_fract:
        num_lts = acc_num_lts_per_species[acc]
        num_lt_from_0fract_species += num_lts

        if num_lts not in number_lt_per_species_FRACT0:
            continue

        number_lt_per_species_FRACT0[num_lts] += 1 / num_species_0fract    

    print("# of LTs from species with fract 0 LT stems:", num_lt_from_0fract_species)

    # We want to create a barplot of the ratio
    x_data, y_data = [], []
    for i in range(1, 11):
        x_data.append(i)

        norm_fract_0 = number_lt_per_species_FRACT0[i]
        norm_fract_all = number_lt_per_species_ALL[i]

        ratio = norm_fract_0 / norm_fract_all

        y_data.append(ratio)
   
    gen_barplot(plt_loc = plot_loc + rrna_type + "-species-numlt-enrichment.png", \
                    title = rrna_type, \
                    xname = "Number of LTs per Species", \
                    yname = "Enrichment Ratio", \
                    x_data = x_data, \
                    y_data = y_data)


    # determine which seqs have a viable LT
    guids_w_viable_lt = []
    with open(rrna_type + "-seq-clusters.txt", "r") as f:
        for row in f:
            prior = row.split("_")[0]
            guid = row.split(prior + "_")[1].replace("\n", "").lower()
            guids_w_viable_lt.append(guid)
    f.close()

    # Next, look at the # of Unique LTs per species to study enrichment ratio at fract = 0
    rrna_lt_cluster = []
    lt_cluster_accs = []
    acc_number_lts_lookup_dict_covariation = {}
    num_missing = 0
    for i in range(len(data["RNACLUST: Cluster ID"])):
        cluster = data["RNACLUST: Cluster ID"][i]
        
        species = data["Species"][i]
        counter = data["Counter"][i]
        guid = species + "-" + str(counter)
        guid = guid.lower().replace(" ","_")
        
        if guid not in guids_w_viable_lt:
            num_missing += 1
            cluster = "NA"
        else:
            rrna_lt_cluster.append(cluster)
            lt_cluster_accs.append(accs[i])

            acc = accs[i]
            if acc not in acc_number_lts_lookup_dict_covariation:
                acc_number_lts_lookup_dict_covariation[acc] = 0
            acc_number_lts_lookup_dict_covariation[acc] += 1

    print("num LTs without sufficient length:", num_missing)
    acc_lt_clusters = {}
    for i in range(len(lt_cluster_accs)):
        acc = lt_cluster_accs[i]
        cluster = rrna_lt_cluster[i]

        if acc not in acc_lt_clusters:
            acc_lt_clusters[acc] = []

        # If no RNAClust cluster, then we add one to the # of unique LT clusters
        if type(cluster) != str and math.isnan(cluster) == True:
            acc_lt_clusters[acc].append(cluster)
        else:
            # but if there is a cluster, only want to count unique clusters
            if cluster not in acc_lt_clusters[acc]:
                acc_lt_clusters[acc].append(cluster)

    acc_num_uniq_clusters = {}
    for acc in acc_lt_clusters:
        clusters = acc_lt_clusters[acc]
        num_clusters = len(clusters)
        acc_num_uniq_clusters[acc] = num_clusters

    # generate colored scatterplot
    gen_colored_scatterplot(acc_num_uniq_clusters, acc_number_lts_lookup_dict_covariation, \
                             acc_stem_fracts, acc_num_lts_per_species, \
                             plt_loc = plot_loc + rrna_type + "-species-uniq-cluster-analysis.png", \
                             title = rrna_type)


    # Now, we want to plot histograms of the z-score = 1 data
    shuffled_mean, shuffled_stdev = [], []
    for mean in data['SINGLE: Mean Scrambled']:
        if math.isnan(mean) == False:
            shuffled_mean.append(float(mean))
    for stdev in data['SINGLE: StDev Scrambled']:
        if math.isnan(stdev) == False:
            shuffled_stdev.append(float(stdev))

    shuffled_z_scores = []
    for i in range(len(shuffled_mean)):
        mean = shuffled_mean[i]
        stdev = shuffled_stdev[i]

        zscore = mean + stdev

        shuffled_z_scores.append(zscore)

    x_min, x_max = 0, 50
    y_max = 1400

    bins = np.linspace(x_min, x_max, num=int(((x_max - x_min) / 1.0)+1))
    gen_histogram(plt_loc = plot_loc + rrna_type + "-shuffled-zscore1-analysis.png", \
                    title = rrna_type, \
                    xname = "Number of SSWs at Z-score = 1", \
                    yname = "Number of LTs", \
                    data = shuffled_z_scores, \
                    bins = bins, \
                    x_lims = [x_min, x_max], \
                    y_lims = [0, y_max], \
                    vert_lines = [15], use_weights = False, use_title = True)