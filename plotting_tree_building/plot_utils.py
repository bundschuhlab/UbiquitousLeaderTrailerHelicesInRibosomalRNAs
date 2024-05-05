import matplotlib.pyplot as plt
import numpy as np
import csv

plt.rcParams["font.family"] = "Arial"

def gen_histogram(plt_loc, title, xname, yname, \
                    data, bins, x_lims, y_lims, \
                    vert_lines, use_weights = True, \
                    use_title = False):
    """
    """
    if use_weights == True:
        weights = np.ones_like(data)/float(len(data))
        plt.hist(data, bins = bins, weights = weights)
    else:
        plt.hist(data, bins = bins)

    plt.xlabel(xname, fontsize=24)
    plt.ylabel(yname, fontsize=24)

    if use_title == True:
        plt.title(title, fontsize=26)

    plt.xlim(x_lims)
    plt.ylim(y_lims)

    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)

    for vert in vert_lines:
        plt.axvline(vert, color="k", linestyle="--")

    plt.savefig(plt_loc, bbox_inches="tight")
    plt.close()

    return


def gen_barplot(plt_loc, title, xname, yname, \
                    x_data, y_data):
    """
    """
    plt.bar(x_data, y_data)

    plt.xlabel(xname, fontsize=24)
    plt.ylabel(yname, fontsize=24)
    plt.title(title, fontsize=22)

    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)

    plt.axhline(y=1, color="k", linestyle="--")

    plt.savefig(plt_loc, bbox_inches="tight")
    plt.close()

    return



def gen_histogram_cut_axes(plt_loc, title, xname, yname, \
                    data, bins, x_lims, y_lims):
    """
    """
    f, (ax1, ax2) = plt.subplots(2, 1, sharex=True)

    ax1.hist(data, bins = bins)
    ax2.hist(data, bins = bins)

    ax1.set_ylim(4100, 4700)  # LT = 1 only
    ax2.set_ylim(0, 400)  # rest of the data

    ax1.spines['bottom'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    #ax1.xaxis.tick_top()
    ax1.tick_params(axis='x', which='both',bottom=False,top=False)
    ax1.tick_params(labeltop=False)  # don't put tick labels at the top
    ax2.tick_params(labeltop=False)
    ax2.xaxis.tick_bottom()

    ax1.tick_params(axis='both', which='major', labelsize=18, width=2.5, length=10)
    ax1.tick_params(axis='both', which='minor', labelsize=18, width=2.5, length=10)
    ax2.tick_params(axis='both', which='major', labelsize=18, width=2.5, length=10)
    ax2.tick_params(axis='both', which='minor', labelsize=18, width=2.5, length=10)

    f.suptitle(title, fontsize=32)
    f.text(0.5, -0.05, xname, ha='center', fontsize=24)
    f.text(-0.05, 0.5, yname, va='center', rotation='vertical', fontsize=24)

    plt.savefig(plt_loc, bbox_inches="tight")
    plt.close()

    return


def gen_barplot_taxonomy(plt_loc, title, xname, yname, \
                            tax_signals_dict, taxa_domain_dict):
    """
    """
    # Order the taxids from most to least abundant
    taxid_counter = []
    tax_fract_dict = {}
    for taxid in tax_signals_dict:
        calls = tax_signals_dict[taxid]
        taxid_counter.append([taxid, len(calls)])

        num_stemtrue = calls.count(True)
        fract = num_stemtrue / len(calls)

        tax_fract_dict[taxid] = fract

    min_num_to_bin = 25
    overflow_taxids = {}

    xtick_colors = []
    ordered_taxids = []
    xaxis_labels = []
    taxid_counter = sorted(taxid_counter, key = lambda x:x[1], reverse=True)
    for taxid, count in taxid_counter:
        if count < min_num_to_bin:
            overflow_taxids[taxid] = count
            continue

        ordered_taxids.append(taxid)

        phyla_name = taxid.split("p__")[1]
        phyla_name = phyla_name[0].upper() + phyla_name[1:]
        xaxis_labels.append(phyla_name + " (" + str(count) + ")")

        domain = taxa_domain_dict[taxid]
        if domain == "d__bacteria":
            xtick_colors.append("red") #red
        elif domain == "d__archaea":
            xtick_colors.append("maroon") #brown
    
    # deal with overflow taxids
    overflow_tot_fract, overflow_tot_count = 0.0, 0
    for taxid in overflow_taxids:
        count = overflow_taxids[taxid]
        fract = tax_fract_dict[taxid]

        overflow_tot_fract += fract * count
        overflow_tot_count += count
    
    avg_overflow_fract = overflow_tot_fract / overflow_tot_count

    ordered_taxids.append("Other Phyla")
    tax_fract_dict["Other Phyla"] = avg_overflow_fract
    xaxis_labels.append("Other Phyla (" + str(overflow_tot_count) + ")")
    xtick_colors.append("k") # black

    # now, write out the barplot
    x_data, y_data = [], []
    counter = 0
    for taxid in ordered_taxids:
        counter += 1
        fract = tax_fract_dict[taxid]

        x_data.append(counter)
        y_data.append(fract)

    fig = plt.figure(figsize = [28, 8])
    ax = fig.add_subplot(111)

    plt.bar(x_data, y_data)

    #plt.xlabel(xname, fontsize=24)
    plt.ylabel(yname, fontsize=32)
    #plt.title(title, fontsize=34)

    plt.ylim([0, 1])
    plt.xticks(x_data, xaxis_labels, rotation=45, ha='right', fontsize = 38)

    for xtick, color in zip(ax.get_xticklabels(), xtick_colors):
        xtick.set_color(color)

    plt.yticks(fontsize=32)

    plt.savefig(plt_loc, bbox_inches="tight")
    plt.close()   

    return


def gen_colored_scatterplot(acc_num_uniq_clusters, acc_number_lts_lookup_dict_covariation, \
                             acc_stem_fracts, acc_num_lts_per_species, plt_loc, title):
    """
    """
    x_pnts, y_pnts = [], []
    tot_num_points, points_outside_window = 0, 0

    for acc in acc_number_lts_lookup_dict_covariation:
        num_lts = acc_number_lts_lookup_dict_covariation[acc]
        num_uniq_clusters = acc_num_uniq_clusters[acc]

        # skip outliers - do this here for plotting boxes
        tot_num_points += 1
        if num_lts > 15 or num_uniq_clusters > 10:
            points_outside_window += 1
            continue

        x_pnts.append(num_lts)
        y_pnts.append(num_uniq_clusters)

    print("percent of data falling outside heatmap:", 100 * (points_outside_window / tot_num_points))

    max_x = max(x_pnts)
    max_y = max(y_pnts)

    x_bins, y_bins = [], []
    x_ticks, y_ticks = [], []
    for i in range(0, max_x):
        x_bins.append(i + 0.5)
        x_ticks.append(i+1)
    for i in range(0, max_y):
        y_bins.append(i + 0.5)
        y_ticks.append(i+1)

    #heatmap, xedges, yedges = np.histogram2d(x_pnts, y_pnts, bins=(max_x, max_y))
    heatmap, xedges, yedges = np.histogram2d(x_pnts, y_pnts, bins=(x_bins, y_bins))
    extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]

    fig, ax = plt.subplots()
    plt.clf()
    im = plt.imshow(heatmap.T, extent=extent, origin='lower')
    im.set_cmap('gist_heat_r')

    cbar = ax.figure.colorbar(im, ax = ax, location='right', shrink=0.6, anchor = (1.5,0.5))
    cbar.ax.set_ylabel("Number of Species", rotation = -90, va = "bottom")

    plt.xlabel("Number of LTs")
    plt.ylabel("Number of Unique LT Clusters")

    plt.xticks(x_ticks)
    plt.yticks(y_ticks)

    plt.title(title)

    plt.savefig(plt_loc, bbox_inches="tight")
    plt.close()

    # write out CSV of #s
    with open(plt_loc.split(".png")[0] + "-data.csv", "w") as f:
        out = csv.writer(f)
        out.writerow(["acc", "Total Number of LTs", \
                        "Number LTs with predicted stem", \
                        "Number LTs without predicted stem", \
                        "Number of LTs with Sufficient nts for Covariation Method", \
                        "Number of Unique LT Clusters from Covariation Method"])

        for acc in acc_num_lts_per_species:
            num_lts = acc_num_lts_per_species[acc]
            num_with, num_without = acc_stem_fracts[acc]

            if acc in acc_number_lts_lookup_dict_covariation:
                num_covariation_lts = acc_number_lts_lookup_dict_covariation[acc]
                num_uniq_clusters = acc_num_uniq_clusters[acc]
            else:
                num_covariation_lts = 0
                num_uniq_clusters = 0

            out.writerow([acc, num_lts, num_with, num_without, num_covariation_lts, num_uniq_clusters])
    f.close()

    return