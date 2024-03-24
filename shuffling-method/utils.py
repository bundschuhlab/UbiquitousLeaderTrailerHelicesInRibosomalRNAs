import os


def check_mk_fldrs(fldrs_in):
    """
    """
    for fldr in fldrs_in:
        if os.path.isdir(fldr) == False:
            os.makedirs(fldr)

    return


def check_species_in_clade(species, cl_of_interest, tree):
    """
    """
    parent_lineage = tree.ascend(species.lower())
    species_in_lineage = False

    for node in parent_lineage:
        node_taxid = node.taxid

        if node_taxid == cl_of_interest.lower():
            species_in_lineage = True

    return species_in_lineage


def load_in_species_accs(assembly_data_summary, cl_of_interest, tree):
    """
    """
    species_genome_dict = {}
    genome_species_dict = {}
    with open(assembly_data_summary, "r") as f:
        for row in f:
            row = row.replace("\n","").split("\t")
            species, accs = row

            species_in_lineage = check_species_in_clade(species, cl_of_interest, tree)
            #if species_in_lineage == False:
            #    continue

            accs = accs.split(";")
            species_genome_dict[species] = accs

            for acc in accs:
                genome_species_dict[acc] = species
    f.close()

    return species_genome_dict, genome_species_dict


def load_fa_to_dict(fa_in):
    """
    """
    seqid_dict = {}
    with open(fa_in, "r") as f:
        subj, seq = "", ""
        for row in f:
            if row.startswith(">"):
                if len(subj) != 0:
                    seqid_dict[subj] = seq
                subj = row[1:].replace("\n","").split()[0]
                seq = ""
            else:
                seq = seq + row.replace("\n","").upper()
        
        if len(subj) != 0:
            seqid_dict[subj] = seq
    f.close()

    return seqid_dict