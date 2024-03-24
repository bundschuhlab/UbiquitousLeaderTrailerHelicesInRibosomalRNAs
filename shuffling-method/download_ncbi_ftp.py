import os
import gtdb_tree


def load_dump_taxdmp(address, loc):
    """
    """
    os.system("wget " + address)
    os.system("mv " + address.split("/")[-1] + " " + loc)
    if ".tar.gz" in address.split("/")[-1]:
        os.system("tar -xvf " + loc + address.split("/")[-1] + " -C " + loc)
    elif ".gz" in address.split("/")[-1]:
        os.system("gunzip " + loc + address.split("/")[-1])		   

    return


def load_ncbi_data(temp_data_root, taxdmp_loc, genbank_assembly_summary_loc, \
                        refseq_assembly_summary_loc, force_redownload):
    """
    """
    if force_redownload == True and os.path.isdir(temp_data_root) == True:
        os.system("rm -rf " + temp_data_root)

    if os.path.isdir(temp_data_root) == False:
        os.system("mkdir " + temp_data_root)
        os.system("mkdir " + taxdmp_loc)

        # download NCBI taxdmp
        gtdb_bact_taxdmp_address = "https://data.gtdb.ecogenomic.org/releases/latest/bac120_taxonomy.tsv.gz"
        load_dump_taxdmp(gtdb_bact_taxdmp_address, taxdmp_loc)
        gtdb_arch_taxdmp_address = "https://data.gtdb.ecogenomic.org/releases/latest/ar53_taxonomy.tsv.gz"
        load_dump_taxdmp(gtdb_arch_taxdmp_address, taxdmp_loc)
        
        # download NCBI assembly summary 
        ncbi_ftp_genbank_address = "https://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt"
        os.system("wget " + ncbi_ftp_genbank_address)
        os.system("mv assembly_summary_genbank.txt " + genbank_assembly_summary_loc)

        ncbi_ftp_refseq_address = "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt"
        os.system("wget " + ncbi_ftp_refseq_address)
        os.system("mv assembly_summary_refseq.txt " + refseq_assembly_summary_loc)

        # download sp_clusters and bact metadata for genome gathering
        gtdb_sp_clusters_address = "https://data.gtdb.ecogenomic.org/releases/latest/auxillary_files/sp_clusters.tsv"
        load_dump_taxdmp(gtdb_sp_clusters_address, temp_data_root)
        gtdb_bact_metadata = "https://data.gtdb.ecogenomic.org/releases/latest/bac120_metadata.tar.gz"
        load_dump_taxdmp(gtdb_bact_metadata, temp_data_root)

    print("loading in taxtree")
    tree = gtdb_tree.TaxonomyTree(taxdmp_loc)
    print("taxtree loaded\n\n\n")	

    return tree