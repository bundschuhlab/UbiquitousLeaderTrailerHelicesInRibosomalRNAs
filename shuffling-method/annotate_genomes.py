import subprocess
import os


def run_barrnap(species_genome_dict, genome_species_dict, \
                    barrnap_exe, assembly_dl_root, \
                    barrnap_annotation_loc, \
                    threads, logging, \
                    cl_of_interest):
    """
    """
    for acc in genome_species_dict:
        fasta_loc = assembly_dl_root + acc + "_genomic.fna"

        # run barrnap
        command = barrnap_exe + \
                    " " + fasta_loc + \
                    " --threads " + str(threads)

        if "archaea" in cl_of_interest.lower():
            command += " --kingdom arc"
        else:
            command += " --kingdom bac"

        command += " > " + barrnap_annotation_loc + acc + "-barrnap.out"
        os.system(command)

    return


def genome_annotation(species_genome_dict, genome_species_dict, \
                        prokka_exe, assembly_dl_root, \
                        prokka_annotation_loc, \
                        threads, logging):
    """
    """
    for acc in genome_species_dict:
        fasta_loc = assembly_dl_root + acc + "_genomic.fna"
        acc_prokka_loc = prokka_annotation_loc + acc + "/"

        # run prokka on the genome
        command = prokka_exe + \
                    " --outdir " + acc_prokka_loc + \
                    " --cpus " + str(threads) + \
                    " --prefix " + acc + \
                    " --force" + \
                    ' "' + fasta_loc + '"'
        print(command)
        os.system(command)
        #run = subprocess.run(command.split(), stdout = subprocess.PIPE)
        #if run.returncode != 0:
        #    logging.write("ERROR\tprokka\t" + acc + "\n")

    return