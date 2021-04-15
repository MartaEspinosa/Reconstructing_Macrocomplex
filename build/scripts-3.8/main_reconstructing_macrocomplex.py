import argparse
import os
import sys
import logging
import gzip
import shutil
from functions_reconstructing_macrocomplex import *
from Bio.PDB import *

parser = argparse.ArgumentParser(description = """This program models macrocomplex structures of biomolecules formed by proteins. The macrocomplex is build from pairwaise interactions.""")

parser.add_argument('-i','--input-directory',
                    required=True,
                    dest="input_path",
                    action="store",
                    help="Directory containig the input structure files.")

parser.add_argument('-s', '--stoichiometry',
                    default=None,
                    dest="stoichiometry",
                    action="store",
                    help="Path containing the stoichiometry of the complex.")

parser.add_argument('-o','--output-directory',
                    dest="output_path",
                    action="store",
                    type=str,
                    default="output",
                    help="Directory where the output will be saved.")

parser.add_argument('-f','--force',
                    default=False,
                    dest="force",
                    action="store",
                    help="Checks the output-directory doesn't exist before the application is executed.")

parser.add_argument('-v','--verbose',
                    default=False,
                    dest="verbose",
                    action="store_true",
                    help="Returns the progression log of the program execution printed in standard error.")

parser.add_argument('-ofn', '--output_filename',
                    dest = "output_filename",
                    action = "store",
                    default = "final_macrocomplex",
                    help = "Optional. Defines the name of the output file where the final model should be stored.")

arguments = parser.parse_args()


if __name__=="__main__":

    absolute_path = os.getcwd()

    try:
        os.mkdir(arguments.output_path)
        subdirs = ['structures', 'analysis']
        for folder in subdirs:
            os.mkdir(os.path.join(arguments.output_path, folder))
    except FileExistsError:
        if arguments.force == False:
            print (f"The {arguments.output_path} already exists.\n")
            sys.exit(-1)
        elif arguments.force == True:
            shutil.rmtree(arguments.output_path)
            os.mkdir(arguments.output_path)
            subdirs = ['structures', 'analysis']
            for folder in subdirs:
                folder = "/" + folder
                os.mkdir(os.path.join(arguments.output_path, folder))

    # Create the logging file
    logging_storage = arguments.output_path + "/analysis/macrocomplex_builder.log"
    logging.basicConfig(filename = logging_storage, format = '%(levelname)s:%(message)s', level = logging.DEBUG, filemode = "w")	#move to output folder
    logging.debug('...STARTING...')

    if arguments.verbose:
        logging.getLogger().addHandler(logging.StreamHandler())

    # Get the list of pdb files in input directory provided
    try:
        list_files = only_pdb_files(arguments.input_path)
    except NotADirectoryError:
        logging.info("Input path provided could not be found as a directory. Aborting program. Please try again.")
        sys.exit(-1)

    # Get the structures of the PDB files
    pdb_structure, molecule_type = read_pdb_files(list_files)

    os.chdir(absolute_path) # get back to main working directory

    # Get the stoichiometry from the file if provided by the user
    if arguments.stoichiometry:
        logging.info("\nStoichiometry file %s was found. This will be used to build the complex.\n" %(arguments.stoichiometry))
        try:
            stoichiometry_file = {}
            with open(arguments.stoichiometry, "r") as file:
                for line in file:
                    stoich_info = line.strip().split(":")
                    stoichiometry_file[stoich_info[0]] = int(stoich_info[1])
                    print(stoichiometry_file)
        except FileNotFoundError:
            logging.info("\nWe could not find the stoichiometry file provided. Aborting program. Please try again.")
            exit()
    else:
        logging.info("\nStoichiometry was not provided. The program will use 1 for each chain as default value.\n")
        stoichiometry_file = {}
        for file in pdb_structure:
            stoichiometry_file[file] = 1

    # Analyse protein-protein interactions to build complex
    if molecule_type == "PROT":
        logging.info("The files provided contain Protein-Protein Interactions.\n")

        # get the first file as reference and put it to the end of the list
        id_list = list(pdb_structure.keys())
        print(id_list)
        reference_id = id_list.pop(0)
        id_list.append(reference_id)
        current_stoichiometry = {reference_id:1}

        # Superimpose C-alphas of those of chains with alignment higher than 0.95
        reference_structure = pdb_structure[reference_id] # get the structure of the reference
        num_chains = 2  # initial number of chains of the complex

        while (num_chains <= sum(list(stoichiometry_file.values()))):    # while the complex has less or same chains as the stoichiometry provided
            sample_id = id_list.pop(0)       # get interaction to compare with reference and superimpose
            if sample_id not in stoichiometry_file:
                sample_id = ""
                continue
            if not sample_id in current_stoichiometry:   # initialise the count for the structure that is currently being analysed
                current_stoichiometry[sample_id] = 0

            sample_structure = pdb_structure[sample_id]

            try:
                superimposition, RMSD = superimpose_chains(reference_structure, sample_structure) # superimpose. We get the superimposition with their best RMSD
            except:
                current_stoichiometry[sample_id] += 1
                num_chains += 1

            # If no superimposition
            if bool(superimposition) == False:
                id_list.append(sample_id)          # if the sample could not be superimposed, add it to the end of the list to analysed it again later
                continue

            # If there are superimpositions made
            for tuple_chains, super_object in superimposition:
                logging.info("Chain %s of the reference structure and chain %s of the sample structure have been superimposed with a RMSD of %.3f.\n" %(tuple_chains[0], tuple_chains[1], RMSD))
                chain_to_add = [chain for chain in sample_structure.get_chains() if chain.id != tuple_chains[1]][0]
                super_object.apply(chain_to_add.get_atoms()) # apply rotation matrix to moving structure
                reference_structure = look_for_clashes (reference_structure, chain_to_add)

                num_chains += 1

            # If the structure is not the same as the stoichoimetry provided, add it to the end of the list to analyse and try to superimpose again later
            if current_stoichiometry[sample_id] != stoichiometry_file[sample_id]:
                id_list.append(sample_id)


        store_output (reference_structure[0], arguments.output_path, arguments.output_filename)
        logging.info("The macrocomplex built has %d chains." %(reference_structure[0].__len__()))
