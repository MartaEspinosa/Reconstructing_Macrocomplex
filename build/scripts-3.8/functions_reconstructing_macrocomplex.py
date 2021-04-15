import argparse
import os
import sys
import logging
import string
import gzip
import shutil
import random
from Bio.PDB import *
from Bio import pairwise2

def only_pdb_files(input_path):
    """Checks the input files are pdb files. No other format is accepted."""

    list_files = []
    error_files = []

    for file in os.listdir(input_path):
        if file.endswith(".pdb") or file.endswith("pdb.gz"):    # compressed files are also accepted
            list_files.append(file)
        else:
            error_files.append(file)

    logging.info("This program is analyzing %i from %i files.\n" %(len(list_files), (len(list_files) + len(error_files))))

    if len(error_files) > 0:
        logging.info("There are some files with wrong format.\n")

    else:
        os.chdir(input_path)
        return list_files


def read_pdb_files(pdb_files):
    """Get the input data of the pdb files provided.

	From these input files, heteroatoms are removed .

	The input files may be compressed (gz) or not. Files in different format as pdb will not be considered.

	The function returns a dictionary with the chain ids and the structure depending on the type of interaction.
    """

    dict_PPI = {}
    dict_DNA_P = {}
    dict_RNA_P = {}

    pdb_parser = PDBParser(PERMISSIVE=1, QUIET=True)

    for file in pdb_files:
        if file.endswith(".gz"):
            with gzip.open(file,"rt") as gzip_file:
                id = file[:-7]
                structure = pdb_parser.get_structure(id, gzip_file)

        if file.endswith(".pdb"):
            id = file[:-4]
            structure = pdb_parser.get_structure(id, file)

        chains_ids = ''.join([chain.id for chain in structure.get_chains()])
        chains = []
        CA_chains = 0

        files_not_well_named = []

        for chain in structure.get_chains():
            if chain.id == id[-1] or chain.id == id[-3]:
                continue
            else:
                files_not_well_named.append(id)	# bad naming structure

        if len(files_not_well_named) != 0:
            logging.info("The naming structure and the chains in the following structure file don't coincide:")
            for wrong_file in files_not_well_named:
                logging.info("%s" %(wrong_file))
            continue

        het_to_remove = []
        for chain in structure.get_chains():
            for residue in chain:
                if residue in chain:
                    if residue.id[0] != ' ':
                        het_to_remove.append(residue.id)

            for residue in het_to_remove:
                chain.detach_child(residue) # remove heteroatoms
            chains.append(chain)

            if len(next(chain.get_residues()).get_resname()) < 3:  # if the structure is too short or the residue name is not the expected, that residue will not be considered
                structure[0].detach_child(chain.id)
            else:
                CA_chains += 1          # number of structures


        key_chain = [x for x in structure.get_chains()][1]  # is this a PPI?

        chain_molecule = get_atoms_and_molecule(key_chain)[1]

        if chain_molecule == "PROT":
            if CA_chains != 2:  # we need binary interactions
                logging.info("The '%s' file does not have the right input format to build a complex." %(file))

                continue
            dict_PPI[id] = structure

        else:
            if chain_molecule == "DNA":
                dict_DNA_P[id] = structure
            elif chain_molecule == "RNA":
                dict_RNA_P[id] = structure

    if bool(dict_DNA_P) == True or bool(dict_RNA_P) == True:
        logging.info("The files to analise contain interactions with nucleic acids. This program do not analise those situations. Further versions will do so.\n")

    return dict_PPI, "PROT"


def get_atoms_and_molecule(chain):
    """This function compares the atom and molecule type of a chain of a PDB structure object.

    The atoms can be CA in case of proteins and C4' in case of nucleic acids and only the same type of molecule will be superimposed later.

    It returns a list with atoms and another one with molecules indicating whether the chain is DNA, RNA or PROTEIN.
    Only superimpositions of the same type of molecule are considered.
    """

    nucleic_acids = ['DA','DT','DC','DG','DI','A','U','C','G','I']
    RNA = ['A','U','C','G','I']
    DNA = ['DA','DT','DC','DG','DI']
    atoms = []
    molecule = ""

    for residue in chain:
        residue_name = residue.get_resname().strip()
        if residue.get_id()[0] == " " and residue_name not in nucleic_acids:
            if 'CA' in residue:
                atoms.append(residue['CA'])
                molecule="PROT"

        elif residue.get_id()[0] == " " and residue_name in nucleic_acids:
            if residue_name in DNA:
                molecule = 'DNA'

            elif residue_name in RNA:
                molecule = 'RNA'
            atoms.append(residue['C4\''])

    return atoms, molecule      # returns the list of alpha carbon atoms and the molecule types

def id_construction(existing_ids):
    """This function returns a new ID for the new chain to be added to the complex.

    It generates a a single character random ID based on uppercase, lowercase or digits.

    The input argument is a list with all the existing ids of the macrocomplex.
    The function returns a new id for the chain to add.
    """

    alphabet = string.ascii_letters + string.digits

    new_id = random.choice(alphabet)
    if new_id not in existing_ids:
        return new_id


def alignment_chains(chain1, chain2):
    """This function aligns two chains.

    The input arguments are the two chains to be aligned.
    The function returns the score of the alignment.
    """

    CA = CaPPBuilder()  # Use CA–CA distance to find polypeptides.

    CA_chain1 = CA.build_peptides(chain1)
    CA_chain1 = CA_chain1[0].get_sequence()

    CA_chain2 = CA.build_peptides(chain2)
    CA_chain2 = CA_chain2[0].get_sequence()

    alignment = pairwise2.align.globalxx(CA_chain1,CA_chain2)
    score = alignment[0][2]/max(len(CA_chain1), len(CA_chain2))

    return score


def superimpose_chains(reference_structure,sample_structure):
    """This function superimposes every combination of pairs of chains and calculates the RMSD.

    Only superimpositions with a RMSD below 2 (RMSD threshold) are considered as good.

	The input arguments are the reference and the sample structure objects.
		- reference : PDB structure object where the macrocomplex is being built.
		- sample : PDB structure object that should be superimposed and whose chains may be added to the complex.

	The function returns a dictionary with a tuple as key with the reference and sample chain that have tried to superimpose and the corresponding RMSD of the superimposition instance.
    """
    threshold = 2   # we accept at least a RMSD of 2 to consider a good superimposition
    superimpositions = {}   # will keep track of the superimpositions performed
    best_RMSD = ""
    reference_chains = [x for x in reference_structure.get_chains()]
    sample_chains = [x for x in sample_structure.get_chains()]
    superimposed = Superimposer()

    for ref_chain in reference_chains:
        for s_chain in sample_chains:
            if alignment_chains(ref_chain,s_chain) > 0.95: # we accept at least a score of 0.95
                reference_atoms, reference_molecule = get_atoms_and_molecule(ref_chain)
                sample_atoms, sample_molecule = get_atoms_and_molecule(s_chain)
                superimposed.set_atoms(reference_atoms, sample_atoms)  # rotation and translation matrix
                RMSD = superimposed.rms # get RMSD

                if RMSD < threshold:
                    if not best_RMSD or RMSD < best_RMSD:
                        best_RMSD = RMSD
                    superimpositions[(ref_chain.id,s_chain.id)] = superimposed

    if bool(superimpositions) == True:  # if there are superimpositions
        superimpositions = sorted(superimpositions.items(), key=lambda x:x[1].rms)    # sort the superimpositions by RMSD

        return superimpositions, best_RMSD


def look_for_clashes(reference_structure, chain_to_add):
    """This function makes Neighbor Search to look for clashes between the atoms of the chain to add and the atoms of the references.

    The input arguments are the reference structure and the chain to add.

    The function returns the reference structure with the chain added if the number of clashes allow so.
    """

    reference_atoms = []
    sample_atoms, sample_molecule = get_atoms_and_molecule(chain_to_add)

    for chain in reference_structure.get_chains():  # all atom positions in the current reference structure
        reference_atoms.extend(get_atoms_and_molecule(chain)[0])

    Neighbor = NeighborSearch(reference_atoms)
    num_clashes = 0
    for atom in sample_atoms:     # possible clashes between the atoms of the chain to add and the atoms in the model
        atoms_clashed = Neighbor.search(atom.get_coord(),5) # produces a Neighbor search that returns all atoms/residues/chains/models/structures that have at least one atom within radius of center

        if len(atoms_clashed) > 0:
            num_clashes += len(atoms_clashed)

    if num_clashes < 30:   # if the number of clashes is low enough (we allow less than 30) the chain will be added to the model
        existing_chain_id = [chain.id for chain in reference_structure.get_chains()]
        prev_id = chain_to_add.id
        if chain_to_add.id in existing_chain_id:
            chain_to_add.id= id_construction(existing_chain_id)  # we cannot have the same id for two chains
        reference_structure[0].add(chain_to_add)

        logging.info("Chain %s has become chain '%s' and has been added to the model.\n" %(prev_id, chain_to_add.id))

    return reference_structure


def store_output (model, output_path, filename):
	"""This function saves the macrocomplex in a pdb file.

	The arguments are:
	 	- model : the macrocomplex that should be saved in a .pdb file.
		- output_path : where it should be saved. It will always be saved inside a subfolder named 'structures' of the directory specified as input argument of the program. By default this main directory is 'output'.
		- filename : name of the pdb file generated. If it is not specified as input argument of the program, by default will be 'final_macrocomplex.pdb'.
	"""

	filename = str(filename.strip())
	io = PDBIO()
	io.set_structure(model)
	output_path = output_path + "/structures/" + filename
	io.save(output_path + ".pdb")
