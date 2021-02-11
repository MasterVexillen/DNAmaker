#!/usr/bin/env python3

"""
DNAmaker.main.py

Author: Neville Yee
Version: 0.1
Date: 11-Feb-2021
"""

import os
import random
import numpy as np

def get_info():
    """
    Get DNA info from user
    """

    rand_seq = input("Generate random sequence (Y/N)? ").lower()
    if rand_seq == 'y':
        rand_seq = True
    else:
        rand_seq = False

    form = input("Conformation of DNA (A/B): ").lower()
    assert (form in ['a', 'b']), 'Input not valid conformation'

    # convert form to nab code
    if form == 'a':
        form = 'adna'
    elif form == 'b':
        form = 'abdna'

    allowed_bases = ['A', 'T', 'C', 'G']
    if rand_seq:
        nbp = int(input("Number of base pairs? "))
        assert (isinstance(nbp, int)), 'Number of base pairs must be an integer'

        seed = int(input("Random seed (-1 if different result every time): "))
        base_rand = []

        for bp in range(nbp):
            if seed != -1:
                random.seed(seed+bp)
            base_rand.append(random.random())

        cg_content = float(input("Projected CG content (%)? ")) / 100

        tc_thres = 1 - cg_content
        at_thres = 0.5 * tc_thres
        cg_thres = tc_thres + 0.5*cg_content

        bins = [at_thres, tc_thres, cg_thres]
        binned_base_rand = np.digitize(base_rand, bins=bins)

        sequence_list = [allowed_bases[index] for index in binned_base_rand]
        final_sequence = ''.join(sequence_list)

    else:
        final_sequence = input("Desired sequence? ").upper()
        # Check input
        for base in final_sequence:
            if not base in allowed_bases:
                raise ValueError("Non DNA base entered. Typo?")

    print("DNA sequence created successfully. Now write nab file... \n")
    return (form, final_sequence)


def write_nab(dna_spec_in=None, nab_name=None, pdb_name=None):
    """
    Subroutine to write nab file
    """

    (form, sequence) = dna_spec_in
    # Define molecule name in nab file - could be kept in the black box WLOG
    mol_name = 'm'

    mol_name_line = "molecule {};".format(mol_name)
    mol_def_line = '{} = fd_helix( "{}", "{}", "dna" );'.format(mol_name, form, sequence)
    pdb_def_line = 'putpdb( "{}", {}, "-wwpdb" );'.format(pdb_name, mol_name)

    with open(nab_file, 'w') as nab:
        nab.write(mol_name_line + "\n")
        nab.write(mol_def_line + "\n")
        nab.write(pdb_def_line)


if __name__ == '__main__':
    # Generate DNA specs
    dna_spec = get_info()

    # Generate NAB file for processing
    nab_file = input('Name (less ".nab") of nab file? ') + '.nab'
    pdb_file = input('Name (less ".pdb") of pdb file? ') + '.pdb'
    write_nab(dna_spec, nab_file, pdb_file)

    # Check AMBERHOME properly set
    amberhome = os.getenv("AMBERHOME")
    assert (len(amberhome) > 0), \
        'Error: length of $AMBERHOME == 0. Environment variable not set?'

    command = "$AMBERHOME/bin/nab {}".format(nab_file)
    os.system(command)
    os.system("./a.out")

    # Clear unwanted junk
    os.system("rm *.c *.out")
