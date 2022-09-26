import sys
import numpy as np
import pandas as pd
import ast
import json
from typing import Tuple

# read and write  all kinds of file types

def read_traj_file(xyz_file: str) -> Tuple[list, list]:
    ''' Read trajectory file.

    Args:
        xyz_file: path to trajectory file from nanoreactor simulation
                  Format: 1 number of atoms
                          2 TIME: time step
                          3 elem x y z
                          4 elem x y z
                                .
                                .
    Returns:
        natoms: number of atoms in the simulation
        atom_map: order of atomic symbols at each time step
        xyz_list: complete trajectory as list of lists
    '''

    index = -1
    atom_map_total = []
    xyz_list = []
    natoms = 0

    traj_file = open(xyz_file, "r")

    for line_number,line in enumerate(traj_file):
        try:
            vals = line.strip().split()
        except:
            raise Exception ("Failed to parse line...")

        if len(vals) == 1:
            natoms = int(vals[0])
        elif len(vals) == 2:
            timestep = vals[1]
            index += 1
            atom_map_total.append([])
            xyz_list.append([])
        elif len(vals) == 4:
            atom_map_total[index].append(vals[0])
            xyz_list[index].append([float(vals[1]), float(vals[2]), float(vals[3])])

    atom_map = atom_map_total[0]
    traj_file.close()

    return atom_map, xyz_list

def read_bo_file(dat_file: str) -> list:
    ''' Read (Wiberg) bond order (wbo) file.

    Args:
        dat_file: path to bond order file from nanoreactor simulation
                  Format: 1 TIME: time step
                          2 wbo(0,1
                          3 wbo(0,2)
                          4 wbo(0,3)
                                .
                                .
                          only upper triangular (without diagonal elements because equal to 0) is stored to reduce file size
    Returns:
        bond_orders = upper half of the bond order matrix stored as list
    '''

    index = -1
    bond_orders = []
    bo_file = open(dat_file,"r")
    for line_number,line in enumerate(bo_file):
        try:
            vals = line.strip().split()
        except:
            raise Exception ("Failed to parse line...")

        if len(vals) == 2:
            timestep = vals[1]
            index += 1
            bond_orders.append([])
        elif len(vals) == 1:
            bond_orders[index].append(float(vals[0]))

    bo_file.close()
    return bond_orders

def read_frag_file(dat_file: str) -> Tuple[list,list,list]:
    ''' Read fragment file containing indices belonging to found molecules in each step.

    Args:
        dat_file: filepath to file containing (on-the-fly) computed fragments stored as lists of atom indices (starting at 1)
                  Format: 1  time step
                          2 at_index1   at_index2
                          3 at_index3   at_index4   at_index5
                          4 at_index6
                                .
                                .
    Returns:
        timesteps: list of stored time steps in fs
        fragment = list of lists storing the fragment indices (starting with 1) found at each time step
        atom_indices_frag = list of lists storing corresponding list of atom indices for each fragment at every time step
        elem_frag = list of lists storing corresponding elements # do we need this? we have atom_map
    '''

    atom_indices_frag = []
    index = -1

    mols_file = open(dat_file, "r")

    for line in mols_file:
        if line.startswith(' '):
            index += 1
            atom_indices_frag.append([])

        else:
            numbers = line.split()
            atom_indices_frag[index].append(numbers)

    mols_file.close()
    return atom_indices_frag

def write_frag_file(name: str, timestep: float, fragments: list):
    ''' Write fragment file containing indices belonging to found molecules in each step.

    Args:
        timestep: time step in fs
        name: name of file to be written
        fragments: found molecules as lists of atom indices
                  Format: 1  time step
                          2 at_index1   at_index2
                          3 at_index3   at_index4   at_index5
                          4 at_index6
                                .
                                .
    Returns:
        -
    '''

    if timestep == 0.0:
        f = open(name,"w")
    else:
        f = open(name,"a")
    string = str("%20.10e\n") % (timestep)
    f.write(string)
    for i in range(len(fragments)):
        for j in fragments[i]:
            string = str("%i\t") % j
            f.write(string)
        string = str("\n")
        f.write(string)
    f.close()

    return

def read_trafo_df(df_file: str) -> pd.DataFrame:
    ''' Restore data frame containing all necessary information from the post-processing from file.

    Args:
        df_file: filepath to .csv file where data frame was stored (sep: ';')
    
    Returns:
        df : pd.DataFrame object with columns 'Time Step', '# Fragment', '# Atom in Fragment', '# Elem in Fragment', 'XYZ', 'SMILES', 'Molecular Formulas'
    '''

    df = pd.read_csv(df_file, sep = ";")

    new_fragments = []
    for ts in range(len(df['# Fragment'])):
        new_fragments.append(ast.literal_eval(df['# Fragment'][ts]))
    
    new_atom_in_fragments = []
    for ts in range(len(df['# Atom in Fragment'])):
        new_atom_in_fragments.append(ast.literal_eval(df['# Atom in Fragment'][ts]))

    new_SMILES = []
    for ts in range(len(df['SMILES'])):
        new_SMILES.append(ast.literal_eval(df['SMILES'][ts]))

    new_molecular_formulas = []
    string = []
    for ts in range(len(df['Molecular Formulas'])):
        string = df['Molecular Formulas'][ts].replace('\'','')
        string = string.replace("[","")
        string = string.replace("]","")
        string = string.replace(" ","")
        new_molecular_formulas.append(string.split(','))     
    
    df['# Fragment'] = new_fragments
    df['# Atom in Fragment'] = new_atom_in_fragments
    df['SMILES'] = new_SMILES
    df['Molecular Formulas'] = new_molecular_formulas
    
    return df

def read_reaction_list(json_file: str) -> list:
    ''' Read JSON file containing reaction data to construct network.

    Args:
        json_file: JSON file path where the reaction list was been stored

    Returns:
        reactions_list: list of reactions
                        Format: [event #, [ts_r, ts_p], [smiles_r...], [smiles_p...]]
    '''
    with open(json_file, "r") as fp:
        reactions_list = json.load(fp)

    return reactions_list
