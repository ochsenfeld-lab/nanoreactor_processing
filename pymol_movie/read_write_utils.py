import sys
import numpy as np
import pandas as pd
import ast
from typing import Tuple

# read and write  all kinds of file types

def read_traj_file(xyz_file: str) -> Tuple[list, list]:
    ''' Read trajectory file.

    Args:
        xyz_file: path to trajectory file from nanoreactor simulation
                  Format: 1 time step
                          2 number of atoms
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
