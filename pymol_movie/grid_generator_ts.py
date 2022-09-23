import pandas as pd
import ast
import os
import sys
import rdkit
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import MCS
import rdkit.Chem.Draw

import read_write_utils

# Read in data frame and manipulate data
df = read_trafo_df(sys.argv[1])

for ts in range(0,len(df['Time step [fs]'])):
    
    abs_prob_ts = {}
    
    for struc in df['SMILES'][ts]:
        for elem in struc:
            if elem in abs_prob_ts:
               abs_prob_ts[elem] += 1
            else:
               abs_prob_ts[elem] = 1
    
    list_mols = []
    list_formulas = []

    for elem in abs_prob_ts:
        if elem != "Revise structure":
            mol = Chem.MolFromSmiles(elem)

            list_formulas.append(CalcMolFormula(mol) + ": " + str(abs_prob_ts[elem]))
            list_mols.append(mol)

    img = Chem.Draw.MolsToGridImage(mols=list_mols, molsPerRow=5, subImgSize = [250,250], legends = list_formulas, maxMols = 999999)
    img.save("grid/grid_ts_"+str(ts)+".png",dpi = (1500,1500))
