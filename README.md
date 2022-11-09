Processing of Computational Nanoreactor Simulations
===================================================

This package implements various tools to automatically post-process and evaluate nanoreactor MD simulations.

## Available tools include
*	Automated molecule recognition [1] from ab initio nanoreactor simulations [2] based on Wiberg Bond Orders [3] 
	
* 	Automated network construction [1]

	* RDKit [4] based event recognition
	
	* stored list of reactions for further refinement
	
* 	Visualization tools (plots, molecular grids, PyMOL movie (available only directly from GitHub repository))

## Install:
To install nanoreactor_processing type:
```shell
$ pip install nanoreactor-processing
```

## Requirements:
* python >= 3.8
* numpy >= 1.23
* scipy >= 1.8
* pandas >= 1.4 
* networkx >= 2.8
* matplotlib >= 3.5
* seaborn >= 0.11
* pyvis >= 0.2
* rdkit >= 2021.03.4

## Basic Usage:
To use the functions implemented in nanoreactor_processing you should be able to generate files with the following format from your MD simulation:

Trajectory: .xyz file

Format: <br />
1   number of atoms <br />
2   TIME: time step <br />
3   elem x y z <br />
4   elem x y z <br />
        . <br />
        . <br />

Bond order file: only the upper triangular matrix needs to be stored

Format:

1   TIME: time step <br />
2   wbo(0,1) <br />
3   wbo(0,2) <br />
4   wbo(0,3) <br />
        . <br />
        . <br />

To start the automated post-processing for your ab initio nanoreactor simulations you have to first create a NanoSim object:
```python
import nanoreactor_processing
from nanoreactor_processing import *
from nanoreactor_processing import NanoSim
from nanoreactor_processing import NanoNetwork

nanoSim = NanoSim(path_to_traj, path_to_bo_file)

# apply functions to your newly created object to generate fragments and compute SMILES:
nanoSim.generate_frag_lists()
nanoSim.generate_mols_smiles()
...
```
If you have already generated the fragment file (mols file), then you can also include it as an argument to speed up the evaluation:
```python
...
nanoSim = NanoSim(path_to_traj, path_to_bo_file, path_to_mols_file)
...
```
To generate the reaction network use the stored data frame:
```python
...
df = read_write_utils.read_trafo_df(path_to_df)
reactions_list = nanoreactor_network.construct_reactions_list(df)
nanoNet = NanoNetwork()
nanoNet.create_network(reactions_list)
...
```
If you encounter any problems with the Draw module in RDKit try to add the following line to your script:
```python
from rdkit.Chem.Draw import IPythonConsole
```
The scripts for generating a PyMOL movie from your trajectory and data frame are available only at GitHub as they require 
PyMOL as an interpreter. A free version of PyMOL can be installed with:
```shell
$ apt-get install pymol
```

## Documentation:
Code documentation can be created with pdoc3:
```shell
$ pip install pdoc3
$ pdoc3 --html nanoreactor_processing -o doc/
```

## References:
This work:
1.  A. Stan et al., J. Chem. Theory Comput. (2022); https://doi.org/10.1021/acs.jctc.2c00754

Other references:

2. L.-P. Wang et al., Nat. Chem. (2014); https://doi.org/10.1038/nchem.2099
3. K. Wiberg, Tetrahedron (1968); https://doi.org/10.1016/0040-4020(68)88057-3
4. G. Landrum, "RDKit: Open-source cheminformatics. https://www.rdkit.org"; https://doi.org/10.5281/zenodo.5085999

