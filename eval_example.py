import sys
from nanoreactor_processing import *
from nanoreactor_processing import NanoSim
from nanoreactor_processing import NanoNetwork

# initialize NanoSim object

# type 1 (w/o 'mols.dat' file)
nanoSim = NanoSim(sys.argv[1], sys.argv[2])

# type 2
#nanoSim = NanoSim(sys.argv[1], sys.argv[2], sys.argv[3])
        
# compute fragments
# comment next row if you use second type of initialization
nanoSim.generate_frag_lists()

# generate SMILES and fill in data frame
nanoSim.generate_mols_smiles() # radicals=False is the default for multiplicity equal 1
nanoSim.generate_df()

# plot mol grid and bar plot
plot_tools.generate_mol_grid(nanoSim.df)
plot_tools.generate_bar_plot(nanoSim.df)

# calc start_ts_index as (equil_time%period_time) / (dt * read_rate) + 475 fs (19th time step) (for the smooth-step periodical potential every 2 ps)                                                              
# here 0 + 19     

# make reactions list and network
reactions_list = nanoreactor_network.construct_reactions_list(nanoSim.df, start_ts_index=19)

if reactions_list==[]:
    print('No events found!')
else:
    nanoNetwork = NanoNetwork()
    nanoNetwork.create_network(reactions_list)

    # plot network graphs
    plot_tools.generate_network_grid(nanoNetwork)
    plot_tools.plot_static_network_kk(nanoNetwork)
    plot_tools.plot_static_network_shell(nanoNetwork)
