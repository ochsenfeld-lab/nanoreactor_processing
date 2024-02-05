#!/bin/bash

# First split the trajectory file
traj_dir="traj_files"

if [ -d "$traj_dir" ]; then
  printf "Trajectory splitting completed.\nProceeding with grid generator...\n"
else
  mkdir traj_files
  traj=${1}
  python3 split_traj.py ${traj}
  printf "Trajectory splitting completed.\n"
fi

if [ -d "$grids" ]; then
  printf "Grids generation completed.\nProceeding with frame generator...\n"
else
  mkdir grids
  df=${2}
  python3 grid_generator_ts.py ${df}
  printf "Grids generation completed.\n" 
fi
