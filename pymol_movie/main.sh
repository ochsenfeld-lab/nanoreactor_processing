#!/bin/bash

# First split the trajectory file
traj_dir="traj_files"

if [ -d "$traj_dir" ]; then
  echo "Proceeding with grid generator..."
else
  mkdir traj_files
  traj=${1}
  python3 split_traj.py ${traj}
fi

echo "Trajectory splitting completed."

if [ -d "$grids" ]; then
  echo "Proceeding with frame generator..."
else
  mkdir grids
  df=${2}
  python3 grid_generator_ts.py ${df}
fi

echo "Grids generation completed."

# Color fragments and generate PNGs for movie
mkdir movie
pymol -cq pymol_frame_script.py

echo "Frames generation completed."

# Make movie and destroy folders
python3 script_make_movie.py
cp movie/pymol_movie.mp4 pymol_movie.mp4
rm -r traj_files/ movie/ grids/ 

echo "Rendering succesful!"

echo "DONE"
