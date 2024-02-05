#!/bin/bash

traj=${1}
df=${2}

bash main_gen_grids.sh $traj $df

# Color fragments and generate PNGs for movie
mkdir movie
pymol -cq pymol_frame_script.py $df

echo "Frames generation completed."

# Make movie and destroy folders
python3 make_movie.py

mv movie/pymol_movie.mp4 pymol_movie.mp4
rm -r traj_files/ movie/ grids/ 

if [ -f "pymol_movie.mp4" ]
then 
	printf "Rendering succesful.\nEnjoy the movie!\n"
else
	printf "Errors encountered. Movie could not be generated.\n"
fi
