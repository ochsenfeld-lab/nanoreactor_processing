from my_helpers import *
from pymol.cgo import *
import colorsys,sys,re
from pymol import cmd
import pandas as pd
import ast
import colorsys
import seaborn as sns
import os
import cv2
from PIL import Image
from utils import *
from read_write_utils import *

# Pymol launching
pymol.finish_launching()

# Set ray tracing on
cmd.set("ray_trace_frames", 1)

# Read in data frame and manipulate data
df = read_trafo_df("dataframe_nanosim.csv")

# Build dictionary of SMILES according to absolute occurence

smiles_prob_dict = {}
for ts in range(len(df['Time step [fs]'])):
    for struc in df['SMILES'][ts]:
        for elem in struc:
            if elem in smiles_prob_dict:
                smiles_prob_dict[elem] += 1
            else:
                smiles_prob_dict[elem] = 1

# Construct color palette and assign colors to fragments

RGB_tuples = sns.color_palette("Spectral", len(smiles_prob_dict))
np.random.shuffle(RGB_tuples)

index = 0
smiles_color_dict = {}
for struc in smiles_prob_dict:
    if smiles_prob_dict[struc] >= 1 and struc not in smiles_color_dict:
        cmd.set_color("color" + str(index), [RGB_tuples[index][0],RGB_tuples[index][1],RGB_tuples[index][2]])
        smiles_color_dict[struc] = "color" + str(index)
        index += 1

# Color fragments

for ts in range(0,len(df['Time step [fs]']),4):
    cmd.load("traj_files/small_traj_" + str(ts) + ".xyz", "obj")
    
    obj = "obj"
    cmd.hide("lines")
    cmd.hide("spheres") 
    cmd.hide("sticks")
    
    cmd.show("lines", obj)
    cmd.set("line_width", 5)
    cmd.util.cbaw("all")
    cmd.set("line_color", "grey50", "elem C")
    
    cmd.show("spheres")
    cmd.set("sphere_scale",0.5)
    cmd.set("sphere_mode",5)
    cmd.set("sphere_transparency", 0.5)
    
    for fragment in range(len(df['# Fragment'][ts])):
        for elem in df['SMILES'][ts][fragment]:
            color_for_stick = smiles_color_dict[elem]
        
        selection = str(df['# Atom in Fragment'][ts][fragment])
        selection = selection.replace(" ","")
        selection = selection.replace(",","+")
        selection = selection.replace("[","")
        selection = selection.replace("]","")
        
        cmd.set("sphere_color", color_for_stick, "index " + selection)

    cmd.bg_color("white")

    cmd.center(obj)
    cmd.zoom("center",25)
    
    path = ""
    dpi_value = 100
    if ts >= 10000:
        path = "movie/state_" + str(ts) + ".png"
    elif ts >= 1000:
        path = "movie/state_0" + str(ts) + ".png"
    elif ts >= 100:
        path = "movie/state_00" + str(ts) + ".png"
    elif ts >= 10:
        path = "movie/state_000" + str(ts) + ".png"
    else:
        path = "movie/state_0000" + str(ts) + ".png"
    
    cmd.png(path,"32cm","25cm", dpi = dpi_value)
    
    cmd.delete(obj)
    
    convertImage_transparent("grids/grid_ts_"+str(ts)+".png")
    background = cv2.imread(path)
    overlay = cv2.imread("grids/grid_ts_"+str(ts)+".png", cv2.IMREAD_UNCHANGED)  # IMREAD_UNCHANGED => open image with the alpha channel
    
    white_bg = np.zeros((background.shape[0], background.shape[1], 3), np.uint8)
    white_bg[:] = (255, 255, 255)

    x_offset = 0
    y_offset = 0
    
    add_transparent_image(white_bg, overlay, x_offset, y_offset)
    
    convertImage_transparent(path)
    mol_img = cv2.imread(path,cv2.IMREAD_UNCHANGED)

    add_transparent_image(white_bg, mol_img, x_offset, y_offset)

    text = f"Time Step: {(ts*0.5*50/1000):>6} ps"
    coordinates = (5,950)
    font = cv2.FONT_HERSHEY_SIMPLEX
    fontScale = 0.8
    color = (0,0,0)
    thickness = 1
    white_bg = cv2.putText(white_bg, text, coordinates, font, fontScale, color, thickness, cv2.LINE_AA)
   
    cv2.imwrite(path,white_bg)

    cv2.destroyAllWindows()

cmd.quit()

