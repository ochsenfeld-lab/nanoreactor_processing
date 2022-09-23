from pymol import cmd
import pandas as pd
import ast
import colorsys
import seaborn as sns
import os
import cv2
from PIL import Image
import rdkit
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
import rdkit.Chem.Draw
from rdkit.Chem.Draw import IPythonConsole
import cairosvg

from read_write_utils import *

# Split trajectory files in N (# time steps) trajectory files
def split_traj(traj_file: filepath):
    ''' Split trajectory file into N (number of simulation steps) small trajectory files.

    Args:
        traj_file: filepath to trajectory file
    Returns:
        -
    '''

    natoms = read_traj_file(traj_file)[1]
    lines_per_file = natoms + 2
    smallfile = None
    with open(traj_file) as bigfile:
        for lineno, line in enumerate(bigfile):
            if lineno % lines_per_file == 0:
                if smallfile:
                    smallfile.close()
                small_filename = 'traj_files/small_traj_{}.xyz'.format(int(lineno/lines_per_file))
                smallfile = open(small_filename, "w")
            smallfile.write(line)
        if smallfile:
            smallfile.close()
    return

def add_transparent_image(background: filepath, foreground: filepath, x_offset=None: float, y_offset=None: float):
    ''' Overlap two images. Foreground image should be transparent.
        Warning: this function modifies the background image!

    Args:
        background: path to image to be used as background image
        foreground: path to image to be used as foreground
        x_offset: pixel on the x axis where the upper corner of the foreground image should be placed
        y_offset: pixel on the y axis where the upper corner of the foreground image should pe placed
    Returns:
        -
    '''

    bg_h, bg_w, bg_channels = background.shape
    fg_h, fg_w, fg_channels = foreground.shape

    assert bg_channels == 3, f'background image should have exactly 3 channels (RGB). found:{bg_channels}'
    assert fg_channels == 4, f'foreground image should have exactly 4 channels (RGBA). found:{fg_channels}'

    # center by default
    if x_offset is None: x_offset = (bg_w - fg_w) // 2
    if y_offset is None: y_offset = (bg_h - fg_h) // 2

    w = min(fg_w, bg_w, fg_w + x_offset, bg_w - x_offset)
    h = min(fg_h, bg_h, fg_h + y_offset, bg_h - y_offset)

    if w < 1 or h < 1: return

    # clip foreground and background images to the overlapping regions
    bg_x = max(0, x_offset)
    bg_y = max(0, y_offset)
    fg_x = max(0, x_offset * -1)
    fg_y = max(0, y_offset * -1)
    foreground = foreground[fg_y:fg_y + h, fg_x:fg_x + w]
    background_subsection = background[bg_y:bg_y + h, bg_x:bg_x + w]

    # separate alpha and color channels from the foreground image
    foreground_colors = foreground[:, :, :3]
    alpha_channel = foreground[:, :, 3] / 255  # 0-255 => 0.0-1.0

    # construct an alpha_mask that matches the image shape
    alpha_mask = np.dstack((alpha_channel, alpha_channel, alpha_channel))

    # combine the background with the overlay image weighted by alpha
    composite = background_subsection * (1-alpha_mask) + foreground_colors * alpha_mask

    # overwrite the section of the background image that has been updated
    background[bg_y:bg_y + h, bg_x:bg_x + w] = composite

def convertImage_transparent(path: filepath):
    ''' Make white background transparent.

    Args:
        path: path to the image to add alpha channel
    Returns:
        -
    '''

    img = Image.open(path)
    img = img.convert("RGBA")

    datas = img.getdata()

    newData = []

    # Make white background transparent
    for item in datas:
        if item[0] == 255 and item[1] == 255 and item[2] == 255:
            newData.append((255,255,255,0))
        else:
            newData.append(item)
    img.putdata(newData)
    img.save(path, dpi = (300,300))

    return

# Video Generating function
def generate_video():
    ''' Function to render the PyMol video.

    Args:
        -
    Returns:
        -
    '''        
    os.chdir(os.getcwd() + "/movie/")
    path = os.getcwd() + "/movie/"
        
    image_folder = '.' # make sure to use your folder

    dirFiles = os.listdir(image_folder) #list of directory files
    dirFiles.sort() #good initial sort but doesnt sort numerically very well
    sorted(dirFiles)

    video_name = 'pymol_movie.mp4'

    images = [img for img in dirFiles
              if img.endswith(".jpg") or
                 img.endswith(".jpeg") or
                 img.endswith("png")]

    # Array images should only consider
    # the image files ignoring others if any
    print(images)

    frame = cv2.imread(os.path.join(image_folder, images[0]),cv2.IMREAD_UNCHANGED)
    # setting the frame width, height width
    # the width, height of first image
    height, width, layers = frame.shape
    fourcc = cv2.VideoWriter_fourcc('m', 'p', '4', 'v')
    video = cv2.VideoWriter(video_name, fourcc, 10, (width, height))

    # Appending the images to the video one by one
    for image in images:
        img = cv2.imread(os.path.join(image_folder, image),cv2.IMREAD_UNCHANGED)
        video.write(img)

    # Deallocating memories taken for window creation
    cv2.destroyAllWindows()
    video.release() # releasing the video generated

    return
