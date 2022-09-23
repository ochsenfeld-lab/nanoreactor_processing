# importing libraries
import os
import cv2
from PIL import Image
from utils import *

# Checking the current directory path
print(os.getcwd())

# Folder which contains all the images from which video is to be generated
os.chdir(os.getcwd() + "/movie/")
path = os.getcwd() + "/movie/"

# Calling the generate_video function
generate_video(path)

