#This python script contains functions used to manage where files are stored (which includes the creation of directories)

import os
from constants import *

def join_path(your_path, new_directory):
    new_path = os.path.join(your_path, new_directory)
    if not os.path.exists(new_path):
        os.makedirs(new_path)
        print("Directory ", new_path, " created.")

    else:
        print("Directory ", new_path, " already exists.")

    return new_path
