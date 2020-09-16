from __future__ import division
import os
import shutil

def setup_folder(folder):
    "Makes simulation directory for later use. Deletes directory before making a new one if it already exists."
    if os.path.exists(folder):
        shutil.rmtree(folder)
    os.makedirs(folder)
    os.chdir(folder)