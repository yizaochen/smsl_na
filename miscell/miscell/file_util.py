from os import path, mkdir
from shutil import copyfile

def check_dir_exist_and_make(file_path):
    if path.exists(file_path):
        print("{0} exists".format(file_path))
    else:
        print("mkdir {0}".format(file_path))
        mkdir(file_path)

def check_file_exist(file_path):
    exists = path.isfile(file_path)
    if exists:
        print("{0} already exist".format(file_path))
        return True
    else:
        print("{0} not exist!!!!!!".format(file_path))
        return False

def copy_verbose(f_old, f_new):
    copyfile(f_old, f_new)
    print(f'cp {f_old} {f_new}')