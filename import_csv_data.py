#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 29 19:54:58 2020

@author: lisa-marianeedham
"""
def fn_import_csv(location):
    
    """
    Imports .txt files and creates a list. 
    Inputs: file location as path
    
    """
    import os
    files=[]
    
    for file in os.listdir(location):
        if file.endswith(".txt"):
            files.append(file)
    for i in range(len(files)):
        files[i]=os.path.join(location, files[i])

    return files


def fn_load_array(file_list):
    
    """
   Loads csv data into array
    
    """
    import numpy as np
    
    for i in range(len(file_list)):
        data = np.genfromtxt(file_list[i], delimiter='\t', dtype=np.str, skip_header=1)
        
        return data