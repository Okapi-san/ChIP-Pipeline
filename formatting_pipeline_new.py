import pandas as pd
import numpy as np
import regex as re
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', required=True, dest="i", help="Specify name of input file.")
args = parser.parse_args()

ID_list=[]

def ID_Check(ID):
    if ID in ID_list:
        return False
    else:
        return True
    
with open(args.i, "r") as file:
    entries = file.readlines()
    for lines in entries:
        gene_ID = "".join(re.findall('ID=(.*?);', lines))
        if ID_Check(gene_ID) == True:
            ID_list.append(gene_ID)            
        else:
            continue
for element in ID_list:
    print(element)
