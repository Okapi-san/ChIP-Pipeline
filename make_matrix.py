import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', required=True, dest="i", help="Specify name of input file.")
args = parser.parse_args()

file = pd.read_csv("matrix_raw", sep = "\t", header = None)
matrix = file.set_index([0, 1])[2].sort_index().unstack()
matrix.to_csv("matrix_final.tsv", index = False, header = False, sep = "\t")
