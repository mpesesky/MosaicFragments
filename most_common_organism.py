import pandas as pd
import argparse

def get_mode(locus, df):
    print(df[df['qseqid'] == locus]['genus'].mode().to_string())

parser = argparse.ArgumentParser(description="Identify the most common origin genus for a particular mosaic fragment")

parser.add_argument("GenusTable", help="Table linking fragments to genera")
parser.add_argument("Fragment", help="The full name of the target fragment")

args = parser.parse_args()

genusDF = pd.read_table(args.GenusTable, sep="\t", index_col=0)

get_mode(args.Fragment, genusDF)