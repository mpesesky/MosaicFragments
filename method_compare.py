import pandas as pd
import argparse

parser = argparse.ArgumentParser(description="Match results from parametric and sequence comparison tests")

parser.add_argument("Para", help="Parametic outlier genes")
parser.add_argument("Mosaic", help="Seq-compare mosaic genes")
parser.add_argument("Outtable", help="Desired name of output file")

args = parser.parse_args()

para = pd.read_table(args.Para, sep="\t", index_col=0)
mosaic = pd.read_table(args.Mosaic, sep="\t", index_col=0)

combo = para.merge(mosaic, how='inner', left_on='protein_id', right_on='geneID')
combo.fillna("None")
combo.to_csv(args.Outtable, sep="\t")
