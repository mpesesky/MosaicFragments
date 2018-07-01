import pandas as pd
import argparse

parser = argparse.ArgumentParser(description="Reduce intersection file to links where gene is a parametric outlier")

parser.add_argument("IntersectionFile", help="File with table of mosaic and parametric genes")
parser.add_argument("Outfile", help="Desired name of output table")

args = parser.parse_args()

df = pd.read_table(args.IntersectionFile, sep="\t", index_col=0)
df = df[(df['locus'] == df['qlocus']) | (df['locus'] == df['slocus'])]
df.reset_index(inplace=True)

df.to_csv(args.Outfile, sep="\t")
