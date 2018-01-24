import pandas as pd
import argparse
import matplotlib.pyplot as plt
import math

parser = argparse.ArgumentParser(description="Get fraction of IS26-containing plasmids by species")

parser.add_argument("TableFile", help="Data table file")
parser.add_argument("-c", "--cutoff", default=10, type=int, help="Minimum number of plasmids per species in results")
parser.add_argument("-o", "--output", default="IS26Taxonomy.tsv", help="Name of output tsv file")
parser.add_argument("-p", "--outplt", default="IS26Taxonomy.png", help="Name of output png file")

args = parser.parse_args()

df = pd.read_csv(args.TableFile, sep="\t")

taxonDF = pd.DataFrame({'count': df.groupby(['Genus', 'IS26']).size()}).reset_index()

taxonPivot = taxonDF.pivot(index='Genus', columns='IS26', values='count')

taxonBaseline = taxonPivot[(taxonPivot.Present + taxonPivot.Absent) >= args.cutoff]
print(taxonBaseline)

orderedOrganisms = ['Citrobacter', 'Edwardsiella', 'Enterobacter', 'Escherichia', 'Klebsiella', 'Pantoea', 'Proteus',
                    'Providencia', 'Salmonella', 'Serratia', 'Shigella', 'Yersinia', 'Acinetobacter', 'Pseudomonas',
                    'Aeromonas', 'Photobacterium', 'Vibrio', 'Rhodobacter', 'Corynebacterium']
colors = ['m', 'm', 'm', 'm', 'm', 'm', 'm', 'm', 'm', 'm', 'm', 'c', 'c', 'c', 'c', 'c', 'c', 'c', 'c']
#mapping = {org: i for i, org in enumerate(orderedOrganisms)}
#key = taxonBaseline.index.map(mapping)
#taxonBaseline = taxonBaseline.iloc[key.argsort()]

taxonBaseline['Proportion with IS26'] = taxonBaseline['Present']/(taxonBaseline['Absent'] + taxonBaseline['Present'])
taxonBaseline['se'] = ((taxonBaseline['Proportion with IS26']*(1-taxonBaseline['Proportion with IS26']))
                       /(taxonBaseline['Present'] + taxonBaseline['Absent'])).map(math.sqrt)

sortedDF = taxonBaseline.reindex(orderedOrganisms)
sortedDF.plot.bar(y='Proportion with IS26', yerr='se', color='m')
ax = plt.gca()
ax.legend_.remove()
ax.set_ylabel('Proportion of plasmids with IS26')
taxonBaseline.to_csv(args.output, sep="\t")

plt.tight_layout()
plt.savefig(args.outplt)