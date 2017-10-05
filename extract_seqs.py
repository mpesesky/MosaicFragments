import argparse
import SharedFragmentHist as sfh
import Fasta_one_line as fol

parser = argparse.ArgumentParser(description="Generate a Fasta from a list of regions")

parser.add_argument("FASTA", help="Initial Fasta file used to generate regions")
parser.add_argument("REGIONS", help="Region file generated from SharedFragmentHist.py")
parser.add_argument("-o", "--outfile", help="Destination fasta file")
parser.add_argument("-l", "--length", default=500, type=int, help="Minimum fragment length (default = 500")

args = parser.parse_args()

regionFile = open(args.REGIONS)
regionDict = {}

for line in regionFile:
    fields = line.split("\t")
    newRegion = sfh.Region(fields[0], int(fields[1]), int(fields[2]), args.length)
    newRegion.numOccur = int(fields[4].rstrip("\n"))
    if fields[0] in regionDict.keys():
        regionDict[fields[0]].append(newRegion)
    else:
        regionDict[fields[0]] = [newRegion]
regionFile.close()

print(len(list(regionDict.keys())))

fastaDict = fol.one_line_d(args.FASTA)

for header in fastaDict.keys():
    uid = header.lstrip(">").split(" ")[0]
    if uid in regionDict:
        for reg in regionDict[uid]:
            reg.seq = fastaDict[header][(reg.start - 1):reg.stop]

outfile = open(args.outfile, 'w')

for plasmid in regionDict.keys():
    for reg in regionDict[plasmid]:
        outfile.write(">{}_{}_{}_numOccur:{}\n{}\n".format(reg.uid, reg.start, reg.stop, reg.numOccur, reg.seq))
outfile.close()
