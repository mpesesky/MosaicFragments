from subprocess import Popen
import argparse
import os

def refinement(regionFile, seqFile, repNumber, wrkDir):
    newSeqFile = "{}mosaic_fragments_v{}.fna".format(wrkDir, repNumber)
    extractCMD = ["python", "extract_seqs.py", "-o", newSeqFile, seqFile, regionFile]
    extract = Popen(extractCMD)
    extract.wait()

    newDB = "{}mosaic_fragments_v{}".format(wrkDir, repNumber)
    blastDBCMD = ["makeblastdb", "-in", newSeqFile, "-dbtype", "nucl", "-out", newDB]
    blastDB = Popen(blastDBCMD)
    blastDB.wait()

    newBlastOut = "{}self_check_v{}.txt".format(wrkDir, repNumber)
    blastCMD = ["blastn", "-query", newSeqFile, "-db", newDB, "-outfmt", "6", "-perc_identity", "100", "-out",
                newBlastOut]
    blast = Popen(blastCMD)
    blast.wait()

    newFiltered = "{}self_check_v{}_filtered.txt".format(wrkDir, repNumber)
    filterCMD = "python blast_filtering.py -n -l 500 {} > {}".format(newBlastOut, newFiltered)
    filtering = Popen(filterCMD, shell=True)
    filtering.wait()

    blastHandle = open(newFiltered)

    if len(list(blastHandle)) <= 2:
        return None

    newRegions = "mosaic_fragments_v{}.txt".format(repNumber + 1)
    reclusterCMD = ["python", "SharedFragmentHist.py", "-l", "500", "-e", regionFile, "-k", "-r", newRegions,
                    newFiltered]
    recluster = Popen(reclusterCMD)
    recluster.wait()

    return newRegions


parser = argparse.ArgumentParser(description="Wrapper script to refine mosaic fragment calls")

parser.add_argument("RegionFile", help="Initial output of SharedFragmentHist.py to be refined")
parser.add_argument("FastaFile", help="Original plasmid FASTA file")

args = parser.parse_args()

currDir = os.getcwd()
newDir = currDir + "/FragmentVersions/"


if not os.path.exists(newDir):
    os.makedirs(newDir)

versionNum = 0
newRegions = args.RegionFile
oldRegions = None

while newRegions is not None:
    versionNum += 1
    oldRegions = newRegions

    newRegions = refinement(oldRegions, args.FastaFile, versionNum, newDir)

os.rename(oldRegions, "mosaic_fragments_final.txt")
