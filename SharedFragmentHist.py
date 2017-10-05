import argparse
import pandas as pd
import blast_filtering as bf
import re

class Region:
    def __init__(self, uid, start, stop, mergeLen):
        self.uid = uid
        start = int(start)
        stop = int(stop)
        if start < stop:
            self.start = start
            self.stop = stop
        else:
            self.stop = start
            self.start = stop
        self.numOccur = 2
        self.identifier = "{}_{}_{}".format(uid,start,stop)
        self.mergeLen= mergeLen

    def length(self):
        return (self.stop - self.start)

    def compare(self, region2):
        if self.uid != region2.uid:
            return False
        if (self.start < region2.start) and (self.stop < region2.start):
            return False
        if (self.start > region2.stop) and (self.stop > region2.stop):
            return False
        if (self.start < region2.start) and ((self.stop - region2.start) < self.mergeLen):
            return False
        if (self.stop > region2.stop) and ((region2.stop - self.start) < self.mergeLen):
            return False
        else:
            return True

    def internal_match(self, region2):
        if ((self.uid == region2.uid) and (self.start == region2.start) and (self.stop == region2.stop)):
            return True
        elif (self.uid == region2.uid) and (self.start <= region2.start) and (
                    self.stop >= region2.stop) and ((self.length() - region2.length()) < self.mergeLen):
            return True
        else:
            return False

    def exact_match(self, region2):
        return (self.uid == region2.uid) and (self.start == region2.start) and (self.stop == region2.stop)

    def merge(self, region2):

        newOccur = self.numOccur + region2.numOccur - 1

        if self.start > region2.start:
            newStart = self.start
            lftstart = region2.start
            lftstop = self.start
            lftOccur = region2.numOccur
        else:
            newStart = region2.start
            lftstart = self.start
            lftstop = region2.start
            lftOccur = self.numOccur


        if self.stop < region2.stop:
            newStop = self.stop
            rgtstart = self.stop
            rgtstop = region2.stop
            rgtOccur = region2.numOccur
        else:
            newStop = region2.stop
            rgtstart = region2.stop
            rgtstop = self.stop
            rgtOccur = self.numOccur

        merged = Region(self.uid, newStart, newStop, self.mergeLen)
        merged.numOccur = newOccur

        lftRegion = Region(self.uid, lftstart, lftstop, self.mergeLen)
        if lftRegion.length() < self.mergeLen:
            lftRegion = None
        else:
            lftRegion.numOccur = lftOccur

        rgtRegion = Region(self.uid, rgtstart, rgtstop, self.mergeLen)
        if rgtRegion.length() < self.mergeLen:
            rgtRegion = None
        else:
            rgtRegion.numOccur = rgtOccur

        return (lftRegion, merged, rgtRegion)

    def exact_match_fast(self, region2):
        return (self.identifier == region2.identifier)

def collect_regions(df, minLength):
    rows = df.iterrows()
    regions = []

    for i, row in rows:
        newRegion = Region(row['sseqid'], row['sstart'], row['send'], minLength)
        merged = False
        if len(regions) == 0:
            regions.append(newRegion)
            continue
        for oldRegion in regions:
            if newRegion.internal_match(oldRegion):
                oldRegion.numOccur += 1
                merged = True
                break
            if newRegion.compare(oldRegion):
                regions.remove(oldRegion)
                (left, center, right) = newRegion.merge(oldRegion)
                if left is not None:
                    regions.append(left)
                if right is not None:
                    regions.append(right)
                regions.append(center)
                merged = True
                break
#If newRegion was not merged
        if not merged:
            regions.append(newRegion)
    return regions


def cleanup(regions):
    existingRegions = list(regions)
    print("Clean up")
    mod = True
    while mod:
        mod = False
        print(len(existingRegions))
        newList = []
        while len(existingRegions) > 1:
            reg = existingRegions.pop()
            comparisonRegions = list(existingRegions)
            added = False
            for newReg in comparisonRegions:
                if reg.internal_match(newReg):
                    existingRegions.remove(newReg)
                    newReg.numOccur += reg.numOccur
                    newList.append(newReg)
                    added = True
                    mod = True
                    break
                elif reg.compare(newReg):
                    existingRegions.remove(newReg)
                    (left, center, right) = reg.merge(newReg)
                    if left is not None:
                        newList.append(left)
                    if right is not None:
                        newList.append(right)
                    newList.append(center)
                    added = True
                    mod = True
                    break
            if not added:
                newList.append(reg)
        else:
            newList.extend(existingRegions)
        existingRegions = list(newList)
    return existingRegions


def make_df(regions, binList):
    binMin = 2
    table = {}

    for limit in binList:
        binMax = limit
        table[limit] = 0
        for reg in regions:
            if (reg.numOccur >= binMin) and (reg.numOccur < binMax):
                table[limit] += 1
 #               regions.remove(reg)
        binMin = limit
    return pd.DataFrame(table, index=['Bin', 'Count'])


def import_region_file(regionFile, minLength):
    regionDict = {}

    with open(regionFile) as rf:
        for line in rf:
            fields = line.rstrip("\n").split("\t")

            newReg = Region(fields[0], fields[1], fields[2], minLength)
            newReg.numOccur = int(fields[4])
            if fields[0] in regionDict.keys():
                regionDict[fields[0]].append(newReg)
            else:
                regionDict[fields[0]] = [newReg]
    return regionDict


def find_matching_region(query, subjectList):
    for reg in subjectList:
        if reg.exact_match_fast(query):
            return reg
    return None


def decide_regions(region1, region2):
    if region1.numOccur > region2.numOccur:
        return region1, region2
#        newOccur = region1.numOccur
    else:
        return region2, region1
#        newOccur = region2.numOccur

#    if region1.length() < region2.length():
#        region1.numOccur = int(newOccur)
#        return region1, region2
#    else:
#        region2.numOccur = int(newOccur)
#        return region2, region1


def check_if_none(someThing, label="something", errMessage="Some problem"):
    if someThing is None:
        print(label)
        print(errMessage)
        exit()
    else:
        return someThing


def evaluate_regions(regions, blastTable, minLength, keepSingles=False):
    rows = blastTable.iterrows()

    newName = re.compile('(?P<uid>\S+)_(?P<start>\d+)_(?P<stop>\d+)')
    newRegions = {}
    dumpRegions = {}

    for i, row in rows:
        query = check_if_none(newName.match(row['qseqid']), "query", row)
        qparts = query.groupdict()
        tmpRegion = Region(qparts['uid'], qparts['start'], qparts['stop'], minLength)
        qRegion = check_if_none(find_matching_region(tmpRegion, regions[qparts['uid']]), "qRegion", row)

        subject = check_if_none(newName.match(row['sseqid']), "subject", row)
        sparts = subject.groupdict()
        tmpRegion = Region(sparts['uid'], sparts['start'], sparts['stop'], minLength)
        sRegion = check_if_none(find_matching_region(tmpRegion, regions[sparts['uid']]), "subject", row)

        keeper, tosser = decide_regions(qRegion, sRegion)

        if tosser.uid in dumpRegions.keys():
            totalToss = find_matching_region(tosser, dumpRegions[tosser.uid])
            if totalToss is None:
                dumpRegions[tosser.uid].append(tosser)
        elif tosser.uid in newRegions.keys():
            retractKeep = find_matching_region(tosser, newRegions[tosser.uid])
            if retractKeep is not None:
                newRegions[tosser.uid].remove(retractKeep)
            dumpRegions[tosser.uid] = [tosser]
        else:
            dumpRegions[tosser.uid] = [tosser]



        if keeper.uid in dumpRegions.keys():
            potentialToss = find_matching_region(keeper, dumpRegions[keeper.uid])
            if potentialToss is not None:
                continue
        print(i)
        if keeper.uid in newRegions:
            oldKeeper = find_matching_region(keeper, newRegions[keeper.uid])
            if oldKeeper is not None:
                newRegions[keeper.uid].remove(oldKeeper)
                ultimate, discard = decide_regions(oldKeeper, keeper)
                newRegions[keeper.uid].append(ultimate)
            else:
                newRegions[keeper.uid].append(keeper)
        else:
            newRegions[keeper.uid] = [keeper]

    if keepSingles:
        for reg in flatten(regions):
            try:
                keepFind = find_matching_region(reg, newRegions[reg.uid])
            except KeyError:
                keepFind = None
            try:
                tossFind = find_matching_region(reg, dumpRegions[reg.uid])
            except KeyError:
                tossFind = None
            if (keepFind is None) and (tossFind is None):
                if reg.uid in newRegions.keys():
                    newRegions[reg.uid].append(reg)
                else:
                    newRegions[reg.uid] = [reg]
    return newRegions

def flatten(dictOfLists):
    outList = []
    for key in dictOfLists.keys():
        for item in dictOfLists[key]:
            outList.append(item)
    return outList


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Make a histogram of region occurances")

    parser.add_argument("BLAST", help="BLAST results input file")
    parser.add_argument("-t", "--total_nuc", action='store_true', help="Chart total nucleotides, not number of segments")
    parser.add_argument("-e", "--eval", default=None, type=str, help="Evaluate existing results")
    parser.add_argument("-r", "--region", default=None, type=str, help="Region File Name")
    parser.add_argument("-k", "--keep_aligned", action='store_true', help="With -e, keep all regions not found in BLAST file")
    parser.add_argument("-l", "--length", default=500, type=int, help="Minimum fragment length (default = 500)")

    args = parser.parse_args()


    upperLimits = [10, 20, 40, 80, 160, 320, 640, 1280, 2560]

    #colNames = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart",
    #                "send", "evalue", "bitscore"]
    df = pd.read_table(args.BLAST, sep="\t")

    if args.eval is not None:
        regions = import_region_file(args.eval, args.length)
        coreRegions = evaluate_regions(regions, df, args.length, args.keep_aligned)
        regions = flatten(coreRegions)
    else:
        regions = collect_regions(df, args.length)
        regions = cleanup(regions)

    #print(len(regions))

    if args.region is not None:
        outfile = open(args.region, 'w')
        for region in regions:
            outfile.write("{}\t{}\t{}\t{}\t{}\n".format(region.uid, region.start, region.stop,
                                                    region.length(), region.numOccur))

    countdf = make_df(regions, upperLimits)

    bf.print_full(countdf)

