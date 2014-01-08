#!/usr/bin/python

from sets import Set		# we're not running version 2.4 yet
import sys, os, re, getopt
sys.path.append(os.path.expanduser("~/lib/Python"))
from RangeFinder import RangeFinder


usage = sys.argv[0]+""" <gtf format inputfile>

Program to merge all overlapping transcripts in a gtf file to genes with multiple trancripts.
The program will automatically delete identical transcripts, unless the -d flag is set
Genes with just CDS features that are identical to the CDS of a UTR-containing genes, are also removed.
Options:         -d  Do not delete anything, just merge
                 -c  Consider genes which match only CDS identical

"""



def splitGtfId(gtfIds):
    """ extracts transcript and gene name from last column in gtf"""
    splitids = gtfIds.split('"')
    return splitids[1], splitids[3]

class FeatureObject(object):
    """ parses gtf line"""
    def __init__(self, gtf_line):
        feat = gtf_line.split('\t')
        self.chr = feat[0]
        self.defn =  feat[1]
        self.descriptor = feat [2]
        self.start= int(feat[3])
        self.stop= int(feat[4])
        self.score= feat[5]
        self.strand= feat[6]
        self.frame = feat[7]
        self.gId, self.tId = splitGtfId(feat[8])
        self.inline = gtf_line

class Gene(object):
    def __init__(self, FeatureObject):
        """ initializes chromosome, strand, and Ids from FeatureObject so we don't have to call these every time"""
        self.features=[]
        self.exons = RangeFinder()
        self.chr = FeatureObject.chr
        self.strand = FeatureObject.strand
        self.tId = FeatureObject.tId
        self.gId = FeatureObject.gId

        self.features.append(FeatureObject)
    def add(self, FeatureObject):
        """adds FeatureObject to Gene object"""
        self.features.append(FeatureObject)
    def makeExons(self):
        """ take features list and create a list of exonstarts and stops. This means merging utr and exon features"""
        self.exonStarts = []
        self.exonEnds= []
        keepStart =''
        keepStop = 0
        for f in self.features:
            if (f.start - 1 > keepStop):         # if it is a new exon
                if(keepStop > 0):
                    self.exonStarts.append(keepStart)
                    self.exonEnds.append(keepStop)
                keepStart = f.start
                    
            keepStop = f.stop
        # now add the last
        self.exonStarts.append(keepStart)
        self.exonEnds.append(keepStop)
    def addExons(self):           # addExons can not be part of add, since the exons must first be created in makeExons under the Gene class
        """ add exons to Gene object """
        for i in xrange(len(self.exonStarts)):
            self.exons.add(self.chr, self.exonStarts[i], self.exonEnds[i], self, self.strand)
    def isIdentical(self, otherGene):
        """ this should only be used after getTranscripts, because it assumes chromosome and strand are identical """
        if not (len(self.features) == len (otherGene.features)):
            return False
        for f in xrange(len(self.features)):
            if not (self.features[f].descriptor == otherGene.features[f].descriptor and self.features[f].start == otherGene.features[f].start and self.features[f].stop == otherGene.features[f].stop):
                return False
        return True    # if the number of exons is the same and the features are all the same
    def justCDS(self):
        """ tests if the gene has UTR features. Returns False if it does"""
        for f in self.features:
            if (f.descriptor == "5UTR" or f.descriptor == "3UTR"):
                return False
        return True
    def cdsIdentical(self, otherGene):
        selfFeatures = []
        otherFeatures = []
        for f in self.features:
            if f.descriptor == "CDS":
                selfFeatures.append(f)
        for f in otherGene.features:
            if f.descriptor == "CDS":
                otherFeatures.append(f)
        # now compare
        if not (len(selfFeatures) == len (otherFeatures)):
            return False
        for f in xrange(len(selfFeatures)):
            if not (selfFeatures[f].descriptor == otherFeatures[f].descriptor and selfFeatures[f].start == otherFeatures[f].start and selfFeatures[f].stop == otherFeatures[f].stop):
                return False
        return True    # if the number of exons is the same and the features are all the same
                


class GeneTable(object):
    """ holds all gene objects """
    def __init__(self):
        self.genes = []
    def add (self, gene):
        """ adds gene object to geneTable"""
        self.genes.append(gene)
    def addFeatToGene(self, FeatureObject):
        """ tests if gene is already in GeneTable object, if not, it is created"""
        for g in self.genes:
            if(g.tId == FeatureObject.tId):     # by definition, all features in a gene have the same tId
                #print FeatureObject.tId, "adding to gene"
                g.add(FeatureObject)
                return True
    def getTranscripts(self, otherGene):
        """ get all overlapping transcripts for inputgene"""
        overlappingTranscripts = Set()             # this is a set because transcripts can be added multiple times
        for i in xrange(len(otherGene.exonStarts)):
            for k in self.genes:
                for g in k.exons.overlapping(otherGene.chr, otherGene.exonStarts[i], otherGene.exonEnds[i], otherGene.strand):
                    overlappingTranscripts.add(g)
        return overlappingTranscripts

def noUTR(overlapTx):
    """ Tests if all genes in the input table have UTRs. If they do not, and their CDS is identical to that of another gene, they are removed"""
    cdsOnly = []
    containsUTR = []
    for gene in overlapTx:
        if gene.justCDS():
            cdsOnly.append(gene)
        else:
            containsUTR.append(gene)

    # now see if a gene from containsUTR has the same CDS as a gene in cdsOnly
    removeCdsOnly = []
    for gene in cdsOnly:
        for otherGene in containsUTR:
            if gene.cdsIdentical(otherGene):
                print >> sys.stderr, gene.tId, "has no UTR and identical CDS to UTR-containing", otherGene.tId
                removeCdsOnly.append(gene)
    return removeCdsOnly


def findDuplicates(overlapTx, cdsOnly):
    """ finds identical transcripts in a list and returns the ones to be removed."""
    duplicates = Set(noUTR(overlapTx))         # populate the set with genes containing only CDS 
    for i in xrange(len(overlapTx)):
        for k in xrange(i+1, len(overlapTx)):
            if overlapTx[k] in duplicates:               # if an earlier test found that this gene is a duplicate of another, do not compare again
                continue
            if cdsOnly:
                if overlapTx[i].cdsIdentical(overlapTx[k]):
                    duplicates.add(overlapTx[k])
            else:
                if overlapTx[i].isIdentical(overlapTx[k]):
                    duplicates.add(overlapTx[k])
    return list(duplicates)

def sortByGeneLength(geneList):
    for x in geneList:
        nlist = [(len(x.features),x) for x in geneList]
        nlist.sort()
        nlist.reverse()
        for key, val in nlist:
            return [val for (key, val) in nlist]
                

def makeClusters(gtfGenes, doNotDelete, cdsOnly):
    """Finds overlapping transcripts and creates clusters with identical gene_name. If a transcript is duplicated, it is removed. Output is gtf format. """
    cleanClusters = []
    skipThese = Set()               # this set keeps track of the genes we have done, so we do not find similar overlaps multiple times
    for gene in gtfGenes.genes:
        if gene.tId in skipThese:
            continue
        overlapTx = gtfGenes.getTranscripts(gene)

        # only worry about removing duplicates and merging genes if there actually are multiple transcripts
        if len(overlapTx)>1:
            if not doNotDelete:
                duplicates = findDuplicates(list(overlapTx), cdsOnly)
                for p in duplicates:
                    print >>sys.stderr, "removing duplicate", p.tId
                    skipThese.add(p.tId)
                    overlapTx.remove(p)
            # gene name should be that of longest gene in cluster
            overlapTx = sortByGeneLength(list(overlapTx))
            replaceGid = overlapTx[0].gId
            for i in overlapTx:
                skipThese.add(i.tId)
                for k in i.features:
                    out = k.inline.replace('gene_id \"'+k.gId,'gene_id \"'+replaceGid)        # rename the geneId for the subsequent transcripts
                    cleanClusters.append(out.strip())
                cleanClusters.append("")
        else:
            for i in list(overlapTx):
                skipThese.add(i.tId)
                for k in i.features:
                    cleanClusters.append(k.inline.strip())
                cleanClusters.append("")
    del cleanClusters[-1]                # remove last newline
    return cleanClusters



# Main

# read in command line and options
try:
    opts, args = getopt.getopt(sys.argv[1:], "dc")
except getopt.GetoptError:
    
        # print help information and exit:
    print usage
    print "ERROR did not recognize input\n"
    sys.exit(2)

doNotDelete = False
cdsOnly = False

for o, a  in opts:
    if o == "-d":
        doNotDelete = True
    if o == "-c":
        cdsOnly = True
    if o == "-h":
        print usage
        sys.exit()

# Read in gtf line, create FeatureObject and if neccessary, a gene object, then add feature to gene

gtfGenes = GeneTable()


if len(args) != 1:
    sys.exit(usage)
    
print >> sys.stderr, "read file"

f = open(args[0],'r')
for gtf_line in f:
    gtf_line = gtf_line.strip()
    if not gtf_line or gtf_line.startswith("#"):            # skip empty lines
        continue

    feats = FeatureObject(gtf_line)
    if not gtfGenes.addFeatToGene(feats):
        # gene does not yet exist in table: create and add first feature
        newGene = Gene(feats)
        gtfGenes.add(newGene)
f.close()


# now that all is read in, create exons and add them

for gene in gtfGenes.genes:
    gene.makeExons()
    gene.addExons()
    
print >> sys.stderr, "making clusters"

outputGtf = makeClusters(gtfGenes, doNotDelete, cdsOnly)


for i in outputGtf:
    print i
