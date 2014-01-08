"""Container index by sequence id, range, and optionally strand that efficiently
searches for overlapping entries."""

from sets import Set

# There's a bin for each 128k (1<<17) segment. The next coarsest is 8x as big
# (1<<13).  That is for each 1M segment, for each 8M segment, for each 64M
# segment, for each 512M segment, and one top level bin for 4Gb.  Note, since
# start and end are int's, the practical limit is up to 2Gb-1, and thus, only
# four result bins on the second level. A range goes into the smallest bin it
# will fit in.

__all__ = ("RangeFinder")

class Entry(object):
    "entry associating a range with a value"
    __slots__ = ("start", "end", "value")

    def __init__(self, start, end, value):
        self.start = start
        self.end = end
        self.value = value

    def overlaps(self, start, end):
        "test if the range is overlapped by the entry"
        return (end > self.start) and (start < self.end)

class RangeBins(object):
    """Range indexed container for a single sequence.  This using a binning
    scheme that implements spacial indexing. Based on UCSC hg browser binRange
    C module.  """

    binOffsets = (4096+512+64+8+1, 512+64+8+1, 64+8+1, 8+1, 1, 0)

    binFirstShift = 17 # How much to shift to get to finest bin.
    binNextShift = 3   # How much to shift to get to next larger bin.

    def __init__(self, seqId, strand):
        self.seqId = seqId
        self.strand = strand
        self.bins = {}  # indexed by bin

    def _getBin(self, start, end):
        "get the bin for a range"

        startBin = start >> self.binFirstShift
        endBin = (end-1) >> self.binFirstShift

        for binOff in self.binOffsets:
            if (startBin == endBin):
                return startBin + binOff;
            startBin >>= self.binNextShift;
            endBin >>= self.binNextShift;

        raise Exception("start %d, end %d out of range in findBin (max is 2Gb)" % (start, end))

    def add(self, start, end, value):
        bin = self._getBin(start, end)
        entries = self.bins.get(bin)
        if (entries == None):
           self.bins[bin] = entries = []
        entries.append(Entry(start, end, value))

    def overlapping(self, start, end):
        "generator over values overlapping the specified range"

        if (start >= end):
            return
        startBin = start >> self.binFirstShift
        endBin = (end-1) >> self.binFirstShift

        for binOff in self.binOffsets:
            for j in xrange(startBin+binOff, endBin+binOff+1):
                bin = self.bins.get(j)
                if (bin != None):
                    for entry in bin:
                        if entry.overlaps(start, end):
                            yield entry.value
            startBin >>= self.binNextShift
            endBin >>= self.binNextShift

    def values(self):
        "generator over all values"
        for bin in self.bins.itervalues():
            for entry in bin:
                yield entry.value

class RangeFinder(object):
    """Container index by sequence id, range, and optionally strand.
    All entries added to the object must either have strand or not
    have strand.  A query without strand will find all overlapping
    entries on either strand if strand was specified when adding entries.
    """

    validStrands = Set((None, "+", "-"))

    def __init__(self):
        self.haveStrand = None
        self.seqBins = {}

    def add(self, seqId, start, end, value, strand=None):
        "add an entry for a sequence and range, and optional strand"
        if self.haveStrand == None:
            self.haveStrand = (strand != None)
        elif self.haveStrand != (strand != None):
            raise Exception("all RangeFinder entries must either have strand or not have strand")
        if strand not in self.validStrands:
            raise Exception("invalid strand: " + str(strand))
        key = (seqId, strand)
        bins = self.seqBins.get(key)
        if bins == None:
           self.seqBins[key] = bins = RangeBins(seqId, strand)
        bins.add(start, end, value)

    def overlapping(self, seqId, start, end, strand=None):
        "generator over values overlaping the specified range on seqId, optional strand"
        if strand not in self.validStrands:
            raise Exception("invalid strand: " + str(strand))
        if self.haveStrand and (strand == None):
            # must try both strands
            bins = self.seqBins.get((seqId, "+"))
            if bins != None:
                for value in bins.overlapping(start, end):
                    yield value
            bins = self.seqBins.get((seqId, "-"))
            if bins != None:
                for value in bins.overlapping(start, end):
                    yield value
        else:
            # must only check a specifc strand, or no strand to check
            if not self.haveStrand:
                strand = None  # no strand to check
            bins = self.seqBins.get((seqId, strand))
            if bins != None:
                for value in bins.overlapping(start, end):
                    yield value
            

    def values(self):
        "generator over all values"
        for bins in self.seqBins:
            for value in bins.values():
                yield value

