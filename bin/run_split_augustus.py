#!/usr/bin/env python
import argparse, sys, os
import operator

from Bio import SeqIO

def make_rundir( run, verbose ):
    if os.path.exists( run ):
        if os.path.isdir( run ):
            if verbose:
                sys.stderr.write(run + ' already exists\n')
            return os.path.abspath(run)
        else:
            raise OSError( run + ' is not a directory' )
    else:
        if verbose:
            sys.stderr.write('Creating ' + run )
        os.makedirs( run )
        return os.path.abspath( run )

def reverse_dict( d ):
    r = {}
    for k in d:
        for v in d[k]:
            if v in r:
                r[v].append( k )
            else:
                r[v] = [k]
    return r

def LPT( jobs, nprocessors, verbose=False ):
    """LPT algorithm (longest processing time) sorts jobs by processing time and assigns to machine with minimum load.
    Upper bound of 4/3 -1/(3m) OPT
    """
    if type(jobs) is dict:
        jobs = sorted(jobs.items(), key=operator.itemgetter(1), reverse=True)
        if verbose:
            sys.stderr.write('Sorted %d jobs by load\n' % (len(jobs),))

    processor_loads = dict([(i, 0) for i in range(nprocessors)])
    processor_jobs = dict([(i, []) for i in range(nprocessors)])
    for job, load in jobs:
        i,min_load  = min(processor_loads.items(), key=operator.itemgetter(1))
        processor_loads[i] += load
        processor_jobs[i].append( job )
    if verbose:
        sys.stderr.write('Assigned %d jobs into %d processors\n' % (len(jobs), nprocessors))

    return processor_jobs

def split_fasta( fasta, chunks, rundir, verbose ):
    """Split fasta into chunks based on size of sequence"""
    records = SeqIO.index( fasta.name, 'fasta' )
    if verbose:
        sys.stderr.write('Indexed %s\n' % fasta.name)

    seqsizes = dict([(r, len(records[r].seq)) for r in records])
    if verbose:
        sys.stderr.write('Extracted %d sequence sizes from index\n' % len(seqsizes))

    seqs_of_chunk = LPT( seqsizes, chunks, verbose )
    out = {}
    for c in seqs_of_chunk:
        if verbose:
            sys.stderr.write('Writing %d seqs to chunk %d\n' % (len(seqs_of_chunk[c]), c) )        
        out[c] = open( chunkfile( rundir, fasta.name, c ), 'w' )
        SeqIO.write([records[r] for r in seqs_of_chunk[c]], 
                    out[c], 'fasta')
        out[c].close()
    return seqs_of_chunk, out

def split_fasta_by_hints( fasta, seqs_of_chunk, chunks_of_seq, rundir, verbose ):
    """Split fasta into chunks based on number of hints per sequence"""
    records = SeqIO.index( fasta.name, 'fasta' )
    if verbose:
        sys.stderr.write('Indexed %s\n' % fasta.name)
    out = {}
    for c in seqs_of_chunk:
        if verbose:
            sys.stderr.write('Writing %d seqs to chunk %d\n' % (len(seqs_of_chunk[c]), c) )        
        out[c] = open( chunkfile( rundir, fasta.name, c ), 'w' )
        SeqIO.write([records[r] for r in seqs_of_chunk[c] if r in records], 
                    out[c], 'fasta')
        out[c].close()
    out[c+1] = open(chunkfile( rundir, fasta.name, c+1 ), 'w')
    nohints = [records[r] for r in records if r not in chunks_of_seq]
    if verbose:
            sys.stderr.write('Writing remaining %d seqs with no hints to chunk %d\n' % (len(nohints), c+1))
    SeqIO.write(nohints, out[c+1], 'fasta')
    out[c+1].close()
    return out


def chunkfile( rundir, f, chunk ):
    return os.path.join( rundir, '%s.%02d' % (os.path.basename( f ), chunk ) )

def split_hints( hints, chunks, rundir, verbose ):
    """Split hints into chunks based on number of hints per sequence"""
    hints_of_seq = {}
    for hint in hints:
        seqid = hint.split('\t')[0].strip()
        if seqid in hints_of_seq:
            hints_of_seq[seqid] += 1
        else:
            hints_of_seq[seqid] = 1
    hints_of_seq = sorted(hints_of_seq.items(), key=operator.itemgetter(1), reverse=True)
    seqs_of_chunk = LPT( hints_of_seq, chunks, verbose )
    chunk_of_seq = reverse_dict( seqs_of_chunk )
    out = {}
    for chunk in seqs_of_chunk:
        out[chunk] = open(chunkfile( rundir, hints.name, chunk ), 'w')
    for hint in open(hints.name, 'r'):
        seqid = hint.split('\t')[0]
        if seqid in chunk_of_seq:
            for chunk in chunk_of_seq[seqid]:
                out[chunk].write( hint )
    for chunk in out:
        out[chunk].close()
    return out, seqs_of_chunk, chunk_of_seq


def prepare_augustus_commands( UTR, gff3, species, uniqueGeneId, genemodel, alternatives, extrinsicCfgFile, hint_files, seq_files, rundir ):
    cmds = []
    for chunk in seq_files:
        augustus = ['augustus']
        augustus.append('--gff3=' + gff3)
        augustus.append('--species=' + species)
        augustus.append('--uniqueGeneId=' + uniqueGeneId)
        augustus.append('--genemodel=' + genemodel)
        augustus.append('--maxtracks=10')
        augustus.append('--UTR=' + UTR)
        augustus.append('--alternatives-from-evidence=' + alternatives)
        if chunk in hint_files:
            augustus.append('--extrinsicCfgFile=' + extrinsicCfgFile )
            augustus.append('--hintsfile=' + hint_files[chunk].name)
        augustus.append(seq_files[chunk].name)
        augustus.append('2> ' + chunkfile( rundir, 'log', chunk ))
        augustus.append('> ' + chunkfile( rundir, 'result', chunk) )
        cmds.append( '%s &' % ' '.join( augustus ) )
    return cmds
        
        
def run_split_augustus( run, fasta, species, chunks, hints, extrinsicCfgFile, cfgPath, gff3, genemodel, UTR, alternatives, out, verbose, uniqueGeneId ):
    rundir = make_rundir( run, verbose )

    if hints and extrinsicCfgFile:
        hint_files, seqs_of_chunk, chunks_of_seq = split_hints( hints, chunks, rundir, verbose )
        seq_files = split_fasta_by_hints( fasta,seqs_of_chunk, chunks_of_seq, rundir, verbose )
    else:
        hint_files = []
        seqs_of_chunk, seq_files = split_fasta( fasta, chunks, rundir, verbose )
    with out:
        cmds = prepare_augustus_commands(UTR, gff3, species, uniqueGeneId, genemodel, alternatives, extrinsicCfgFile, hint_files, seq_files, rundir )
        out.write( '\n'.join( cmds ) )

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='split the genome and the hints into chunks. Run augustus on each chunk')
    parser.add_argument('run', help='name of run directory')
    parser.add_argument('-f','--fasta', help="FASTA Genome (default %(default)s)", type=argparse.FileType('r'), default=os.getenv('GENOME_PATH'))
    parser.add_argument('-s', '--species', help="Augustus species parameter (default %(default)s)", default=os.getenv('SPECIES'))
    parser.add_argument('-n', '--chunks', help='Number of chunks to split Genome into (default %(default)s)' , type=int, default=os.getenv('LOCAL_CPUS'))
    parser.add_argument('-t', '--hints', help="Extrinsic evidence to be used as hints to augustus (Optional, but requires extrinsicCfgFile to be specified)", type=argparse.FileType('r'))
    parser.add_argument('-c', '--extrinsicCfgFile', help='Configuration file to tell Augustus how to weight different lines of evidence')
    parser.add_argument('-p', '--cfgPath', help='Location of species configuration files (default %(default)s)', default=os.getenv('AUGUSTUS_CONFIG_PATH'))
    parser.add_argument('-g', '--gff3', help='GFF3 output (default on)', default='on', choices=['on', 'off'])
    parser.add_argument('-m', '--genemodel', default='complete', help='Gene model (default %(default)s)', choices=['complete', 'partial'] )
    parser.add_argument('-u', '--UTR', default='off', help='Predict UTR (default %(default)s)', choices=['on', 'off'])
    parser.add_argument('-a', '--alternatives', default='false', help='Use alternatives from evidence (default %(default)s)', choices=['true', 'false'] )
    parser.add_argument('-o', '--out', default=sys.stdout, type=argparse.FileType('w') , help='command file (default STDOUT)')
    parser.add_argument('-q', '--uniqueGeneId', default='true', choices=['true', 'false'])
    parser.add_argument('-v', '--verbose', action='store_true', default=False)
    args = parser.parse_args()
    if args.hints and not args.extrinsicCfgFile:
        sys.stderr.write('Please specify an extrinsicCfgFile in order to use Hints\n')
        parser.print_help()
    else:
        run_split_augustus( **vars(args) )
