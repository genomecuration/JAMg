#!/usr/bin/env python3

import argparse
import sys
import os
from Bio import SeqIO
import operator


def make_rundir(run, verbose):
    """Create or validate a run directory."""
    if os.path.exists(run):
        if os.path.isdir(run):
            if verbose:
                sys.stderr.write(run + ' already exists\n')
            return os.path.abspath(run)
        else:
            raise OSError(run + ' is not a directory')
    else:
        if verbose:
            sys.stderr.write('Creating ' + run + '\n')
        os.makedirs(run)
        return os.path.abspath(run)


def reverse_dict(d):
    """Reverse a dictionary."""
    r = {}
    for k in d:
        for v in d[k]:
            if v in r:
                r[v].append(k)
            else:
                r[v] = [k]
    return r


def LPT(jobs, nprocessors, verbose=False):
    """Longest Processing Time (LPT) algorithm for job scheduling."""
    if type(jobs) is dict:
        jobs = sorted(jobs.items(), key=operator.itemgetter(1), reverse=True)
    if verbose:
        sys.stderr.write('Sorted %d jobs by load\n' % len(jobs))
    processor_loads = dict([(i, 0) for i in range(nprocessors)])
    processor_jobs = dict([(i, []) for i in range(nprocessors)])
    for job, load in jobs:
        i, min_load = min(processor_loads.items(), key=operator.itemgetter(1))
        processor_loads[i] += load
        processor_jobs[i].append(job)
    if verbose:
        sys.stderr.write('Assigned %d jobs into %d processors\n' % (len(jobs), nprocessors))
    return processor_jobs


def validate_file(file_path, description):
    """Validate that a file exists and is readable."""
    if not file_path:
        raise ValueError(f"The {description} is required but was not provided.")
    if not os.path.isfile(file_path):
        raise FileNotFoundError(f"The {description} '{file_path}' does not exist.")
    if not os.access(file_path, os.R_OK):
        raise PermissionError(f"The {description} '{file_path}' is not readable.")


def validate_output_files(file_paths, verbose=False):
    """Validate that all output files exist and are non-empty."""
    for file_path in file_paths:
        if not os.path.isfile(file_path):
            raise ValueError(f"Output file '{file_path}' does not exist.")
        if os.path.getsize(file_path) == 0:
            raise ValueError(f"Output file '{file_path}' is empty.")
        if verbose:
            print(f"Validated: {file_path} (size: {os.path.getsize(file_path)} bytes)")


def chunkfile(rundir, f, chunk):
    """Generate a filename for a specific chunk."""
    return os.path.join(rundir, '%s.%02d' % (os.path.basename(f), chunk))

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

def split_fasta(fasta, chunks, rundir, verbose):
    """Split fasta into chunks based on size of sequence."""
    records = SeqIO.index(fasta.name, 'fasta')
    if verbose:
        sys.stderr.write(f'Indexed {fasta.name}\n')
    seqsizes = dict([(r, len(records[r].seq)) for r in records])
    seqs_of_chunk = LPT(seqsizes, chunks, verbose)
    out = {}
    output_files = []
    for c in seqs_of_chunk:
        if verbose:
            sys.stderr.write(f'Writing {len(seqs_of_chunk[c])} seqs to chunk {c}\n')
        chunk_file = chunkfile(rundir, fasta.name, c)
        out[c] = open(chunk_file, 'w')
        SeqIO.write([records[r] for r in seqs_of_chunk[c]], out[c], 'fasta')
        out[c].close()
        output_files.append(chunk_file)
    validate_output_files(output_files, verbose)
    return seqs_of_chunk, out


def split_fasta_by_hints(fasta, seqs_of_chunk, chunks_of_seq, rundir, verbose):
    """Split fasta into chunks based on number of hints per sequence."""
    records = SeqIO.index(fasta.name, 'fasta')
    if verbose:
        sys.stderr.write(f'Indexed {fasta.name}\n')
    out = {}
    output_files = []
    
    for c in seqs_of_chunk:
        if verbose:
            sys.stderr.write(f'Writing {len(seqs_of_chunk[c])} seqs to chunk {c}\n')
        chunk_file = chunkfile(rundir, fasta.name, c)
        out[c] = open(chunk_file, 'w')
        SeqIO.write([records[r] for r in seqs_of_chunk[c] if r in records], out[c], 'fasta')
        out[c].close()
        output_files.append(chunk_file)

    # Handle sequences with no hints
    nohints = [records[r] for r in records if r not in chunks_of_seq]
    
    if len(nohints) > 0:  # Only write this chunk if there are sequences with no hints
        chunk_file_no_hints = chunkfile(rundir, fasta.name, len(seqs_of_chunk))
        out[len(seqs_of_chunk)] = open(chunk_file_no_hints, 'w')
        
        if verbose:
            sys.stderr.write(f'Writing remaining {len(nohints)} seqs with no hints to chunk {len(seqs_of_chunk)}\n')
        
        SeqIO.write(nohints, out[len(seqs_of_chunk)], 'fasta')
        out[len(seqs_of_chunk)].close()
        output_files.append(chunk_file_no_hints)
    
    else:  # Log skipping empty "no hints" chunk
        if verbose:
            sys.stderr.write(f'Skipping creation of chunk {len(seqs_of_chunk)} as there are no sequences with no hints.\n')

    validate_output_files(output_files, verbose)
    
    return out


def prepare_augustus_commands(UTR, gff3, species, uniqueGeneId,
                              genemodel, alternatives,
                              extrinsicCfgFile,
                              hint_files,
                              seq_files,
                              rundir,
                              softmasking):
    """Prepare Augustus commands for each chunk."""
    cmds = []
    
    # Validate species argument
    if not species:
        raise ValueError("The --species argument or SPECIES environment variable must be specified.")

    for chunk in seq_files:
        augustus = ['augustus']
        
        augustus.append('--min_intron_len=30')
        augustus.append('--softmasking=' + softmasking)
        augustus.append('--gff3=' + gff3)
        augustus.append('--species=' + species)
        augustus.append('--uniqueGeneId=' + uniqueGeneId)
        augustus.append('--genemodel=' + genemodel)
        
        augustus.append('--maxtracks=10')
        
        augustus.append('--UTR=' + UTR)
        
        augustus.append('--alternatives-from-evidence=' + alternatives)

        if chunk in hint_files:
            augustus.append('--extrinsicCfgFile=' + extrinsicCfgFile)
            augustus.append('--hintsfile=' + hint_files[chunk].name)

        augustus.append(seq_files[chunk].name)
        
        augustus.append('2> ' + chunkfile(rundir, 'log', chunk))
        
        augustus.append('> ' + chunkfile(rundir, 'result', chunk))

        cmds.append(' '.join(augustus))
    
    return cmds


def run_split_augustus(run,
                       fasta,
                       species,
                       chunks,
                       hints=None,
                       extrinsicCfgFile=None,
                       cfgPath=None,
                       gff3='on',
                       genemodel='complete',
                       UTR='off',
                       alternatives='false',
                       out=sys.stdout,
                       verbose=False,
                       uniqueGeneId='true',
                       softmasking='1'):
    rundir = make_rundir(run, verbose)

    # Validate input files
    validate_file(fasta.name, "FASTA genome")
    
    if hints and extrinsicCfgFile:
        validate_file(hints.name, "Hints file")
    
    if extrinsicCfgFile:
        validate_file(extrinsicCfgFile, "Extrinsic configuration file")

    # Split data into chunks
    if hints and extrinsicCfgFile:
        hint_files, seqs_of_chunk, chunks_of_seq = split_hints(hints, chunks, rundir, verbose)
        seq_files = split_fasta_by_hints(fasta, seqs_of_chunk, chunks_of_seq, rundir, verbose)
    else:
        hint_files = []
        seqs_of_chunk, seq_files = split_fasta(fasta, chunks, rundir, verbose)

    # Prepare Augustus commands and write to output
    with out:
        cmds = prepare_augustus_commands(UTR=UTR,
                                         gff3=gff3,
                                         species=species,
                                         uniqueGeneId=uniqueGeneId,
                                         genemodel=genemodel,
                                         alternatives=alternatives,
                                         extrinsicCfgFile=extrinsicCfgFile,
                                         hint_files=hint_files,
                                         seq_files=seq_files,
                                         rundir=rundir,
                                         softmasking=softmasking)
        out.write('\n'.join(cmds) + '\n')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Split the genome and the hints into chunks. Run Augustus on each chunk.')

    parser.add_argument('run', help='Name of run directory')
    parser.add_argument('-f', '--fasta', help="FASTA Genome (default %(default)s)", type=argparse.FileType('r'), default=os.getenv('GENOME_PATH'))
    parser.add_argument('-s', '--species', help="Augustus species parameter (default %(default)s)", default=os.getenv('SPECIES'))
    parser.add_argument('-n', '--chunks', help='Number of chunks to split Genome into (default %(default)s)', type=int, default=os.getenv('LOCAL_CPUS'))
    parser.add_argument('-t', '--hints', help="Extrinsic evidence to be used as hints to Augustus (Optional, but requires extrinsicCfgFile to be specified)", type=argparse.FileType('r'))
    parser.add_argument('-c', '--extrinsicCfgFile', help='Configuration file to tell Augustus how to weight different lines of evidence')
    parser.add_argument('-p', '--cfgPath', help='Location of species configuration files (default %(default)s)', default=os.getenv('AUGUSTUS_CONFIG_PATH'))
    parser.add_argument('-g', '--gff3', help='GFF3 output (default on)', default='on', choices=['on', 'off'])
    parser.add_argument('-m', '--genemodel', default='complete', help='Gene model (default %(default)s)', choices=['complete', 'partial', 'intronless'])
    parser.add_argument('-u', '--UTR', default='off', help='Predict UTR (default %(default)s)', choices=['on', 'off'])
    parser.add_argument('-a', '--alternatives', default='false', help='Use alternatives from evidence (default %(default)s)', choices=['true', 'false'])
    parser.add_argument('-o', '--out', default=sys.stdout, type=argparse.FileType('w'), help='Command file (default STDOUT)')
    parser.add_argument('-q', '--uniqueGeneId', default='true', choices=['true', 'false'])
    parser.add_argument('-x', '--softmasking', default='1', help='Use softmasked regions as repeats', choices=['1', '0'])
    parser.add_argument('-v', '--verbose', action='store_true', default=False)

    args = parser.parse_args()

    # Check for missing required arguments
    if args.hints and not args.extrinsicCfgFile:
        sys.stderr.write('Please specify an extrinsicCfgFile in order to use Hints\n')
        parser.print_help()
        sys.exit(1)

    # Run the main function
    run_split_augustus(**vars(args))

