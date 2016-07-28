import argparse, re, sys, os
from gff_utils import read_gff, by_key, add_ID
import urllib

def read_interval( infasta, lineRE=re.compile(r'>([^:]+):([0-9]+)-([0-9]+)\(([+-])\)') ):
    intervals = {}
    for line in infasta:
        m  = lineRE.match( line )
        if m:
            seqid = m.group(1)
            start = int(m.group(2))
            end = int(m.group(3))
            strand = m.group(4)
            intervals[(seqid, str(start + 1), str(end))] = ''
        else:
            intervals[(seqid, str(start+1), str(end))] += line
    return intervals

def get_name_from_intervals( infasta, ingff, outfasta, field, filt ):
    intervals = read_interval( infasta )
    fields = ['seqid','source','type','start','end','score','strand','phase','attributes']
    template = '\t'.join(['%%(%s)s' % f for f in fields])  + '\n'
    with outfasta as out:
        for feature in read_gff( ingff, fields ):
            if type(feature) is dict  and (feature['seqid'], feature['start'], feature['end']) in intervals:
                    if field in feature['atts']:
                        out.write( '>%s %s:%s(%s)\n%s' % (feature['atts'][field], feature['start'], feature['end'], feature['strand'], intervals[(feature['seqid'], feature['start'], feature['end'])]))
                    elif field in feature:
                        out.write( '>%s\n%s' % (feature[field], intervals[(feature['seqid'], feature['start'], feature['end'])]))


def combine_field_from_intervals( infasta, ingff, outfasta, field, filt ):
    intervals = read_interval( infasta )
    fields = ['seqid','source','type','start','end','score','strand','phase','attributes']
    template = '\t'.join(['%%(%s)s' % f for f in fields])  + '\n'
    features = read_gff( ingff, fields )
    #features = add_ID( features )
    parents = by_key( features, 'ID' )
    fff = by_key( features, field )
    with outfasta as out:
        for ff in sorted(fff.keys()):
            seqs = ''
            note = ''
            name = ''
            for f in parents[ff]:
                if f['type'] in ['mRNA', 'transcript', 'gene']:
                    if 'Note' in f['atts']:
                        note = urllib.unquote(f['atts']['Note'])
                    elif 'note' in f:
                        note = urllib.unquote(f['atts']['note'])
                    if 'Name' in f['atts'] and f['atts']['Name'] != ff:
                        name = f['atts']['Name']
                    elif 'name' in f['atts'] and f['atts']['name'] != ff:
                        name = f['atts']['name']
            for f in sorted(fff[ff], lambda x,y: cmp(int(x['start']), int(y['start']))):
                if f['type'] in [filt] and (f['seqid'], f['start'], f['end']) in intervals:
                    if 'name' in f['atts'] and f['atts']['name'] != ff:
                        name = f['atts']['name']
                    if f['strand'] == '-':
                        seqs = intervals[(f['seqid'], f['start'], f['end'])] + seqs
                    else:
                        seqs += intervals[(f['seqid'], f['start'], f['end'])]
            if seqs != '':
                out.write( '>%s %s|%s %s(%s)\n%s' % (name, ff, note,  ','.join(['%s:%s' %(f['start'], f['end']) for f in fff[ff]]), f['strand'], seqs ) )

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='convert GFF interval to ID')
    parser.add_argument('infasta', nargs='?',type=argparse.FileType('r'), default=sys.stdin)
    parser.add_argument('ingff', nargs='?',type=argparse.FileType('r'), default=sys.stdin)
    parser.add_argument('outfasta',  nargs='?',type=argparse.FileType('w'), default=sys.stdout)
    parser.add_argument('--field', dest='field', default='ID')
    parser.add_argument('--filter', dest='filt', default='CDS')
    parser.add_argument('--combine', dest='get_or_combine', action='store_const', const=combine_field_from_intervals, default=get_name_from_intervals)
    args = parser.parse_args()
    args.get_or_combine(args.infasta, args.ingff, args.outfasta, args.field, args.filt)

