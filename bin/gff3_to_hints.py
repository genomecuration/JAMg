#!/usr/bin/env python
import gff_utils
import argparse, sys, os
def gff3_to_hints( ingff3, outhints ):
    gff = gff_utils.read_gff( ingff3 )
    hints  = []
    for f in gff:
        if type(f) is str:
            hints.append( f )
        elif f['type'] == 'CDS':
            f['attributes'] = 'src=M;pri=100;grp=%s' % f['atts']['Parent']
            f['type'] = 'CDSpart'
            hints.append( f )
    gff_utils.write_gff( outhints, hints )

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='convert canonical genes to hints')
    parser.add_argument( '--ingff3', type=argparse.FileType('r'), default=sys.stdin)
    parser.add_argument( '--outhints', type=argparse.FileType('w'), default=sys.stdout)
    args = parser.parse_args()
    gff3_to_hints( args.ingff3, args.outhints ) 