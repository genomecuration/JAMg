#!/usr/bin/env python

import argparse, sys

def mips_to_repeatmasker_format( infa, outfa ):
    with outfa as out:
        for line in infa:
            if line[0] == '>':
            #three-letter_classification_REdat_ID | name | repeatmasker_classification | REcat_key | ncbi_tax_id | Genus
                REdat_ID, name, rmclass, REcat_key, ncbi_tax_id, genus = line[1:-1].split('|')
                out.write('>%s#%s (%s)\n' % (REdat_ID, rmclass, line[1:-1]))
            else:
                out.write(line)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Converts MIPS format ftp://ftpmips.helmholtz-muenchen.de/plants/REdat/mipsREdat_9.3p.READme.txt to RepeatMasker format https://helix.nih.gov/Applications/repeatmasker.help')
    parser.add_argument('fain',nargs='?', type=argparse.FileType('r'), default=sys.stdin)
    parser.add_argument('faout',nargs='?', type=argparse.FileType('w'), default=sys.stdout)
    args = parser.parse_args()
    mips_to_repeatmasker_format( args.fain, args.faout )
