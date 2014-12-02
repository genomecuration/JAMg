from gff_utils import order_rule, clean_type_rule, add_exons_rule, sort_gff_rule, add_uniquid, remove_bad_parents_rule, encompass_children_rule, remove_big_genes, remove_orphans_rule, clean_gff2, MAX_INTRON, sanitize_notes_rule, merge_adjacent_features_rule, add_header_rule, add_UTR_rule
import re, argparse, sys



rules = [#('add header', add_header_rule),
         ('start <= end',order_rule),
         ('clean types that GBrowse does not accept', clean_type_rule), 
         #('sanitize Notes field', sanitize_notes_rule),
         ('Add Exons for each CDS', add_exons_rule),
         ('sort gff', sort_gff_rule),
         ('add uniqueid', add_uniquid),
         ('Add UTRs for each mRNA', add_UTR_rule),
         ('remove parents whose children have a different seqid or strand', remove_bad_parents_rule),
         ('parent encompasses children', encompass_children_rule),
         ('merge adjacent features', merge_adjacent_features_rule),
         ('remove genes larger than %d bp' % MAX_INTRON, remove_big_genes),
         ('remove orphans', remove_orphans_rule)
         ]


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='clean gff')
    parser.add_argument('gffin',nargs='?', type=argparse.FileType('r'), default=sys.stdin)
    parser.add_argument('gffout',nargs='?', type=argparse.FileType('w'), default=sys.stdout)
    args = parser.parse_args()
    clean_gff2( args.gffin, args.gffout, rules=rules )

