import re, argparse, sys, os, urllib
import copy

MAX_INTRON = 2e5

def alignment2genes( ingff, outgff ):
    alignment = read_gff( ingff )
    alignIds = by_key( alignment, 'ID' )
    features = []
    for a in alignIds:
         gene = copy.deepcopy( alignIds[a][0] )
         gene['start'] = str(min([int(aa['start']) for aa in alignIds[a] if aa['strand'] == gene['strand'] and aa['seqid'] == gene['seqid']]))
         gene['end'] = str(max([int(aa['end']) for aa in alignIds[a] if aa['strand'] == gene['strand'] and aa['seqid'] == gene['seqid']]))
         gene['type'] = 'gene'
         mRNA = copy.deepcopy( gene )
         mRNA['type'] = 'mRNA'
         mRNA['atts']['Parent'] = mRNA['atts']['ID'] + '.gene'
         mRNA['attributes'] = ';'.join(['%s=%s' % (aaa, mRNA['atts'][aaa]) for aaa in mRNA['atts']])
         gene['atts']['ID'] = gene['atts']['ID'] + '.gene'
         gene['attributes'] = ';'.join(['%s=%s' % (aaa, gene['atts'][aaa]) for aaa in gene['atts']])
         features.append( gene )
         features.append( mRNA )
         count = 0
         for cds in alignIds[a]:
              cds['atts']['Parent'] = cds['atts']['ID']
              cds['atts']['ID'] = cds['atts']['ID'] + '.cds.' + str(count)                
              cds['attributes'] = ';'.join(['%s=%s' % (aaa, cds['atts'][aaa]) for aaa in cds['atts']])
              exon = copy.deepcopy( cds )
              exon['atts']['ID'] = exon['atts']['ID'].replace('cds', 'exon')
              exon['attributes'] = exon['attributes'].replace('cds','exon')
              features.append( exon )
              features.append( cds )
              count += 1

    with outgff as out:
        out.write('##gff-version 3\n')
        features.sort(cmp=cmpgff)
        for feature in features:
            out.write( template % feature )

def get_distance( r1, r2 ):
     # sort the two ranges such that the range with smaller first element
     # is assigned to x and the bigger one is assigned to y
     x, y = sorted((r1, r2))

     #now if x[1] lies between x[0] and y[0](x[1] != y[0] but can be equal to x[0])
     #then the ranges are not overlapping and return the difference of y[0] and x[1]
     #otherwise return 0 
     if x[0] <= x[1] < y[0] and all( y[0] <= y[1] for y in (r1,r2)):
        return y[0] - x[1]
     return 0

def cmpgff( x, y ):
    if type(x) is str and type(y) is str:
        return cmp(x,y)
    elif type(x) is str:
        return -1
    elif type(y) is str:
        return 1

    elif x['seqid'] != y['seqid']:
        return cmp(x['seqid'], y['seqid'])
    z = int(x['start']) - int(y['start'])
    if z != 0:
        return z
    z = int(y['end']) - int(x['end'])
    if z != 0:
        return z
    else:
        typeorder = ['transcription_end_site','transcription_start_site','gene','mRNA', 'exon', 'three_prime_utr', 'five_prime_utr', 'five_prime_cis_splice_site', 'three_prime_cis_splice_site', 'five_prime_UTR', 'three_prime_UTR', 'CDS', 'start_codon', 'stop_codon','intron', 'cDNA_match', 'EST_match']
        return typeorder.index(x['type']) - typeorder.index(y['type'])

def cmpsib( x, y ):
    if x['seqid'] != y['seqid']:
        return cmp(x['seqid'], y['seqid'])
    if x['type'] != y['type']:
        typeorder = ['transcription_end_site','transcription_start_site','gene','mRNA', 'exon', 'three_prime_utr', 'five_prime_utr', 'five_prime_cis_splice_site', 'three_prime_cis_splice_site', 'five_prime_UTR', 'three_prime_UTR', 'CDS', 'start_codon', 'stop_codon','intron', 'cDNA_match', 'EST_match']
        return typeorder.index(x['type']) - typeorder.index(y['type'])
    z = int(x['start']) - int( y['start'])
    if z != 0:
        return z
    return int(x['end']) - int(y['end'])

def merge_adjacent_features_rule( features, children, uniquid ):
    newfeatures = []
    for parent in children:
        if len(children[parent]) > 1:
            siblings = sorted(children[parent], cmp=cmpsib)
            conjoin = copy.deepcopy(siblings[0])
            for sib in siblings[1:]:
                if sib['type'] == conjoin['type'] and (int(sib['start']) - int(conjoin['end'])) == 1:
                    conjoin['end'] = sib['end']
                    del uniquid[sib['atts']['ID']]
                else:
                    uniquid[conjoin['atts']['ID']] = [ copy.deepcopy(conjoin ) ]
                    conjoin = copy.deepcopy( sib )
    for feature in features:
        if type(feature) is str:
            newfeatures.append( feature )
        elif feature['atts']['ID'] in uniquid:
            newfeatures.append( uniquid[feature['atts']['ID']][0] )
    return newfeatures, by_key( newfeatures, 'Parent' ), by_key( newfeatures, 'ID' )
def sanitize_notes_rule( features, children, uniquid ):
    urlRE = re.compile(r"[^A-Za-z0-9%._/]")
    for feature in features:
        if type(feature) is dict and 'Note' in feature['atts'] and urlRE.search( features['atts']['Note'] ):
            feature['atts']['Note'] = urllib.quote( feature['atts']['Note'] )
            feature['attributes'] = ';'.join(['%s=%s' % (a, feature['atts'][a]) for a in feature['atts']])
    return features, children, uniquid


def order_rule( features, children, uniquid ):
    for feature in features:
        if type(feature) is dict:
            if int(feature['start']) > int(feature['end']):
                feature['start'], feature['end'] = feature['end'], feature['start']
    return features,children, uniquid

def sort_gff_rule( features, children, uniquid ):
    features.sort(cmp=cmpgff)
    return features, children, uniquid

def encompass_children_rule( features, children, uniquid ):
    for feature in features:
        if type(feature) is dict:
            if feature['type'] in ['gene', 'mRNA', 'transcript'] and feature['atts']['ID'] in children:
                feature['start'] = str(min([int(c['start']) for c in children[feature['atts']['ID']] if c['seqid'] == feature['seqid']]))
                feature['end'] = str(max([int(c['end']) for c in children[feature['atts']['ID']] if c['seqid'] == feature['seqid']] ) )
    return features, children, uniquid

def remove_comments_rule( features, children, uniquid ):
    newfeatures = ['##gff-version 3\n']
    for feature in features:
        if type(feature) is dict:
            newfeatures.append( feature )
    return newfeatures, children, uniquid 


def clean_type_rule( features, children, uniquid ):
    newfeatures = []
    for feature in features:
        if type(feature) is dict: 
            if feature['type'] in ['singleexon', 'firstexon', 'lastexon']:
                feature['type'] = 'exon'
                newfeatures.append( feature )
            elif feature['type'] in ['transcript','.']:
                feature['type'] = 'mRNA'
                newfeatures.append( feature )
            else: # feature['type'] not in ['intron']:
                newfeatures.append( feature )
        else:
            newfeatures.append( feature )
    return newfeatures, by_key(newfeatures, 'Parent'), by_key( newfeatures, 'ID' )

def remove_bad_parents_rule(features, children, uniquid ):
    newfeatures = []
    for feature in features:
        if type( feature) is str or 'ID' not in feature['atts'] or feature['atts']['ID'] not in children:
            newfeatures.append( feature )
        elif all([child['strand'] == feature['strand'] and child['seqid'] == feature['seqid'] 
                  for child in children[feature['atts']['ID']]]) and any(child['type'] in ['mRNA', 'exon', 'CDS'] for child in children[feature['atts']['ID']]): 
            newfeatures.append( feature )
    return newfeatures, by_key( newfeatures, 'Parent'), by_key( newfeatures, 'ID' )

def remove_orphans_rule ( features, children, uniquid ):
    newfeatures = []
    for feature in features:
        if type( feature ) is str or not orphan_p( feature, features, children, uniquid ):
            newfeatures.append( feature )
    return newfeatures, by_key( newfeatures, 'Parent'), by_key( newfeatures, 'ID' )

def orphan_p(feature, features, children, uniquid ):
    if 'Parent' not in feature['atts']:
        return False
    elif any([parentID not in uniquid 
              for parentID in feature['atts']['Parent'].split(',')]):
        return True
    else:
        return any([orphan_p(parent, features, children, uniquid) 
                    for parentID in feature['atts']['Parent'].split(',') 
                    for parent in uniquid[parentID]])

def remove_big_genes(features, children, uniquid, MAX_INTRON=MAX_INTRON ):
    newfeatures = []
    for feature in features:
        if type(feature) is str or int(feature['end']) - int(feature['start']) < MAX_INTRON:
            newfeatures.append( feature )
    return newfeatures, by_key(newfeatures, 'Parent'), by_key(newfeatures, 'ID' )

def read_mapping( mappingfile ):
    mapping = {}
    if type(mappingfile) is str:
	mappingfile = open(mappingfile)
    for line in mappingfile:
        cols = [c.strip() for c in line.split(' ')]
        if len(cols) == 3:
            mapping[cols[1]] = cols[2]
    return mapping

def add_header_rule( features, children, uniquid ):
    gff_version = '##gff-version 3\n'
    if type(features[0]) is str and features[0] == gff_version:
        return features, children, uniquid
    print "feature[0]=" + str(features[0])
    features.insert(0, gff_version )
    return features, children, uniquid
        
def add_names(features, children, uniquid, MAPPING_FILE='mapping_old2new_names.txt'):
    mapping = read_mapping( MAPPING_FILE )
    for feature in features:
        if type(feature) is dict and feature['atts']['ID'] in mapping:
            feature['atts']['Name'] = mapping[feature['atts']['ID']]
            feature['attributes'] = ';'.join(['%s=%s' % (a, feature['atts'][a]) for a in sorted( feature['atts'])])
    return features, children, uniquid 

def add_uniquid( features, children, uniquid ):
    for feature in features:
        if type(feature) is dict and 'ID' not in feature['atts']:
            feature['atts']['ID'] = feature['atts']['Parent'] + ':' + feature['type'] 
            if feature['atts']['ID'] in uniquid:
                uniquid[feature['atts']['ID']].append( feature )
            else:
                uniquid[feature['atts']['ID']] = [feature]
    for u in uniquid:
        if len(uniquid[u]) > 1:
            count = 0
            if u in children:
                for child in children[u]:
                    child['atts']['Parent'].replace(u, ','.join(['%s:%d' % (u, c) for c in range(len(uniquid[u]))]))
                    child['attributes'] = ';'.join(['%s=%s' % (a, child['atts'][a]) for a in sorted(child['atts'])])
            for uu in uniquid[u]:
                uu['atts']['ID'] += ':%d' % count
                uu['attributes'] = ';'.join(['%s=%s' % (a, uu['atts'][a]) for a in sorted(uu['atts'])])
                count += 1
    return features, by_key(features, 'Parent'), by_key(features, 'ID')


def clean_gff2( gffin, gffout, rules ):
    features = read_gff( gffin )
    children = by_key( features,  'Parent' )
    uniquid = by_key( features, 'ID' )
    nfeatures = 1e100
    count = 0
    while nfeatures != len(features):
        print "Applying round %d to %d features:" % (count, len(features))
        count +=1 
        nfeatures = len(features)
        for desc, apply_rule in rules:
            lastfeatures = len(features)
            features, children, uniquid = apply_rule(features, children, uniquid )
            if len(features) > lastfeatures:
                print "\t%s added %d out of %d features (%0.2f)%%" % (desc, len(features)-lastfeatures,lastfeatures, len(features)/float(lastfeatures)*100.0 - 100)
            elif len(features) < lastfeatures:
                print "\t%s removed %d out of %d features (%0.2f)%%" % (desc, lastfeatures - len(features),lastfeatures, 100 - len(features)/float(lastfeatures)*100.0)
            else:
                print "\t%s did not add or remove features" %desc
    with gffout as out:
        for feature in features:
            if type( feature ) is dict:
                out.write( template % feature )
            else:
                out.write( feature )

def make_utr( feature, parents ):
    left_cds = min( [cds for cds in parents[feature['atts']['ID']] if cds['type'] == 'CDS'], key=lambda x: int(x['start']) )
    right_cds = max( [cds for cds in parents[feature['atts']['ID']] if cds['type'] == 'CDS'], key=lambda x: int(x['end']) )    
    utrs = []
    if feature['strand'] == '+':
        left_utr = 'five_prime_UTR'
        right_utr = 'three_prime_UTR'
    else:
        left_utr = 'three_prime_UTR'
        right_utr = 'five_prime_UTR'

    for exon in parents[feature['atts']['ID']]:
        if exon['type'] == 'exon':
            if int(exon['start']) < int(left_cds['start']):
                utr = copy.deepcopy( exon )
                utr['type'] = left_utr
                utr['atts']['ID'] = exon['atts']['ID'].replace('exon', left_utr)
                utr['attributes'] = ';'.join(['%s=%s' % (a, utr['atts'][a]) for a in utr['atts']])
                utr['end'] = str(min(int(left_cds['start']) -1, int(exon['end'])))
                utrs.append( utr )
            if int(exon['end']) > int(right_cds['end']):
                utr = copy.deepcopy( exon )
                utr['type'] = right_utr
                utr['atts']['ID'] = exon['atts']['ID'].replace('exon', right_utr )
                utr['attributes'] = ';'.join(['%s=%s' % (a, utr['atts'][a]) for a in utr['atts']])
                utr['start'] = str(max(int(right_cds['end']) +1, int(exon['start'])))
                utrs.append( utr )
    return utrs

def add_UTR_rule( features, children, uniquid ):
    newfeatures = []
    parents = by_key( features, 'Parent' )
    featuretypes = by_key( features, 'type' )
    if 'three_prime_UTR' not in featuretypes and 'five_prime_UTR' not in featuretypes and 'three_prime_utr' not in featuretypes and 'five_prime_utr' not in featuretypes:
        for feature in features:
            if type(feature) is dict and feature['type'] == 'mRNA':
                newfeatures.extend( make_utr( feature, parents ) )
            newfeatures.append( feature )
        return newfeatures, by_key( newfeatures, 'Parent'), by_key( newfeatures, 'ID' )
    else:
        return features, children, uniquid
            
def add_exons_rule( features, children, uniquid ):
    newfeatures = []
    featuretypes = by_key( features, 'type' )
    if 'exon' not in featuretypes:
        for feature in features:
            if type( feature ) is dict:
                if feature['type'] in ['CDS', 'five_prime_UTR', 'three_prime_UTR']:
                    exon = copy.deepcopy(feature)
                    exon['type'] = 'exon'
                    exon['atts']['ID'] = feature['atts']['ID'].replace('cds', 'exon')
                    exon['attributes'] = ';'.join(['%s=%s' % (a, exon['atts'][a]) for a in exon['atts']])
                    newfeatures.append( exon )
            newfeatures.append( feature )
        return newfeatures, by_key( newfeatures, 'Parent'), by_key( newfeatures, 'ID' )
    else:
        return features, children, uniquid


fields=['seqid','source','type','start','end','score','strand','phase','attributes']
template = '\t'.join(['%%(%s)s' % f for f in fields])  + '\n'

def write_gff( filename, features ):
    if type(filename) is str:
        filename = open(filename, 'w')
    with filename as out:
        for f in features:
            if type(f) is dict:
                out.write(template % f)
            else:
                out.write( f )


def validate_gff( ingff, outgff ):
    features = read_gff( ingff )
    children = by_key( features, 'Parent' )
    uniquid = by_key( features, 'ID' )
    MAX_INTRON = 2e5
    valid_types =  ['CDS', 'match_part', 'similarity', 'exon', 'mRNA', 'gene', 'five_prime_UTR', 'three_prime_UTR', 'stop_codon', 'start_codon']
    seqidRE = re.compile(r'[^a-zA-Z0-9\.:^\*\$@!\+\?-\|]+')
    with outgff as out:
	    for u in uniquid:
	        if len(uniquid[u]) != 1:
	            out.write( "Non-unique ids for %s\n" % u)
	    for uu in features:
	        if type(uu) is str:
	            if uu.strip()[0] != '#':
	                out.write( "uncommented illegal line: %s\n" % uu )
	        elif 'ID' not in uu['atts']:
	                out.write( "feature has no ID attribute: %s\n" % (template % uu) )
	        else:
	            u = uu['atts']['ID']
	            if seqidRE.search(uu['seqid']):
	                out.write( "Seqid for %s has unescaped illegal characters not in %s\n" % (u, seqidRE.pattern))
	            if int(uu['start']) >= int(uu['end']):
	                out.write( "%s has start %s >= end %s\n" % (u, uu['start'], uu['end']))
	            if uu['strand'] not in ['+', '-','?']:
	                out.write( "%s has illegal strand characters: %s\n" % (u, uu['strand']))
	            if uu['type'] not in valid_types:
	                out.write( "%s has an invalid sequence type %s. It is not one of %s.\n" % (u, uu['type'], ','.join(valid_types)))
	            if uu['type'] == 'CDS' and uu['phase'] not in ['0', '1','2']:
	                out.write( "%s has an illegal phase\n" % (u, uu['phase']))
	            
	            if int(uu['end']) - int(uu['start']) > MAX_INTRON:
	                out.write( "%s has an unusually large size (%d bp)\n" % (u, int(uu['end']) - int(uu['start'])))
                        
	    for p in children:
	
	        pstart = int(uniquid[p][0]['start'])
	        pend = int(uniquid[p][0]['end'])
	        pstrand = uniquid[p][0]['strand']
	        pseqid = uniquid[p][0]['seqid']
	        badseqids = [c['seqid'] for c in children[p] if c['seqid'] != pseqid]
	        badstrands = [c['strand'] for c in children[p] if c['strand'] != pstrand]
	        badstart = [int(c['start']) for c in children[p] if int(c['start']) < pstart]
	        badend = [int(c['end']) for c in children[p] if int(c['end']) > pend]
	        nchildren = len(children[p])        
	                     
	        if len(badseqids) > 0:
	            out.write( "%d of %d children mapped to different seqids than parent %s (%s): %s\n" % (len(badseqids), nchildren, p, pseqid, ','.join(set(badseqids))))
	        if len(badstrands) > 0:
	            out.write( "%d of %d children mapped to different strand than parent %s (%s): %s\n" % (len(badstrands), nchildren, p, pstrand, ','.join(set(badstrands))))
	        if len(badstart) > 0:
	            out.write( "%d of %d children start before parent %s (%d): %s\n" % (len(badstart), nchildren, p, pstart, ','.join(set([str(b) for b in badstart]))))
	        if len(badend) > 0:
	            out.write( "%d of %d children end after parent %s (%d): %s\n" % (len(badend), nchildren, p, pend, ','.join(set([str(b) for b in badend]))))
	        
	        for c in children[p]:
	            if c['seqid'] != pseqid:
	                out.write( "\t%s different seqid (%s) than parent %s (%s)\n" % (c['atts']['ID'], c['seqid'], p, pseqid))
	            if int(c['start']) < pstart:
	                out.write( "\t%s (%s) starts %d bases before parent %s (%d)\n" % (c['atts']['ID'],c['start'], pstart - int(c['start']),  p, pstart))
	            if int(c['end']) > pend:
	                out.write( "\t%s (%s) ends %d bases after parent %s (%d)\n" % (c['atts']['ID'], c['end'], int(c['end']) - pend, p, pend ))
	            if c['strand'] != pstrand:
	                out.write( "\t%s different strand (%s) from parent %s (%s)\n" % (c['atts']['ID'], c['strand'], p, pstrand ))
	        
    
def diff_gff( ff1, ff2, outfile ):
    gff1_features = read_gff( ff1 )
    gff2_features = read_gff( ff2 )
    gff2 = by_key( gff2_features, 'ID' )
    gff1 = by_key( gff1_features, 'ID' )
    difflist = ['features in gff2 not in gff1', 'features in gff1 not in gff2','duplicate_features', 'seqid','start', 'end']
    diffs = dict([(f, 0) for f in difflist])
    with outfile as out:

        for feature in gff1_features:
            if type(feature) is dict:
                featID = feature['atts']['ID']
                if featID not in gff2:
                    diffs['features in gff1 not in gff2'] +=1
                    out.write("%s\tFeature in gff1 missing in gff2: %s\n" % (featID, template % feature ))
                elif len( gff2[featID] ) > 1:
                    diffs['duplicate_features'] +=1
                    out.write("%s\tMultiple features in gff2 contain same ID as gff1:\n\t%s" % (featID, '\n\t'.join([template % g for g in gff2[featID]])))
                elif gff2[featID][0]['seqid'] != feature['seqid']:
                    diffs['seqid'] += 1
                    out.write("%s\t%s seqid in gff1 does not match %s seqid in gff2\n" % (featID, gff2[featID][0]['seqid'], feature['seqid'] ))
                elif gff2[featID][0]['start'] != feature['start']:
                    diffs['start'] += 1
                    out.write("%s\t%s start in gff1 does not match %s start in gff2\n" % (featID, gff2[featID][0]['start'], feature['start'] ))
                elif gff2[featID][0]['end'] != feature['end']:
                    diffs['end'] += 1
                    out.write("%s\t%s end in gff1 does not match %s end in gff2\n" % (featID, gff2[featID][0]['end'], feature['end'] ) )
        for features in gff2_features:
            if type(feature) is dict:
                featID = feature['atts']['ID']
                if featID not in gff1:
                    diffs['features in gff2 not in gff1'] += 1
        out.write('%s\t%s\n' % (ff1.name, ff2.name))
        for t in ['gene', 'mRNA', 'exon', 'CDS', 'five_prime_UTR', 'three_prime_UTR']:
            out.write('%s\t%d\t%d\n' % (t,len([g for g in gff1 if type(g) is dict and g['type'] == t]), len([gg for g in gff2 for gg in gff2[g] if gg['type'] == t]) ) )
        out.write('type\tdifferences\n')
        for d in difflist:
            out.write('%s\t%d\n' % (d, diffs[d]))
def by_key( features, k, empty=False ):
    d = {}
    for f in features:
        if type(f) is dict:
            if k in f:
                if f[k] in d:
                    d[f[k]].append( f )
                else:
                    d[f[k]] = [f]
            elif 'atts' in f and k in f['atts']:
                for kk in f['atts'][k].split(','):
                    if kk in d:
                        d[kk].append( f )
                    else:
                        d[kk] = [f]
            elif empty and k in d:
                d[k].append(f)
            elif empty:
                d[k] = [f]
                
            
    return d

def exonerategff2_to_gff3( exonerategff, outgff3 ):
    semicolonRE = re.compile(r'''((?:[^;"']|"[^"]*"|'[^']*')+)''')
    equalsRE =  re.compile(r'''((?:[^="']|"[^"]*"|'[^']*')+)''')
    commaRE = re.compile(r'''((?:[^,"']|"[^"]*"|'[^']*')+)''')
    ID = {}
    with outgff3 as out:
        out.write('##gff-version 3\n')
        for line in exonerategff:
            col = [c.strip() for c in line.split('\t')]
            if len(col) == len(fields):
                feature = {}
                for i, field in enumerate(fields):
                    feature[field] = col[i].strip()
                    if field == 'attributes':
                        feature['atts'] = {}
                        for attval in semicolonRE.split(feature[field])[1::2]:
                            att, vals = attval.strip().split(' ', 1)
                            if att in feature['atts']:
                                feature['atts'][att] += ',' + vals
                            else:
                                feature['atts'][att] = vals
                if 'Query' in feature['atts']:
                    feature['attributes'] = 'ID=%s' % feature['atts']['Query']
                out.write(template % feature )

def parse_classification(classfile):
    t = {}
    mips = {}
    if type(classfile) is str:
        classfile = open(classfile)
    for line in classfile:
        if line[0] == '>':
            if line.find('#') != -1:
                target, rest = line[1:-1].split('#', 1)
                rmclass, notes = rest.split(' ',1)
                notes = notes.strip('()')
            else:
                notes = line[1:-1]
                target = notes.split('|')[0]
                rmclass = notes.split('|')[2]
            t[target]  = {'Name': rmclass, 'Note': urllib.quote(rmclass + ' ' + target + ' ' + notes) }
            c = notes.split('|')
            if len(c) == 4:
                mips[target] = {'MIPS': c[1], 'Name': c[0], 'Note': urllib.quote(rmclass + ' ' + target + ' ' + notes) }
    return t, mips

def parse_rmgff( gff, classfile={}, fields=['seqid','source','type','start','end','score','strand','phase','attributes']):
    ssr = ['##gff-version 3\n']
    retro = ['##gff-version 3\n']
    semicolonRE = re.compile(r'''((?:[^;"']|"[^"]*"|'[^']*')+)''')
    equalsRE =  re.compile(r'''((?:[^="']|"[^"]*"|'[^']*')+)''')
    commaRE = re.compile(r'''((?:[^,"']|"[^"]*"|'[^']*')+)''')
    spaceRE = re.compile(r'''((?:[^ "']|"[^"]*"|'[^']*')+)''')
    if type(gff) is str:
        gff = open(gff)
    uniquids = []
    count = 1
    for line in gff:
        col = [c.strip() for c in line.split('\t')]
        if len(col) == len(fields):
            feature = {}
            for i, field in enumerate(fields):
                feature[field] = col[i].strip()
                if field == 'attributes':
                    feature['atts'] = {}
                    uniquid = 'retroelement:%s:%d-%d' % (feature['seqid'], int(feature['start']), int(feature['end']))
                    if uniquid in uniquids:
                        uniquid += ':%d' % count
                        count += 1
                        uniquids.append( uniquid )
                    feature['atts']['ID'] = uniquid
                    for attval in semicolonRE.split(feature[field])[1::2]:
                        if attval.strip() != '':
                            att, val = attval.split(' ', 1)
                            motif, start, stop = spaceRE.split(val)[1::2]
                            motif, target = motif.strip('"').split(':')
                            feature['atts'][att] = target
                            if target in classfile:
                                for k in classfile[target]:
                                    feature['atts'][k] = classfile[target][k]
                    feature[field] = ';'.join(['%s=%s' % (a,feature['atts'][a]) for a in feature['atts']])
            if feature['atts']['Target'].split('|')[0] in classfile:
                retro.append(feature)
            else:
                ssr.append(feature)
    return retro, ssr

def add_ID( features ):
    count = 0
    for f in features:
        if type(f) is dict and 'ID' not in f['atts']:
            f['atts']['ID'] = ':'.join([f['atts']['Parent'], f['type'], str(count)])
            f['attributes'] = ';'.join(['%s=%s' % (a, f['atts'][a]) for a in f['atts']])
            count += 1
    return features
            
def read_gff( gff, fields=['seqid','source','type','start','end','score','strand','phase','attributes']):
    features = []
    semicolonRE = re.compile(r'''((?:[^;"']|"[^"]*"|'[^']*')+)''')
    equalsRE =  re.compile(r'''((?:[^="']|"[^"]*"|'[^']*')+)''')
    commaRE = re.compile(r'''((?:[^,"']|"[^"]*"|'[^']*')+)''')
    noteRE = re.compile(r'''Note=([^;]+)(;?)''')
    urlRE = re.compile(r"[^A-Za-z0-9%._/]")
    attsRE = re.compile(r'''Note=|ID=|Name=|Alias=|Parent=|Target=|Gap=|Derives_from=|Dbxref=|Ontology_term=|Is_circular=|exontype=|orf_start=|orf_end=|full_length_transcript=|_AED=|_eAED=|_QI=|src=|grp=|mult=|pri=|aaloc=''')
    for line in gff:
        col = [c.strip() for c in line.split('\t')]
        if len(col) == len(fields) and col[0][0] != '#':
            feature = {}
            for i, field in enumerate(fields):
                feature[field] = col[i].strip()
                if field == 'attributes':
                    atts = attsRE.findall( feature[field] )
                    feature['atts'] = {}
                    for i, val in enumerate(attsRE.split(feature[field])):
                        if i > 0:
                            val = val.strip(';')
                            if atts[i-1] in [ 'Note=', 'Name='] and urlRE.search( val ):
                                val = urllib.quote( val )
                            att = atts[i-1].strip('=')
                            feature['atts'][att] = val
                    feature[field] = ';'.join(['%s=%s' % (a, feature['atts'][a]) for a in feature['atts']])
            features.append( feature )
        else:
            features.append(line)
    return features

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='gffutils')
    parser.add_argument('--validate', dest='gff', type=argparse.FileType('r'), default=None)
    args = parser.parse_args()
    validate_gff( args.gff)
        
