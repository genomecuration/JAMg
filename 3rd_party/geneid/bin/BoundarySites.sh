#!/bin/sh
####################################################
#
#  NAME:
#    BoundarySites.sh (sh)
#  AUTHOR:
#    Moises Burset
#  HISTORY:
#    Version 1.0 -> 11-08-98
#    Last modification -> 27-09-98
#  DESCRIPTION:
#    This program obtain small pieces from sequence in sites
#    boundary
#  INPUT:
#    variables:
#    parameters:
#      $1 -> Pattern length
#      $2 -> Lateral extension from pattern limits in both directions
#    files:
#      <STDIN> -> label seq pos1 pos2 ... posn
#      where each pos is one site position
#  OUTPUT:
#    files:
#      <STDOUT> -> label pos1 seqsite1 pos2 seqsite2 ... posn seqsiten
#      where each seqsite is the sequence bounding the site
#  RETURN:
#  DEPENDENCES:
#  BUGS:
#  COMMENTS:
#    If you try to obtain nucleotides bounding some site without
#    sufficient nucleotides bounding this site, the program add the
#    needed n's
####################################################

gawk 'BEGIN{
  len=ARGV[1]
  ext=ARGV[2]
  ARGV[1]=ARGV[2]=""
}
{
  a = ""
  for (i=3;i<=NF;i++) {
    ini = $i-ext
    fin = len+(ext*2)
    a = a $i" "
    if (ini < 1) {
      ini = 1
      for (j=1;j<=(ext-$i+1);j++) {
        a = a "n"
        fin--
      }
    }
    if ($i+len+ext >= length($2)) {
     # for (j=1;j<=(ext-(length($2)-$i-1));j++) {
	for (j=1;j<=(ext-(length($2)-$i-len+1));j++) {
        fin--
        b = b "n"
      }
    }
    a = a substr($2,ini,fin)
    a = a b " "
    b = ""
  }
  print $1,a
}' "$@"
