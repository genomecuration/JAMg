<?php

$sequence = '';
$family = '';

if (!empty($_REQUEST)){
    if (!empty($_REQUEST['sequence'])){
      $sequence = $_REQUEST['sequence'];
    }
    if (!empty($_REQUEST['family'])){
      $family = $_REQUEST['family'];
   }
}
if (empty($family)){$family = 'Harm_OR'; }
if (empty($sequence)){
 $sequence = '>TEST
MPPNTRSLLRKICRYVYYAGAGNCWYEDIYHETKPYKLYSIVTFSIYTIMIFLENVAALF
GKFPEVEKNSAVMFSAIHDIVLTKMFLLLYHKQSIRKLNYDMSTVGSSFEEDHVMRKQYL
KTTVGIWLYVISVYLSLGAGCTTLLLHSIVEGTPFYTVVTYLPFYDDNSVVALIFRIFFY
ITWLYMMLPMMSADCMPITHLITMTYKFITLCHHFERIRTEFDEDTKIMNRREAIDKLRA
GCLEGIRMHQKLLWLADEIHRVFGIIMSLQVCESSAVAVLLLLRLALSPHLDLTNAFMTY
TFVCSLFLLLALNLWNAGEVTYQASLLSNAMFQCGWHLCELEKHNHRDIRRLVLIGCTQA
QKPLILKAFGIQDLSYETFVSVARMTYSIFAVFYQRGE
';
}

$return = align_sequences($sequence,$family);
echo $return;
/////////////////////////////////////////////////////////////////////////////////////
function align_sequences($sequence,$family){

  $family_file = 'profiles/'.$family.'.mafft';
  if (empty($family)){return FALSE;}
  if (!file_exists($family_file)){return FALSE; }
  if (empty($sequence)){ return FALSE; }


   $sequence_file = tempnam('tmp','seq_');
   $sequence_filename = basename($sequence_file);
   $sequence_file_handle = fopen($sequence_file,'w');
   fwrite($sequence_file_handle, $sequence."\n");
   fclose($sequence_file_handle);
   $cmd = "mafft_bin/mafft --amino --anysymbol --globalpair --maxiterate 1000 --seed $family_file $sequence_file > $sequence_file.result";
   system($cmd." 2>$sequence_file.err",$rtn);
   $response = '';
   if (file_exists("$sequence_file.result") && filesize("$sequence_file.result")>0){
     system("sed -i -r '~s/^>_seed__[a-z0-9]+_[a-z0-9]+_[a-z0-9]+_/>/' $sequence_file.result");
     system("sed -i -r '~s/^>_[a-z0-9]+_[a-z0-9]+_[a-z0-9]+_/>/' $sequence_file.result");
     system("mafft_bin/f2cl < $sequence_file.result > $sequence_file.result.clu");
//     system("fprotpars -sequence $sequence_file.result -outfile $sequence_file.result.parsimony -auto -noprogress -outtreefile $sequence_file.result.parsimony.newick -njumble 10 ");
//<h2>Parsimony tree</h2>NB: Numbers on tree nodes are not bootstrap, just unique identifiers for each node<pre>".file_get_contents("$sequence_file.result.parsimony")."</pre>
   $response = "
<h2>Alignment using a profile created from the $family gene family</h2>
<h3>View with JalView</h3>
<p>first time it needs browser's JAVA plugin to be activated (<a target='_blank' href='/ws/pointout.png'>see here for example on firefox</a>)</p>
<applet name='jal'
width='150' height='35'
code='jalview.bin.JalviewLite' archive='/ws/jalviewApplet.jar'>
<PARAM NAME='file' value='/ws/tmp/$sequence_filename.result'>
<PARAM NAME='userDefinedColour' value='K,R,H=00eeaa; L,V,I,M=ff2b00; Y,F,W=ff1bff; A,G,S,T,P=cccccc; C=yellow; D,E,N,Q,B,Z=00aaee; k,r,h=00eeaa; l,v,i,m=ff2b00; y,f,w=ff1bff; a,g,s,t,P=cccccc; c=yellow; d,e,n,q,b,z=00aaee;'>
<PARAM NAME='showFullId' value='false'>
<param name='APPLICATION_URL' value='http://www.jalview.org/services/launchApp'>
</applet>
<h2>Alignment raw data</h2>
<h3>CLUSTAL format</h3>
<pre>".file_get_contents("$sequence_file.result.clu")."</pre>
<h3>FASTA format</h3>
<pre>" .file_get_contents("$sequence_file.result")."</pre>
<h2>Target sequence aligned</h2>
<pre>" .file_get_contents("$sequence_file")."</pre>
   ";
}else{
     $response = "<pre>".file_get_contents("$sequence_file.err")."</pre>";
}
  unlink($sequence_file);
  return $response;

}

?>
