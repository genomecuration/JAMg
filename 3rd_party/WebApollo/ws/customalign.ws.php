<?php
$sequence = '';
$custom_profile = '';

if (!empty($_REQUEST)){
    if (!empty($_REQUEST['sequence'])){
      $sequence = $_REQUEST['sequence'];
    }
    if (!empty($_REQUEST['custom_profile'])){
      $custom_profile = $_REQUEST['custom_profile'];
   }
}
if (empty($custom_profile)){
 $custom_profile = ">profile1\nMPPNTRSLLRKICRYVYYAGAGNCWYEDIYHETKPYKLYSIVTFSIYTIMIFLENVAALFGKFPEVEKNSAVMFSAIHDIVLTKMFLLLYHKQSIRKLNYDMSTVGSSFEEDHVMRKQYLKTTVGIWLYVISVYLSLGAGCTTLLLHSIVEGTPFYTVVTYLPFYDDNSVVALIFRIFFY\n"
 .">profile2\nMPPNTRSL------YVYYAGAGNCWYEDIYHETKPYKLYSIVTFSIYTIMIFLENVAALFGKFPEVEAGGGGGYSAIHDIVLTKMFLLLYHKQSIRKLNYDMSTVGSSFEEDHVMRKQYLKTTVGIWLYVISVYLSLGAGCTTLLLHSIVEGTPFYTVVTYLPFYDDNSVVALIFRIFFY";
}
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

$method = 'unaligned';
if (!empty($_REQUEST['prealigned'])) {
  $method = 'prealigned';
}else{
   $array = explode("\n",$custom_profile);
   foreach ($array as &$ln){
        if (strpos($ln,'>') === FALSE){
		$ln = preg_replace('/\-+/','',$ln);
	}
   }
   $custom_profile = implode("\n",$array);
}

$return = align_sequences($sequence,$custom_profile,$method);
echo $return;
/////////////////////////////////////////////////////////////////////////////////////
function align_sequences($sequence,$custom_profile,$method){
  if (empty($custom_profile)){return FALSE;}
  if (empty($sequence)){ return FALSE; }
  if (empty($method)){ $method = 'unaligned';}
  $cmd = '';
  $response = '';
   $sequence_file = tempnam('tmp','seq_');
   $sequence_filename = basename($sequence_file);
   $sequence_file_handle = fopen($sequence_file,'w');
   fwrite($sequence_file_handle, $sequence."\n");
   fclose($sequence_file_handle);
   $profile_file = $sequence_file.'.profile';
   $profile_file_handle = fopen($profile_file,'w');
   fwrite($profile_file_handle, $custom_profile."\n");
   fclose($profile_file_handle);
  if ($method == 'prealigned'){
   $response ="<h2>Global alignment using a custom profile provided by user</h2>";
   $cmd = "mafft_bin/mafft --amino --anysymbol --globalpair --maxiterate 1000 --seed $profile_file $sequence_file > $sequence_file.result";
  }else{
   $response="<h2>Global alignment using unaligned sequences provided by user</h2>";
   $cmd = "cat $sequence_file $profile_file | mafft_bin/mafft --amino --anysymbol --globalpair --maxiterate 1000 -  > $sequence_file.result";
  }

   system($cmd." 2>$sequence_file.err",$rtn);
   if (file_exists("$sequence_file.result") && filesize("$sequence_file.result")>0){
     if ($method == 'prealigned'){
       system("sed -i -r '~s/^>_seed__[a-z0-9]+_[a-z0-9]+_[a-z0-9]+_/>/' $sequence_file.result");
     }
     system("sed -i -r '~s/^>_[a-z0-9]+_[a-z0-9]+_[a-z0-9]+_/>/' $sequence_file.result");
     system("mafft_bin/f2cl < $sequence_file.result > $sequence_file.result.clu");
   $response .= "
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
