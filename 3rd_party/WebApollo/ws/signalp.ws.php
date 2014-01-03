<?php

$sequence = '';

if (!empty($_REQUEST)){
    if (!empty($_REQUEST['sequence'])){
      $sequence = $_REQUEST['sequence'];
    }
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

$return = blast_sequences($sequence);
echo $return;
/////////////////////////////////////////////////////////////////////////////////////
function blast_sequences($sequence){

  if (empty($sequence)){ return FALSE; }
   $sequence_file = tempnam('tmp','seq_');
   $sequence_filename = basename($sequence_file);
   $sequence_file_handle = fopen($sequence_file,'w');
   fwrite($sequence_file_handle, $sequence."\n");
   fclose($sequence_file_handle);
   $cmd = "signalp-4.1/signalp -f summary -g png -t euk -w $sequence_file > $sequence_file.result ";
   system($cmd." 2>$sequence_file.err",$rtn);
   $response = '';
   if (file_exists("$sequence_file.result") && filesize("$sequence_file.result")>0){
     $response = file_get_contents("$sequence_file.result");
   }else{
     $response = "<pre>".file_get_contents("$sequence_file.err")."</pre>";
   }
  unlink($sequence_file);
  return $response;

}

?>
