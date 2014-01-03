<?php

$id = '';

if (!empty($_REQUEST)){
    if (empty($_REQUEST['id'])){
      $id = filter_var($$_REQUEST['id'], FILTER_SANITIZE_STRING, FILTER_FLAG_STRIP_LOW | FILTER_FLAG_STRIP_HIGH);
    }
}
if (empty($id)){ $id = 'UniRef90_P23443'; }
print "ID is $id\n";
// this ought to happen once only
$dbconn = connect_annotationdb();

if (!empty($dbconn)){
  print "Connection succeeded\n";
}else{
  print "Connection failed\n";
}

$return = show_annotations($dbconn,$id);
print_r($return);
//echo $return;


//////////////////////////////////////////
function show_annotations($dbconn,$id){
  $arr = array();

  if (empty($id)){return FALSE;}
  // these ought to be prepared once only outside ws
  $ids = array();
  $protein_cluster_ids_sql = "SELECT distinct uniprot_id as id from uniprot_assoc where xref_db='UniRef90' and xref_accession='$id'";
  $res = pg_query($dbconn, $protein_cluster_ids_sql);
  while ($row = pg_fetch_array($res)){
    $ids[]="'".$row['id']."'";
  }

  $ec_annotations_sql = 'SELECT distinct ec_id as id FROM enzyme_assoc where uniprot_id IN ';
  get_terms($dbconn,$ec_annotations_sql,'ec',$ids,$arr);
  $go_annotations_sql = 'SELECT distinct go_id as id from go_assoc where uniprot_id IN ';
  get_terms($dbconn,$go_annotations_sql,'go',$ids,$arr);
  $kegg_annotations_sql = 'SELECT distinct pathway_id as id from kegg_pathway_assoc where uniprot_id IN ';
  get_terms($dbconn,$kegg_annotations_sql,'kegg',$ids,$arr);

  return $arr;
}

function connect_annotationdb(){
   //all connections go via pgbouncer 
   $pass_file = '/var/lib/tomcat7/conf/db_connection.sql';
   include("$pass_file");
   $dbconn = pg_connect("host=$host port=$port dbname=$dbname user=$user password=$pass");
   return $dbconn;
}

function get_terms($dbconn,$sql,$cv,$ids,&$arr){

  $res = pg_query($dbconn, $sql.'('.implode(',',$ids).')');
  while ($row = pg_fetch_array($res)){
     $arr[$cv][] = $row['id'];
  }

}

?>
