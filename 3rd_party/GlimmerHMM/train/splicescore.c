//Copyright (c) 2003  by  Mihaela Pertea.

// This is a GHMM of GlimmerM 
// mpertea@tigr.org

#include  "delcher.h"
#include  "gene.h"

//#define MESSAGE

const unsigned int  OK = 0x0;

// ******************** FUNCTION DECLARATIONS ***********************

int  Is_Acceptor  (const int *, double *, double, int);
int  Is_Donor  (const int *, double *, double, int);
int  Is_Atg  (const int *, double *,double);
int  Is_Stop  (const int *, double *,double);

// ********************* VARIABLE DECLARATIONS **********************

char  Data[1000000];
long int Data_Len;
char TRAIN_DIR[500]="";
double Acc_Thr = -99;
double Don_Thr = -99;
double Atg_Thr = -99;
double Stop_Thr = -99;

// ************************** Main PROGRAM ***************************

int main  (int argc, char * argv [])
{
  FILE  * fp;
  char  File_Name [MAX_LINE], Name [MAX_LINE];
  int ret, istacc, istdon;
  int i,j,k,Input_Size,len;
  int B[200];
  double score;


  // ---------------------------- OPTIONS ------------------------------

  if  (argc < 9) {
    fprintf (stderr,
	     "USAGE:  %s <multi-fasta-file> <acc_thr> <istacc> <don_thr> <istdon> <atg_thr> <stop_thr> <train_dir> \n",
	     argv [0]);
    exit (-1);
  }

  Acc_Thr = strtod (argv[2],NULL);
  istacc = atoi(argv[3]);
  Don_Thr = strtod (argv[4],NULL);
  istdon = atoi(argv[5]);
  Atg_Thr = strtod (argv[6],NULL);
  Stop_Thr = strtod (argv[7],NULL);
  strcpy(TRAIN_DIR,argv[8]);
  
  //printf("%f %d %f %d %f %f %s\n", Acc_Thr,istacc ,Don_Thr,istdon,Atg_Thr,Stop_Thr,TRAIN_DIR);
  

  // --------------------------- READ DATA ------------------------------


  fp = File_Open (argv [1], "r");
  
  while(fgets(Data,999999,fp)) {

    i=0;
    while(Data[i]!=' ') { i++;}
    strncpy(Name,Data,i); Name[i]='\0';
    
    strcpy(Data,Data+i);

    //printf("Name=%s\n",Name);
    //printf("Data=%s\n",Data);exit(0);
    
    Data_Len = strlen (Data + 1);

    for  (i = 1;  i <= Data_Len;  i ++) {
      // Converts all characters to  acgt
      Data [i] = Filter (tolower (Data [i]));
    }

    for(i=1;i<Data_Len+1;i++) {
      
      /* atg's : 0 */

      if(i<Data_Len-1 && Data[i]=='a' && Data[i+1]=='t' && Data[i+2]=='g') { // Deal w/ start sites
      
	if(i>12 && i<=Data_Len-6) {
	  k=0;
	  for(j=i-12;j<i+7;j++) {
	    switch (Data[j]){
	    case 'a': B[k]=0;break;
	    case 'c': B[k]=1;break;
	    case 'g': B[k]=2;break;
	    case 't': B[k]=3;break;
	    }
	    k++;
	  }
	  ret=Is_Atg(B,&score,Atg_Thr);
	  if(ret) printf("%s atg %d %f\n",Name,i,score);
	  
	}
      }

      /* gt's : 1 */

      if(i<Data_Len && Data[i]=='g' && Data[i+1]=='t') { // Deal with donors
	
	if(i>80 && i<=Data_Len-81) {
	  k=0;
	  for(j=i-80;j<i+82;j++){
	    switch (Data[j]){
	    case 'a': B[k]=0;break;
	    case 'c': B[k]=1;break;
	    case 'g': B[k]=2;break;
	    case 't': B[k]=3;break;
	    }
	    k++;
	  }

	  ret=Is_Donor(B,&score, Don_Thr, istdon);
	  if(ret) printf("%s gt %d %f\n",Name,i,score);
	}
      }
      
      /* ag's : 2 */

      if(i<Data_Len && Data[i]=='a' && Data[i+1]=='g') { // Deal with acceptors

	if(i>80 && i<=Data_Len-81) {
	  k=0;
	  for(j=i-80;j<i+82;j++){
	    switch (Data[j]){
	    case 'a': B[k]=0;break;
	    case 'c': B[k]=1;break;
	    case 'g': B[k]=2;break;
	    case 't': B[k]=3;break;
	    }
	    k++;
	  }
	  
	  ret=Is_Acceptor(B,&score, Acc_Thr, istacc);
	  if(ret) printf("%s ag %d %f\n",Name,i,score);
	}
      }

      /* stop codons : 3  only if they are in frame 0 */
    
      if( i < Data_Len-1 && (i-1)%3==0 && 
	  ((Data[i]=='t' && Data[i+1]=='a' && Data[i+2]=='a') ||
	   (Data[i]=='t' && Data[i+1]=='g' && Data[i+2]=='a') ||
	   (Data[i]=='t' && Data[i+1]=='a' && Data[i+2]=='g'))) {

	score=-99;
	
	if(i>4 && i<=Data_Len-14) {
	  k=0;
	  for(j=i-4;j<i+15;j++) {
	    switch (Data[j]){
	    case 'a': B[k]=0;break;
	    case 'c': B[k]=1;break;
	    case 'g': B[k]=2;break;
	    case 't': B[k]=3;break;
	    }
	    k++;
	  }
	  
	  
	  ret=Is_Stop(B,&score,Stop_Thr);
	}
	printf("%s stop %d %f\n",Name,i,score);
      }
    }

  }
  fclose (fp);


}

