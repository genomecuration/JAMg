/******************************************************************\
* This program evaluate gene predictions in command-line or in    *
* cgi script, you can use your database of real sequences in      *
* format ==> <LocusName> <LocusLength> <Real Positions>           *
* to define where is the database you must specify the enviroment *
* variable "EvaluationDB" or you can modify the default value of  *
* variable in this source code. (begin of main function)          *
*                                                                 *
* BUGS: There are some problems with last genes without           *
*       prediction due to ReadWord routine probably               *
\*****************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <stdarg.h>
#include <string.h>

#define LOCUSLENGTH  100     /*****  Max length of locus name  *****/
#define MAXEXONS    80000*2    /*****  Max number of exons       *****/
#define MAXGENES    5     /*****  Max number of genes       *****/
#define MAXWORD     1000     /*****  Max characters in a word  *****/

typedef struct s_gene {
  char LocusName[LOCUSLENGTH];
  long Reals[MAXEXONS];
  long Predicts[MAXEXONS];
  long LongSeq;
} gene;

void error (char *, ...);
void warning (char *, ...);
char ReadWord (char *,FILE *,long);
void OrderPositions (long *);
char LoadGenes(char *,long *,long,FILE *);
void CheckErrors(gene *,long);
long CalculateTp(long *,long *);
void CalculateExons(long *,long *,int *,int *,int *,int *,int *);
double sqrt(double);            /* defined in math libs */

int main (int argc,char *argv[]) { 

  gene Genes[MAXGENES];
  long nContent,i,j,nGenes=0;
  char *Method, *CharTmp, Separations[MAXWORD];
  char Status=' ',Vector[MAXWORD];
  int Length;

  FILE *Database;
  char tmpdata[MAXWORD];
  char data[MAXWORD]="/usersgenome2/mburset/public_html/evaluation",*pdata;
  long Tp,Tn,Fp,Fn,CDSr,CDSp;          /* variables for evaluation */
  int Tpe,Fpe,Fne,OverR,OverP,RealEx,PredEx;
  double sTp=0,sTn=0,sFp=0,sFn=0,sCDSr=0,sCDSp=0,sLongSeq=0;
  double sTpe=0,sFpe=0,sFne=0,sRealEx=0,sPredEx=0;
  double Cc,Angle,SnSp,Sne,Spe,S1,S2,S3,S4,ME,WE;
  double sCc=0,sAngle=0,sSnSp=0,sSne=0,sSpe=0,sS1=0,sS2=0,sS3=0,sS4=0;
  double sME=0,sWE=0;
  long cCc=0,cSnSp=0,cSne=0,cSpe=0,cS1=0,cS2=0,cS3=0,cS4=0,k=0;
  double T,RC,RN,PC,PN,TP,TN,FP,FN,SEN1,SEN2,CC,ANGLE,SNE,SPE;

  for (i=0;i<MAXGENES;i++) {     /* Initialize genes */
    Genes[i].LocusName[0]='\0';
    Genes[i].LongSeq=0;
    for (j=0;j<MAXEXONS;j++) {
      Genes[i].Reals[j]=-1;
      Genes[i].Predicts[j]=-1;
    }
  }
  Method = getenv("REQUEST_METHOD");   /* Get METHOD */
  if (Method == NULL)
    Method = "nothing";

  if (strcmp(Method,"POST")) {        /* if METHOD != POST */
    FILE *R;

    if (argc!=2)
      error("Incorrect number of arguments\n");
    if ((R=fopen(argv[1],"rb"))==NULL)
      error("I cannot open reals file %s\n",argv[1]);
    
    while (Status!='f') {            /* Load reals */
      Status=ReadWord(Genes[nGenes].LocusName,R,LONG_MAX);
      Status=ReadWord(Vector,R,LONG_MAX);
      Genes[nGenes].LongSeq=atol(Vector);
      if (Genes[nGenes].LongSeq==0)
        error("Incorrect length for %s gene in reals file\n",Genes[nGenes].LocusName);
      if (Status!='n' && Status!='f')
	Status=LoadGenes(Genes[nGenes].LocusName,Genes[nGenes].Reals,LONG_MAX,R);
      else
	Genes[nGenes].Reals[0]=0;
      nGenes++;
    }
    fclose(R);

    Status=' ';
    while (Status!='f') {            /* Load predicts */
      Status=ReadWord(Vector,stdin,LONG_MAX);
      for (i=0;i<=nGenes;i++)
        if (strcmp(Vector,Genes[i].LocusName) == 0) {
	  if (Status!='f' && Status!='n')
	    Status=LoadGenes(Genes[i].LocusName,Genes[i].Predicts,LONG_MAX,stdin);
	  else
	    Genes[i].Predicts[0]=0;
          break;
        }
        else
	  warning("More than one gene for %s locusname!!!",Vector);
      if ((i-1) == nGenes) {
        strcpy(Genes[nGenes].LocusName,Vector);
	if (Status!='n' && Status!='f')
	  Status=LoadGenes(Genes[nGenes].LocusName,Genes[nGenes].Predicts,LONG_MAX,stdin);
	else
	  Genes[i].Predicts[0]=0;
	nGenes++;
      }      
    }
  }
  else {          /* if METHOD == POST */
    printf ("Content-type: text/html\n\n");          /* Header of HTML */
    printf ("<HTML><HEAD><TITLE>Evaluation of gene prediction programs</TITLE></HEAD><BODY>\n");
    printf ("<center>\n");
    printf ("<h2>THE RESULTS OF YOUR QUERY ARE:</h2>\n");
    printf ("</center>\n");
    printf ("<pre>\n");

    pdata = getenv("Evaluation");
    if (pdata != NULL)
      strcpy(data,pdata);
    sprintf(tmpdata,"%s/Reals",data);
    if ((Database=fopen(tmpdata,"rb"))==NULL) {   /* Load database */
      printf ("Our database of real positions aren't disponible, \n");
      printf ("please contact us for more information\n");
    }
    else {
      while (Status!='f') {
        ReadWord(Genes[nGenes].LocusName,Database,LONG_MAX);
        ReadWord(Vector,Database,LONG_MAX);
        Genes[nGenes].LongSeq=atol(Vector);
        if (Genes[nGenes].LongSeq==0)
          error("Reals database internal error in %s\n",Genes[nGenes].LocusName);
        Status=LoadGenes(Genes[nGenes].LocusName,Genes[nGenes].Reals,LONG_MAX,Database);
        nGenes++;
      }
    fclose(Database);
    }

    Status=' ';
    ReadWord("Initialize",stdin,1);   /* Initialize Positions */
    CharTmp=getenv("CONTENT_LENGTH");
    nContent=atol(CharTmp);
    ReadWord(Separations,stdin,nContent);
    Length=strlen(Separations);

    while ((Status=ReadWord(Vector,stdin,nContent)) != 'f') {   /* Main post loop */

      if (!strcmp(Vector,"name=\"p\"") || !strcmp(Vector,"name=\"p\";")) {   /* Read predicts */
	while(Status!='n')
	  Status=ReadWord(Vector,stdin,nContent);
        Status=ReadWord(Vector,stdin,nContent);
        while (strncmp(Separations,Vector,Length)) {
          for (i=0;i<=nGenes;i++)
            if (!strcmp(Vector,Genes[i].LocusName)) {
              Status=LoadGenes(Genes[i].LocusName,Genes[i].Predicts,nContent,stdin);
              break;
            }
          if ((i-1) == nGenes) {
            strcpy(Genes[nGenes].LocusName,Vector);
            Status=LoadGenes(Genes[nGenes].LocusName,Genes[nGenes].Predicts,nContent,stdin);
            nGenes++;
          }
        if (Status != 'f')
          Status=ReadWord(Vector,stdin,nContent);
        }      /* close the while of load genes */
      }        /* close the while of read predicts */

      if (!strcmp(Vector,"name=\"r\"") || !strcmp(Vector,"name=\"r\";")) {   /* Read reals */
	while(Status!='n')
	  Status=ReadWord(Vector,stdin,nContent);
        Status=ReadWord(Vector,stdin,nContent);
        while (strncmp(Separations,Vector,Length)) {
          for (i=0;i<=nGenes;i++)
            if (!strcmp(Vector,Genes[i].LocusName)) {
              Status=LoadGenes(Genes[i].LocusName,Genes[i].Reals,nContent,stdin);
              break;
            }
          if ((i-1) == nGenes) {
            strcpy(Genes[nGenes].LocusName,Vector);
            Status=LoadGenes(Genes[nGenes].LocusName,Genes[nGenes].Reals,nContent,stdin);
            nGenes++;
          }
        if (Status != 'f')
          Status=ReadWord(Vector,stdin,nContent);
        }      /* close the while of load genes */
      }        /* close the while of read reals */


      if (!strcmp(Vector,"name=\"l\"") || !strcmp(Vector,"name=\"l\";")) {   /* Read LongSeq */
        while (Status!='n')
          Status=ReadWord(Vector,stdin,nContent);
        Status=ReadWord(Vector,stdin,nContent);
        while (strncmp(Separations,Vector,Length)) {
          for (i=0;i<=nGenes;i++)
            if (!strcmp(Vector,Genes[i].LocusName)) {
              Status=ReadWord(Vector,stdin,nContent);
              Genes[i].LongSeq=atol(Vector);
              if (Genes[i].LongSeq == 0)
                error("Incorrect locus length for %s gene\n",Genes[i].LocusName);
              if (Status != 'n')
                error("More than sequence length in %s gene\n",Genes[i].LocusName);
              break;
            }
          if ((i-1) == nGenes) {
            strcpy(Genes[nGenes].LocusName,Vector);
            Status=ReadWord(Vector,stdin,nContent);
            Genes[i].LongSeq=atol(Vector);
            if (Genes[i].LongSeq == 0)
	      error("Incorrect locus length for %s gene\n",Genes[i].LocusName);
            if (Status != 'n')
	      error("More than sequence length in %s gene\n",Genes[i].LocusName);
            nGenes++;
          }
        if (Status != 'f')
          Status=ReadWord(Vector,stdin,nContent);
        }      /* close the while of load genes */
      }        /* close the while of read reals */
    }          /* close the main while loop */
  }            /* Close the METHOD=POST */

  CheckErrors(Genes,nGenes);     /* Check errors in Genes */

  /* Evaluation of Genes */
  printf ("Locus        LongSeq  CDSRea  CDSPred    TP       TN     ");
  printf ("FP     FN   ExR  ExP  TPe  Sn    Sp    AC    CC    ");
  printf ("Sne   Spe  SnSp    ME    WE\n\n");

  for (i=0;i<nGenes;i++) {
    if (Genes[i].LocusName[0]) {
      k++;
      CDSr=CDSp=0;
      Tpe=Fpe=Fne=OverR=OverP=RealEx=PredEx=0;
      Tp=CalculateTp(Genes[i].Reals,Genes[i].Predicts);
      for (j=0;Genes[i].Reals[j]!=-1;j+=2)
        CDSr += Genes[i].Reals[j+1] - Genes[i].Reals[j] + 1;
      for (j=0;Genes[i].Predicts[j]!=-1;j+=2)
        CDSp += Genes[i].Predicts[j+1] - Genes[i].Predicts[j] + 1;
      Fn = CDSr-Tp;
      Fp = CDSp-Tp;
      Tn = Genes[i].LongSeq-Tp-Fp-Fn;
      CalculateExons(Genes[i].Reals,Genes[i].Predicts,
                     &Tpe,&Fpe,&Fne,&OverR,&OverP);
      RealEx = Tpe+Fne+OverR;
      PredEx = Tpe+Fpe+OverP;
      j=4;
      Angle=0;
      if ((Tp+Fn)==0) {S1=-1;j--;}
      else {S1=(double)Tp/(Tp+Fn);Angle+=S1;sS1+=S1;cS1++;}
      if ((Tp+Fp)==0) {S2=-1;j--;}
      else {S2=(double)Tp/(Tp+Fp);Angle+=S2;sS2+=S2;cS2++;}
      if ((Tn+Fn)==0) {S3=-1;j--;}
      else {S3=(double)Tn/(Tn+Fn);Angle+=S3;sS3+=S3;cS3++;}
      if ((Tn+Fp)==0) {S4=-1;j--;}
      else {S4=(double)Tn/(Tn+Fp);Angle+=S4;sS4+=S4;cS4++;}
      Angle=((Angle/j)-0.5)*2;
      if (j!=4)
        Cc=-2;
      else {
	Cc=(((double)Tp*Tn)-((double)Fp*Fn))/sqrt(((double)(Tn+Fn)*(Tn+Fp)*(Tp+Fn)*(Tp+Fp)));
        sCc+=Cc;
        cCc++;
      }
      j=2;
      if (RealEx==0) {
        Sne=-1;
        j--;
      }
      else {
        Sne=((double)Tpe/RealEx);
        sSne+=Sne;
        cSne++;
      }
      if (PredEx==0) {
        Spe=-1;
        j--;
      }
      else {
        Spe=((double)Tpe/PredEx);
        sSpe+=Spe;
        cSpe++;
      }
      if (j!=2)
        SnSp=-1;
      else {
        SnSp=(Sne+Spe)/2;
        sSnSp+=SnSp;
        cSnSp++;
      }
      if (RealEx!=0)
	ME=(double)Fne/RealEx;
      else
	ME=0;
      if (PredEx!=0)
        WE=(double)Fpe/PredEx;
      else
        WE=0;        
      sLongSeq+=Genes[i].LongSeq;
      sCDSr += CDSr;
      sCDSp += CDSp;
      sTp += Tp;
      sTn += Tn;
      sFp += Fp;
      sFn += Fn;
      sRealEx += RealEx;
      sPredEx += PredEx;
      sTpe += Tpe;
      sFne += Fne;
      sFpe += Fpe;
      sME += ME;
      sWE += WE;
      sAngle += Angle;
      
      printf ("%-10s %8ld ",Genes[i].LocusName,Genes[i].LongSeq);
      printf ("%7ld %7ld %7ld %8ld %6ld %6ld ",CDSr,CDSp,Tp,Tn,Fp,Fn);
      printf ("%4d %4d %4d ",RealEx,PredEx,Tpe);
      if (S1==-1)
        printf ("  In. ");
      else
        printf ("%5.2f ",S1);
      if (S2==-1)
        printf ("  In. ");
      else
        printf ("%5.2f ",S2);
      printf ("%5.2f ",Angle);
      if (Cc==-2)
        printf ("  In. ");
      else
        printf ("%5.2f ",Cc);
      if (Sne==-1)
        printf ("  In. ");
      else
	printf ("%5.2f ",Sne);
      if (Spe==-1)
        printf ("  In. ");
      else
	printf ("%5.2f ",Spe);
      if (SnSp==-1)
        printf ("  In. ");
      else
        printf ("%5.2f ",SnSp);
      printf ("%5.2f %5.2f\n",ME,WE);
    }
  }
  if (k==0) {
    if (!strcmp(Method,"POST"))
      printf ("<center><h2>\n");
    printf ("you MUST have any correct data to evaluate!!!\n");
    if (!strcmp(Method,"POST"))
      printf ("</h2></center>\n");
    exit(1);
  }
  printf ("\nAverage    %8.0f %7.0f %7.0f ",sLongSeq/k,sCDSr/k,sCDSp/k);
  printf ("%7.0f %8.0f %6.0f %6.0f ",sTp/k,sTn/k,sFp/k,sFn/k);
  printf ("%4.1f %4.1f %4.1f ",sRealEx/k,sPredEx/k,sTpe/k);
  printf ("%5.2f %5.2f %5.2f %5.2f ",sS1/cS1,sS2/cS2,sAngle/k,sCc/cCc);
  printf ("%5.2f %5.2f %5.2f ",sSne/cSne,sSpe/cSpe,sSnSp/cSnSp);
  printf ("%5.2f %5.2f\n",sME/k,sWE/k);

  T=sLongSeq;                             /* Calculate TOTAL */
  RC=sCDSr;
  PC=sCDSp;
  TP=sTp;
  FP=PC-TP;
  FN=RC-TP;
  RN=T-RC;
  PN=T-PC;
  TN=RN-FP;
  CC=((TP*TN)-(FN*FP))/sqrt(RC*RN*PC*PN);
  SEN1=TP/RC;                                    /* RC=TP+FN */
  SEN2=TP/PC;                                    /* PC=TP+FP */
  ANGLE=((TP/RC+TP/PC+TN/RN+TN/PN)/4-0.5)*2;
  SNE=sTpe/sRealEx;
  SPE=sTpe/sPredEx;
  printf ("\nTOTAL      %8.0f %7.0f %7.0f ",sLongSeq,sCDSr,sCDSp);
  printf ("%7.0f %8.0f %6.0f %6.0f ",sTp,sTn,sFp,sFn);
  printf ("%4.0f %4.0f %4.0f ",sRealEx,sPredEx,sTpe);
  printf ("%5.2f %5.2f %5.2f %5.2f ",SEN1,SEN2,ANGLE,CC);
  printf ("%5.2f %5.2f %5.2f %5.2f %5.2f\n",SNE,SPE,(SNE+SPE)/2,sFne/sRealEx,sFpe/sPredEx);

  
/* Tail of HTML (if necessari) */
  if (!strcmp(Method,"POST")) {
    printf ("</pre>\n");

    /* form to add results to file */
    printf ("<form action=\"/~mburset/cgi-bin/add.cgi\" method=\"post\" enctype=\"multipart/form-data\">\n");
    printf ("<p><br><hr><p><br>\n");
    printf ("<center>\n");
    printf ("<h2>If you want to add your results in the table ");
    printf ("please fill the form</h2><br><p><br>\n");
    printf ("</center>\n");
    printf ("<h3><nobr>YOUR NAME:\n");
    printf ("<input type=\"text\" name=\"name\"></nobr><p>\n");
    printf ("<nobr>YOUR E-MAIL:\n");
    printf ("<input type=\"text\" name=\"email\"></nobr><p>\n");
    printf ("THE NAME OF PROGRAM USED:<br>\n");
    printf ("<nobr>(Please, use only one word) \n");
    printf ("<input type=\"text\" name=\"program\"></nobr><p>\n");
    printf ("COMMENTS:<br>\n");
    printf ("<textarea name=\"comment\" rows=6 cols=60>\n");
    printf ("</textarea><br>\n");
    printf ("<input type=\"hidden\" name=\"results\" value=\"");
    printf ("%ld %.2f %.2f %.2f %.2f %.2f ",k,sS1/cS1,sS2/cS2,sAngle/k,sCc/cCc,sSne/cSne);
    printf ("%.2f %.2f %.2f %.2f \">\n",sSpe/cSpe,sSnSp/cSnSp,sME/k,sWE/k);
    printf ("<center><p><br><p><br>\n");
    printf ("<input type=\"submit\" value=\"send data\">\n");
    printf ("<input type=\"reset\" value=\"clear data\">\n");
    printf ("</center>\n");
    printf ("</form>\n");

    printf ("</BODY></HTML>\n");
  }
  exit(0);
}

void error(char *s, ...) {
  char format[]="evaluation <LengthSeq+CDSrea> < <CDSpred>";
  va_list args;

  va_start(args,s);
  fprintf(stderr, "error: ");
  vfprintf(stderr, s, args);
  fprintf(stderr,"\nTo run this program try with the format:\n");
  fprintf (stderr,"\t%s\n",format);
  va_end(args);
  exit(1);
}

void warning(char *s, ...) {
  va_list args;

  va_start(args,s);
  fprintf(stderr, "warning: ");
  vfprintf(stderr, s, args);
  fprintf(stderr,"\n");
  va_end(args);
}

char ReadWord (char *Word, FILE *F, long Limit) {
  char Status=' ';
  int a, Cont=0, Bit=0;
  static long Position=0;

  if (!strcmp(Word,"Initialize")) {
    Position=0;
    return(' ');
  }
  while (1) {
    a=fgetc(F);
    Position++;
    if (Position>=Limit)
      return('f');
    else if (a==9) {
      if (Cont!=0)
	Bit=1;
    }
    else if (a==10) {
      Status='n';
      if (Cont!=0)
	Bit=1;
    }
    else if (a==13) ;
    else if (a==32) {
      if (Cont!=0)
	Bit=1;
    }
    else if (a<0 || a>256) {
      *(Word+Cont)='\0';
      return ('f');
    }
    else if (a<32 || (a>126 && a<=256))
      error("Strange ASCII code %u present in your files\n",a);
    else {
      if (Bit==1) {
	Position--;
        ungetc(a,F);
	*(Word+Cont)='\0';
	return(Status);
      }
      else {
	*(Word+Cont)=a;
	Cont++;
      }
    }
  }
}

void OrderPositions (long *vector) {
  int i,j;

  for (i=0;vector[i]!=-1;i+=2)
    for (j=i+2;vector[j]!=-1;j+=2)
      if (vector[i]>vector[j]) {
	vector[i]=vector[i]+vector[j];
	vector[j]=vector[i]-vector[j];
	vector[i]=vector[i]-vector[j];
	vector[i+1]=vector[i+1]+vector[j+1];
	vector[j+1]=vector[i+1]-vector[j+1];
	vector[i+1]=vector[i+1]-vector[j+1];
      }
}

char LoadGenes(char *name, long *gene, long Limit, FILE *F) {
  char Status=' ',TmpVector[MAXWORD];
  long i=0;

  while (Status!='f' && Status!='n') {
    Status=ReadWord(TmpVector,F,Limit);
    *(gene+i)=atol(TmpVector);
    if (*(gene+i)==0)
      error("Position error in %s position in %s gene\n",TmpVector,name);
    if (i>0)
      if (*(gene+i-1)>*(gene+i))
        warning("Non increasing position values in %s gene",name);
    i++;
  }
  if (i%2 !=0)
    error("Odd positions for %s gene\n",name);
  OrderPositions(gene);
  return(Status);
}

void CheckErrors(gene *Genes,long nGenes) {
  long i,j,MaxReal,MaxPredict;

  for (i=0;i<nGenes;i++) {
    if (Genes[i].Reals[0] == -1) {
      warning ("You haven't got the reality for %s gene, ignored",Genes[i].LocusName);
      Genes[i].LocusName[0]='\0';
    }
    else if (Genes[i].Predicts[0]==-1)
      Genes[i].LocusName[0]='\0';
    else if (Genes[i].Reals[0]<0 || Genes[i].Predicts[0]<0)
      error ("Negative positions in %s gene\n",Genes[i].LocusName);
    else {
      j=0;
      while (Genes[i].Reals[j] != -1)
        j++;    
      MaxReal=Genes[i].Reals[--j];
      j=0;
      while (Genes[i].Predicts[j] != -1)
        j++;
      MaxPredict=Genes[i].Predicts[--j];
      if (MaxReal>Genes[i].LongSeq || MaxPredict>Genes[i].LongSeq)
        error ("Positions bigger than sequence length in %s gene\n",Genes[i].LocusName);
    }
  }
}
    
long CalculateTp(long *Reals,long *Predicts) { 
  register int cR,cP;
  int long Tp=0,Max=0;

  for (cR=0;Reals[cR]!=-1;cR+=2) {
    cP=0;
    while (Predicts[cP]<Reals[cR] && Predicts[cP]!=-1)
      cP++;
    while (Predicts[cP]<=Max && Predicts[cP]!=-1)
      cP++;
    if (Predicts[cP] == -1)
      break;
    if (cP%2==0) {
      if (Predicts[cP]<=Reals[cR+1]) {
        if (Predicts[cP+1]>Reals[cR+1])
          Tp+=Reals[cR+1]-Predicts[cP]+1;
        else
          Tp+=Predicts[cP+1]-Predicts[cP]+1;
      }
      if (Predicts[cP+2]<=Reals[cR+1] && Predicts[cP+2]!=-1) {
        cR-=2;
        Max=Predicts[cP+1];
      }
    }
    else {
      if (Predicts[cP]<Reals[cR+1])
        Tp+=Predicts[cP]-Reals[cR]+1;
      else
        Tp+=Reals[cR+1]-Reals[cR]+1;
      if (Predicts[cP+1]<=Reals[cR+1] && Predicts[cP+1]!=-1) {
        cR-=2;
        Max=Predicts[cP];
      }
    }
  }
return(Tp);
}

void CalculateExons(long *Reals, long *Predicts, int *TPe, int *FPe,
                    int *FNe, int *OverR, int *OverP) {
  register int cR=0,cP=0;

  /* Check if reals or predicts have no prediction */
  if (Reals[0] == 0) {
    while(Predicts[cP*2+1] != -1)
      cP++;
    (*FPe)=cP;
    return;
  }
  else if (Predicts[0] == 0) {
    while(Reals[cP*2+1] != -1)
      cP++;
    (*FNe)=cP;
    return;
  }

  /* If not, calculate measures */
  for(cR=0;Reals[cR]!=-1;cR+=2) {
    cP=0;
    while(Predicts[cP]<Reals[cR] && Predicts[cP]!=-1)
      cP++;
    if (Predicts[cP] == -1) {
      (*FNe)++;
      continue;
    }
    if (cP%2==0) {
      if (Predicts[cP]==Reals[cR])
        if (Predicts[cP+1]==Reals[cR+1])
          (*TPe)++;
        else
          (*OverR)++;
      else
        if (Predicts[cP]>Reals[cR+1])
          (*FNe)++;
        else (*OverR)++;
    }
    else
      (*OverR)++; 
  }
  for(cR=0;Predicts[cR]!=-1;cR+=2) {
    cP=0;
    while(Reals[cP]<Predicts[cR] && Reals[cP]!=-1)
      cP++;
    if (Reals[cP] == -1) {
      (*FPe)++;
      continue;
    }
    if (cP%2==0) {
      if (Reals[cP]==Predicts[cR])
        if (Reals[cP+1]==Predicts[cR+1])
          (*TPe)++;
        else
          (*OverP)++;
      else
        if (Reals[cP]>Predicts[cR+1])
          (*FPe)++;
        else 
          (*OverP)++;
    }
    else
      (*OverP)++;
  }
(*TPe)/=2;
}
