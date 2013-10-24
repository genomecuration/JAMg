
/***********************************************************************************/
/*                                                                                 */
/* pictogram.c                                     by Chris Burge : February, 2000 */
/*                                                                                 */
/* Creates a compositional profile from a library of precisely aligned sequences.  */
/*                                                                                 */
/***********************************************************************************/

#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "info.h"

#include "ps.h"

#define LOG2(x) (log(x)/0.693147181)

#define MAX_SEQ_LEN	1000

int SMALL_SAMPLE = 0;		/* if 1, use small sample correction	*/

#define PSEUDOCOUNTS	0

#define SS_SD		2.0	/* number of standard deviations above mean	*/
				/* for significance of positional information	*/

int PRINT_MONO	= 0;	
int PRINT_BITS	= 0;
int PRINT_RF	= 0;
int PRINT_YF	= 0;
int PRINT_COUNT = 1;
int PRINT_DI	= 0;
int PRINT_TRI	= 0;
int PRINT_TETRA	= 0;
int PRINT_BIT_PLOT = 0;

int PRINT_PS	= 0;

int PRINT_MAX_Y_RUN = 0;

int OFFSET = 0;

char dna_alphabet[] = "ACGT";
char rna_alphabet[] = "ACGU";

char *alphabet = dna_alphabet;

char a_color[] = "0 1 0.5 setrgbcolor";
char c_color[] = "1 0 0 setrgbcolor";
char g_color[] = "0 0 1 setrgbcolor";
char t_color[] = "1 0.7 0 setrgbcolor";
char black[] = "0 0 0 setrgbcolor";

struct bf
  {
  float f;
  int b;
  };

bfcomp(pi,pj)
struct bf *pi,*pj;
{
if (pi->f > pj->f) return(1);
if (pi->f < pj->f) return(-1);
if (pi->f == pj->f) return(0);
}


main(argc,argv)
int argc;
char *argv[];

{
double freq,r;
int s;
char infname[100],outfname[100],psfname[100];

/* char C,*cseq[MAX_NUM_SEQ],line[MAX_SEQ_LEN*2],test[10]; */

char C,cseq[MAX_SEQ_LEN];
char line[MAX_SEQ_LEN*2],test[10];

int landscape=0;

int profc[MAX_SEQ_LEN][4];
double EI,OI,WF;		/* expected information, observed information	*/
double profI[MAX_SEQ_LEN];
int profsum[MAX_SEQ_LEN];
int prof2sum[MAX_SEQ_LEN];
int prof3sum[MAX_SEQ_LEN];
int prof4sum[MAX_SEQ_LEN];
int prof2c[MAX_SEQ_LEN][4][4];
int prof3c[MAX_SEQ_LEN][4][4][4];
int prof4c[MAX_SEQ_LEN][4][4][4][4];
double proff[MAX_SEQ_LEN][4],f;
double N,D,H,I,Isum;

int seq[MAX_SEQ_LEN];

int snum=0,len,nhits,linelen;
int maxYrun,Yrun;
int maxY1Erun,Y1Erun,nE;
int Ybegin,Y1begin;
float F,cumF,SF;
int i,j,k,h,a,b,c,d,n,lnum;

double bgfreq[4],tf[4];

struct bf *tmpf;

FILE *infp,*outfp,*psfp;

if (argc<3) usage(argv[0]);

strcpy(infname,argv[1]);
strcpy(outfname,argv[2]);

infp = fopen(infname,"r");

outfp = fopen(outfname,"w");

for (i=0;i<4;++i) bgfreq[i] = 0.25;

for (i=3;i<argc;++i)
  {
  if ( strncmp(argv[i],"-bg",3)==0 && i+4 < argc)
    {
    bgfreq[0]=atof(argv[++i])/100.0;
    bgfreq[1]=atof(argv[++i])/100.0;
    bgfreq[2]=atof(argv[++i])/100.0;
    bgfreq[3]=atof(argv[++i])/100.0;
    } 
  if (strncmp(argv[i],"-bits",5)==0)    { PRINT_BITS  = 1; 		}
  if (strncmp(argv[i],"-ss",3)==0)      { SMALL_SAMPLE = 1; 		}
  if (1)				{ PRINT_COUNT = 1; 		}
  if (1)				{ PRINT_MONO = 1; 		}
  if (1)				{ PRINT_PS    = 1; 		}
  if (strncmp(argv[i],"-land",5)==0)    { landscape = 1;		}
  if (strncmp(argv[i],"-offset",7)==0)  { OFFSET = atoi(argv[++i]);	 }
  if (strncmp(argv[i],"-bitplot",8)==0) { PRINT_BIT_PLOT = 1;		}
  }


while ( fgets(line,150,infp) != NULL) 
  {
  if (line[0] == '#') continue;
  if (line[0] == '>') continue;
  if (line[0] == '/') continue;
  if (line[0] == ' ') continue;
  if (line[0] == '\n') continue;
  if (!(isalpha(line[0]))) continue;
  strcpy(cseq,line);
  ++snum;
  if (snum==1)
    {
    len=0;
    linelen = strlen(cseq); 
    for (j=0;j<linelen;++j) { if ( isalpha(cseq[j]) ) ++len; }
    for (j=0;j<len;++j)
      {
      profsum[j]=0;
      prof2sum[j]=0;
      prof3sum[j]=0;
      prof4sum[j]=0;
    
      for (a=0;a<4;++a) { profc[j][a]=0;
        for (b=0;b<4;++b) { prof2c[j][a][b]=0;
          for (c=0;c<4;++c) { prof3c[j][a][b][c]=0;
            for (d=0;d<4;++d) { prof4c[j][a][b][c][d]=0;
              }
            }
          }
        }
    
      }
    }

  k=0;
  for (j=0;j<linelen;++j) {
    a = 4;
    c = cseq[j];
    if (isalpha(c))
      { 
      C = toupper(c);
      if (C == 'A') a = 0;
      if (C == 'C') a = 1;
      if (C == 'G') a = 2;
      if (C == 'T' || C == 'U') a = 3;
      seq[k] = a;
      ++k;
      }
    }

  maxYrun=0; Yrun=0; Ybegin=1;
  maxY1Erun=0; Y1Erun=0; nE=0; Y1begin=1;
  for (j=0;j<len;++j) {
    a = seq[j];
    if (a < 4) { profc[j][a]+=1; profsum[j]+=1; }
    if (a==1 || a==3) {
      Yrun++; Y1Erun++;
      if (Yrun>maxYrun) { maxYrun=Yrun; Ybegin=j+2-Yrun; }
      if (Y1Erun>maxY1Erun) { maxY1Erun=Y1Erun; Y1begin=j+1-Y1Erun; }
      }
    else { Y1Erun=Yrun; Yrun=0; }
     
    if (j<len-1) {
      b = seq[j+1];
      if (a<4 && b<4) {
        prof2c[j][a][b] += 1;
        prof2sum[j] += 1;
        }
      }
    if (j<len-2) {
      c = seq[j+2];
      if (a<4 && b<4 && c<4) {
        prof3c[j][a][b][c] += 1;
        prof3sum[j]+=1;
        }
      }
    if (j<len-3) {
      d = seq[j+3];
      if (a<4 && b<4 && c<4 && d<4) {
        prof4c[j][a][b][c][d] += 1;
        prof4sum[j]+=1;
        }      
      }
    }
  if (PRINT_MAX_Y_RUN) {
    fprintf(stdout,"%d\t%d\t%d\t%d\n",maxYrun,maxY1Erun,Ybegin,Y1begin);
    }

  } 

fprintf(stderr,"\n%d sequences of length %d read...",snum,len);


if (PRINT_COUNT)
  {

  for (j=0;j<len;++j) {
    for (i=0;i<4;++i) { 
      fprintf(outfp,"%d",profc[j][i]);
      if (i<3) fprintf(outfp,"\t");
      else fprintf(outfp,"\n");
      }
    }

  }

if (PRINT_MONO || PRINT_RF || PRINT_YF || PRINT_BITS)
  {
  for (lnum=0;lnum<=(len-1)/8;++lnum)
    {

    fprintf(outfp,"\nPos");
    for (j=lnum*8;j<len && j-(lnum*8)<8;++j)
      {
      fprintf(outfp,"\t%4d",j+OFFSET);
      }
    fprintf(outfp,"\n\n");

    if (PRINT_MONO) 
      {

#ifdef HORIZONTAL
      for (i=0;i<4;++i)
        {
        fprintf(outfp,"%c",alphabet[i]);
        for (j=lnum*8;j<len && j-(lnum*8)<8;++j)
          {
          proff[j][i] = (double)profc[j][i]/(double)profsum[j];
          if (PSEUDOCOUNTS)
            proff[j][i] = ((double)profc[j][i]+0.25)/((double)profsum[j]+1.0);
         
#ifdef INTEGER
          fprintf(outfp,"\t%4d",(int)(100.0*(proff[j][i]+0.005)));
#endif
#ifndef INTEGER
          fprintf(outfp,"\t%5.3f",proff[j][i]);
#endif
          }
        fprintf(outfp,"\n");
        }
#endif

#ifndef HORIZONTAL
      for (j=lnum*8;j<len && j-(lnum*8)<8;++j)
        {
        fprintf(outfp,"%c",alphabet[i]);
        for (i=0;i<4;++i)
          {
          proff[j][i] = (double)profc[j][i]/(double)profsum[j];
          if (PSEUDOCOUNTS)
            proff[j][i] = ((double)profc[j][i]+0.25)/((double)profsum[j]+1.0);
#ifdef INTEGER
          fprintf(outfp,"\t%4d",(int)(100.0*(proff[j][i]+0.005)));
#endif
#ifndef INTEGER
          fprintf(outfp,"\t%5.3f",proff[j][i]);
#endif
          }
        fprintf(outfp,"\n");
        }
#endif

      }

    for (j=lnum*8;j<len && j-(lnum*8)<8;++j)
      {
      for (k=0;k<4;++k)
        {
        freq = ((double)profc[j][k]+0.25)/((double)profsum[j]+1.0);
        r = freq / bgfreq[k];
        s = (r > 0.0) ? (int)(10.0*LOG2(r)) : -100;
        fprintf(stdout,"%d",s);
        if (k<3) fprintf(stdout,"\t");
        else fprintf(stdout,"\n");
        }
      }


    if (PRINT_RF)
      {
      for (j=lnum*8;j<len && j-(lnum*8)<8;++j)
        {
        fprintf(outfp,"\t%4d",(int)(100.0*(proff[j][0]+proff[j][2]+0.005))); 
        }
      fprintf(outfp,"\n");
      }

    if (PRINT_YF)
      {
      for (j=lnum*8;j<len && j-(lnum*8)<8;++j)
        {
        fprintf(outfp,"\t%4d",(int)(100.0*(proff[j][1]+proff[j][3]+0.005))); 
        }
      fprintf(outfp,"\n");
      }

    if (PRINT_BITS)
      {
      EH[0] = 2.0 - ( 3.0 / (2.0 * log(2.0) * snum) );
      fprintf(outfp,"\nH:");
      for (j=lnum*8;j<len && j-(lnum*8)<8;++j)
        {

        /* H = entropy(proff[j],4); */

        I = rel_entropy(proff[j],bgfreq,4);

        H = 2.0 - I; 
				/* subtract from expected sample entropy */


        if (SMALL_SAMPLE) {	
          if (snum <= 25) EI = (double) ( 2.0 - EH[snum] );
          else            EI = (double) ( 2.0 - EH[0] );
          OI = 2.0 - H;
          I = ( OI - (double) SS_SD * EI > 0.0 ) ? OI - EI : 0.0;
          }

        profI[j] = I;

        fprintf(outfp,"\t%4.2f",I);
        }

      fprintf(outfp,"\n");
      }

    }

  }

tf[0]=0.97;
tf[1]=0.01;
tf[2]=0.01;
tf[3]=0.01;

I = rel_entropy(tf,bgfreq,4);

if (PRINT_PS == 1)
  {
  strcpy(psfname,outfname);
  strcat(psfname,".ps");
 
  psfp = fopen(psfname,"w");

  for (j=0;j<sizeof(PS_HEAD)/100;++j) fprintf(psfp,"%s",PS_HEAD[j]); 

  fprintf(psfp,"/Times-Roman findfont 40 scalefont setfont\n");

  fflush(psfp);

	/* Scale of usable image: 2000 vertical x 1500 horizontal	*/

  if (landscape)
    {
    fprintf(psfp,"90 rotate\n");
    fprintf(psfp,"300 -2700 T\n");
    }

  fprintf(psfp,"/dY { 400 } def\n");	/* height of 1 bit	*/

  if (landscape) {
    fprintf(psfp,"/dX { 120 15 mul %d div } def\n",len);
    }

  else
    {

    if (len == 15)
      fprintf(psfp,"/dX { 100 } def\n");	/* width of 1 base	*/

    else
      fprintf(psfp,"/dX { 100 15 mul %d div } def\n",len);

    }

  fprintf(psfp,"200 200 T\n");

  fprintf(psfp,"0 1500 T\n");

  fprintf(psfp,"/Times-Bold findfont 60 scalefont setfont\n");

  fprintf(psfp,"0 dY 200 add  M\n");
  fprintf(psfp,"(Compositional profile of %s) show\n",infname);

  fprintf(psfp,"/Times-Roman findfont 40 scalefont setfont\n");

  fprintf(psfp,"0 0 M\n");

  fprintf(psfp,"0 dY RL\n");

  if (len <= 20) fprintf(psfp,"dX %d mul 5 %d 15 sub mul add 0 RL\n",len,len);

  if (len > 20 && len <= 30) 
    fprintf(psfp,"dX %d 0.5 sub mul 5 %d 0.5 sub 15 sub mul add 0 RL\n",len,len);

  if (len > 30 && len <= 40)
    fprintf(psfp,"dX %d 1 sub mul 5 %d 1 sub 15 sub mul add 0 RL\n",len,len);

  if (len > 40)
    fprintf(psfp,"dX %d 1.5 sub mul 5 %d 1.5 sub 15 sub mul add 0 RL\n",len,len);

  fprintf(psfp,"0 dY -1 mul RL\n");
  fprintf(psfp,"closepath S\n");

  fprintf(psfp,"-125 10 T\n");

/*
  fprintf(psfp,"0 0 M\n");
  fprintf(psfp,"(0.00) Lshow\n");
  fprintf(psfp,"0 dY 0.25 mul M\n");
  fprintf(psfp,"(0.25) Lshow\n");
  fprintf(psfp,"0 dY 0.5 mul M\n");
  fprintf(psfp,"(0.50) Lshow\n");
  fprintf(psfp,"0 dY 0.75 mul M\n");
  fprintf(psfp,"(0.75) Lshow\n");
  fprintf(psfp,"0 dY 1.00 mul M\n");
  fprintf(psfp,"(1.00) Lshow\n");
*/

  fprintf(psfp,"125 -10 T\n");

/*
  fprintf(psfp,"0 dY 0.00 mul M\n-25 0 RL\n");
  fprintf(psfp,"0 dY 0.25 mul M\n-25 0 RL\n");
  fprintf(psfp,"0 dY 0.50 mul M\n-25 0 RL\n");
  fprintf(psfp,"0 dY 0.75 mul M\n-25 0 RL\n");
  fprintf(psfp,"0 dY 1.00 mul M\n-25 0 RL\n");
*/

  fprintf(psfp,"0 1 T\n");

  fprintf(psfp,"/Helvetica-Bold findfont 140 scalefont setfont\n");

  for (j=0;j<len;++j)
    {
    tmpf = (struct bf *) malloc(4*sizeof(struct bf));
    for (k=0;k<4;++k) { tmpf[k].b = k; tmpf[k].f = (float)proff[j][k]; }
    qsort(tmpf,4,sizeof(struct bf),bfcomp);

    cumF = 0.0;

    for (k=0;k<4;++k)
      {
      fprintf(psfp,"%d dX mul %f dY mul 20 add M\n",j,0.99*cumF);
      F = (float)tmpf[k].f;
      SF = (float)(F * 4.0 * 0.99);
      cumF += F;
      C = alphabet[tmpf[k].b];
      if (F > 0.01)
        {
        fprintf(psfp,"0 -20 RM\n"); 
        fprintf(psfp,"0 1 RM\n");
        if (C == 'T') fprintf(psfp,"9 1 RM\n");
        if (C == 'G') fprintf(psfp,"-2 0 RM\n");
        fprintf(psfp,"1 %1.6f scale\n",SF);
        WF = 1.0 + 0.03 * (double) ( 15 - len );
        if (C == 'A') 
          {
          fprintf(psfp,"%s\n",a_color);
          fprintf(psfp,"%4.2f 0.97 scale\n",WF);
          fprintf(psfp,"(%c) show\n",C);
          fprintf(psfp,"1 %4.2f div 1 0.97 div scale\n",WF);
          fprintf(psfp,"%s\n",black);
          }
        if (C == 'T') 
          {
          fprintf(psfp,"%s\n",t_color);
          fprintf(psfp,"%4.2f 0.97 scale\n",WF);
          fprintf(psfp,"(%c) show\n",C);
          fprintf(psfp,"1 %4.2f div 1 0.97 div scale\n",WF);
          fprintf(psfp,"%s\n",black);
          }

        if (C == 'U')
          {
          fprintf(psfp,"%s\n",t_color);
          fprintf(psfp,"%4.2f 0.98 scale\n",WF);
          fprintf(psfp,"0 2 RM\n");
          fprintf(psfp,"(%c) show\n",C);
          fprintf(psfp,"1 %4.2f div 1 0.98 div scale\n",WF);
          fprintf(psfp,"%s\n",black);
          }
        if (C == 'C' || C == 'G')
          {
          if (C == 'C') fprintf(psfp,"%s\n",c_color);
          if (C == 'G') fprintf(psfp,"%s\n",g_color);
          fprintf(psfp,"%4.2f 0.93 scale\n",WF);
          fprintf(psfp,"0 3 RM\n");
          fprintf(psfp,"(%c) show\n",C);
          fprintf(psfp,"1 %4.2f div 1 0.93 div scale\n",WF);
          fprintf(psfp,"%s\n",black);
          }
        fprintf(psfp,"1 %f scale\n",1.0/SF);
        }
      }

    }

  fprintf(psfp,"0 -1 T\n");

  fprintf(psfp,"0 dY 50 add T\n");

  fprintf(psfp,"/Times-Roman findfont 40 scalefont setfont\n");

  fprintf(psfp,"-125 0 M\n");
  fprintf(psfp,"(Pos:) show\n");
  
  fprintf(psfp,"25 0 T\n");

  if (landscape) fprintf(psfp,"-15 0 T\n");

  for (j=0;j<len;++j)
    {
    fprintf(psfp,"dX %d mul 0 M\n",j);
    k = (j-OFFSET<0) ? j-OFFSET : j+1-OFFSET;
    if (k>0 && k < 10) fprintf(psfp,"15 0 RM\n");
    if (k >= 10) fprintf(psfp,"5 0 RM\n");
    fprintf(psfp,"%2d prt-n\n",k);
    }

  if (landscape) fprintf(psfp,"20 0 T\n");

  fprintf(psfp,"-25 0 T\n");

  fprintf(psfp,"0 dY 50 add -1 mul T\n");

  if (PRINT_BIT_PLOT)
    {
    fprintf(psfp,"0 -225 T\n");

    fprintf(psfp,"0 0 M\n"); 
    fprintf(psfp,"0 225 RL\n");
    fprintf(psfp,"dX %d mul 0 RL\n",len);
    fprintf(psfp,"0 -225 RL\n");
    fprintf(psfp,"closepath S\n");
    
    fprintf(psfp,"0 0   M\n-20 0 RL\n");
    fprintf(psfp,"0 50  M\n-20 0 RL\n");
    fprintf(psfp,"0 100 M\n-20 0 RL\n");
    fprintf(psfp,"0 150 M\n-20 0 RL\n");
    fprintf(psfp,"0 200 M\n-20 0 RL\n");

    fprintf(psfp,"0 0 M\n");
    fprintf(psfp,"dX 2 div 0 T\n");
    for (j=0;j<len;++j) 
      {
      /* fprintf(psfp,"dX %d mul 100 %5.3f mul L\n",j,profI[j]); */
      fprintf(psfp,"dX %d mul 0 M\n",j);
      fprintf(psfp,"0 100 %5.3f mul RL\n",profI[j]);
      }
    fprintf(psfp,"dX 2 div -1 mul 0 T\n");   
    /* fprintf(psfp,"dX %d mul 0 L\n",len); */

    fprintf(psfp,"0 225 T\n");
    } 

  if (PRINT_BITS)
    {
    if (!PRINT_BIT_PLOT) fprintf(psfp,"0 75 T\n");
    else fprintf(psfp,"0 -150 T\n");

    fprintf(psfp,"0 -75 T\n");

    fprintf(psfp,"dX %d mul 0 M\n",len);
    fprintf(psfp,"50 -25 RM\n");
    fprintf(psfp,"/Times-Bold findfont 40 scalefont setfont\n");
    fprintf(psfp,"(Total) show\n");

    fprintf(psfp,"0 75 T\n");

    fprintf(psfp,"0 -150 T\n");

    if (len <= 15) fprintf(psfp,"/Times-Bold findfont 50 scalefont setfont\n");

    else
      {
      if (len <= 20) fprintf(psfp,"/Times-Bold findfont 46 scalefont setfont\n");
      else 
        {
        if (len <= 25) fprintf(psfp,"/Times-Bold findfont 42 scalefont setfont\n");
        else 
          {
          if (len <= 30) fprintf(psfp,"/Times-Bold findfont 38 scalefont setfont\n");
          else fprintf(psfp,"/Times-Bold findfont 34 scalefont setfont\n");
          }
        }
      }
  
    if (landscape) fprintf(psfp,"-125 20 M\n");

    else fprintf(psfp,"-125 -25 M\n");

    fprintf(psfp,"(Bits:) show\n");

    Isum = 0.0;

    fprintf(psfp,"10 0 T\n");

    if (landscape) fprintf(psfp,"-20 0 T\n");

    for (j=0;j<len;++j)
      {
      fprintf(psfp,"dX %d mul 0 M\n",j);
      I = profI[j];
      Isum += I;
      /* if (j%2==1) fprintf(psfp,"0 -50 RM\n"); */
      fprintf(psfp,"10 20 RM\n");
      fprintf(psfp,"%3.1f prt-n\n",I); 
      }  

    if (landscape) fprintf(psfp,"20 0 T\n");

    fprintf(psfp,"dX %d mul 50 add 0 M\n",len);
    fprintf(psfp,"0 -30 RM\n");
    fprintf(psfp,"/Times-Bold findfont 40 scalefont setfont\n");
    fprintf(psfp,"%3.1f prt-n\nS\n",Isum); 

	/* Draw box around bit total	*/

/*
    fprintf(psfp,"S\nlwidth 2 div setlinewidth\n");
    fprintf(psfp,"dX %d mul 25 add 30 M\n",len);
    fprintf(psfp,"200 0 RL\n");
    fprintf(psfp,"0 -80 RL\n");
    fprintf(psfp,"-200 0 RL\n");
    fprintf(psfp,"closepath\n");
*/

    fprintf(psfp,"-10 0 T\n");

    fprintf(psfp,"0 150 T\n");

    if (!PRINT_BIT_PLOT) fprintf(psfp,"0 -75 T\n");
    else fprintf(psfp,"0 150 T\n");
    }

  for (j=0;j<sizeof(PS_TAIL)/100;++j) fprintf(psfp,"%s",PS_TAIL[j]); 

  fflush(psfp);

  fclose(psfp);
  }

if (PRINT_DI)
  {
  for (j=0;j<len-1;++j)
    {
    for (a=0;a<4;++a) { 
      for (b=0;b<4;++b) {
        N = (double)prof2c[j][a][b]/(double)prof2sum[j];
        D = (double)profc[j][a]/(double)prof2sum[j];
        D *= (double)profc[j+1][b]/(double)prof2sum[j];
        fprintf(outfp,"%6.3f\t",N/D);
        }
      }
    fprintf(outfp,"\n");
    }

  }

if (PRINT_TRI)
  {
  fprintf(outfp,"WORD");
  for (j=0;j<len-3;++j) {
    fprintf(outfp,"\t%d",j+1);
    }
  fprintf(outfp,"\n\n");
  
  for (a=0;a<4;++a) {
    for (b=0;b<4;++b) {
      for (c=0;c<4;++c) {
        fprintf(outfp,"%c%c%c",alphabet[a],alphabet[b],alphabet[c]);
        for (j=0;j<len-3;++j) {
          N = (double)prof3c[j][a][b][c]/(double)prof3sum[j];
          D = (double)prof2c[j][a][b]/(double)prof3sum[j];
          D *= (double)prof2c[j+1][b][c]/(double)prof3sum[j];
          D /= (double)profc[j+1][b]/(double)prof3sum[j];
          f = N/D;
          fprintf(outfp,"\t%6.3f",f);
          }
        fprintf(outfp,"\n");
        } 
      }
    }
  }

if (PRINT_TETRA)
  {
  fprintf(outfp,"WORD");
  for (j=0;j<len-3;++j) {
    fprintf(outfp,"\t%d",j+1);
    }
  fprintf(outfp,"\n\n");
  
  for (a=0;a<4;++a) {
    for (b=0;b<4;++b) {
      for (c=0;c<4;++c) {
        for (d=0;d<4;++d) {
          fprintf(outfp,"%c%c%c%c",alphabet[a],alphabet[b],alphabet[c],alphabet[d]); 
          for (j=0;j<len-3;++j) {
            n = prof4c[j][a][b][c][d];
            /* f = (double)prof4c[j][a][b][c][d]/(double)prof4sum[j]; */
            fprintf(outfp,"\t%d",n);
            }
          fprintf(outfp,"\n");
          }
        } 
      }
    }
  }

fprintf(stderr,"\n\ndone\n\n");

}

/*****************************************************************/

usage(outname)
char *outname;
{
printf("\nusage: %s infile outfile [-bg fA fC fG fT] [other options]\n",outname);
printf("\n       infile   : input file containing one seq per line\n");
printf("\n       outfile  : file for program output\n");
printf("\n       -bg      : background percent frequencies of A, C, G, T (default: 25,25,25,25)\n");
printf("\n       -bits    : print information content at each position\n");
printf("\n       -bitplot : plot information content at each position\n");
printf("\n       -offset c: start coordinates at 1 + c (default: c=0)\n");
printf("\n       -ss      : do small sample size correction for info content (default: don't)\n");
printf("\n       -land    : print PostScript in landscape orientation (default: portrait)\n\n");
exit(0);
}

/*****************************************************************/

