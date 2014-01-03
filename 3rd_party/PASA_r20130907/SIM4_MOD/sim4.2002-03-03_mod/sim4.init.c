/*
* sim4 - Align a cDNA sequence with a genomic sequence for it.
*
* The basic command syntax is
*
*       sim4 seq1 seq2 [[AXWRKCDHEPNBS]=]
*
* where sequence1 and sequence2 name files containing DNA sequences; a name
* like "foo[300,400]" refers to sequence entries 300-400 in file "foo".  The
* files are to be in FASTA format.  Thus a typical sequence file might begin:
*
*       >BOVHBPP3I  Bovine beta-globin psi-3 pseudogene, 5' end.
*       GGAGAATAAAGTTTCTGAGTCTAGACACACTGGATCAGCCAATCACAGATGAAGGGCACT
*       GAGGAACAGGAGTGCATCTTACATTCCCCCAAACCAATGAACTTGTATTATGCCCTGGGC
*
* Typically, sequence1 name file contains one DNA sequence, and sequence2
* name file contains a database of sequences in FASTA format. Sim4 can
* compare a genomic sequence with a database of cDNAs or a cDNA sequence
* with a genomic database. If highly accurate sequences are being compared,
* specify option N=1.
*
* The command permits optional additional arguments, e.g. "A=1" and "W=8",
* that reset certain internal parameters.  The available parameters are:
*
*       W gives the word size.
*       X gives the value for terminating word extensions.
*       K gives the MSP score threshold for the first pass.
*       C gives the MSP score threshold for the second pass. 
*       R direction of search; 0 - search the '+' strand only;
*         1 - search the '-' strand only; 2 - search both strands and
*         report the best match. (R=2)
*       D adjusts the range of diagonals in building the exons.
*       H adjusts the re-linking weight factor
*       A specifies the output format: exon endpoints only (A=0),
*         alignment text (A=1), alignment in lav format (A=2) or both
*         exon endpoints and alignment text (A=3, A=4). For A=3, positions
*         in sequence 1 are given in the original sequence, and for A=4 in
*         its reverse complement. A=5 prints the exon and CDS coordinates 
*         (the latter, if known) in the `exon file' format required by PipMaker.
*       N if !=0, highly accurate exon detection is expected, for highly
*         accurate sequence data. 
*       P remove polyA tails; if match on complementary strand, change
*         coordinates in sequence 1 according to the '+' strand and print
*         #lav alignment header for all alignment options.
*       B control the presence of ambiguity codes in sequence data. If
*         1 (default), allow ambiguity codes (ABCDGKMNRCTVWXY); if 0,  
*         restrict the set of accepted characters to ACGTNX.
*       S specifies a known coding region (to be used only with A=5, and for
*         comparing full-mRNA sequences).
*/      


#include "psublast.h"

#include "sim4.h" 
#include "align.h"
#include "poly.h"
               
#ifndef __lint     
/*@unused@*/       
static const char rcsid[] =
"$Id: sim4.init.c,v 1.1.1.1 2006-07-07 18:15:06 bhaas Exp $";
#endif         
           

static void init_stats(sim4_stats_t *);
static void sim4_argvals(sim4_args_t *);
static void cds_range(char *,int *,int *);
static char *extract_tok(char *);

static void add_offset_exons(Exon *,int);
static void add_offset_aligns(edit_script_list *,int);
static void print_align_blk(uchar *,uchar *,int,int,edit_script_list **,int,int);
static void print_align_lat(uchar *,uchar *,int,int,edit_script_list **,Exon *,int,int);

static const char Usage[] =
"%s seq1 seq2_db [[WXKCRDHAPNBS]=]\n\n\
       W  -  word size. (W=12)\n\
       X  -  value for terminating word extensions. (X=12)\n\
       K  -  MSP score threshold for the first pass. (e.g., K=16)\n\
       C  -  MSP score threshold for the second pass. (e.g., C=12)\n\
       R  -  direction of search; 0 - search the '+' (direct) strand only; \n\
             1 - search the '-' strand only; 2 - search both strands and \n\
             report the best match. (R=2)\n\
       D  -  bound for the range of diagonals within consecutive msps in an\n\
             exon. (D=10)\n\
       H  -  weight factor for MSP scores in relinking. (H=500)\n\
       A  -  output format: exon endpoints only (A=0), alignment text (A=1),\n\
             alignment in lav (block) format (A=2), or both exon endpoints\n\
             and alignment text (A=3, A=4). If complement match, A=0,1,2,3\n\
             give direct positions in the long sequence and complement \n\
             positions in the short sequence. A=4 gives direct positions in \n\
             the first sequence, regardless of the relative lengths.\n\
             A=5 prints the exon and CDS coordinates (the latter, if known)\n\
             in the `exon file' format required by PipMaker. To be used\n\
             with full-length mRNA sequences.\n\
       P  -  if not 0, remove poly-A tails; report coordinates in the \n\
             '+' (direct) strand for complement matches; use lav alignment \n\
             headers in all display options. (P=0) \n\
       N  -  accuracy of sequences (non-zero for highly accurate). (N=0)\n\
       B  -  if 0, dis-allow ambiguity codes (other than N and X) in the\n\
             sequence data. (B=1)\n\
       S  -  coding region specification (available only with A=5);\n\
             format: S=n1..n2 \n";
       
#ifdef _STATS
static void print_stats(sim4_stats_t, char *, int);
static void print_exon_stats(Exon *, char *, int);
#endif

#ifdef DEBUG
static void polystats(int,int,int,int,int,int,int,int);
#endif

       
#include "sim4b1.h"
       
int main(int argc, char *argv[])
{
        uchar  *revseq1=NULL; 
        int    len1, len2, count, dist, match_ori, in_K, in_C, in_H;
        int pA, pT, xpT, xpA, rev_xpT, rev_xpA;
        int cds_from, cds_to;
        uchar *seq1, *seq2;
        char  *h2, *h1, *cds_gene=NULL, *line, *tok1;
        SEQ   *sf1, *sf2, *rf1=NULL;
        argv_scores_t ds;
        ss_t ss;
       
        Exon *Exons=NULL, *rev_Exons=NULL; 
        edit_script_list *Aligns=NULL, *rev_Aligns=NULL;
       
        if (argc < 3) fatalf(Usage, argv[0]);
        ckargs("AXWRKCDHEPNBS", argc, argv, 2);

        sim4_argvals(&rs);
        DNA_scores(&ds, ss);
       
        /* read seq1 */ 
        sf1 = seq_open(argv[1], 0, (rs.B ? SEQ_ALLOW_AMB : 0));
        if (!seq_read(sf1)) fatalf("Cannot read sequence from %s.", argv[1]);
        
        seq1 = SEQ_CHARS(sf1);
        len1 = SEQ_LEN(sf1);
        h1 = SEQ_HEAD(sf1);
        tok1 = extract_tok(h1);
        if (tok1==NULL) {
            tok1 = ckalloc(strlen("(no header)")+1);
            strcpy(tok1, "(no header)");
        }

        if (!is_DNA(seq1, len1))
                fatal("The first sequence is not a DNA sequence.");
        seq_toupper(seq1, len1, argv[1]);
        
        /* read seq2 */
        sf2 = seq_open(argv[2], 0, (rs.B ? SEQ_ALLOW_AMB : 0));
        if (!seq_read(sf2)) fatalf("Cannot read sequence from %s.", argv[2]);
       
        seq2 = SEQ_CHARS(sf2);
        len2 = SEQ_LEN(sf2);
        h2 = SEQ_HEAD(sf2);
        if (!is_DNA(seq2, len2))
                fatal("The first sequence is not a DNA sequence.");
        seq_toupper(seq2, len2, argv[2]);
        
        /* determine the type of comparison */
        file_type = (len2<=len1) ? GEN_EST : EST_GEN;
        if (file_type== EST_GEN) {
            rf1 = seq_copy(sf1);
            rf1 = seq_revcomp_inplace(rf1);
            revseq1 = SEQ_CHARS(rf1);

            if (rs.ali_flag==5) {
                if (rs.CDS_to>len1) 
                   fatal("Command line CDS endpoint exceeds sequence length.");
                cds_gene = extract_tok(h1);
                if (cds_gene==NULL) {  /* no FastaA header */
                    cds_from = rs.CDS_from; cds_to = rs.CDS_to;
                } else {
                    line = strstr(h1, "CDS="); 
                    if (line && rs.S) {
                       fprintf(stderr, "Warning: Command line CDS specification overrides header CDS specification."); 
                       cds_from = rs.CDS_from; cds_to = rs.CDS_to;
                    } else if (line) {
                       cds_range(line+4, &cds_from, &cds_to); 
                    } else if (rs.S) {
                       cds_from = rs.CDS_from; cds_to = rs.CDS_to;
                    } else {
                       cds_from = cds_to = 0;
                    }
                }
                if (cds_to>len1) 
                   fatal("CDS endpoints exceed sequence length.");
            }
        }
        
        if (rs.poly_flag && file_type==EST_GEN)  {
            get_polyAT(seq1,len1,&pT,&pA,BOTH_AT);
        } else pT = pA = 0;

        bld_table(seq1-1+pT, len1-pA-pT, rs.W, INIT);
        
        count = 0; 
        while (!count || (seq_read(sf2)>0)) {
           sim4_stats_t st, rev_st; char *tok;
        
           if (count) { /* skip the first seq2, already in memory */
              h2 = SEQ_HEAD(sf2);
              seq2 = SEQ_CHARS(sf2);
              len2 = SEQ_LEN(sf2);
              tok = extract_tok(h2);
        
              if (!is_DNA(seq2, len2)) {
                char tmp[200];
                (void)sprintf(tmp,"%s sequence is not a DNA sequence.", tok);
                perror(tmp); continue;
              } 
              seq_toupper(seq2, len2, argv[2]);
           } else {  /* first sequence in the file, seq2 is already in memory */
               tok = extract_tok(h2);
               if (tok==NULL) { 
                   tok = ckalloc(strlen("(no header)")+1); 
                   strcpy(tok, "(no header)"); 
               } 
           }

           if ((rs.ali_flag==5) && (file_type==GEN_EST)) {
               cds_gene = tok;
               if (!cds_gene && !count) {
                   cds_from = rs.CDS_from; cds_to = rs.CDS_to; 
               } else if (!count) {
                   line = strstr(h2, "CDS=");
                   if (rs.S) {
                       if (line) fprintf(stderr, "Warning: Command line CDS specification overrides header CDS specification.");
                       cds_from = rs.CDS_from; cds_to = rs.CDS_to;
                   } else if (line) {
                       cds_range(line+4, &cds_from, &cds_to); 
                   }
               } else if (count) {
                   line = strstr(h2, "CDS=");
                   if (line) {
                       cds_range(line+4, &cds_from, &cds_to);
                    } else {
                       cds_from = cds_to = 0;
                    }
               }
               if (cds_to>len2) fatal("CDS endpoints exceed sequence length.");
           }

           if (rs.poly_flag && file_type==GEN_EST)  {
               get_polyAT(seq2, len2, &pT, &pA, BOTH_AT);
           }

           ++count; 
           init_stats(&st); init_stats(&rev_st);
           in_K = (rs.set_K==TRUE) ? rs.K:-1;
           in_C = (rs.set_C==TRUE) ? rs.C:-1;
           in_H = (rs.set_H==TRUE) ? rs.weight:-1;

           switch (rs.reverse) {
             case  0: Aligns = (file_type==EST_GEN) ?
                      SIM4(seq2,seq1+pT,len2,len1-pT-pA,rs.W,rs.X,in_K,in_C,
                           in_H,&dist,&xpT,&xpA,&Exons,&st):
                      SIM4(seq1,seq2+pT,len1,len2-pT-pA,rs.W,rs.X,in_K,in_C,
                           in_H,&dist,&xpT,&xpA,&Exons,&st);
                      break;
              
             case  1: sf2 = seq_revcomp_inplace(sf2);
                      seq2 = SEQ_CHARS(sf2);
                      rev_Aligns = (file_type==EST_GEN) ?
                      SIM4(seq2,seq1+pT,len2,len1-pT-pA,rs.W,rs.X,in_K,in_C,
                           in_H,&dist,&rev_xpT,&rev_xpA,&rev_Exons,&rev_st):
                      SIM4(seq1,seq2+pA,len1,len2-pT-pA,rs.W,rs.X,in_K,in_C,
                           in_H,&dist,&rev_xpT,&rev_xpA,&rev_Exons,&rev_st);
                      break; 
              
             case  2: Aligns = (file_type==EST_GEN) ?
                      SIM4(seq2,seq1+pT,len2,len1-pT-pA,rs.W,rs.X,in_K,in_C,
                           in_H,&dist,&xpT,&xpA,&Exons,&st):
                      SIM4(seq1,seq2+pT,len1,len2-pT-pA,rs.W,rs.X,in_K,in_C,
                           in_H,&dist,&xpT,&xpA,&Exons,&st);

                      sf2 = seq_revcomp_inplace(sf2);
                      seq2 = SEQ_CHARS(sf2);
                      rev_Aligns = (file_type==EST_GEN) ?
                      SIM4(seq2,seq1+pT,len2,len1-pT-pA,rs.W,rs.X,in_K,in_C,
                           in_H,&dist,&rev_xpT,&rev_xpA,&rev_Exons,&rev_st):
                      SIM4(seq1,seq2+pA,len1,len2-pT-pA,rs.W,rs.X,in_K,in_C,
                           in_H,&dist,&rev_xpT,&rev_xpA,&rev_Exons,&rev_st);
                      break;
                               
              default: fatal ("Unrecognized request for EST orientation.");
           }
              
           if (st.nmatches>=rev_st.nmatches) {
               /* forward ('+') strand match */
               match_ori = FWD;
               if (rs.reverse && rs.ali_flag) {
                   /* reverse-complement back seq2 for alignment */
                   sf2 = seq_revcomp_inplace(sf2); seq2 = SEQ_CHARS(sf2);
               }
               if (rev_Exons) { free_list(rev_Exons); rev_Exons = NULL; }
               if (rev_Aligns) { free_align(rev_Aligns); rev_Aligns = NULL; }
           } else {            
               match_ori = BWD;
               if (Exons) { free_list(Exons);  Exons = NULL; }
               if (Aligns) { free_align(Aligns); Aligns = NULL; }
           }          
                      
           if (rs.poly_flag) {
             if (match_ori==FWD) {
                 add_offset_exons(Exons, pT);  
                 add_offset_aligns(Aligns, pT);
             } else {
                 add_offset_exons(rev_Exons,(file_type==EST_GEN)?pT:pA); 
                 add_offset_aligns(rev_Aligns, (file_type==EST_GEN)?pT:pA);
             }
           } 

#ifdef DEBUG
           polystats(pT,pA,xpT,xpA,rev_xpT,rev_xpA,file_type,match_ori);
#endif

           switch (rs.ali_flag) {
              case 0: (void)printf("\nseq1 = %s (%s), %d bp\n", argv[1],tok1,len1);
                      (void)printf("seq2 = %s (%s), %d bp\n",argv[2],tok,len2);
                      if (match_ori==FWD) {
                          if (Exons) (void)printf("\n");
                          print_exons(Exons);
                      } else { 
                          (void)printf("\n(complement)");
                          if (file_type==EST_GEN) 
                              complement_exons(&rev_Exons, len2, len1);
                          if (rev_Exons) (void)printf("\n"); 
                          print_exons(rev_Exons);
                      }
                      break;
               
              case 1: (void)printf("\nseq1 = %s (%s), %d bp\n", argv[1],tok1,len1);
                      (void)printf("seq2 = %s (%s), %d bp\n",argv[2],tok,len2);
                      (void)printf("\n");
                      if (match_ori==FWD) {
                          (void)printf("\n");
                          print_align_lat(seq1, seq2, len1, len2,
                                          &Aligns, Exons, file_type, FWD);
                      } else {
                          (void)printf("\n(complement)\n");

                          if (file_type==EST_GEN) {
                              complement_exons(&rev_Exons, len2, len1);
                              sf2 = seq_revcomp_inplace(sf2); 
                              seq2 = SEQ_CHARS(sf2);
                          }
                          print_align_lat((file_type==EST_GEN) ? revseq1:seq1, 
                          seq2,len1,len2,&rev_Aligns,rev_Exons,file_type,BWD);
                      }
                      break;

              case 2: print_align_header(sf1, sf2, &ds);
                      if (match_ori==FWD) {
                          (void)printf("\n");
                          print_align_blk(seq1,seq2,len1,len2,&Aligns,
                                          file_type,FWD);
                      } else {
                          if (file_type==EST_GEN) {
                              complement_exons(&rev_Exons, len2, len1);
                              sf2 = seq_revcomp_inplace(sf2); 
                              seq2 = SEQ_CHARS(sf2);
                          }
                          (void)printf("\n(complement)\n");
                          print_align_blk((file_type==EST_GEN) ? revseq1:seq1,
                          seq2,len1,len2,&rev_Aligns,file_type,BWD);
                      }
                      break;

              case 3: (void)printf("\nseq1 = %s (%s), %d bp\n", argv[1],tok1,len1);
                      (void)printf("seq2 = %s (%s), %d bp\n",argv[2],tok,len2);
                      (void)printf("\n");
                      if (match_ori==FWD) {
                          if (Exons) (void)printf("\n");
                          print_exons(Exons); (void)printf("\n");
                          print_align_lat(seq1, seq2, len1, len2, &Aligns,
                                          Exons, file_type, FWD);
                      } else {
                          if (file_type==EST_GEN) {
                              complement_exons(&rev_Exons, len2, len1);
                              sf2 = seq_revcomp_inplace(sf2);
                              seq2 = SEQ_CHARS(sf2);
                          }
                          (void)printf("\n(complement)");
                          if (rev_Exons) (void)printf("\n"); 
                          print_exons(rev_Exons); (void)printf("\n");
                          print_align_lat((file_type==EST_GEN) ? revseq1:seq1, 
                          seq2,len1,len2,&rev_Aligns,rev_Exons,file_type,BWD);
                      }
                      break;

              case 4: (void)printf("\nseq1 = %s (%s), %d bp\n",argv[1],tok1,len1);
                      (void)printf("seq2 = %s (%s), %d bp\n",argv[2],tok,len2);
                      (void)printf("\n%s\n", SEQ_HEAD(sf1));
                      (void)printf("%s", SEQ_HEAD(sf2));
                      (void)printf("\n");
                      if (match_ori==FWD) {
                          if (Exons) (void)printf("\n");
                          print_exons(Exons); (void)printf("\n");
                          print_align_lat(seq1,seq2,len1,len2,&Aligns,Exons,file_type,FWD);
                      } else {
                          /* seq2 is already reversed, need only interchange
                             from1<->from2, etc, in Exons; the alignment should
                             be read in the forward orientation */ 
                          (void)printf("\n(complement)\n");
                          if (rev_Exons) (void)printf("\n"); 
                          print_exons(rev_Exons); (void)printf("\n");
                          print_align_lat(seq1,seq2,len1,len2,&rev_Aligns,rev_Exons,file_type,FWD);
                      }

                      break;

              case 5: if (match_ori==FWD) {
                         (void)printf("\n");
                         if (file_type==EST_GEN) {
                           print_pipmaker_exons(Exons,Aligns,cds_gene, 
                                 cds_from,cds_to,len2,len1,seq2,seq1,FWD); 
                         } else { 
                           print_pipmaker_exons(Exons, Aligns, cds_gene,
                                 cds_from,cds_to,len1,len2,seq1,seq2,FWD);
                         }
                      } else {
                         (void)printf("\n");
                         /* give it the "real" sequences */
                         sf2 = seq_revcomp_inplace(sf2);
                         seq2 = SEQ_CHARS(sf2);
                         if (file_type==EST_GEN) {
                            complement_exons(&rev_Exons, len2, len1);
                            print_pipmaker_exons(rev_Exons,rev_Aligns,cds_gene,
                                  cds_from,cds_to,len2,len1,seq2,seq1,BWD);
                         } else {
                            print_pipmaker_exons(rev_Exons,rev_Aligns,cds_gene,
                                  cds_from,cds_to,len1,len2,seq1,seq2,BWD);
                         }
                      }
                      break;
 
              default:fatal("Unrecognized option for alignment output.");
           }
#ifdef _STATS
           print_exon_stats((match_ori==FWD) ? Exons:rev_Exons,
                            (file_type==EST_GEN) ? argv[1]:tok+1,
                            (file_type==EST_GEN) ? len1-pA-pT:len2-pA-pT);
           print_stats((match_ori==FWD) ? st:rev_st,
                          (file_type==EST_GEN) ? argv[1]:tok+1,
                          (file_type==EST_GEN) ? len1-pA-pT:len2-pA-pT);
#endif                

           (void)printf("\n");
           if (Aligns) { free_align(Aligns); Aligns = NULL; }
           if (rev_Aligns) { free_align(rev_Aligns); rev_Aligns = NULL; }
           if (Exons) { free_list(Exons); Exons = NULL; }
           if (rev_Exons) { free_list(rev_Exons); rev_Exons = NULL; }
           if (tok) { free(tok); tok = NULL; }
       }

       if ((count==1) && (file_type==GEN_EST))
           (void)fprintf(stderr,"Try shorter sequence first for better performance.\n");  

       if (tok1) free(tok1);
       free_table();
       if (file_type==EST_GEN) seq_close(rf1);
       seq_close(sf1);
       seq_close(sf2);

       return 0;
}

static void print_align_blk(uchar *seq1, uchar *seq2, int len1, int len2,
                            edit_script_list **Aligns, 
                            int file_type, int match_ori)
{
    int *S;
    edit_script_list *head, *aligns;

    if (*Aligns==NULL) return;

    aligns = *Aligns;
    while (aligns!=NULL) {
       head = aligns;
       aligns = aligns->next_script;

       S = (int *)ckalloc((2*head->len2+1+1)*sizeof(int));
       S++;
       S2A(head->script, S, (file_type==EST_GEN) ? 1:0);
       Free_script(head->script);

       if (file_type==EST_GEN) {
           if (match_ori==FWD) {
               print_align(head->score,seq1,seq2,
                           head->offset2,   
                           head->offset2+head->len2-1,   
                           head->offset1,    
                           head->offset1+head->len1-1,S);   
           } else {
               align_reverse(S);
               print_align(head->score,seq1,seq2,
                           len1+1-(head->offset2+head->len2-1),
                           len1+1-head->offset2,
                           len2+1-(head->offset1+head->len1-1),
                           len2+1-head->offset1,S);
           } 
       } else {    /* file_type==GEN_EST */
           print_align(head->score, seq1, seq2,
                       head->offset1, 
                       head->offset1+head->len1-1,
                       head->offset2,
                       head->offset2+head->len2-1,S);
       }
       free(S-1);
       free(head);
    }
    *Aligns = NULL;
    return;
}

static void print_align_lat(uchar *seq1, uchar *seq2, int len1, int len2, 
                            edit_script_list **Aligns, Exon *Exons, 
                            int file_type, int match_ori) 
{
    int *S;
    edit_script_list *head, *aligns;
                           
    if (*Aligns==NULL) return;
           
    aligns = *Aligns;
    while (aligns!=NULL) {
       head = aligns;
       aligns = aligns->next_script; 
                       
       S = (int *)ckalloc((2*head->len2+1+1)*sizeof(int));
       S++;            
       S2A(head->script, S, (file_type==1) ? 1:0);
       Free_script(head->script);
       
       if (file_type==EST_GEN) {
           if (match_ori==FWD) {
               IDISPLAY(seq1+ head->offset2-1-1, seq2+ head->offset1-1-1,
                        head->len2, head->len1, S,
                        head->offset2, head->offset1, 1, Exons);
           } else {        
               align_reverse(S);
               IDISPLAY(seq1+len1+1-(head->offset2+head->len2-1)-1-1,
                        seq2+len2+1-(head->offset1+head->len1-1)-1-1,
                        head->len2, head->len1, S,
                        len1+1-(head->offset2+head->len2-1),
                        len2+1-(head->offset1+head->len1-1), 1, Exons);
           }
       } else {    /* file_type==GEN_EST */
           IDISPLAY(seq1+ head->offset1-1-1, seq2+ head->offset2-1-1,
                       head->len1, head->len2, S,
                       head->offset1, head->offset2, 2, Exons);
       }
       free(S-1);
       free(head);
    }                      
    *Aligns = NULL;
    return;                
}


static void sim4_argvals(sim4_args_t *args)
{
        
        if (!get_argval('A', &(args->ali_flag)))
             args->ali_flag = 0; 
        if ((args->ali_flag>5) || (args->ali_flag<0))
             fatal("A options: 0, 1, 2, 3, 4, 5.");
        
        if (!get_argval('P', &(args->poly_flag)))
             args->poly_flag = 0;
        
        if (get_argval('R',&(args->reverse))) {
                if ((args->reverse<0) || (args->reverse>2))
                        fatal("Direction R must be 0, 1, or 2.");
        } else
                args->reverse = 2;
        
        if (get_argval('E', &(args->cutoff))) {
                if ((args->cutoff < 3) || (args->cutoff > 10))
                        fatal("Cutoff must be within [3,10].");
        } else  
                args->cutoff = DIST_CUTOFF;
        
        if (get_argval('D', &(args->DRANGE))) {
                if (args->DRANGE < 0)
                        fatal("Positive number required for D.");
        } else
                args->DRANGE = DEFAULT_DRANGE;
        
        if (get_argval('H', &(args->weight))) {
                if (args->weight<0)
                        fatal("Positive number required for H.");
                args->set_H = TRUE;
        } else {
            /*  args->weight = DEFAULT_WEIGHT; */
                args->set_H = FALSE;
        }
         
                        
/*      
        if (get_fargval('v',&V)) {
                if ((V<.7) || (V>1.0))
                        fatal("V must be in the range [0.7,1.0].");
        } else  
                V = DEFAULT_MIN_COV;
*/      
                
        if (get_argval('W', &(args->W))) {
                if (args->W < 1)
                        fatal("W must be positive.");
                if (args->W > 15)
                        fatal("W must be <= 15.");
        } else  
                args->W = DEFAULT_W;
        
        if (get_argval('X', &(args->X))) {
                if (args->X < 1)
                        fatal("X must be positive.");
        } else  
                args->X = DEFAULT_X;
        
        if (get_argval('K',&(args->K))) {
                if (args->K<0) fatal("K must be positive."); 
                args->set_K = TRUE;
        } else {
                args->K = DEFAULT_K;
                args->set_K = FALSE;
        }

        if (get_argval('C',&(args->C))) { 
                if (args->C<0) fatal("C must be positive.");                    
                args->set_C = TRUE;
        } else {
                args->C = DEFAULT_C;
                args->set_C = FALSE;
        }
        
        if (!get_argval('N', &(args->acc_flag)))
             args->acc_flag = 0;
                
        if (get_argval('B', &(args->B))) {
                 if (args->B && (args->B!=1))
                         fatal("B must be either 0 or 1.");
        } else  
                args->B = 1;
        
        if (get_strargval('S', &(args->S))) {
            cds_range(args->S, &(args->CDS_from), &(args->CDS_to));
            if ((args->CDS_from<=0) || (args->CDS_to<=0) || 
                (args->CDS_from>args->CDS_to))
                fatal("Illegal endpoints for the CDS region.");
        } else 
                args->S = NULL;

        if (args->S && (args->ali_flag!=5))
           fatal ("A=5 must accompany CDS specification.");

        return;
}

/* extract the CDS endpoints from the command line specification <n1>..<n2> */
static void cds_range(char *line, int *from, int *to)
{
     char *s = line;

     if (line == NULL) fatal ("NULL CDS specification.");

     if (!isdigit((int)(*s))) 
         fatal("Non-numerical value in the CDS specification."); 
     while (*s && isdigit((int)(*s))) s++;
     if (*s!='.') fatal ("Illegal CDS specification."); s++;
     if (*s!='.') fatal ("Illegal CDS specification."); s++;
     if (!isdigit((int)(*s))) 
         fatal ("Non-numerical value in the CDS specification.");
     while (*s && isdigit((int)(*s))) s++;   
     if (*s && !isspace((int)(*s))) 
        fatal ("Garbage at the end of the CDS numerical specification."); 
     
     /* now extract the CDS elements */
     if (sscanf(line, "%d..%d", from, to)!=2) 
         fatal ("Error when reading the CDS endpoints.");

     return;
}

static void add_offset_exons(Exon *exons, int offset)
{
    Exon *t;

    if (!offset || !(exons)) return;
 
    t = exons;
    while (t) {
       if (t->to1) { t->from2 += offset; t->to2 += offset; }
       t = t->next_exon;
    }
}

static void add_offset_aligns(edit_script_list *aligns, int offset)
{
    edit_script_list *head;
                 
    if (!offset || !aligns) return;
             
    head = aligns;
    while (head) { head->offset2 += offset; head = head->next_script; }
    
    return;
}

static char *extract_tok(char *h2)
{
    char *s, *tmp, *q;

    if ((h2==NULL) || (*h2=='\0')) return NULL;

    if (*h2!='>') fatal("Not a FASTA header.");
    s = h2+1; while (isspace((int)(*s)) && *s!='\n') s++;
    if (*s=='\n') return NULL;
    q = s; while (*s && !isspace((int)(*s))) s++; 
    tmp = (char *) ckalloc((unsigned int)(s-q+1));
    strncpy(tmp, q, (int)(s-q));
    tmp[s-q] = '\0'; 

    return tmp;
}


/*  ----------  utilities for collecting and reporting statistics ---------- */

static void init_stats(sim4_stats_t *st)
{
       (void)memset(st,0,sizeof(sim4_stats_t));
}

#ifdef _STATS
static void print_exon_stats(Exon *exons, char *seq_name, int len)
{
       FILE *xp;
       Exon *t;

       xp = ckopen("EXONS","a");
       (void)fprintf(xp,"%s:   ", seq_name);
       (void)fprintf(xp,"{\n  %3d\n",len);
                   
       t = exons;
       while (t!=NULL) {
          if (t->to1) (void)fprintf(xp,"  %6d %6d %3d %3d\n",
                                        t->from1,t->to1,t->from2,t->to2);
          t = t->next_exon;
       } 
       (void)fprintf(xp,"}\n");
       if (fclose(xp)==EOF) perror("sim4.init.c: fclose failed.");

       return;
}
#endif

#ifdef _STATS
static void print_stats(sim4_stats_t st, char *seq_name, int len)
{      
       FILE *fp;

       fp = ckopen("SCORES","a"); 
       (void)fprintf(fp,"%s:\t   %5.4f   (%1d / %1d)  %d %f %d",
                         seq_name, st.fcoverage, st.icoverage, len,
                         st.internal, st.marginals, st.nmatches);
                   
       if (st.mult>1) (void)fprintf(fp,"   (%d)\n",st.mult);
       else if (!st.mult) (void)fprintf(fp,"\t\t*\n");
       else (void)fprintf(fp,"\n");
       if (fclose(fp)==EOF)  perror("sim4.init.c: fclose failed.");

       return;
}
#endif

#ifdef DEBUG
static void polystats(int pT,int pA,int xpT,int xpA,int rev_xpT,int rev_xpA, 
                      int file_type,int match_ori)
{
      int tmp;

      if (file_type==EST_GEN) {
          printf("Poly(EG): %d %d (%d) %d %d (%d)\n",
                  pT, tmp=((match_ori==FWD) ? xpT:rev_xpT), pT+tmp,
                  pA, tmp=((match_ori==FWD) ? xpA:rev_xpA), pA+tmp);
      } else {
          printf("Poly(GE): %d %d (%d) %d %d (%d)\n",
                  pT, tmp=((match_ori==FWD) ? xpT:rev_xpA), pT+tmp,
                  pA, tmp=((match_ori==FWD) ? xpA:rev_xpT), pA+tmp);
      }
}
#endif

