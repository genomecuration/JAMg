/****************************************************************\
*                                                                *
*  Analysis module for exonerate                                 *
*                                                                *
*  Guy St.C. Slater..   mailto:guy@ebi.ac.uk                     *
*  Copyright (C) 2000-2008.  All Rights Reserved.                *
*                                                                *
*  This source code is distributed under the terms of the        *
*  GNU General Public License, version 3. See the file COPYING   *
*  or http://www.gnu.org/licenses/gpl.txt for details            *
*                                                                *
*  If you use this code, please keep this notice intact.         *
*                                                                *
\****************************************************************/

#include "analysis.h"
#include "ungapped.h"
#include "compoundfile.h"
#include "hspset.h"

#include <stdlib.h> /* For atoi() */
#include <unistd.h> /* For sleep() */
#include <string.h> /* For strcmp() */
#include <ctype.h>  /* For isdigit() */

Analysis_ArgumentSet *Analysis_ArgumentSet_create(Argument *arg){
    register ArgumentSet *as;
    static Analysis_ArgumentSet aas;
    if(arg){
        as = ArgumentSet_create("Analysis Options");
        ArgumentSet_add_option(as, 'E', "exhaustive", NULL,
            "Perform exhaustive alignment (slow)", "FALSE",
            Argument_parse_boolean, &aas.use_exhaustive);
        ArgumentSet_add_option(as, 'B', "bigseq", NULL,
            "Allow rapid comparison between big sequences", "FALSE",
            Argument_parse_boolean, &aas.use_bigseq);
        ArgumentSet_add_option(as, '\0', "forcescan", "[q|t]",
            "Force FSM scan on query or target sequences", "none",
            Argument_parse_string, &aas.force_scan);
        /**/
        ArgumentSet_add_option(as, 0, "saturatethreshold", "int",
            "Word saturation threshold", "0",
            Argument_parse_int, &aas.saturate_threshold);
        /**/
        ArgumentSet_add_option(as, 0, "customserver", "command",
            "Custom command to send non-standard server", "NULL",
            Argument_parse_string, &aas.custom_server_command);
        Argument_absorb_ArgumentSet(arg, as);
        }
    return &aas;
    }

/**/

static void Analysis_report_func(Comparison *comparison,
                                 gpointer user_data){
    register Analysis *analysis = user_data;
    register GAM_Result *gam_result;
    g_assert(Comparison_has_hsps(comparison));
    if(analysis->scan_query){
        /* Swap back query and target after a query scan */
        Comparison_swap(comparison);
    } else {
        /* Revert target scan alignments to +ve query strand when possible */
        if((comparison->query->alphabet->type == Alphabet_Type_DNA)
        && (comparison->target->alphabet->type == Alphabet_Type_DNA)
        && (comparison->query->strand == Sequence_Strand_REVCOMP)
        && (comparison->target->strand == Sequence_Strand_FORWARD)
        && (!analysis->gam->translate_both))
            Comparison_revcomp(comparison);
        }
    if(Model_Type_is_gapped(analysis->gam->gas->type)){
        gam_result = GAM_Result_heuristic_create(analysis->gam,
                                                 comparison);
    } else {
        gam_result = GAM_Result_ungapped_create(analysis->gam,
                                                comparison);
        }
    if(gam_result){
        GAM_Result_submit(gam_result);
        GAM_Result_destroy(gam_result);
        }
    return;
    }

/**/

static void Analysis_FastaPipe_Pair_init_func(gpointer user_data){
    register Analysis *analysis = user_data;
    g_assert(!analysis->curr_query);
    return;
    }
/* Called before query pipeline loading */

static void Analysis_FastaPipe_Pair_prep_func(gpointer user_data){
    register Analysis *analysis = user_data;
    g_assert(analysis->curr_query);
    return;
    }
/* Called after query pipeline loading */

static void Analysis_FastaPipe_Pair_term_func(gpointer user_data){
    register Analysis *analysis = user_data;
    FastaDB_Seq_destroy(analysis->curr_query);
    analysis->curr_query = NULL;
    return;
    }
/* Called after query pipeline analysis */

static gboolean Analysis_FastaPipe_Pair_query_func(FastaDB_Seq *fdbs,
                                                   gpointer user_data){
    register Analysis *analysis = user_data;
    g_assert(!analysis->curr_query);
    /*
    if(analysis->aas->use_exhaustive)
        g_strup(fdbs->seq->seq);
    */
    if(analysis->verbosity > 1)
        g_message("Load query for pairwise comparision [%s] (%d)",
                  fdbs->seq->id, fdbs->seq->len);
    analysis->curr_query = FastaDB_Seq_share(fdbs);
    return TRUE; /* take queries one at a time */
    }
/* Called on query loading */

static void Analysis_BSAM_compare(Analysis *analysis,
                                  FastaDB_Seq *query,
                                  FastaDB_Seq *target){
    register Comparison *comparison = BSAM_compare(analysis->bsam,
                                                   query->seq,
                                                   target->seq);
    if(comparison){
        if(Comparison_has_hsps(comparison))
            Analysis_report_func(comparison, analysis);
        Comparison_destroy(comparison);
        }
    return;
    }

/**/

static void Analysis_Pair_compare(Analysis *analysis,
                                  FastaDB_Seq *fdbs){
    register GAM_Result *gam_result;
    if(analysis->aas->use_exhaustive){
        gam_result = GAM_Result_exhaustive_create(analysis->gam,
                         analysis->curr_query->seq, fdbs->seq);
        if(gam_result){
            GAM_Result_submit(gam_result);
            GAM_Result_destroy(gam_result);
            }
    } else {
        Analysis_BSAM_compare(analysis, analysis->curr_query, fdbs);
        }
    return;
    }

static gboolean Analysis_FastaPipe_Pair_target_func(FastaDB_Seq *fdbs,
                                                    gpointer user_data){
    register Analysis *analysis = user_data;
    g_assert(analysis->gam);
    g_assert(analysis->curr_query);
    /*
    if(analysis->aas->use_exhaustive)
        g_strup(fdbs->seq->seq);
    */
    if(analysis->verbosity > 1)
        g_message("Load target for pairwise comparison [%s] (%d)",
                fdbs->seq->id, fdbs->seq->len);
    /**/
    Analysis_Pair_compare(analysis, fdbs);
    return FALSE;
    }
/* Called on target loading */

/**/

static void Analysis_FastaPipe_Seeder_init_func(gpointer user_data){
    register Analysis *analysis = user_data;
    g_assert(!analysis->curr_seeder);
    register Comparison_Param *comparison_param
           = Comparison_Param_share(analysis->comparison_param);
    if(analysis->scan_query)
        comparison_param = Comparison_Param_swap(comparison_param);
    analysis->curr_seeder = Seeder_create(analysis->verbosity,
            comparison_param, analysis->aas->saturate_threshold,
            Analysis_report_func, analysis);
    Comparison_Param_destroy(comparison_param);
    return;
    }
/* Called before query pipeline loading */

static void Analysis_FastaPipe_Seeder_prep_func(gpointer user_data){
    /* does nothing */
    return;
    }
/* Called after query pipeline loading */

static void Analysis_FastaPipe_Seeder_term_func(gpointer user_data){
    register Analysis *analysis = user_data;
    Seeder_destroy(analysis->curr_seeder);
    if(analysis->verbosity > 2){
        g_message("### Seeder destroyed ###");
        RecycleBin_profile();
        }
    analysis->curr_seeder = NULL;
    return;
    }
/* Called after query pipeline analysis */

static gboolean Analysis_FastaPipe_Seeder_query_func(FastaDB_Seq *fdbs,
                                                    gpointer user_data){
    register Analysis *analysis = user_data;
    if(analysis->verbosity > 1)
        g_message("Load query for Seeder [%s] (%d)",
                fdbs->seq->id, fdbs->seq->len);
    return Seeder_add_query(analysis->curr_seeder, fdbs->seq);
    }
/* Called on query loading */

static gboolean Analysis_FastaPipe_Seeder_target_func(FastaDB_Seq *fdbs,
                                                    gpointer user_data){
    register Analysis *analysis = user_data;
    if(analysis->verbosity > 1)
        g_message("Load target for Seeder [%s] (%d)",
                fdbs->seq->id, fdbs->seq->len);
    Seeder_add_target(analysis->curr_seeder, fdbs->seq);
    return FALSE;
    }
/* Called on target loading */

/**/

static gboolean Analysis_decide_scan_query(FastaDB *query_fdb,
                                           FastaDB *target_fdb,
                                           gchar *force_scan){
    register CompoundFile_Pos query_size, target_size;
    if(!g_strcasecmp(force_scan, "none")){
        query_size = CompoundFile_get_length(query_fdb->cf);
        target_size = CompoundFile_get_length(target_fdb->cf);
        if((query_size >> 4) < target_size)
            return FALSE;
        else
            return TRUE;
    } else if((!g_strcasecmp(force_scan, "query"))
           || (!g_strcasecmp(force_scan, "q"))){
        return TRUE;
    } else if((!g_strcasecmp(force_scan, "target"))
           || (!g_strcasecmp(force_scan, "t"))){
        return FALSE;
    } else {
        g_error("Unknown force_scan command [%s]", force_scan);
        }
    return FALSE; /* not reached */
    }

static void Analysis_find_matches(Analysis *analysis,
                       Match **dna_match, Match **protein_match,
                       Match **codon_match){
    register GPtrArray *transition_list
        = C4_Model_select_transitions(analysis->gam->model,
                                      C4_Label_MATCH);
    register gint i;
    register C4_Transition *transition;
    register Match *match;
    for(i = 0; i < transition_list->len; i++){
        transition = transition_list->pdata[i];
        g_assert(transition);
        g_assert(transition->label == C4_Label_MATCH);
        match = transition->label_data;
        if(match){
            switch(match->type){
                case Match_Type_DNA2DNA:
                    if((*dna_match) && ((*dna_match) != match))
                        g_error("Multiple DNA matches not implemented");
                    (*dna_match) = match;
                    break;
                case Match_Type_PROTEIN2PROTEIN:
                case Match_Type_DNA2PROTEIN:
                case Match_Type_PROTEIN2DNA:
                    if((*protein_match) && ((*protein_match) != match))
                        g_error("Multiple protein matches"
                                " not supported");
                    (*protein_match) = match;
                    break;
                case Match_Type_CODON2CODON:
                    if((*codon_match) && ((*codon_match) != match))
                        g_error("Multiple codon matches"
                                " not implemented");
                    (*codon_match) = match;
                    break;
                default:
                    break;
                }
            }
        }
    g_ptr_array_free(transition_list, TRUE);
    return;
    }

/**/

static SocketClient *Analysis_Client_connect(gchar *path){
    register SocketClient *sc;
    register gint i, divider = 0, port;
    register gchar *server;
    register gint connection_attempts = 10;
    for(i = 0; path[i]; i++)
        if(path[i] == ':'){
            divider = i;
            break;
            }
    if(divider){
        port = atoi(path+divider+1);
        server = g_strndup(path, divider);
        for(i = connection_attempts; i >= 0; i--){
            sc = SocketClient_create(server, port);
            if(sc)
                break;
            g_warning("Failed connection to server (retrys left: %d", i);
            sleep(1);
            }
        g_free(server);
        return sc;
        }
    return NULL;
    }

static gchar *Analysis_Client_send(Analysis_Client *aclient, gchar *msg,
                                  gchar *expect, gboolean multi_line_reply){
    register gchar *reply = SocketClient_send(aclient->sc, msg);
    register gchar *line, *p = reply, *processed_reply;
    register gint line_count = 0;
    register GString *str = g_string_sized_new(64);
    if(aclient->verbosity >= 3){
        g_print("Message: client sent message [%s]\n", msg);
        if((aclient->verbosity == 3) && (strlen(reply) >= 80))
            g_print("Message: client received reply [%.*s<truncated>] \n",
                               80, reply);
        else
            g_print("Message: client received reply [%s]\n", reply);
        }
    do {
        while(*p){
            if(!isspace(*p)) /* Strip blank lines */
                break;
            p++;
            }
        line = p;
        if(!strncmp(p, "error:", 6))
            g_error("Error from server: [%s]", p);
        if(!strncmp(p, "warning:", 8)){
            g_warning("Warning from server: [%s]", p);
        } else if(!strncmp(p, "linecount:", 8)){
            while(*p){ /* Skip to end of line */
                if(*p == '\n')
                    break;
                p++;
                }
        } else if(!strncmp(p, expect, strlen(expect))){
            while(*p){ /* Skip to end of line */
                if(*p == '\n')
                    break;
                p++;
                }
            line_count++;
            *p = '\0';
            g_string_append(str, line);
            g_string_append_c(str, '\n');
            *p = '\n';
        } else {
            if(!*p)
                break;
            g_error("Unexpected line from server [%s]", p);
            }
    } while(TRUE);
    g_free(reply);
    if(!line_count)
        g_error("No reply received from server msg=[%s]", msg);
    if((!multi_line_reply) && (line_count > 1))
        g_error("Unexpected multi-line reply from server");
    processed_reply = str->str;
    g_string_free(str, FALSE);
    return processed_reply;
    }

static void Analysis_Client_set_param(Analysis_Client *aclient, GAM *gam){
    register HSPset_ArgumentSet *has = HSPset_ArgumentSet_create(NULL);
    register Analysis_ArgumentSet *aas = Analysis_ArgumentSet_create(NULL);
    register gchar *msg, *reply;
    /**/
    if(aas->custom_server_command){
        msg = g_strdup_printf("%s\n", aas->custom_server_command);
        reply = Analysis_Client_send(aclient, msg, "ok:", FALSE);
        g_free(msg);
        g_free(reply);
        }
    /**/
    msg = g_strdup_printf("set param seedrepeat %d",
                          has->seed_repeat);
    reply = Analysis_Client_send(aclient, msg, "ok:", FALSE);
    g_free(msg);
    g_free(reply);
    /**/
    msg = g_strdup_printf("set param dnahspthreshold %d",
                          has->dna_hsp_threshold);
    reply = Analysis_Client_send(aclient, msg, "ok:", FALSE);
    g_free(msg);
    g_free(reply);
    /**/
    msg = g_strdup_printf("set param proteinhspthreshold %d",
                          has->protein_hsp_threshold);
    reply = Analysis_Client_send(aclient, msg, "ok:", FALSE);
    g_free(msg);
    g_free(reply);
    /**/
    msg = g_strdup_printf("set param codonhspthreshold %d",
                          has->codon_hsp_threshold);
    reply = Analysis_Client_send(aclient, msg, "ok:", FALSE);
    g_free(msg);
    g_free(reply);
    /**/
    msg = g_strdup_printf("set param dnahspdropoff %d",
                          has->dna_hsp_dropoff);
    reply = Analysis_Client_send(aclient, msg, "ok:", FALSE);
    g_free(msg);
    g_free(reply);
    /**/
    msg = g_strdup_printf("set param proteinhspdropoff %d",
                          has->protein_hsp_dropoff);
    reply = Analysis_Client_send(aclient, msg, "ok:", FALSE);
    g_free(msg);
    g_free(reply);
    /**/
    msg = g_strdup_printf("set param codonhspdropoff %d",
                          has->codon_hsp_dropoff);
    reply = Analysis_Client_send(aclient, msg, "ok:", FALSE);
    g_free(msg);
    g_free(reply);
    /**/
    msg = g_strdup_printf("set param geneseedthreshold %d",
                          has->geneseed_threshold);
    reply = Analysis_Client_send(aclient, msg, "ok:", FALSE);
    g_free(msg);
    g_free(reply);
    /**/
    msg = g_strdup_printf("set param geneseedrepeat %d",
                          has->geneseed_repeat);
    reply = Analysis_Client_send(aclient, msg, "ok:", FALSE);
    g_free(msg);
    g_free(reply);
    /**/
    msg = g_strdup_printf("set param maxqueryspan %d",
                          gam->max_query_span);
    reply = Analysis_Client_send(aclient, msg, "ok:", FALSE);
    g_free(msg);
    g_free(reply);
    /**/
    msg = g_strdup_printf("set param maxtargetspan %d",
                          gam->max_target_span);
    reply = Analysis_Client_send(aclient, msg, "ok:", FALSE);
    g_free(msg);
    g_free(reply);
    /**/
    return;
    }

static Analysis_Client *Analysis_Client_create(gchar *path, gint verbosity){
    register Analysis_Client *aclient;
    register SocketClient *sc = Analysis_Client_connect(path);
    register Alphabet_Type alphabet_type;
    register gchar *dbinfo, **dbinfo_word;
    if(!sc)
        return NULL;
    aclient = g_new(Analysis_Client, 1);
    aclient->sc = sc;
    aclient->verbosity = verbosity;
    aclient->probe_fdb = NULL;
    dbinfo = Analysis_Client_send(aclient, "dbinfo", "dbinfo:", FALSE);
    dbinfo_word = g_strsplit(dbinfo, " ", 8);
    /**/
    g_assert(dbinfo_word[0]);
    g_assert(dbinfo_word[1]);
    g_assert(dbinfo_word[2]);
    g_assert(dbinfo_word[3]);
    g_assert(dbinfo_word[4]);
    /**/
    alphabet_type = Alphabet_Type_UNKNOWN;
    if(!strcmp(dbinfo_word[1], "dna"))
        alphabet_type = Alphabet_Type_DNA;
    if(!strcmp(dbinfo_word[1], "protein"))
        alphabet_type = Alphabet_Type_PROTEIN;
    /**/
    aclient->is_masked = FALSE;
    if(!strcmp(dbinfo_word[2], "softmasked"))
        aclient->is_masked = TRUE;
    else if(!strcmp(dbinfo_word[2], "unmasked"))
        aclient->is_masked = TRUE;
    else
        g_error("Unrecognised dbinfo masking state [%s]", dbinfo_word[2]);
    /**/
    aclient->server_alphabet = Alphabet_create(alphabet_type, aclient->is_masked);
    aclient->num_seqs      = atoll(dbinfo_word[3]);
    aclient->max_seq_len   = atoll(dbinfo_word[4]);
    aclient->total_seq_len = atoll(dbinfo_word[5]);
    /**/
    g_strfreev(dbinfo_word);
    g_free(dbinfo);
    aclient->curr_query = NULL;
    aclient->seq_cache = g_new0(Sequence*, aclient->num_seqs);
    return aclient;
    }

static void Analysis_Client_destroy(Analysis_Client *aclient){
    register gint i;
    register Sequence *seq;
    if(aclient->curr_query)
        Sequence_destroy(aclient->curr_query);
    for(i = 0; i < aclient->num_seqs; i++){
        seq = aclient->seq_cache[i];
        if(seq)
            Sequence_destroy(seq);
        }
    g_free(aclient->seq_cache);
    if(aclient->probe_fdb)
        FastaDB_close(aclient->probe_fdb);
    SocketClient_destroy(aclient->sc);
    Alphabet_destroy(aclient->server_alphabet);
    g_free(aclient);
    return;
    }

static void Analysis_Client_set_probe_fdb(Analysis_Client *aclient,
                                          FastaDB *probe_fdb){
    g_assert(!aclient->probe_fdb);
    aclient->probe_fdb = FastaDB_share(probe_fdb);
    return;
    }

static void Analysis_Client_set_query(Analysis_Client *aclient, Sequence *seq){
    register gchar *seq_str = Sequence_get_str(seq);
    register gchar *msg = g_strdup_printf("set query %s", seq_str);
    register gchar *reply = Analysis_Client_send(aclient, msg, "ok:", FALSE);
    register gchar **word;
    register gint len, checksum;
    if(strncmp(reply, "ok:", 3))
        g_error("Could not set query [%s] on server", seq->id);
    word = g_strsplit(reply+4, " ", 4);
    len = atoi(word[0]);
    checksum = atoi(word[1]);
    if(seq->len != len)
        g_error("Query length mismatch on server %d %d", seq->len, len);
    if(Sequence_checksum(seq) != checksum)
        g_error("Query checksum mismatch on server %d %d",
                Sequence_checksum(seq), checksum);
    g_strfreev(word);
    g_free(reply);
    g_free(msg);
    g_free(seq_str);
    if(aclient->curr_query)
        Sequence_destroy(aclient->curr_query);
    aclient->curr_query = Sequence_share(seq);
    return;
    }

static void Analysis_Client_revcomp_query(Analysis_Client *aclient){
    register gchar *reply = Analysis_Client_send(aclient,
            "revcomp query", "ok: query strand revcomp", FALSE);
    register Sequence *curr_query;
    curr_query = aclient->curr_query;
    aclient->curr_query = Sequence_revcomp(aclient->curr_query);
    Sequence_destroy(curr_query);
    g_free(reply);
    return;
    }

static void Analysis_Client_revcomp_target(Analysis_Client *aclient){
    register gchar *reply = Analysis_Client_send(aclient,
            "revcomp target", "ok: target strand", FALSE);
    g_free(reply);
    return;
    }

/**/

typedef struct {
    Analysis_Client *aclient;
               gint  target_id;
               gint  seq_len;
} Analysis_Client_Key;

static Analysis_Client_Key *Analysis_Client_Key_create(Analysis_Client *aclient,
                                                   gint target_id, gint seq_len){
    register Analysis_Client_Key *key = g_new(Analysis_Client_Key, 1);
    key->aclient = aclient;
    key->target_id = target_id;
    key->seq_len = seq_len;
    return key;
    }

static void Analysis_Client_Key_destroy(Analysis_Client_Key *key){
    g_free(key);
    return;
    }

static gpointer Analysis_Client_SparseCache_get_func(gint pos,
                                                     gpointer page_data,
                                                     gpointer user_data){
    return GINT_TO_POINTER((gint)((gchar*)page_data)[pos]);
    }

static SparseCache_Page *Analysis_Client_SparseCache_fill_func(gint start,
                                                     gpointer user_data){
    register Analysis_Client_Key *key = user_data;
    register SparseCache_Page *page = g_new(SparseCache_Page, 1);
    register gint len = MIN(SparseCache_PAGE_SIZE, key->seq_len-start),
                  page_len;
    register gchar *msg = g_strdup_printf("get subseq %d %d %d",
            key->target_id, start, len);
    register gchar *reply = Analysis_Client_send(key->aclient, msg,
                                                "subseq:", FALSE);
    if(strncmp(reply, "subseq:", 7))
        g_error("Failed to get subseq for target (%d,%d,%d) [%s]",
                key->target_id, start, len, reply);
    page->get_func = Analysis_Client_SparseCache_get_func;
    page->copy_func = NULL;
    page_len = strlen(reply+8)-1;
    page->data = g_strndup(reply+8, page_len);
    page->data_size = sizeof(gchar)*page_len;
    g_free(msg);
    g_free(reply);
    FastaDB_SparseCache_compress(page, page_len);
    return page;
    }
/* FIXME: move compression stuff to SeqPage in Sequence */

static void Analysis_Client_SparseCache_free_func(gpointer user_data){
    register Analysis_Client_Key *key = user_data;
    Analysis_Client_Key_destroy(key);
    return;
    }

static SparseCache *Analysis_Client_get_SparseCache(
                    Analysis_Client *aclient,
                    gint sequence_id, gint len){
    register Analysis_Client_Key *key
        = Analysis_Client_Key_create(aclient, sequence_id, len);
    return SparseCache_create(len, Analysis_Client_SparseCache_fill_func,
                         NULL, Analysis_Client_SparseCache_free_func, key);
    }

static Sequence *Analysis_Client_get_Sequence(Analysis_Client *aclient,
                                              gint sequence_id,
                                              gboolean revcomp_target){
    register gchar *msg, *reply, *id, *def;
    register SparseCache *cache;
    register Sequence *seq = aclient->seq_cache[sequence_id];
    register gint len, checksum;
    register gchar **seqinfo_word;
    if(seq){
        if(revcomp_target)
            return Sequence_revcomp(seq);
        return Sequence_share(seq);
        }
    msg = g_strdup_printf("get info %d", sequence_id);
    reply = Analysis_Client_send(aclient, msg, "seqinfo:", FALSE);
    /**/
    if(strncmp(reply, "seqinfo:", 8))
        g_error("Failed to set info for target [%d]", sequence_id);
    /* parse seqinfo for <len> <checksum> <id> and [<def>] */
    seqinfo_word = g_strsplit(reply+9, " ", 4);
    len = atoi(seqinfo_word[0]);
    checksum = atoi(seqinfo_word[1]);
    id = seqinfo_word[2];
    def = seqinfo_word[3];
    /* Strip any trailing newlines */
    if(id[strlen(id)-1] == '\n')
        id[strlen(id)-1] = '\0';
    if(def && (def[strlen(def)-1] == '\n'))
        def[strlen(def)-1] = '\0';
    /**/
    cache = Analysis_Client_get_SparseCache(aclient, sequence_id, len);
    seq = Sequence_create_extmem(id, def, len,
                       (aclient->server_alphabet->type == Alphabet_Type_DNA)
                       ?Sequence_Strand_FORWARD:Sequence_Strand_UNKNOWN,
                       aclient->server_alphabet, cache);
    g_assert(!aclient->seq_cache[sequence_id]);
    aclient->seq_cache[sequence_id] = seq;
    g_strfreev(seqinfo_word);
    SparseCache_destroy(cache);
    g_free(reply);
    g_free(msg);
    if(revcomp_target)
        return Sequence_revcomp(seq);
    return Sequence_share(seq);
    }

typedef enum {
   Analysis_Client_HSP_TOKEN_BEGIN_SET,
   Analysis_Client_HSP_TOKEN_END_SET,
   Analysis_Client_HSP_TOKEN_INT,
   Analysis_Client_HSP_TOKEN_FINISH
} Analysis_Client_HSP_TOKEN;

static Analysis_Client_HSP_TOKEN Analysis_Client_get_hsp_token(gchar *str,
                                                   gint *pos, gint *intval){
    register gint ch;
    gchar *endptr;
    g_assert(str);
    while((ch = str[(*pos)])){
        switch(ch){
            case 'h':
                if(strncmp(str+(*pos), "hspset:", 7))
                    g_error("Unexpected string in HSPset list");
                (*pos) += 7;
                if(strncmp(str+(*pos), " empty\n", 7)){
                    return Analysis_Client_HSP_TOKEN_BEGIN_SET;
                } else {
                    (*pos) += 7;
                    break;
                    }
            case ' ':
                (*pos)++;
                break;
            case '\n':
                (*pos)++;
                return Analysis_Client_HSP_TOKEN_END_SET;
            default:
                if(isdigit(ch)){
                    (*intval) = strtol(str+(*pos), &endptr, 10);
                    (*pos) = endptr-str;
                    return Analysis_Client_HSP_TOKEN_INT;
                } else {
                    g_error("Unexpected character [%d] in HSPset list", ch);
                    }
                break;
            }
        }
    return Analysis_Client_HSP_TOKEN_FINISH;
    }

static void Analysis_Client_get_hsp_sets(Analysis_Client *aclient,
                                         Analysis *analysis,
                                         gboolean swap_chains,
                                         gboolean revcomp_target){
    register gchar *reply = Analysis_Client_send(aclient,
                                                "get hsps", "hspset:", TRUE);
    gint pos = 0, intval = 0;
    register gboolean ok = TRUE;
    register Analysis_Client_HSP_TOKEN token;
    register gint target_id = -1, query_pos = -1, target_pos = -1, length;
    register Comparison *comparison = NULL;
    register Sequence *target = NULL;
    register Match_Type match_type
           = Match_Type_find(aclient->curr_query->alphabet->type,
                             aclient->server_alphabet->type, FALSE);
    /* FIXME: use Match_Type_find with translate_both for codon alignments */
    register HSPset *hsp_set = NULL;
    do {
        token = Analysis_Client_get_hsp_token(reply, &pos, &intval);
        switch(token){
            case Analysis_Client_HSP_TOKEN_BEGIN_SET:
                break;
            case Analysis_Client_HSP_TOKEN_INT:
                if(target_id == -1){
                    target_id = intval;
                    target = Analysis_Client_get_Sequence(aclient,
                                               target_id, revcomp_target);
                    g_assert(!comparison);
                    g_assert(aclient->curr_query);
                    /* FIXME: temp : make work with other HSP types */
                    /* FIXME: should take necessary HSP params from server */
                    if(swap_chains)
                        comparison = Comparison_create(analysis->comparison_param,
                                                   target, aclient->curr_query);
                    else
                        comparison = Comparison_create(analysis->comparison_param,
                                                   aclient->curr_query, target);
                    /* FIXME: should ensure that the HSPset is created
                     * without a horizon
                     */
                    Sequence_destroy(target);
                    target = NULL;
                } else if(query_pos == -1){
                    query_pos = intval;
                } else if(target_pos == -1){
                    target_pos = intval;
                } else {
                    length = intval;
                    g_assert(comparison);
                    /*
                    g_message("adding one [%d,%d,%d]", query_pos, target_pos,
                            length);
                    */
                    /* FIXME: need fix to work with for other match types */
                    switch(match_type){
                        case Match_Type_DNA2DNA:
                            hsp_set = comparison->dna_hspset;
                            break;
                        case Match_Type_PROTEIN2PROTEIN:
                        case Match_Type_PROTEIN2DNA:
                        case Match_Type_DNA2PROTEIN:
                            hsp_set = comparison->protein_hspset;
                            break;
                        default:
                            g_error("Match_Type not supported [%s]",
                                    Match_Type_get_name(match_type));
                        }
                    g_assert(hsp_set);
                    if(swap_chains)
                        HSPset_add_known_hsp(hsp_set, target_pos, query_pos,
                                             length);
                    else
                        HSPset_add_known_hsp(hsp_set, query_pos, target_pos,
                                             length);
                    query_pos = target_pos = -1;
                    }
                break;
            case Analysis_Client_HSP_TOKEN_END_SET:
                /* FIXME: needs to work for other hsp_set types */
                Comparison_finalise(comparison);
                if(Comparison_has_hsps(comparison)){
#if 0
                    /* FIXME: move to use scan_query swap in report */
                    if(swap_chains)
                        Comparison_swap(comparison);
#endif /* 0 */
                    Analysis_report_func(comparison, analysis);
                    }
                Comparison_destroy(comparison);
                comparison = NULL;
                target_id = -1;
                break;
            case Analysis_Client_HSP_TOKEN_FINISH:
                ok = FALSE;
                break;
            }
    } while(ok);
    g_assert(target_id == -1);
    /* format: <HSPSET> <TARGETID> { <QSTART TSTART LEN> } */
    /* tokens <BEGIN_HSPSET> <INT> <ENDHSPSET> */
    g_free(reply);
    return;
    }
/* FIXME: only working for single hspset comparisons */

static void Analysis_Client_process_query(Analysis_Client *aclient,
                                          Analysis *analysis,
                                          Sequence *query,
                                          gboolean swap_chains,
                                          gboolean revcomp_target){
    Analysis_Client_set_query(aclient, query);
    Analysis_Client_get_hsp_sets(aclient, analysis, swap_chains, revcomp_target);
    /* Revcomp query if DNA */
    if(aclient->curr_query->alphabet->type == Alphabet_Type_DNA){
        Analysis_Client_revcomp_query(aclient);
        Analysis_Client_get_hsp_sets(aclient, analysis,
                                     swap_chains, revcomp_target);
        }
    return;
    }

static void Analysis_Client_process(Analysis_Client *aclient, Analysis *analysis,
                                    gboolean swap_chains){
    register FastaDB_Seq *fdbs;
    /* FIXME: need to check for appropriate database type */
    while((fdbs = FastaDB_next(aclient->probe_fdb, FastaDB_Mask_ALL))){
        Analysis_Client_process_query(aclient, analysis, fdbs->seq,
                                      swap_chains, FALSE);
        /* Revcomp target if protein vs DNA or translate_both */
        if(((aclient->curr_query->alphabet->type == Alphabet_Type_PROTEIN)
          && (aclient->server_alphabet->type == Alphabet_Type_DNA))
           || analysis->gam->translate_both){
            Analysis_Client_revcomp_target(aclient);
            Analysis_Client_process_query(aclient, analysis,
                                          fdbs->seq, swap_chains, TRUE);
            Analysis_Client_revcomp_target(aclient);
            }
        FastaDB_Seq_destroy(fdbs);
        }
    return;
    }

/**/

Analysis *Analysis_create(
              GPtrArray *query_path_list, Alphabet_Type query_type,
              gint query_chunk_id, gint query_chunk_total,
              GPtrArray *target_path_list, Alphabet_Type target_type,
              gint target_chunk_id, gint target_chunk_total,
              gint verbosity){
    register Analysis *analysis = g_new0(Analysis, 1);
    register FastaDB *query_fdb = NULL, *target_fdb = NULL,
                     *seeder_query_fdb, *seeder_target_fdb;
    register Match *match;
    Match *dna_match = NULL, *protein_match = NULL,
          *codon_match = NULL;
    register HSP_Param *dna_hsp_param, *protein_hsp_param,
                       *codon_hsp_param;
    register Match_ArgumentSet *mas = Match_ArgumentSet_create(NULL);
    register gboolean use_horizon;
    g_assert(query_path_list);
    g_assert(target_path_list);
    g_assert(query_path_list->len);
    g_assert(target_path_list->len);
    analysis->aas = Analysis_ArgumentSet_create(NULL);

    analysis->verbosity = verbosity;
    /**/
    if(query_path_list->len == 1)
        analysis->query_ac
                = Analysis_Client_create((gchar*)query_path_list->pdata[0],
                                         verbosity);
    if(target_path_list->len == 1)
        analysis->target_ac
                = Analysis_Client_create((gchar*)target_path_list->pdata[0],
                                       verbosity);
    /**/
    if(query_type == Alphabet_Type_UNKNOWN){
        if(analysis->query_ac){
            query_type = analysis->query_ac->server_alphabet->type;
        } else {
            query_type = FastaDB_guess_type(
                              (gchar*)query_path_list->pdata[0]);
            if(verbosity > 1)
                g_message("Guessed query type [%s]",
                        Alphabet_Type_get_name(query_type));
            }
        }
    if(target_type == Alphabet_Type_UNKNOWN){
        if(analysis->target_ac){
            target_type = analysis->target_ac->server_alphabet->type;
        } else {
            target_type = FastaDB_guess_type(
                              (gchar*)target_path_list->pdata[0]);
            if(verbosity > 1)
                g_message("Guessed target type [%s]",
                        Alphabet_Type_get_name(target_type));
            }
        }
    g_assert((query_type == Alphabet_Type_DNA)
           ||(query_type == Alphabet_Type_PROTEIN));
    g_assert((target_type == Alphabet_Type_DNA)
           ||(target_type == Alphabet_Type_PROTEIN));
    if(verbosity > 1)
        g_message("Creating analysis with query[%s] target[%s]",
                Alphabet_Type_get_name(query_type),
                Alphabet_Type_get_name(target_type));
    analysis->gam = GAM_create(query_type, target_type,
                               mas->dna_submat,
                               mas->protein_submat,
                               mas->translate,
                               analysis->aas->use_exhaustive,
                               verbosity);
    /**/
    if(analysis->query_ac)
        Analysis_Client_set_param(analysis->query_ac, analysis->gam);
    if(analysis->target_ac)
        Analysis_Client_set_param(analysis->target_ac, analysis->gam);
    /**/
    Analysis_find_matches(analysis, &dna_match, &protein_match,
                                    &codon_match);
    match = dna_match;
    if(!match)
        match = protein_match;
    if(!match)
        match = codon_match;
    g_assert(match);
    if(!analysis->query_ac)
        query_fdb = FastaDB_open_list_with_limit(query_path_list,
                match->query->alphabet, query_chunk_id, query_chunk_total);
    if(!analysis->target_ac)
       target_fdb = FastaDB_open_list_with_limit(target_path_list,
               match->target->alphabet, target_chunk_id, target_chunk_total);
    if(analysis->aas->use_exhaustive){
        if(analysis->query_ac || analysis->target_ac) /* FIXME: ni */
            g_error("Exhaustive alignment against server not implemented");
        analysis->fasta_pipe = FastaPipe_create(
                                  query_fdb, target_fdb,
                                  Analysis_FastaPipe_Pair_init_func,
                                  Analysis_FastaPipe_Pair_prep_func,
                                  Analysis_FastaPipe_Pair_term_func,
                                  Analysis_FastaPipe_Pair_query_func,
                                  Analysis_FastaPipe_Pair_target_func,
                                  FastaDB_Mask_ALL,
                                  analysis->gam->translate_both);
        analysis->curr_query = NULL;
    } else { /* Not exhaustive */
        use_horizon = (analysis->aas->use_bigseq
                    || analysis->query_ac
                    || analysis->target_ac)?FALSE:TRUE;
        dna_hsp_param = dna_match
                      ? HSP_Param_create(dna_match, use_horizon)
                      : NULL;
        protein_hsp_param = protein_match
                          ? HSP_Param_create(protein_match, use_horizon)
                          : NULL;
        codon_hsp_param = codon_match
                        ? HSP_Param_create(codon_match, use_horizon)
                        : NULL;
        analysis->comparison_param = Comparison_Param_create(
                query_type, target_type,
                dna_hsp_param, protein_hsp_param, codon_hsp_param);
        if(dna_hsp_param)
            HSP_Param_destroy(dna_hsp_param);
        if(protein_hsp_param)
            HSP_Param_destroy(protein_hsp_param);
        if(codon_hsp_param)
            HSP_Param_destroy(codon_hsp_param);
        /* Raise HSP thresholds to score if ungapped */
        if(!Model_Type_is_gapped(analysis->gam->gas->type)){
            if(analysis->comparison_param->dna_hsp_param
            && (analysis->comparison_param->dna_hsp_param->threshold
              < analysis->gam->gas->threshold))
                analysis->comparison_param->dna_hsp_param->threshold
                    = analysis->gam->gas->threshold;
            if(analysis->comparison_param->protein_hsp_param
            && (analysis->comparison_param->protein_hsp_param->threshold
              < analysis->gam->gas->threshold))
                analysis->comparison_param->protein_hsp_param->threshold
                    = analysis->gam->gas->threshold;
            if(analysis->comparison_param->codon_hsp_param
            && (analysis->comparison_param->codon_hsp_param->threshold
              < analysis->gam->gas->threshold))
                analysis->comparison_param->codon_hsp_param->threshold
                    = analysis->gam->gas->threshold;
            }
        /* Don't need HSP horizon for bigseq comparison */
        if(analysis->query_ac || analysis->target_ac){
            if(analysis->query_ac && analysis->target_ac)
                g_error("Server vs server comparison not impelemented");
            analysis->fasta_pipe = NULL;
            if(analysis->query_ac){
                Analysis_Client_set_probe_fdb(analysis->query_ac, target_fdb);
            } else {
                g_assert(analysis->target_ac);
                Analysis_Client_set_probe_fdb(analysis->target_ac, query_fdb);
                }
        } else {
            if(analysis->aas->use_bigseq){
                analysis->bsam = BSAM_create(analysis->comparison_param,
                        analysis->aas->saturate_threshold,
                        verbosity);
                analysis->fasta_pipe = FastaPipe_create(
                                      query_fdb, target_fdb,
                                      Analysis_FastaPipe_Pair_init_func,
                                      Analysis_FastaPipe_Pair_prep_func,
                                      Analysis_FastaPipe_Pair_term_func,
                                      Analysis_FastaPipe_Pair_query_func,
                                      Analysis_FastaPipe_Pair_target_func,
                                      FastaDB_Mask_ALL,
                                      analysis->gam->translate_both);
                analysis->curr_query = NULL;
            } else { /* Use Seeder */
                analysis->scan_query = Analysis_decide_scan_query(query_fdb,
                                                        target_fdb,
                                                 analysis->aas->force_scan);
                if(verbosity > 1)
                    g_message("Applying FSM scan to [%s]",
                              analysis->scan_query?"query":"target");
                /* Swap paths and types
                 * for query and target when scan on query
                 */
                if(analysis->scan_query){
                    seeder_query_fdb = target_fdb;
                    seeder_target_fdb = query_fdb;
                } else {
                    seeder_query_fdb = query_fdb;
                    seeder_target_fdb = target_fdb;
                    }
                analysis->curr_seeder = NULL;
                analysis->fasta_pipe = FastaPipe_create(
                    seeder_query_fdb, seeder_target_fdb,
                    Analysis_FastaPipe_Seeder_init_func,
                    Analysis_FastaPipe_Seeder_prep_func,
                    Analysis_FastaPipe_Seeder_term_func,
                    Analysis_FastaPipe_Seeder_query_func,
                    Analysis_FastaPipe_Seeder_target_func,
                    FastaDB_Mask_ALL, analysis->gam->translate_both);
                }
            }
        }
    if(query_fdb)
        FastaDB_close(query_fdb);
    if(target_fdb)
        FastaDB_close(target_fdb);
    /**/
    return analysis;
    }

void Analysis_destroy(Analysis *analysis){
    if(analysis->fasta_pipe)
        FastaPipe_destroy(analysis->fasta_pipe);
    if(analysis->curr_query)
        FastaDB_Seq_destroy(analysis->curr_query);
    if(analysis->curr_seeder)
        Seeder_destroy(analysis->curr_seeder);
    if(analysis->bsam)
        BSAM_destroy(analysis->bsam);
    if(analysis->comparison_param)
        Comparison_Param_destroy(analysis->comparison_param);
    if(analysis->query_ac)
        Analysis_Client_destroy(analysis->query_ac);
    if(analysis->target_ac)
        Analysis_Client_destroy(analysis->target_ac);
    GAM_destroy(analysis->gam);
    g_free(analysis);
    return;
    }

void Analysis_process(Analysis *analysis){
    if(analysis->query_ac){
        Analysis_Client_process(analysis->query_ac, analysis, TRUE);
    } else if(analysis->target_ac){
        Analysis_Client_process(analysis->target_ac, analysis, FALSE);
    } else {
        while(FastaPipe_process(analysis->fasta_pipe, analysis));
        }
    GAM_report(analysis->gam);
    return;
    }

/**/

