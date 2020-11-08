
#include <stdio.h>
// #include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <stdint.h>
#include <limits.h>
#include <map>
#include <set>
#include <utility>
#include <vector>

#include "htslib/sam.h"
#include "htslib/bgzf.h"


enum test_op {
    READ_COMPRESSED    = 1,
    WRITE_COMPRESSED   = 2,
    READ_CRAM          = 4,
    WRITE_CRAM         = 8,
    WRITE_UNCOMPRESSED = 16,
};


void insert_or_increment(std::map<int32_t, int32_t> & pos_map, int32_t rpos) {
    auto it = pos_map.find(rpos);
    if (it != pos_map.end()) {
        ++(it->second);
    }
    else {
        pos_map.insert({rpos, 1});
    }
}


void help() {
	fprintf(stderr, "Usage: bamsifter [-c max_coverage] [-i max_identical_cigar_pos] [-o out.bam] [--FLAGS] <in.bam>\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "-c: Max coverage value.\n");
    fprintf(stderr, "-o: Output file name. Default is to stdout.\n");
    fprintf(stderr, "-i: Max number of reads with an identical cigar starting at the some position to keep.\n");
    fprintf(stderr, "--keep_unmapped: keep unmapped reads (0x4 flag). \n");
    fprintf(stderr, "--keep_secondary: keep alignments flagged as secondary (0x100 flag).\n");
    fprintf(stderr, "--keep_supplementary: keep alignments flagged as supplementary (0x800 flag).\n");
    fprintf(stderr, "File to process.\n");
    fprintf(stderr, "\n");
}

int main(int argc, char *argv[])
{

    samFile *in;  // open input alignment file
    int flag = 0;
    int clevel = -1;  // compression level
    bam_hdr_t *input_header; // alignment header
    htsFile *out;
    char modew[800];
    int exit_code = 0;
    int coverage_limit = 100;
    int similar_cigar_limit = INT_MAX;
    int keep_unmapped = 0;
    int keep_supplementary = 0;
    int keep_secondary = 0;
    const char *out_name = "-";

    int c;  // for parsing input arguments

    // add option to keep only "proper pairs" (both reads mapped)
    // add option to run strand specificlly 
    // add option for BAM index loading (+ generation if needed) to be able to run each chromosome on a seperate thread

    // while ((c = getopt(argc, argv, "DSIt:i:bCul:o:N:BZ:@:M")) >= 0) {
    while (1) {

        static struct option long_options[] = {
            /* These options set a flag. */
            {"keep_unmapped",           no_argument,        &keep_unmapped,         1},
            {"keep_supplementary",      no_argument,        &keep_supplementary,    1},
            {"keep_secondary",          no_argument,        &keep_secondary,        1},

            /* No flags set. */
            {"coverage_limit",          required_argument,                          0, 'c'},
            {"out_name",                required_argument,                          0, 'o'},
            {"similar_cigar_limit",     required_argument,                          0, 'i'},

            {NULL, 0, NULL, 0}
        };

        int option_index = 0;
        c = getopt_long(argc, argv, "c:o:i:", long_options, &option_index);

        if (c == -1) {
            break;
        }

        switch (c)
        {
            case 0:
                if (long_options[option_index].flag != 0)
                    break;
                printf("option %s", long_options[option_index].name);
                if (optarg)
                    printf (" with arg %s", optarg);
                printf ("\n");
                break;
            case 'c': coverage_limit = atoi(optarg); break;
            case 'o': out_name = optarg; break;
            case 'i': similar_cigar_limit = atoi(optarg); break;
            default:
				help();
                exit(2);
        }
    }

    if (argc == optind) {  // missing input file, print help
		help();
		return (1);
    }

    in = sam_open(argv[optind], "r");


    const htsFormat *in_format = hts_get_format(in);

    // Enable multi-threading (only effective when the library was compiled with -DBGZF_MT)
    // bgzf_mt(in, n_threads, 0);

    if (in == NULL) {
        fprintf(stderr, "Error opening \"%s\"\n", argv[optind]);
        return EXIT_FAILURE;
    }

    input_header = sam_hdr_read(in);
    if (input_header == NULL) {
        fprintf(stderr, "Couldn't read header for \"%s\"\n", argv[optind]);
        return EXIT_FAILURE;
    }

    // check that sort order (SO) flag is set
    std::string tmp_text(input_header->text, input_header->l_text);
    size_t start_pos = tmp_text.find("@HD");
    if (start_pos == std::string::npos) {
        fprintf(stderr, "Error, missing @HD field in header.\n");
        return EXIT_FAILURE;
    }
    size_t end_pos = tmp_text.find("\n", start_pos);
    if (end_pos == std::string::npos) {
        end_pos = tmp_text.length();
    }
    tmp_text = tmp_text.substr(start_pos, end_pos - start_pos);
    start_pos = tmp_text.find("SO:");
    if (start_pos == std::string::npos) {
        fprintf(stderr, "Error, missing SO field in @HD header line.\n");
        return EXIT_FAILURE;
    }
    tmp_text.find("coordinate", start_pos);
    if (start_pos == std::string::npos) {
        fprintf(stderr, "Error, file is not coordinate sorted.\n");
        return EXIT_FAILURE;
    }

    strcpy(modew, "w");
    if (clevel >= 0 && clevel <= 9) sprintf(modew + 1, "%d", clevel);
    if (flag & WRITE_CRAM) strcat(modew, "c");
    else if (flag & WRITE_COMPRESSED) strcat(modew, "b");
    else if (flag & WRITE_UNCOMPRESSED) strcat(modew, "bu");

    // out = hts_open(out_name, modew);
    out = hts_open_format(out_name, modew, in_format);
    if (out == NULL) {
        fprintf(stderr, "Error opening standard output\n");
        return EXIT_FAILURE;
    }

    if (sam_hdr_write(out, input_header) < 0) {  // write header from input to output
        fprintf(stderr, "Error writing output header.\n");
        exit_code = 1;
    }

    // map for start/end positions of alignments that are already selected
    std::map<int32_t, int32_t> starts;
    std::map<int32_t, int32_t> ends;

    // keep a set of mate reads we decided to keep when encountering the first read
    std::set<std::string> mates_to_keep;
    std::map<std::vector<uint32_t>, int> current_cigar_counts;

    bam1_t *aln = bam_init1(); //initialize an alignment
    int32_t current_rname_index = 0; // index compared to header: input_header->target_name[current_rname_index]
    int32_t current_coverage = 0;

    int32_t current_pos = 0;

    while(sam_read1(in, input_header, aln) > 0) {

        if (current_rname_index != aln->core.tid) {
            // should have finished writing reads from current_rname_index contig, so can just reset vars
            current_coverage = 0;
            starts.clear();
            ends.clear();
            // mates_to_keep[current_rname_index].clear();
            fprintf(stdout, "Done with chr %s.\n", input_header->target_name[current_rname_index]);

            current_rname_index = aln->core.tid;
        }

        // make sure the read is mapped
        if (!keep_unmapped && ((aln->core.flag & BAM_FUNMAP) != 0))
            continue;

        // make sure the alignment is not a secondary alignment or a supplementary alignment
        if ((!keep_secondary && (aln->core.flag & BAM_FSECONDARY) != 0) || (!keep_supplementary && (aln->core.flag & BAM_FSUPPLEMENTARY) != 0))
            continue;

        if (current_pos != aln->core.pos) { // left most position, does NOT need adjustment for reverse strand if summing their coverage

            // add the range we want and then erase all the matching entries at once rather than 1 by 1
            auto it = starts.begin();
            if (it->first <= aln->core.pos) {
                for (; it != starts.end(); ++it) {
                    if (it->first <= aln->core.pos) { // or equal because already selected reads take priority in coverage
                        current_coverage += it->second;
                    }
                    else break;
                }
                starts.erase(starts.begin(), it);
            }

            it = ends.begin();
            if (it->first <= aln->core.pos) {
                for (; it != ends.end(); ++it) {
                    if (it->first <= aln->core.pos) { // or equal because already selected reads take priority in coverage
                        current_coverage -= it->second;
                    }
                    else break;
                }
                ends.erase(ends.begin(), it);
            }

            current_cigar_counts.clear();
            current_pos = aln->core.pos;
        }

        // if we are below the max coverage or the read has already been selected to keep through its pair
        if ((current_coverage < coverage_limit) || 
            (mates_to_keep.find(bam_get_qname(aln)) != mates_to_keep.end())) {
            // get cigar
            uint32_t *cigar = bam_get_cigar(aln);

            if (similar_cigar_limit < coverage_limit) {
                std::vector<uint32_t> tmp_cigar(aln->core.n_cigar);
                std::copy(cigar, cigar + aln->core.n_cigar, tmp_cigar.begin());
                auto it = current_cigar_counts.find(tmp_cigar);
                if (it != current_cigar_counts.end()) {  // found
                    if (it->second < similar_cigar_limit) {
                        it->second++;
                    }
                    else {
                        continue;
                    }
                }
                else {
                    current_cigar_counts.emplace(std::move(tmp_cigar), 1);
                }
            }

            int32_t rpos = aln->core.pos;  // update position on the ref with cigar
            for (uint32_t k = 0; k < aln->core.n_cigar; ++k) {

                if ((bam_cigar_type(bam_cigar_op(cigar[k]))&2)) {  // consumes reference
                    if (bam_cigar_op(cigar[k]) == BAM_CREF_SKIP) {
                        insert_or_increment(ends, rpos);

                        rpos += bam_cigar_oplen(cigar[k]);

                        insert_or_increment(starts, rpos);

                    }
                    else {
                        rpos += bam_cigar_oplen(cigar[k]);
                    }
                }
            }

            insert_or_increment(ends, rpos);
            ++current_coverage;

            // save pair mate in their target reference id if not already passed (and cleared)
            if (aln->core.mtid >= current_rname_index) {
                mates_to_keep.insert(bam_get_qname(aln));
            }
            
            // output the alignment
            if (sam_write1(out, input_header, aln) == -1) {
                fprintf(stderr, "Could not write selected record \"%s\"\n", bam_get_qname(aln));
                return EXIT_FAILURE;
            }
        }
    }

    int ret;
    
    ret = hts_close(in);
    if (ret < 0) {
        fprintf(stderr, "Error closing input.\n");
        exit_code = EXIT_FAILURE;
    }

    ret = hts_close(out);
    if (ret < 0) {
        fprintf(stderr, "Error closing output.\n");
        exit_code = EXIT_FAILURE;
    }

    return exit_code;
}
