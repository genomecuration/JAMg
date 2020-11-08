#include <stdlib.h>
#include <map>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iterator>
#include <math.h>
#include <algorithm>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_max_threads() 1
#define omp_get_num_threads() 1
#define omp_get_thread_num() 0
#endif

#include "IRKE.hpp"
#include "Fasta_reader.hpp"
#include "sequenceUtil.hpp"
#include "KmerCounter.hpp"
#include "stacktrace.hpp"
#include "irke_common.hpp"
#include "argProcessor.hpp"

unsigned int IRKE_COMMON::MONITOR = 0;
size_t KMER_SIZE = 25;
int MAX_THREADS = 6;

// various devel params
bool IRKE_COMMON::__DEVEL_no_kmer_sort = false;

void populate_kmer_counter_from_kmers (KmerCounter& kcounter, string& kmers_fasta_file);
void populate_kmer_counter_from_reads (KmerCounter& kcounter, string& reads_fasta_file);

vector<unsigned int> compute_kmer_coverage (string& sequence, KmerCounter& kcounter);
unsigned int median_coverage (vector<unsigned int>);
string usage(ArgProcessor args);
long sum(vector<unsigned int>& vals);
float mean(vector<unsigned int>& vals);
float stDev(vector<unsigned int>& vals);
int runMe(int argc, char* argv[]);

int main(int argc, char* argv[]) {

    try {
        return(runMe(argc, argv));
    }
    catch (string errmsg) {
        
        cerr << "ERROR encountered: " << errmsg;
        exit(1);
        
    }

    return(0);
}

int runMe(int argc, char* argv[]) {

    ArgProcessor args(argc, argv);
    if(args.isArgSet("--help") ||
       (!(args.isArgSet("--reads")
          &&
          ( args.isArgSet("--kmers") || args.isArgSet("--kmers_from_reads") )
          )) ) {
        cerr << usage(args) << endl << endl;
        exit(1);
    }

    string reads_fasta_file = args.getStringVal("--reads");
    
    bool is_DS = (! args.isArgSet("--SS"));
    if(args.isArgSet("--kmer_size")) {
        KMER_SIZE = args.getIntVal("--kmer_size");
        if(KMER_SIZE < 20) {
            cerr << "Error, min kmer size is 20";
            exit(2);
        }
    }
    if(args.isArgSet("--monitor")) {
        IRKE_COMMON::MONITOR = args.getIntVal("--monitor");
    }
    if (args.isArgSet("--num_threads")) {
        int num_threads = args.getIntVal("--num_threads");
        if (num_threads < MAX_THREADS) {
            omp_set_num_threads(num_threads);
        }
        else {
            // set to max
            omp_set_num_threads(MAX_THREADS);
        }
    }
    
    if(omp_get_max_threads() > MAX_THREADS) {
        omp_set_num_threads(MAX_THREADS);
    }
    KmerCounter kcounter (KMER_SIZE, is_DS);

    if (args.isArgSet("--kmers")) {
        string kmers_fasta_file = args.getStringVal("--kmers");
        populate_kmer_counter_from_kmers(kcounter, kmers_fasta_file);
    }
    else {
        string kmer_read_source_fasta_file = args.getStringVal("--kmers_from_reads");
        populate_kmer_counter_from_reads(kcounter, kmer_read_source_fasta_file);
    }
    
    Fasta_reader fasta_reader(reads_fasta_file);
    bool write_coverage_info = args.isArgSet("--capture_coverage_info");
    
    int start_time = time(NULL);

    cout << "acc" << "\t"
         << "median_cov" << "\t"
         << "mean_cov" << "\t"
         << "stdev" << "\t"
         << "tid" << endl;
    
    #pragma omp parallel
    while (true) {

        if (! fasta_reader.hasNext())
            break;
        
        int myTid = omp_get_thread_num();
        
        Fasta_entry fe = fasta_reader.getNext();
        string sequence = fe.get_sequence();
        if(sequence == "")
            continue;

        string header = fe.get_header();
        vector<unsigned int> kmer_coverage = compute_kmer_coverage(sequence, kcounter);
        unsigned int median_cov = median_coverage(kmer_coverage);
        float mean_cov = mean(kmer_coverage);
        float stdev = stDev(kmer_coverage);
        stringstream stats_text;
                
        stats_text << fe.get_accession() << "\t"
                   << median_cov << "\t"
                   << mean_cov << "\t"
                   << stdev;
        

        stats_text << "\tthread:" << myTid;
        
        if(write_coverage_info) {
            // add the coverage info
            stats_text << "\t";
            for (size_t i = 0; i < kmer_coverage.size(); i++) {
                stats_text<< kmer_coverage[i];
                if(i != kmer_coverage.size() - 1) {
                    stats_text<< ",";
                }
            }
        }
        stats_text << endl;
        
        #pragma omp critical 
        {
            cout << stats_text.str();
        }
        
        if (mean_cov < 0) {
            cerr << "ERROR, cannot have negative coverage!!" << endl;
            exit(1);
        }
        
    }

    int end_time = time(NULL);

    cerr << "STATS_GENERATION_TIME: " << (end_time - start_time) << " seconds." << endl;
    
    return(0);
}

void populate_kmer_counter_from_kmers(KmerCounter& kcounter, string& kmers_fasta_file) {
    // code largely copied from IRKE.cpp
    int i, myTid;
    unsigned long sum,
        *record_counter = new unsigned long[omp_get_max_threads()];
    unsigned long start, end;
    // init record counter
    for (int i = 0; i < omp_get_max_threads(); i++) {
        record_counter[i] = 0;
    }
    cerr << "-reading Kmer occurrences..." << endl;
    start = time(NULL);
    Fasta_reader fasta_reader(kmers_fasta_file);
    #pragma omp parallel private (myTid)
    {
        myTid = omp_get_thread_num();
        record_counter[myTid] = 0;
        while (true) {
            Fasta_entry fe = fasta_reader.getNext();
            if(fe.get_sequence() == "") break;
            record_counter[myTid]++;
            if(IRKE_COMMON::MONITOR) {
                if(myTid == 0 && record_counter[myTid] % 100000 == 0)
                    {
                        sum = record_counter[0];
                        for (i=1; i<omp_get_num_threads(); i++)
                            sum+= record_counter[i];
                        cerr << "\r [" << sum/1000000 << "M] Kmers parsed.     ";
                    }
            }
            string seq = fe.get_sequence();
            if(seq.length() != KMER_SIZE) {
                cerr << "ERROR: kmer " << seq << " is not of length: " << KMER_SIZE << endl;
                continue;
            }
            kmer_int_type_t kmer = kcounter.get_kmer_intval(seq);
            unsigned int count = atoi(fe.get_header().c_str());
            kcounter.add_kmer(kmer, count);
        }
    }
    end = time(NULL);
    sum = record_counter[0];
    for (i=1; i<omp_get_max_threads(); i++)
        sum+= record_counter[i];
    delete [] record_counter;
    cerr << endl << " done parsing " << sum << " Kmers, " << kcounter.size() << " added, taking " << (end-start) << " seconds." << endl;
    return;
}

void populate_kmer_counter_from_reads (KmerCounter& kcounter, string& fasta_filename) {
    unsigned int kmer_length = kcounter.get_kmer_length();
    int i, myTid;
    unsigned long sum,
        *record_counter = new unsigned long[omp_get_max_threads()];
    unsigned long start, end;

    // init record counter
    for (int i = 0; i < omp_get_max_threads(); i++) {
        record_counter[i] = 0;
    }


    cerr << "-storing Kmers..." << "\n";
    start = time(NULL);

    Fasta_reader fasta_reader(fasta_filename);

    unsigned int entry_num = 0;

#pragma omp parallel private (myTid)
    {
        myTid = omp_get_thread_num();
        record_counter[myTid] = 0;

        while (fasta_reader.hasNext()) {
            Fasta_entry fe = fasta_reader.getNext();
            string accession = fe.get_accession();

#pragma omp atomic
            entry_num++;
            record_counter[myTid]++;
            
            if (IRKE_COMMON::MONITOR >= 4) {
                cerr << "[" << entry_num << "] acc: " << accession << ", by thread no: " << myTid << "\n";;
            }
            else if (IRKE_COMMON::MONITOR) {
                if (myTid == 0 && record_counter[myTid] % 1000 == 0)
                    {
                        sum = record_counter[0];
                        for (i=1; i<omp_get_num_threads(); i++)
                            sum+= record_counter[i];
                        cerr << "\r [" << sum << "] sequences parsed.     ";
                    }
            }
            
            string seq = fe.get_sequence();

            if (seq.length() < KMER_SIZE + 1) {
                continue;
            }
            kcounter.add_sequence(seq);

        }
        
        cerr << "\n" << " done parsing " << sum << " sequences, extracted " << kcounter.size() << " kmers, taking " << (end-start) << " seconds." << "\n";
        
        
    }


    return;
    
}






vector<unsigned int> compute_kmer_coverage(string& sequence, KmerCounter& kcounter) {
    if(IRKE_COMMON::MONITOR) {
        cerr << "processing sequence: " << sequence << endl;
    }
    if (sequence.length() < KMER_SIZE)
    {
        // Can't rely on length() - KMER_SIZE for this as length is unsigned
        cerr << "Sequence: " << sequence << "is smaller than " << KMER_SIZE << " base pairs, skipping" << endl;
	return vector<unsigned int>();
    }

    vector<unsigned int> coverage;
    for (size_t i = 0; i <= sequence.length() - KMER_SIZE; i++) {
        // cerr << "i: " << i << ", <= " << sequence.length() - KMER_SIZE << endl;
        string kmer = sequence.substr(i, KMER_SIZE);
        if(IRKE_COMMON::MONITOR >= 2) {
            for (size_t j = 0; j <= i; j++) {
                cerr << " ";
            }
            cerr << kmer << endl;
        }
        unsigned int kmer_count = 0;
        if(!contains_non_gatc(kmer)) {
            kmer_count = kcounter.get_kmer_count(kmer);
        }

        // Note, in the jellyfish run, we restrain it to min kmer coverage of 2. 
        // If we don't find a kmer catalogued, it must have a kmer count of 1.
        if (kmer_count < 1) {
            kmer_count = 1;
        }
        
        coverage.push_back(kmer_count);
    }
    return(coverage);
}

unsigned int median_coverage(vector<unsigned int> coverage) {
    int num_pts = coverage.size();
    if(num_pts == 0) {
        return(0);
    }
    sort(coverage.begin(), coverage.end());
    if(num_pts % 2 == 1){
      return(coverage[(num_pts)/2]);
    }
    return((coverage[(num_pts - 1)/2] + coverage[(num_pts)/2]) / 2);
}


string usage(ArgProcessor) {
    stringstream usage_info;
    usage_info << endl << endl;
    usage_info << "Usage: " << endl
               << "  --reads  <str>             " << ":fasta file containing target reads for kmer coverage stats" << endl
               << endl 
               << " and source of kmers: " << endl
               << "  --kmers  <str>             " << ":fasta file containing kmers" << endl
               << "      or " << endl 
               << "  --kmers_from_reads <str>   " << ":fasta file containing reads as source of kmers" << endl 
               << endl 
               << "* optional:" << endl
               << "  --kmer_size <int>          " << ":default = 25" << endl
               << "  --DS                             " << ":double-stranded RNA-Seq mode (not strand-specific)" << endl
               << "  --capture_coverage_info    " << ":writes coverage info file." << endl
               << "  --monitor <int>            " << ":verbose level for debugging" << endl
               << "  --num_threads <int>        " << ":number of threads" << endl
               << endl;
    return(usage_info.str());
}

long sum(vector<unsigned int>& vals) {
    long sum_vals = 0;
    for (int i = 0; i < (int) vals.size(); i++) {
        sum_vals += vals[i];
    }

    return(sum_vals);
}

float mean(vector<unsigned int>& vals) {
  if(vals.size() == 0){
    return(0);
  }
  long sum_vals = sum(vals);
  float avg = (float) sum_vals / vals.size();
  return(avg);
}

float stDev( vector<unsigned int>& vals) {
    float avg = mean(vals);
    int num_vals = vals.size();

    float sum_avg_diffs_sqr = 0;
    for (int i = 0; i < num_vals; i++) {
        float delta = vals[i] - avg;
        sum_avg_diffs_sqr += (delta * delta);
    }

    float stdev = sqrt( sum_avg_diffs_sqr/(num_vals-1) );

    return(stdev);
}
