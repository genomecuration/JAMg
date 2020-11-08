#include <string>
#include <cstring>
#include <sstream>
#include <algorithm>
#include <map>
#include <fstream>

#include "base/CommandLineParser.h"
#include "analysis/DNAVector.h"
#include "aligns/KmerAlignCore.h"

#include <math.h>
#include "analysis/KmerTable.h"
#include "analysis/NonRedKmerTable.h"
#include "util/mutil.h"
#include <omp.h>
#include "analysis/DeBruijnGraph.h"
#include "analysis/sequenceUtil.h"
#include "analysis/Pool.h"


static bool DEBUG = false;
static int MAX_CLUSTER_SIZE = 25;
static int MIN_CONTIG_LENGTH = 24;


// print nucleotide sequence 80 chars per line
void PrintSeq(const DNAVector & d) {
    int i;
    for (i=0; i<d.isize(); i++) {
        cout << d[i];
        if ((i+1) % 80 == 0)
            cout << "\n";
    }
    cout << "\n";
}


void Add(vecDNAVector & all, DNAVector & add, int & counter) 
{
    #pragma omp critical
    {
        all.push_back(add);
        counter++;
    }
    
}


void describe_poolings (vector<Pool>& pool_vec) {
    
    cerr << "\n" << "Pools described as follows:" << "\n";
    
    for (unsigned int i=0; i<pool_vec.size(); i++) {
        int oldpool_id = i;
        
        stringstream old_pool_info;
        
        Pool& oldpool = pool_vec[oldpool_id];

        cerr << oldpool.str() << endl;
        
    }
    
}

void describe_bubblings (vector<Pool>& pool_vec, map<int,Pool>& pool_idx_to_containment, int round) {
    
    for (unsigned int i=0; i<pool_vec.size(); i++) {

        stringstream old_pool_info;
        
        Pool& oldpool = pool_vec[i];
        
        old_pool_info << "Bubbling[" << round << "], Pool [" << oldpool.get_id() << "] size: " << oldpool.size() << ", links: "; 
        for (int j=0; j<oldpool.size(); j++) {
            int iworm_id = oldpool.get(j);
            
            old_pool_info << iworm_id << " ";
        }
        if (pool_idx_to_containment.find(oldpool.get_id()) != pool_idx_to_containment.end()) {
            old_pool_info << "; Containments: ";
            Pool& containment_pool = pool_idx_to_containment[oldpool.get_id()];
            for (int j = 0; j < containment_pool.size(); j++) {
                old_pool_info << containment_pool.get(j) << " ";
            }
        }
        
        old_pool_info << "\n";
        
        cerr << old_pool_info.str();
    }
    
}


bool sort_pool_sizes_descendingly(const Pool& a, const Pool& b) {
    return (a.size() > b.size());
    
}

bool sort_pool_sizes_ascendingly(const Pool& a, const Pool& b) {

    // cerr << "size comparison: " << a.size() << " vs. " << b.size() << endl;
    return (a.size() < b.size());
}




svec<Pool> grow_prioritized_clusters(string& weld_graph_file, map<int,Pool>& weld_reinforced_iworm_clusters) {

    cerr << "Parsing weld graph file: " << weld_graph_file << endl;
    
    ifstream in (weld_graph_file.c_str());
    if (! in.is_open()) {
        cerr << "Error, cannot open file: " << weld_graph_file << endl;
        exit(3);
    }

    // lines look like this:
	// 716 -> 725 weldmers: 9 scaff_pairs: 69 total: 78
    // 234 -> 717 weldmers: 61 scaff_pairs: 0 total: 61
    // 717 -> 234 weldmers: 61 scaff_pairs: 0 total: 61
    // 457 -> 466 weldmers: 0 scaff_pairs: 58 total: 58
    // 2 -> 459 weldmers: 57 scaff_pairs: 0 total: 57
    // 459 -> 2 weldmers: 57 scaff_pairs: 0 total: 57
    // 2 -> 714 weldmers: 0 scaff_pairs: 43 total: 43
    // 232 -> 493 weldmers: 7 scaff_pairs: 34 total: 41


    svec<Pool> clustered_contigs;
    map<unsigned int,unsigned int> id_to_cluster;
        
    while (! in.eof()) {
        string line;
        getline(in, line);
        //cerr << "input_line: " << line << "\n";
        
        istringstream token (line);

        string tok;
        token >> tok;

        unsigned int node_id = atoi(tok.c_str());
                
        token >> tok; // pointer '->' separator string.

        token >> tok;
        unsigned int linked_node_id = atoi(tok.c_str());

        // check both assigned.
                
        if (id_to_cluster.find(node_id) != id_to_cluster.end()
            &&
            id_to_cluster.find(linked_node_id) != id_to_cluster.end()) {

            // try to merge the clusters
            unsigned int node_id_cluster_id = id_to_cluster[node_id];
            unsigned int linked_node_id_cluster_id = id_to_cluster[linked_node_id];

            if (node_id_cluster_id == linked_node_id_cluster_id) {
                // already in the same cluster. Must have been linked transitively
                continue;
            }
            else if (clustered_contigs[node_id_cluster_id].size() + clustered_contigs[linked_node_id_cluster_id].size() + 2 <= MAX_CLUSTER_SIZE) {
                // different clusters.
                // merge clusters.
                for (int i = 0; i < clustered_contigs[linked_node_id_cluster_id].size(); i++) {
                    unsigned int node_to_migrate_id = clustered_contigs[linked_node_id_cluster_id][i];
                    clustered_contigs[node_id_cluster_id].add(node_to_migrate_id);
                    id_to_cluster[node_to_migrate_id] = node_id_cluster_id;
                }
                clustered_contigs[linked_node_id_cluster_id].clear();
            }
            else {
                // ignoring this linkage between iworm contigs, as it'll over-aggregate the cluster sizes
            }
        }
        else if (id_to_cluster.find(node_id) != id_to_cluster.end()
                 || 
                 id_to_cluster.find(linked_node_id) != id_to_cluster.end()) {

            // one of them has a cluster assignment.
            // add the other
            if (id_to_cluster.find(node_id) != id_to_cluster.end()) {
                unsigned int node_id_cluster_id = id_to_cluster[node_id];
                if (clustered_contigs[node_id_cluster_id].size() < MAX_CLUSTER_SIZE) {
                    clustered_contigs[node_id_cluster_id].add(linked_node_id);
                    id_to_cluster[linked_node_id] = node_id_cluster_id;
                }
            }
            else {
                // linked_node_id found in cluster.
                unsigned int linked_node_id_cluster_id = id_to_cluster[linked_node_id];
                if (clustered_contigs[linked_node_id_cluster_id].size() < MAX_CLUSTER_SIZE) {
                    clustered_contigs[linked_node_id_cluster_id].add(node_id);
                    id_to_cluster[node_id] = linked_node_id_cluster_id;
                }
            }
        }
        else {
            // neither is in a cluster.
            // make a new cluster, store both entries.
            Pool tmp_pool;
            tmp_pool.add(node_id);
            tmp_pool.add(linked_node_id);
            unsigned int new_cluster_idx = clustered_contigs.size();
            clustered_contigs.push_back(tmp_pool);
            id_to_cluster[node_id] = new_cluster_idx;
            id_to_cluster[linked_node_id] = new_cluster_idx;
        }
    }


    // remove the emtpy entries
    svec<Pool> ret_clustered_contigs;
    for (size_t i = 0; i < clustered_contigs.size(); i++) {
        if (clustered_contigs[i].size() > 0) {
            ret_clustered_contigs.push_back(clustered_contigs[i]);
        }
    }

    return(ret_clustered_contigs);
    
            
}
    
// Single-linkage clustering of iworm contigs in pools
svec<Pool> sl_cluster_pools(map<int,Pool>& pool, map<int,bool>& ignore) {

    // Just pull out the ordered pools, initial order doesn't matter... they get trickled upward.
    svec<Pool> pool_vec;
    
    for (map<int,Pool>::iterator it = pool.begin(); it != pool.end(); it++) {
        pool_vec.push_back(it->second);
    }
        
    // init entries to loweset pool index they're found in.
    map<int,int> mapped;
    
    for (size_t i = 0; i < pool_vec.size(); i++) {
        
        Pool p = pool_vec[i];
        
        for (int j = 0; j < p.size(); j++) {
            int ele_id = p.get(j);
            
            if (mapped.find(ele_id) == mapped.end()) {
                
                mapped[ele_id] = i; // first time seeing it, so lowest pool index value
                
                // cerr << "Mapping: " << " ele: " << ele_id << " to pool " << i << "\n";
                
            }
            
        }
    }
    
    int round = 0;
    bool remapped_flag = true;
    while (remapped_flag) {
        


        round++;
        cerr << "Transitive closure, round: " << round << "\n";
        remapped_flag = false;
        
        
        if (DEBUG) {
            cerr << "Before remapping: " << "\n";
            
            describe_poolings(pool_vec);
        }
        

        // do transitive closure on the poolings
        map<int,int> remapped = mapped;
        
        // trickle pool mappings upward
        
        for (unsigned int i=0; i<pool_vec.size(); i++) {
            int oldpool_id = i;
            
            Pool& oldpool = pool_vec[oldpool_id];
            
            int lowest_pool_map_id = oldpool_id;
            bool pool_reassignment_required = false;
            for (int j=0; j<oldpool.size(); j++) {
                int iworm_id = oldpool.get(j);
                int pool_mapping = remapped[iworm_id];
                
                
                if (pool_mapping != lowest_pool_map_id) {
                    pool_reassignment_required = true;
                }
                if (pool_mapping >= 0 && pool_mapping < lowest_pool_map_id) {
                    lowest_pool_map_id = pool_mapping;
                }
            }
            
            // reassign all members to lowest_pool_map_id (if same as pool_id, some entries may still be higher even if not less)
            if (pool_reassignment_required) {
                for (int j=0; j<oldpool.size(); j++) {
                    int iworm_id = oldpool.get(j);
                    remapped[iworm_id] = lowest_pool_map_id;
                    // cerr << "-remapping: iworm(" << iworm_id << ") to pool " << lowest_pool_map_id << "\n"; 
                }
            }
            
            if (lowest_pool_map_id < oldpool_id) {
                // reassign to lower pool
                Pool tmp;
                for (int j=0; j<oldpool.size(); j++) {
                    int iworm_id = oldpool.get(j);
                    
                    if (! pool_vec[lowest_pool_map_id].contains(iworm_id)) {

                        pool_vec[lowest_pool_map_id].add(iworm_id);
                    }
                    //#pragma omp critical
                    //cout << "RELOCATING: " << iworm_id << " from old pool " << oldpool_id << " to new pool " << lowest_pool_map_id << "\n";
                }
                //cerr << "-clearing out old pool: " << oldpool_id << "\n";
                pool_vec[oldpool_id] = tmp;
                remapped_flag = true;
            }    
        }
        mapped = remapped; // new assignments
        
        if (DEBUG) {
            cerr << "After remapping: " << "\n";
            describe_poolings(pool_vec);
        }
        

    }
    
         
    
    //cerr << "...done (" << end_time - start_time << " seconds)" << "\n";

    //------------------------------------------
    // Clustered inchworm contigs now defined.
    // Generate final output
    
    svec<Pool> nr_pools;
    
     for (int i=0; i < (int) pool_vec.size(); i++) {
        Pool & p = pool_vec[i];
        
        if (p.size() == 0) { continue; }
        
        p.sortvec();

        Pool nr_entries; // remove the redundant entries in p
        for (int j=0; j<p.size(); j++) {
            int z = p.get(j);
            
            if (ignore[z]) {
                continue;
            }
            

            if (j > 0 && p.get(j-1) == z) {
                // same entry ended up on the pool vec, already reported it.
                continue;
            }
            
            nr_entries.add(z);
        }
        nr_pools.push_back(nr_entries);
    }
    

    return(nr_pools);
    

}


void add_unclustered_iworm_contigs (svec<Pool>& clustered_pools, vecDNAVector& dna) {

    if (DEBUG) {
        cerr << "Adding unclustered iworm contigs." << "\n";
    }

    map<int,bool> found;

    for (svec<Pool>::iterator it = clustered_pools.begin(); it != clustered_pools.end(); it++) {

        Pool& p = *it;

        for (int i = 0; i < p.size(); i++) {
            int id = p[i];
            found[id] = true;
            if (DEBUG)
                cerr << ":found iworm: " << id << " in cluster." << "\n";
        }
    }


    // add in the missing entries
    for (size_t i = 0; i < dna.size(); i++) {
        if (found.find(i) == found.end()) {
            Pool p;
            p.add(i);
            clustered_pools.push_back(p);
            if (DEBUG)
                cerr << "-iworm: " << i << " added as unclustered." << "\n";
        }
    }

}



bool Exists(const string & s) 
{
    FILE * p = fopen(s.c_str(), "r");
    if (p != NULL) {
        fclose(p);
        return true;
    }
    // cout << "FATAL ERROR: Could not open file for read: " << s << "\n";
    // cout << "Please make sure to enter the correct file name(s). Exiting now." << "\n";
    
    return false;
}


void populate_weld_reinforced_iworm_clusters(string& weld_graph_file, map<int,Pool>& weld_reinforced_iworm_clusters) {

    cerr << "Parsing weld graph file: " << weld_graph_file << endl;
    
    ifstream in (weld_graph_file.c_str());
    if (! in.is_open()) {
        cerr << "Error, cannot open file: " << weld_graph_file << endl;
        exit(3);
    }

    /*  old way (pre-June 2018)
    
    // lines look like this:
    //    210 -> 735
    //    907 -> 242 506
    //    735 -> 515 13 807 779 85 92 210 372 397 567 819    
    
    while (! in.eof()) {
        string line;
        getline(in, line);
        //cerr << "input_line: " << line << "\n";
        
        istringstream token (line);

        string tok;
        token >> tok;

        int node_id = atoi(tok.c_str());
        Pool p(node_id);
        
        token >> tok; // pointer '->' separator string.

        
        while (! token.eof()) {
            token >> tok;
            int adjacent_node = atoi(tok.c_str());
            p.add(adjacent_node);
            //cerr << "Tok: [" << tok << "]" << "\n";
        }
        
        if (node_id >= 0 && p.size() > 0) {
            weld_reinforced_iworm_clusters[ node_id ] = p;
            if (DEBUG) {
                cerr << "Assigning node_id [" << node_id << "] to pool: " << p.str() << endl;
            }
        }
    }
    
    */
    

    //------------------------------------------------------
    // new way:

    // lines look like this:
	// 716 -> 725 weldmers: 9 scaff_pairs: 69 total: 78
    // 234 -> 717 weldmers: 61 scaff_pairs: 0 total: 61
    // 717 -> 234 weldmers: 61 scaff_pairs: 0 total: 61
    // 457 -> 466 weldmers: 0 scaff_pairs: 58 total: 58
    // 2 -> 459 weldmers: 57 scaff_pairs: 0 total: 57
    // 459 -> 2 weldmers: 57 scaff_pairs: 0 total: 57
    // 2 -> 714 weldmers: 0 scaff_pairs: 43 total: 43
    // 232 -> 493 weldmers: 7 scaff_pairs: 34 total: 41

        
    while (! in.eof()) {
        string line;
        getline(in, line);
        //cerr << "input_line: " << line << "\n";
        
        istringstream token (line);

        string tok;
        token >> tok;

        int node_id = atoi(tok.c_str());
                
        token >> tok; // pointer '->' separator string.

        token >> tok;
        int linked_node_id = atoi(tok.c_str());

        if (weld_reinforced_iworm_clusters.find(node_id) != weld_reinforced_iworm_clusters.end()) {
            Pool& p = weld_reinforced_iworm_clusters[node_id];
            p.add(linked_node_id);
            if (DEBUG) {
                cerr << "Adding node_id [" << linked_node_id << "] to pool of " << node_id << " yielding " << p.str() << endl;
            }
        }
        else {
            Pool p(node_id);
            p.add(linked_node_id);
            weld_reinforced_iworm_clusters[ node_id ] = p;
            if (DEBUG) {
                cerr << "Initializing node_id [" << node_id << "] to pool: " << p.str() << endl;
            }
        }
        
    }
        
    return;
}



int main(int argc,char** argv)
{
    
    
    commandArg<string> aStringCmmd("-i","input fasta");
    commandArg<string> weldGraphCmmd("-weld_graph", "iworm index weld graph");
    commandArg<bool> debugCmmd("-debug", "verbosely describes operations", false);
    commandArg<int>  minContigLengthCmmd("-min_contig_length", "min sum cluster contig length", MIN_CONTIG_LENGTH);
    commandArg<bool>  debugWeldAllCmmd("-debug_weld_all", "creates a single cluster of all contigs, for debugging only", false);
    commandArg<int> max_cluster_size_Cmmd("-max_cluster_size", "max num if iworms in a single cluster", MAX_CLUSTER_SIZE);
    
    commandLineParser P(argc,argv);
    P.SetDescription("Makes a graph out of a fasta");
    P.registerArg(aStringCmmd);
    P.registerArg(weldGraphCmmd);
    P.registerArg(debugCmmd);
    P.registerArg(minContigLengthCmmd);
    P.registerArg(debugWeldAllCmmd);
    P.registerArg(max_cluster_size_Cmmd);
    
    P.parse();
    
    cerr << "-------------------------------------------------------" << "\n"
         << "--- Chrysalis: BubbleUpClustering ---------------------" << "\n"
         << "-- (generating the final clusters of iworm contigs) ---" << "\n"
         << "-------------------------------------------------------" << "\n" << "\n";
    
    string iworm_contigs_filename = P.GetStringValueFor(aStringCmmd); //inchworm contigs file
    DEBUG = P.GetBoolValueFor(debugCmmd);
    string weld_graph_file = P.GetStringValueFor(weldGraphCmmd);
    MIN_CONTIG_LENGTH = P.GetIntValueFor(minContigLengthCmmd);
    
    bool DEBUG_WELD_ALL = P.GetBoolValueFor(debugWeldAllCmmd);
    MAX_CLUSTER_SIZE = P.GetIntValueFor(max_cluster_size_Cmmd);
    if (! Exists(iworm_contigs_filename)) {
        cerr << "ERROR, cannot open iworm contigs file: " << iworm_contigs_filename << "\n";
        exit(2);
    }
    if (! Exists(weld_graph_file)) {
        cerr << "ERROR, cannot open weld graph file: " << weld_graph_file << "\n";
        exit(3);
    }
    
    // read inchworm contigs into memory
    vecDNAVector dna;
    cerr << "BubbleUpClustering: Reading file: " << iworm_contigs_filename << "\n";
    dna.Read(iworm_contigs_filename, false, false, true, 1000000);
    cerr << "done!" << "\n";
    
        
    map<int,Pool> weld_reinforced_iworm_clusters;
    svec<Pool> clustered_pools;
    
    if (DEBUG_WELD_ALL) {

        // put all entries into a single pool.
        Pool p;
        for (size_t i = 0; i < dna.size(); i++) {
            p.add(i);
        }
        clustered_pools.push_back(p);
    }
    else {

        clustered_pools = grow_prioritized_clusters(weld_graph_file, weld_reinforced_iworm_clusters);
        
        if (DEBUG) {
            cerr << "Final pool description: " << "\n";
            describe_poolings(clustered_pools);
        }
        
        
        add_unclustered_iworm_contigs(clustered_pools, dna);
    }
        
    //-----------------------------------------------------------------------------------
    // Generate final output
    
    
        
    int component_count = 0;
    
    for (int i=0; i<clustered_pools.isize(); i++) {
        Pool & p = clustered_pools[i];
        
        if (p.size() == 0) { continue; }

        p.sortvec();
    
        int sum_iworm_length = 0;
        for (size_t j = 0; j < p.size(); j++) {
            int z = p.get(j);
            sum_iworm_length += dna[z].isize();
        }

        if (sum_iworm_length < MIN_CONTIG_LENGTH) {
            continue;
        }
        
        stringstream pool_info;                
        pool_info << "#POOL_INFO\t" << component_count << ":" << "\t";

        cout << "COMPONENT " << component_count << "\t" << p.size() << "\n";
        for (size_t j = 0; j < p.size(); j++) {
            int z = p.get(j);
            
            pool_info << z << " ";
            cout << ">Component_" << component_count << " " << p.size() << " " << z << " [iworm" << dna.Name(z) << "]" << "\n";
            PrintSeq(dna[z]);
        }
        pool_info << "\n";
        cout << pool_info.str();
        
        cout << "END" << "\n";
    
        component_count++;

    }
    
    return 0;
    
}
  
