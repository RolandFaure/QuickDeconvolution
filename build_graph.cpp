#include "build_graph.h"
#include "cluster_graph.h"

#include <random> //for the iota function
#include <map>

using namespace std::chrono;

using std::cout;
using std::endl;
using std::vector;
using std::list;
using std::min;
using robin_hood::unordered_map;
using robin_hood::unordered_flat_map;
using std::map;
using std::string;
using std::array;
using std::unordered_set;
using std::set;

void thread_deconvolve(short minCommonKmers, unordered_map <string, long int> &tagIDs, const vector <vector<long long int>> &readClouds, std::vector <Read> &reads, const vector<vector<vector<long int>>> &kmers, int thread_id, int num_thread, string folderOut){

    int count = 0;
    for (robin_hood::pair<string, long int> p : tagIDs){

        //a condition because we want each thread to work separately
        if (p.second % num_thread == thread_id){
            vector <int> clusters (readClouds[p.second].size(), -1);
            build_graph(minCommonKmers, p.first, p.second, readClouds, reads, kmers, clusters, folderOut);
            //cout << "Deconvolved " << p.second << endl;

            if (count % 100 == 0) cout << "thread " << thread_id << " deconvolved " << count << " tags over " << readClouds.size() << " in total" << endl;
            count ++;
         }

    }
}

//the function takes as an input the list of all reads having the same tag

void build_graph(short minCommonKmers, string tag, long int tagCloud, const vector <vector<long long int>> &readClouds, std::vector <Read> &reads, const vector<vector<vector<long int>>> &kmers, vector<int> &clusters, string folderOut){
	
	auto t0 = high_resolution_clock::now();

    int adjMatrixSize = clusters.size();
    vector<int> zeros (adjMatrixSize, 0);
    vector<vector<int>> adjMatrix (adjMatrixSize, zeros);
		
    unordered_map<long int, unordered_set<int>> matching_tags; //maps to each tagID the index of the reads of the cloud aligning with this tag
    build_adj_matrix(minCommonKmers, tagCloud, readClouds, reads, kmers, adjMatrix, matching_tags);

    auto t1 = high_resolution_clock::now();

    cluster_graph_chinese_whispers(adjMatrix, clusters, tag);
    //now store the result in the reads
    for (int r = 0 ; r < readClouds[tagCloud].size() ; r++ ){

        reads[readClouds[tagCloud][r]].barcode_extension = clusters[r];

    }

    auto t2 = high_resolution_clock::now();

 //   fast_clustering(tagCloud, readClouds, reads, kmers, clusters);

    auto t3 = high_resolution_clock::now();

//    cout << "Building adjacency matrix : " << duration_cast<microseconds>(t1-t0).count()/1000 << "ms, clustering the matrix : " << duration_cast<microseconds>(t2-t1).count()/1000 << "ms, fast clustering : " << duration_cast<microseconds>(t3-t2).count()/1000 << "ms" << endl;
	

    if (adjMatrix.size() > 3000){

        if (folderOut[folderOut.size()-1] != '/'){
            folderOut += '/';
        }

        string f = folderOut + "cluster_"+tag+"_adj.csv";
        string f2 = folderOut + "cluster_"+tag+"_nodes.csv";
        cout << "exporting..."  << adjMatrix.size()  << " "<< clusters.size()<< endl;

        export_as_CSV(adjMatrix, f, f2, clusters);
        f = folderOut + "cluster_"+tag+"_matching-tag.csv";
        export_as_CSV(matching_tags, f);
    }
}

void build_adj_matrix(short minCommonKmers, long int tagCloud, const vector <vector<long long int>> &readClouds, std::vector <Read> &reads, const vector<vector<vector<long int>>> &kmers, vector<vector<int>> &adjMatrix, unordered_map<long int, unordered_set<int>>& matching_tags){

    auto t0 = high_resolution_clock::now();

    double d = 0;
    double d2 = 0;

    //int r = 0;
    for(int r = 0, sizer = readClouds[tagCloud].size(); r<sizer ; r++){


        unordered_map<long int, int> alreadySeen; //a map to keep track of how many times the read has already been attached to that tag: you need to have at least minCommonKmers common minimizer to attach a read to a tags
        long long int name = readClouds[tagCloud][r];

        vector<vector<long int>> &minis = reads[name].get_minis();

//        if (tagCloud == 24246 && r==368){
//            cout << "Here are the overlappings of read 368 ("<< name << ") on thread 0 and 1 :" << endl;
//             for (long int m : minis[0]){
//                 cout << "[";
//                 for (long int tag : kmers[0][m]){
//                     cout << tag << ",";
//                 }
//                 cout << "] (" << m << ")" << endl;
//             }
//             cout << endl;
//             for (long int m : minis[1]){
//                 cout << "[";
//                 for (long int tag : kmers[1][m]){
//                     cout << tag << ",";
//                 }
//                 cout << "] (" << m << ")" << endl;
//             }
//        }

//        if (tagCloud == 1755 && r==69){
//            cout << "Here are the overlappings of read 69 ("<< name << ") on thread 1 :" << endl;
//             for (long int m : minis[1]){
//                 cout << "[";
//                 for (long int tag : kmers[1][m]){
//                     cout << tag << ",";
//                 }
//                 cout << "] (" << m << ")" << endl;
//             }
//        }

        for (short t = 0 ; t < kmers.size() ; t++){ //here, we iterate through all threads that built the index

            for (long int m : minis[t]){

                for (long int tag : kmers[t][m]){

                        auto tt0 = high_resolution_clock::now();
                        alreadySeen[tag] += 1;

//                        if (tagCloud == 212564 && tag == 10212){

//                            if (r == 59){
//                                cout << "59 got there because of kmer " << t << "," << m << endl;
//                            }
//                            if (r == 61){
//                                cout << "61 got there because of kmer " << t << "," << m << endl;
//                            }
//                            if (r == 63){
//                                cout << "63 got there because of kmer " << t << "," << m << endl;
//                            }
//                        }

                        if (alreadySeen[tag] == minCommonKmers){ // it would be equivalent to put >= here, but a bit slower
                            matching_tags[tag].emplace(r);
                        }
                        auto tt1 = high_resolution_clock::now();
                        d += duration_cast<nanoseconds>(tt1-tt0).count();
                }
            }
        }

    }
    matching_tags[tagCloud] = {}; //this line to avoid self-loops

    auto t1 = high_resolution_clock::now();

    //now build the interaction matrix between reads
    for (robin_hood::pair<const long int, unordered_set<int>> matchs : matching_tags){

//		cout << "\nLet's look at what matched with tag " << matchs.first << endl;
        int i = 0;
        for (int a  : matchs.second) {

            int j = 0;
            for (int b : matchs.second){

//                if (tagCloud == 1755 && a == 10 && b == 42){
//                    cout << "With 42, common tag: " << matchs.first << endl;
//                }
//                if (tagCloud == 1755 && a == 10 && b == 43){
//                    cout << "With 43, common tag: " << matchs.first << endl;
//                }
                if (j>i){
                    adjMatrix[a][b] += 1;
                    adjMatrix[b][a] += 1;
                }
                j++;
            }
            i++;
        }
    }
    auto t2 = high_resolution_clock::now();


    // TagID of tag AGATCTGCAAATTCCG-1 is 1755
//    if (tagCloud == 1755){
//        cout << adjMatrix[10][42] << " " << adjMatrix[10][43] << endl;
//    }

    //cout << "While building adjMat, took me " << duration_cast<microseconds>(t1-t0).count() << "us to create matching tags, among it " << int(d/1000) << "us handling maps (against potentially "<< int(d2/1000) <<"us) " << duration_cast<microseconds>(t2-t1).count() << "us to build the adjMat " << endl;
}

void fast_clustering(long int tagCloud, const std::vector <std::vector<long long int>> &readClouds, std::vector <Read> &reads, const vector<vector<vector<long int>>> &kmers, vector<int> &clusters){

    double time = 0;
    int overlapLimit = 50; //limit of how many tags we look at for each kmer
    vector<int> clusterReps;

    auto t0 = high_resolution_clock::now();
    unordered_map <long int, int> alreadySeenTags; //mapping all the already seen tag to the reps where they were seen

    find_reps(tagCloud, readClouds, reads, kmers, clusterReps, alreadySeenTags);

    auto t1 = high_resolution_clock::now();
    //cout << "Finding the reps took : " << duration_cast<milliseconds>(t1-t0).count() << "ms" << endl;

    //now that we have all the representants of the clusters, map each read to the best cluster

    short newTagsToAdd = 2; //because rep may not always be 100% representative, allow this number of new tags to be indexed at each read to enrich the reference

    //go through all the reads by going through the vector order (we do that instead of a for loop in range(0,clusters.size()) because we'd like to re-cluster the reads we are the less sure about :
    vector<int> order (clusters.size());
    std::iota(order.begin(), order.end(), 0); //fill this vector with increasing int starting from 0

    std::map <std::pair<long int, long int>, int> fusionProposal; //a map of pairs to see if two clusters should be merged
    vector<int> fusionCluster(clusterReps.size()); //a vector indicating which clusters should be merged, for example [0,1,2,0,4] Indicates that cluster 0 and 3 should be merged
    std::iota(fusionCluster.begin(), fusionCluster.end(), 0);

    for(int n = 0 ; n < order.size() ; n++){

        int r = order[n];

        set <long int> newTags;

        vector<int> clusterScores (clusterReps.size());

        long long int name = readClouds[tagCloud][r];

        vector<vector<long int>> &minis = reads[name].get_minis();
        for (short t = 0 ; t < kmers.size() ; t++){
             for (long int m : /*reads[name].*/minis[t]){

                 int index = 0;
                 for (auto tag : kmers[t][m]){

                     if (index < overlapLimit){
                         if (alreadySeenTags.find(tag) != alreadySeenTags.end()){
                             clusterScores[alreadySeenTags[tag]] += 1;
                         }
                         else if (newTags.size() < newTagsToAdd && tag != tagCloud){
                             newTags.emplace(tag);
                         }
                     }
                     index++;
                 }

             }
        }

         //now find the cluster the read seems closest to
         int max = clusterScores[0];
         int idxmax = 0;
         int max2 = -1; //second max
         int idxMax2 = -1;
         for (int i = 1 ; i < clusterScores.size() ; i++){
             if (clusterScores[i]>max){
                 max2 = max;
                 idxMax2 = idxmax;
                 max = clusterScores[i];
                 idxmax = i;
             }
             else if (clusterScores[i] > max2){
                 max2 = clusterScores[i];
                 idxMax2 = i;
             }
         }



         //merge clusters that are too close
         if (clusterScores.size()>1 && max<max2*2 && n<2*clusters.size()){ //if the clustering is not fully convincing
             clusters[r] = idxmax;
             //order.push_back(r);
             fusionProposal[std::make_pair(min(idxmax, idxMax2), std::max(idxmax, idxMax2))] += 1;

             //look if there is enough evidence that idxmax and idxmax2 should in fact be merged in one cluster: if not, then r remains mysterious and you may want to re-inspect it
             int fusionScore = fusionProposal[std::make_pair(min(idxmax, idxMax2), std::max(idxmax, idxMax2))];

             //first check if the fusion has already been decided anyway
             if (fusionCluster[idxmax] != fusionCluster[idxMax2]) {

                 //then check if the fusion should be decided
                 if (fusionScore > 30 && fusionScore > fusionProposal[std::make_pair(min(idxmax, idxmax), std::max(idxmax, idxmax))]+fusionProposal[std::make_pair(min(idxMax2, idxMax2), std::max(idxMax2, idxMax2))] ){

                     //then merge :
                     fusionCluster[std::max(idxmax, idxMax2)] = min(idxmax, idxMax2);
                 }
                 else{ //it is then still unclear what cluster r is from
                    order.push_back(r);
                 }
             }

         }
         else{
             clusters[r] = idxmax;
             fusionProposal[std::make_pair(idxmax,idxmax)] += 1;
             //enrich the reference with new tags
             for (long int nt : newTags){
                 alreadySeenTags[nt] = idxmax;
            }

         }

    }

    for (auto potential : fusionProposal){
            if (potential.first.first != potential.first.second){
                    if (potential.second > fusionProposal[std::make_pair(potential.first.first, potential.first.first)]+fusionProposal[std::make_pair(potential.first.second, potential.first.second)]){
                            fusionCluster[potential.first.second] = potential.first.first;
                    }
            }
    }


    //merge all clusters that should be merged:


    bool cont = true;
    while(cont){

        cont = false;
        for (int c = 0 ; c < clusters.size() ; c++){


            if (fusionCluster[clusters[c]] != clusters[c]){
                cont = true;
                clusters[c] = fusionCluster[clusters[c]];
                //cout << "New value : " << clusters[c] << endl;
            }
        }
    }

    //now store the result in the reads
    for (int r = 0 ; r < readClouds[tagCloud].size() ; r++ ){

        reads[readClouds[tagCloud][r]].barcode_extension = clusters[r];

    }

}

void find_reps(long int tagCloud, const std::vector <std::vector<long long int>> &readClouds, std::vector <Read> &reads, const vector<vector<vector<long int>>> &kmers, vector<int> &clusterReps, unordered_map <long int, int> &alreadySeenTags){

    int limit = 5; //how many tags two reads can share without considering they are linked

    for(int r = 0, sizer = readClouds[tagCloud].size(); r<sizer ; r++){

        short known = 0; //bool keeping trace of with how many already seen barcodes the new read overlaps

        long long int name = readClouds[tagCloud][r];
        vector<vector<long int>> &minis = reads[name].get_minis();

        for (short t = 0 ; t < kmers.size() ; t++){ //different kmers were indexed by different threads
            for (long int m : /*reads[name].*/minis[t]){
                if (known < limit){

                    known --; //that is because tagCloud will always be there as a common tag
                    for (long int tag : kmers[t][m]){

                        if (known < limit){
                            if (alreadySeenTags.find(tag) != alreadySeenTags.end()){
                                known++;
                            }
                        }
                        else{
                            break;
                        }
                    }
                }
                else{
                    break;
                }
            }
        }

        if (known < limit){ //then add the new tags to alreadySeenTags

            clusterReps.push_back(r);
            int clustIdx = clusterReps.size()-1;

            for (short t = 0 ; t < kmers.size() ; t++){
                for (long int m : /*reads[name].*/minis[t]){
                    for (long int tag : kmers[t][m]){
                          alreadySeenTags[tag] = clustIdx;
                    }
                }
            }
        }
    }
//    cout << "Fast clustering, the reps are : ";
//    for (int i : clusterReps){
//        cout << i << " ";
//    }
//    cout << endl;
}
