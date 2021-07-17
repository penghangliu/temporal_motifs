//
//  main.cpp
//  tmc
//
//  Created by Penghang Liu on 11/24/19.
//  Copyright © 2019 Penghang Liu. All rights reserved.
//

#include <iostream>
#include "tmc.hpp"

int main(int argc, char * argv[]) {
    const auto t1 = chrono::steady_clock::now();
    
    
    string tmp (argv[1]);
    string gname = tmp.substr (tmp.find_last_of("/") + 1);
    string out_file;
    string tmp_s, method;
    int d_c, d_w, N_event, N_vtx;
    bool multi;
    N_vtx=stoi(argv[2]);    //doesn't matter
    N_event=stoi(argv[3]);
    out_file = "out_static_" + gname + "_" + argv[2] + "_" + argv[3];
    multi = false;
    

    FILE* fp = fopen (out_file.c_str(), "w");
    string oofile = out_file + "_motifs";
    FILE* oofp = fopen (oofile.c_str(), "w");
    
    cout << "delta W: " << d_w << endl;
    cout << "delta C: " << d_c << endl;
    
//Read file and create a sorted list of temporal events
    WGraph graph;
    IPC ipc;
    adj_edges AE;
    createGraph(tmp, graph, AE);
    
    TGraph graph_s;
    SGraph g;
    adj_edges BE;
    if (method == "v1") {
        createGraph(tmp_s, g);
    }
    if (method == "v2") {
        createUndirectedGraph(tmp_s, graph_s);
//        createGraph(tmp_s, graph_s, BE);
    }
        
    const auto t2 = chrono::steady_clock::now();
    print_time (fp, "Read data time: ", t2 - t1);
    
    map<string, int> motif_count;
    Graph2motif(graph, AE, N_vtx, N_event, motif_count, oofp);
    
    const auto t3 = chrono::steady_clock::now();
    print_time (fp, "Count motifs time: ", t3 - t2);
    print_time (fp, "End-to-end Time: ", t3 - t1);
    
    for (auto it=motif_count.begin(); it!=motif_count.end(); ++it) {
        fprintf(fp, "%s,%d \n", (*it).first.c_str(), (*it).second);
    }
//    for (auto it=graph.begin(); it!=graph.end(); ++it) {
//        cout << it->first.first << " " << it->first.second << " : ";
//        for (int i=0; i<it->second.size(); i++) {
//            cout << it->second[i] << " ";
//        }
//        cout << endl;
//    }
//    if (method=="v1") {
//        N_vtx=stoi(argv[5]);    //number of vertices in the motif
//        N_event=stoi(argv[6]);  //number of events in the motif
//        //Enumerate all instances that satisfy the given constrains(delta C, delta W, and motif size)
//        instancemap instances;
//        set<vector<event>> keys;
//        for (size_t i=0; i<events.size(); i++) {
//            countInstance(events[i], instances, keys, N_vtx, N_event, d_c, d_w);
//        }
//        const auto t3 = chrono::steady_clock::now();
//        print_time (fp, "Count instances time: ", t3 - t2);
//        cout << "# of instances: " << instances.size() << endl;
//
//        //Classify instances to corresponding type of motif
//        map<string, int> motif_count;
//        for (auto it=instances.begin(); it!=instances.end(); ++it) {
//            if (it->first.size()==N_event && it->second.second.size()==N_vtx) {
//                string motif = encodeMotif(it->first);
//                motif_count[motif] += it->second.first;
//            }
//        }
//
//        const auto t4 = chrono::steady_clock::now();
//        print_time (fp, "Count motifs time: ", t4 - t3);
//        print_time (fp, "End-to-end Time: ", t4 - t1);
//
//        for (auto it=motif_count.begin(); it!=motif_count.end(); ++it) {
//            fprintf(fp, "%s \t %d \n", (*it).first.c_str(), (*it).second);
//        }
//        fclose (fp);
//    }
    
//    if (method=="v1") {
//        N_vtx=stoi(argv[5]);    //number of vertices in the motif
//        N_event=stoi(argv[6]);  //number of events in the motif
//        //Directly count motif from events stream
//        map<string, int> motif_count;
//        set<key> pre;
//        for (size_t i=0; i<events.size(); i++) {
//            countMotif(events[i], pre, motif_count, N_vtx, N_event, d_c, d_w);
//        }
//        removeIsomorphic(motif_count);
//        
//        const auto t4 = chrono::steady_clock::now();
//        print_time (fp, "Count motifs time: ", t4 - t2);
//        print_time (fp, "End-to-end Time: ", t4 - t1);
//        
//        for (auto it=motif_count.begin(); it!=motif_count.end(); ++it) {
//            fprintf(fp, "%s \t %d \n", (*it).first.c_str(), (*it).second);
//        }
//        fclose (fp);
//    }
    
//    if (method=="v2") {
//        string motif (argv[5]);
//        string edges;
//        for (int i=0; i<motif.size(); i++) {
//            if (motif[i]>='0'&&motif[i]<='9') {
//                edges.push_back(motif[i]);
//            }
//        }
//        int l = edges.length();
//        if (l%2!=0) {
//            printf ("Invalid input, input example: 011202\n");
//            exit(1);
//        }
//        N_event = l/2;
//        N_vtx = 0;
//        for (int i=0; i<edges.length(); i++) {
//            int a = edges[i] - '0';
//            N_vtx = max(a, N_vtx);
//        }
//        N_vtx++;
//        int motif_count;
//        set<key> pre;
//        for (size_t i=0; i<events.size(); i++) {
//            countSpecificmotif (events[i], pre, motif_count, motif, N_vtx, N_event, d_c, d_w);
//        }
//        const auto t3 = chrono::steady_clock::now();
//        print_time (fp, "Count motif time: ", t3 - t2);
//        print_time (fp, "End-to-end Time: ", t3 - t1);
//        cout << "result: " << motif << "\t" << motif_count << endl;
//        fprintf(fp, "%s \t %d \n", motif.c_str(), motif_count);
//    }
    fclose (fp);
    return 0;
}
