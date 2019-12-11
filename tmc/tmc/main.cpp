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
    int N_vtx=stoi(argv[2]);    //number of vertices in the motif
    int N_event=stoi(argv[3]);  //number of events in the motif
    int d_c=stoi(argv[4]);      //delta C
    int d_w=stoi(argv[5]);      //delta W
    string out_file = "out_" + gname;
    FILE* fp = fopen (out_file.c_str(), "w");
    
    cout << "delta W: " << d_w << endl;
    cout << "delta C: " << d_c << endl;
    
//Read file and create a sorted list of temporal events
    vector<event> events;
    createEvents(gname, events);
    const auto t2 = chrono::steady_clock::now();
    print_time (fp, "Read data time: ", t2 - t1);
    cout << "# of events: " << events.size() << endl;
    
//Enumerate all instances that satisfy the given constrains(delta C, delta W, and motif size)
//    instancemap instances;
//    set<vector<event>> keys;
//    for (size_t i=0; i<events.size(); i++) {
//        countInstance(events[i], instances, keys, N_vtx, N_event, d_c, d_w);
//    }
//    const auto t3 = chrono::steady_clock::now();
//    print_time (fp, "Count instances time: ", t3 - t2);
//    cout << "# of instances: " << instances.size() << endl;

//Classify instances to corresponding type of motif
//    map<string, int> motif_count;
//    for (auto it=instances.begin(); it!=instances.end(); ++it) {
//        if (it->first.size()==N_event && it->second.second.size()==N_vtx) {
//            string motif = encodeMotif(it->first);
//            motif_count[motif] += it->second.first;
//        }
//    }
    
//Directly count motif from events stream
    map<string, int> motif_count;
    vector<key> pre;
    for (size_t i=0; i<events.size(); i++) {
        countMotif(events[i], pre, motif_count, N_vtx, N_event, d_c, d_w);
    }
    
    const auto t4 = chrono::steady_clock::now();
    print_time (fp, "Count motifs time: ", t4 - t2);
//    print_time (fp, "Count motifs time: ", t4 - t3);
    print_time (fp, "End-to-end Time: ", t4 - t1);
    
    for (auto it=motif_count.begin(); it!=motif_count.end(); ++it) {
        fprintf(fp, "%s \t %d \n", (*it).first.c_str(), (*it).second);
    }
    fclose (fp);
    return 0;
}
