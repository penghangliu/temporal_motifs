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
	string gname = tmp.substr (tmp.find_last_of("/") + 1, tmp.find_last_of(".") - tmp.find_last_of("/") -1);
    string target (argv[2]);
    string tname = target.substr (target.find_last_of("/") + 1, target.find_last_of(".") - target.find_last_of("/") -1);
	int d_c=stoi(argv[3]);      //delta C
	int N_vtx=stoi(argv[4]);    //number of vertices in the motif
    int N_event=stoi(argv[5]);  //number of events in the motif
	string out_file = "out_" + gname + "_" + tname + "_" + argv[3] + "_" + argv[4] + "_" + argv[5];
	FILE* fp = fopen (out_file.c_str(), "w");

	cout << "delta C: " << d_c << endl;

	//Read file and create a sorted list of temporal events
	vector<event> events;
    vector<event> target_events;
	createEvents(tmp, events);
    createEvents(target, target_events);
	const auto t2 = chrono::steady_clock::now();
	print_time (fp, "Read data time: ", t2 - t1);
	cout << "# of events: " << events.size() << endl;
    cout << "# of target events: " << target_events.size() << endl;
//
//	if (method=="v1") {
//		N_vtx=stoi(argv[5]);    //number of vertices in the motif
//		N_event=stoi(argv[6]);  //number of events in the motif
//        string consecutive (argv[7]);
//		//Enumerate all instances that satisfy the given constrains(delta C, delta W, and motif size)
//		instancemap instances;
//		set<vector<event>> keys;
//		for (size_t i=0; i<events.size(); i++) {
//			countInstance(events[i], instances, keys, N_vtx, N_event, d_c, d_w, consecutive);
//		}
//
//		const auto t3 = chrono::steady_clock::now();
//		print_time (fp, "Count instances time: ", t3 - t2);
//		cout << "# of instances: " << instances.size() << endl;
//
//		//Classify instances to corresponding type of motif
//		map<string, int> motif_count;
//        string out_file1 = out_file + "_list";
//        FILE* fp1 = fopen(out_file1.c_str(), "w");
//		for (auto it=instances.begin(); it!=instances.end(); ++it) {
//			if (it->first.size()==N_event && it->second.second.size()==N_vtx) {
//				string motif = encodeMotif(it->first);
//				vector<event> motif_events = it->first;
//
//				// Only consider delta_C
//				if (N_event == 3) {
//					int second_event_percentile = 100 * ((double) (motif_events[1].first - motif_events[0].first) / (motif_events[2].first - motif_events[0].first));
//					fprintf (fp1, "%s %d %d %d perc dc: %d\n", motif.c_str(), 0, second_event_percentile, 100, d_c);
//					fprintf (fp1, "%s %d span dc: %d\n", motif.c_str(), (motif_events[2].first - motif_events[0].first), d_c);
//				}
//				if (N_event == 4) {
//					int second_event_percentile = 100 * ((double) (motif_events[1].first - motif_events[0].first) / (motif_events[3].first - motif_events[0].first));
//					int third_event_percentile = 100 * ((double) (motif_events[2].first - motif_events[0].first) / (motif_events[3].first - motif_events[0].first));
//					fprintf (fp1, "%s %d %d %d %d perc dc: %d\n", motif.c_str(), 0, second_event_percentile, third_event_percentile, 100, d_c);
//					fprintf (fp1, "%s %d span dc: %d\n", motif.c_str(), (motif_events[3].first - motif_events[0].first), d_c);
//				}
//				motif_count[motif] += it->second.first;
//			}
//		}
//
//		const auto t4 = chrono::steady_clock::now();
//		print_time (fp, "Count motifs time: ", t4 - t3);
//		print_time (fp, "End-to-end Time: ", t4 - t1);
//
//		for (auto it=motif_count.begin(); it!=motif_count.end(); ++it) {
//			fprintf(fp, "%s \t %d \n", (*it).first.c_str(), (*it).second);
//		}
//		fclose (fp);
//	}
    
        //Directly count motif from events stream
        map<event, map<string, int>> motif_count;
        set<pair<event, key>> pre;
        for (size_t i=0; i<events.size(); i++) {
            if(i%100==0) cout << i << " of " << events.size() << endl;
            bool is_target {false};
            if(find(target_events.begin(),target_events.end(),events[i])!=target_events.end()) is_target = true;
            countMotif(events[i], pre, motif_count, N_vtx, N_event, d_c, is_target);
        }
        // removeIsomorphic(motif_count);

        const auto t4 = chrono::steady_clock::now();
        print_time (fp, "Count motifs time: ", t4 - t2);
        print_time (fp, "End-to-end Time: ", t4 - t1);

        for (auto it=motif_count.begin(); it!=motif_count.end(); ++it) {
            for (auto itt=it->second.begin(); itt!=it->second.end(); ++itt)
            fprintf(fp, "%s,%s,%d,%s,%d\n", (*it).first.second.first.c_str(), (*it).first.second.second.c_str(), (*it).first.first, (*itt).first.c_str(), (*itt).second);
        }
        fclose (fp);

	return 0;
}
