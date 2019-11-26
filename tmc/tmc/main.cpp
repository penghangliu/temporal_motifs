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
    // insert code here...
    const auto t1 = chrono::steady_clock::now();
    
    string tmp (argv[1]);
    string gname = tmp.substr (tmp.find_last_of("/") + 1);
    int N_vtx=stoi(argv[2]);
    int N_event=stoi(argv[3]);
    int d_c=stoi(argv[3]);
    int d_w=stoi(argv[4]);
    string out_file = "out_" + gname;
    FILE* fp = fopen (out_file.c_str(), "w");
    
    vector<event> events;
    createEvents(gname, events);
    cout << "# of events: " << events.size() << endl;
    
    instancemap instances;
    set<vector<event>> keys;
    for (size_t i=0; i<events.size(); i++) {
        countInstance(events[i], instances, keys, N_vtx, N_event, d_c, d_w);
    }
    cout << "# of instances: " << instances.size() << endl;

    map<string, int> motif_count;
    for (auto it=instances.begin(); it!=instances.end(); ++it) {
        if (it->first.size()==N_event && it->second.second.size()==N_vtx) {
            string motif = encodeMotif(it->first);
            motif_count[motif] += it->second.first;
        }
    }
    for (auto it=motif_count.begin(); it!=motif_count.end(); ++it) {
        fprintf(fp, "%s \t %d \n", (*it).first.c_str(), (*it).second);
    }
    
    std::cout << "Hello, World!\n";
    const auto t2 = chrono::steady_clock::now();
    print_time (fp, "End-to-end Time: ", t2 - t1);
    
//    for (int i=0; i<events.size(); i++) {
//        fprintf(fp, "%d \t %d \t %d \n", events[i].second.first, events[i].second.second, events[i].first);
//    }
//    unordered_map<string, int> mapOfWords;
//    mapOfWords.insert(make_pair("earth", 1));
//    mapOfWords.insert(make_pair("moon", 2));
//    mapOfWords["sun"] = 3;
//    cout << "earth: " << mapOfWords["earth"] << endl;
//    cout << "moon: " << mapOfWords.at("moon") << endl;
//    cout << "sun: " << mapOfWords.find("sun")->second << endl;
    
//    set<vector<pair<int, int>>> x;
//    vector<pair<int, int>> a, b;
//    a.push_back(make_pair(0, 1));
//    a.push_back(make_pair(1, 2));
//    a.push_back(make_pair(2, 3));
//    b.push_back(make_pair(0, 1));
//    b.push_back(make_pair(1, 2));
//    b.push_back(make_pair(2, 3));
//    x.insert(a);
//    x.insert(b);
//
//    for (vector<pair<int,int>> const &myvec: x) {
//        for (pair<int, int> mypair: myvec) {
//            cout << mypair.first << mypair.second << endl;
//        }
//    }
//    vector<event> j, k, l;
//    j.push_back(make_pair(1, make_pair(3, 5)));
//    j.push_back(make_pair(2, make_pair(4, 6)));
//    k.push_back(make_pair(1, make_pair(3, 5)));
//    k.push_back(make_pair(2, make_pair(4, 6)));
//    l.push_back(make_pair(2, make_pair(4, 6)));
//    l.push_back(make_pair(1, make_pair(3, 5)));
//
//    if (j==k) {
//        cout << "great!, j=k" << endl;
//    }
//    if (j!=l) {
//        cout << "great!, j!=l" << endl;
//    }
    
//    set<int> a;
//    a.insert(1);
//    set<int> b = a;
//    a.insert(2);
//    cout<< "a:" << a.size() << endl;
//    cout<< "b:" << b.size() << endl;
//
//    unordered_map<string, pair<int, set<vertex>>> mapOfWords;
//    mapOfWords["moon"].second.insert(4);
//    mapOfWords["earth"].first += 1;
//    cout << "earth: " << mapOfWords["earth"].first << mapOfWords["earth"].second.size() << endl;
//    cout << "moon: " << mapOfWords["moon"].first << mapOfWords.at("moon").second.size() << endl;
//    cout << "sun: " << mapOfWords.find("sun")->second.size() << endl;
    
//    a.insert(3);
//    a.insert(4);
//    a.insert(5);
//    for (auto it = a.begin(); it != a.end();) {
//        cout << *it << endl;
//        if (*it!=2) {
//            if (*it!=4) {
//                ++it;
//            } else {
//                it = a.erase(it);
//            }
//        } else {
//            it = a.erase(it);
//        }
//    }
//    cout << "end size" << a.size() << endl;
//    unordered_map<int, int> t;
//    t[5] = 15;
//    t[2] = 22;
//    for (auto it=t.begin(); it!=t.end(); ++it) {
//        cout << "star" << it->first << it->second << endl;
//    }
    
    fclose (fp);
    return 0;
}
