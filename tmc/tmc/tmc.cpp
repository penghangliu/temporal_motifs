//
//  tmc.cpp
//  tmc
//
//  Created by Penghang Liu on 11/25/19.
//  Copyright © 2019 Penghang Liu. All rights reserved.
//

#include "tmc.hpp"

void createEvents (string filename, vector<event>& events){
    ifstream in(filename);
    string line;
    
    while (getline(in, line)) {
        if (line[0] != '%' && line[0] != '#'){
            stringstream ss (line);
            vertex u, v;
            timestamp t;
            edge e;
            ss >> u >> v >> t;
            if (u != v) {
                e = make_pair(u, v);
                events.push_back(make_pair(t, e));
            }
        }
    }
    sort(events.begin(), events.end());
    return;
}

void countInstance (event e, instancemap& imap, set<vector<event>>& keys, int N_vtx, int N_event, int d_c, int d_w){
    vertex u = e.second.first;
    vertex v = e.second.second;
    vector<vector<event>> new_motif;
//    int a = 1;
//    cout << "check event " << u << " " << v << ": " << endl;
    for (auto it = keys.begin(); it != keys.end();) {
//        cout << "create" << endl;
        vector<event> key = *it;
//        cout << "check key " << a << " size " << key.size() << ": ";
        
//        cout << key.size() << endl;
//        cout << "event time: " << e.first << endl;
//        cout << "d_w: " << e.first - key.front().first << endl;
//        cout << "d_c: " << e.first - key.back().first << endl;
        if (e.first - key.front().first <= d_w && e.first - key.back().first <= d_c) {
//            a++;
            if (key.size() < N_event) {
                set<vertex> nodes = imap[key].second;
//                cout << "node size " << nodes.size() << endl;
                nodes.insert(u);
                nodes.insert(v);
                if (nodes.size() <= N_vtx) {
                    if (imap[key].second.find(u)!=imap[key].second.end() || imap[key].second.find(v)!=imap[key].second.end()) {
//                        cout << "set: ";
//                        for (auto tt = imap[key].second.begin(); tt!=imap[key].second.end(); ++tt) {
//                            cout << *tt << " ";
//                        }
//                        cout << endl;
                        vector<event> motif = key;
                        motif.push_back(e);
                        new_motif.push_back(motif);
                        imap[motif].first += imap[key].first;
                        imap[motif].second = nodes;
//                        cout << "time " << e.first << " "<< key.front().first << endl;
                    }
                }
                ++it;
            } else {
                it = keys.erase(it);
            }
        } else {
            it = keys.erase(it);
        }
    }
    if (!new_motif.empty()) {
        for (vector<event> const &mt: new_motif) {
            keys.insert(mt);
        }
    }
//    cout << "a " << a << endl;
//    cout << "map size " << imap.size() << endl;
    vector<event> E;
    E.push_back(e);
    imap[E].first += 1;
    imap[E].second.insert(u);
    imap[E].second.insert(v);
    keys.insert(E);
//    cout << "add edge " << u << " " << v << endl;
    return;
}

string encodeMotif(vector<event> instance){
    string motif;
    map<vertex, string> code;
    int i=0;
    for (auto it=instance.begin(); it!=instance.end(); ++it) {
        vertex u = it->second.first;
        vertex v = it->second.second;
        if (code.find(u)==code.end()) {
            code[u] = to_string(i);
            i++;
        }
        motif.append(code[u]);
        if (code.find(v)==code.end()) {
            code[v] = to_string(i);
            i++;
        }
        motif.append(code[v]);
    }
    return motif;
}
