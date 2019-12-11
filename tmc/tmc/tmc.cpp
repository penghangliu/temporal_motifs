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
    events.erase(unique(events.begin(), events.end()),events.end());
    return;
}

void countInstance (event e, instancemap& imap, set<vector<event>>& keys, int N_vtx, int N_event, int d_c, int d_w){
    vertex u = e.second.first;
    vertex v = e.second.second;
    vector<vector<event>> new_motif;    //used to store the new motifs
    for (auto it = keys.begin(); it != keys.end();) {   //for each current prefix
        vector<event> key = *it;
        if (e.first - key.front().first <= d_w && e.first - key.back().first <= d_c) {  //check delta C and delta W
            if (key.size() < N_event) { //check the number of events
                set<vertex> nodes = imap[key].second;
                nodes.insert(u);
                nodes.insert(v);
                if (nodes.size() <= N_vtx) {    //check the number of vertices
                    if (imap[key].second.find(u)!=imap[key].second.end() || imap[key].second.find(v)!=imap[key].second.end()) {
                        if (key.back().first!=e.first) {
                            vector<event> motif = key;
                            motif.push_back(e);
                            new_motif.push_back(motif);
                            imap[motif].first += imap[key].first;
                            imap[motif].second = nodes;
                        }
                    }
                }
                ++it;
            } else {
                it = keys.erase(it);    //remove prefix if it exceeds the size constrain
            }
        } else {
            it = keys.erase(it);    //remove prefix if it exceeds the delat constrain
        }
    }
    //add the new motifs to the current prefix list
    if (!new_motif.empty()) {
        for (vector<event> const &mt: new_motif) {
            keys.insert(mt);
        }
    }
    vector<event> E;
    E.push_back(e);
    imap[E].first += 1;
    imap[E].second.insert(u);
    imap[E].second.insert(v);
    keys.insert(E); // add the new event to the current prefix list
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
