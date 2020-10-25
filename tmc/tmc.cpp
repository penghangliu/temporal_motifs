//
//  tmc.cpp
//  tmc
//
//  Created by Penghang Liu on 11/25/19.
//  Copyright © 2019 Penghang Liu. All rights reserved.
//

#include "tmc.hpp"

void Graph2motif(TGraph graph, adj_edges AE, int d_c, int d_w, int N_vtx, int N_event, map<string, int>&  motif_count){
    for (auto it=graph.begin(); it!=graph.end(); ++it) { // second edge
        edge e = it->first;
        vector<timestamp> Tm = it->second;
        vertex u = e.first;
        vertex v = e.second;
        set<edge> edges(AE[u]);
        edges.insert(AE[v].begin(), AE[v].end());

        for (auto ap=edges.begin(); ap!=edges.end(); ++ap) { // first edge
            edge a = *ap;
            set<vertex> nodes = get_Nodes(a, e);
            set<edge> third = nodes.size()==2 ? edges : get_third(AE, nodes);
//            cout << a.first << a.second << " " << e.first << e.second << " " << third.size() << " ";
            for (auto bp=third.begin(); bp!=third.end(); ++bp) { // third edge
                edge b = *bp;
//                cout << b.first << b.second << ": ";
//                if(!checkConnect(a, b, e)) continue;
                if(graph.find(b)==graph.end()) continue;;
                vector<timestamp> Ta = graph[a];
                vector<timestamp> Tb = graph[b];
//                string s1 = easyEncode(a, e, b);
                string s2 = easyEncode(b, e, a);
//                cout << a.first << a.second << " " << e.first << e.second << " " << b.first << b.second << " ";
                for (int j=0; j<Tm.size(); j++) {
                    for (int i=0; i<Ta.size(); i++) {
                        if(abs(Ta[i]-Tm[j])>d_c) continue;
                        //binary search
                        if (Ta[i] > Tm[j]) {
                            timestamp upper = Tm[j];
                            timestamp lower = max(Ta[i]-d_w, Tm[j]-d_c);
                            int temp_count = n_larger_eq(Tb, lower);
                            temp_count -= n_larger_eq(Tb, upper);
//                            cout << s2 << " " << temp_count << " ";
                            if(temp_count>0) motif_count[s2] += temp_count;
                        }
//                        else if (Ta[i] < Tm[j] && (a.first!=b.first || a.second!=b.second)){
//                            timestamp lower = Tm[j];
//                            timestamp upper = min(Ta[i]+d_w, Tm[j]+d_c);
//                            int temp_count = n_less_eq(Tb, upper);
//                            temp_count -= n_less_eq(Tb, lower);
//                            cout << s1 << " " << temp_count << " ";
//                            if(temp_count>0) motif_count[s1] += temp_count;
//                        }
                        //original
//                        for (int k=0; k<Tb.size(); k++) {
//                            if(abs(Tb[k]-Tm[j])>d_c) continue;
//                            if(abs(Ta[i]-Tb[k])>d_w) continue;
//                            if(Ta[i] < Tm[j] && Tm[j] < Tb[k]){
//                                motif_count[s1] += 1;
////                                cout << "+1" << " ";
////                                cout << Ta[i] << Tm[j] << Tb[k] << " ";
//                            } else if (Ta[i] > Tm[j] && Tm[j] > Tb[k] && (a.first!=b.first || a.second!=b.second)) {
//                                motif_count[s2] += 1;
////                                cout << "+2" << " ";
////                                cout << Ta[i] << Tm[j] << Tb[k] << " ";
//                            }
//                        }
                    }
                }
//                cout << endl;
            }
//            cout << endl;
        }
    }
    return;
}

set<edge> get_third(adj_edges AE, set<vertex> nodes){
    set<edge> output;
    for (auto it=nodes.begin(); it!=nodes.end(); ++it) {
        auto itt=it;
        itt++;
        while (itt!=nodes.end()) {
            edge a = make_pair(*it, *itt);
            edge b = make_pair(*itt, *it);
//            if(AE[*it].find(a)!=AE[*it].end()) output.insert(a);
//            if(AE[*it].find(b)!=AE[*it].end()) output.insert(b);
            output.insert(a);
            output.insert(b);
            itt++;
        }
    }
    return output;
}

set<vertex> get_Nodes(edge a, edge b){
    set<vertex> nodes;
    nodes.insert(a.first);
    nodes.insert(a.second);
    nodes.insert(b.first);
    nodes.insert(b.second);
    return nodes;
}

int n_larger_eq(vector<timestamp> T, timestamp t){
    int l=0;
    int n=T.size();
    int r=n-1;
    int output = n;
    while (l<=r) {
        int m = (l + r) / 2;
        if (T[m]>t) {
            output = r;
            r = m - 1;
        } else if (T[m]<t){
            l = m + 1;
        } else{
            return n-m;
        }
    }
    return n-output;
}

int n_less_eq(vector<timestamp> T, timestamp t){
    int l=0;
    int n=T.size();
    int r=n-1;
    int output = 0;
    while (l<=r) {
        int m = (l + r) / 2;
        if (T[m]>t) {
            r = m - 1;
        } else if (T[m]<t){
            output = l+1;
            l = m + 1;
        } else{
            return m+1;
        }
    }
    return output;
}

string easyEncode(edge a, edge b, edge c){
    string motif;
    map<vertex, string> code;
    code[a.first] = "0";
    code[a.second] = "1";
    motif.append("01");
    vector<vertex> temp;
    temp.push_back(b.first);
    temp.push_back(b.second);
    temp.push_back(c.first);
    temp.push_back(c.second);
    for (int i=0; i<temp.size(); i++) {
        if (code.find(temp[i])==code.end()){
            code[temp[i]] = "2";
        }
        motif.append(code[temp[i]]);
    }
    return motif;
}

bool checkConnect(edge a, edge b, edge e){
    set<vertex> V;
    V.insert(a.first);
    V.insert(a.second);
    V.insert(b.first);
    V.insert(b.second);
    V.insert(e.first);
    V.insert(e.second);
    if (V.size()>3) {
        return false;
    }
    return true;
}

void createGraph (string filename, TGraph& graph, adj_edges& AE){
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
                graph[e].push_back(t);
                AE[u].insert(e);
                AE[v].insert(e);
            }
        }
    }
    for (auto it=graph.begin(); it!=graph.end(); ++it) {
        sort(it->second.begin(), it->second.end());
        it->second.erase(unique(it->second.begin(), it->second.end()),it->second.end());
    }
    return;
}

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
    events.erase(unique(events.begin(), events.end()),events.end()); //remove duplicates
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
                        if (key.back().first!=e.first) { //check synchronous events
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
            it = keys.erase(it);    //remove prefix if it exceeds the delta constrain
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
    string temp;
    bool concurrent {false};
    timestamp t0 = instance[0].first - 1;
    map<vertex, string> code;
    int i=0;
    for (auto it=instance.begin(); it!=instance.end(); ++it) {
        vertex u = it->second.first;
        vertex v = it->second.second;
        timestamp t = it->first;
        if (t==t0) {
            motif.append("{");
            concurrent = true;
        }
        motif.append(temp);
        temp.clear();
        if (t!=t0&&concurrent){
            motif.append("}");
            concurrent = false;
        }
        t0 = t;
        if (code.find(u)==code.end()) {
            code[u] = to_string(i);
            i++;
        }
        temp.append(code[u]);
        if (code.find(v)==code.end()) {
            code[v] = to_string(i);
            i++;
        }
        temp.append(code[v]);
    }
    motif.append(temp);
    if (concurrent) {
        motif.append("}");
    }
    return motif;
}

//string encodeMotif(vector<event> instance){
//    string motif;
//    map<vertex, string> code;
//    int i=0;
//    unordered_map<timestamp, int> concurrent_count;
//    for (auto it=instance.begin(); it!=instance.end(); ++it) {
//        vertex u = it->second.first;
//        vertex v = it->second.second;
//        timestamp t = it->first;
//        concurrent_count[t] += 1;
//        if (concurrent_count[t]==2) {
//            motif.append("a");
//        }
//        if (code.find(u)==code.end()) {
//            code[u] = to_string(i);
//            i++;
//        }
//        motif.append(code[u]);
//        if (code.find(v)==code.end()) {
//            code[v] = to_string(i);
//            i++;
//        }
//        motif.append(code[v]);
//        if (concurrent_count[t]>1) {
//            char c = sconvert(concurrent_count[t]);
//            motif.push_back(c);
//        }
//    }
//    return motif;
//}

char sconvert (int i) {
    string s("abcdefghijklmnopqrstuvwxyz");
    return s.at(i-1);
}

set<vertex> getNodes(vector<event> key){
    set<vertex> nodes;
    for (int i=0; i<key.size(); i++) {
        nodes.insert(key[i].second.first);
        nodes.insert(key[i].second.second);
    }
    return nodes;
}

void countMotif (event e, set<key>& pre, map<string, int>& motif_count, int N_vtx, int N_event, int d_c, int d_w){
    vertex u = e.second.first;
    vertex v = e.second.second;
    vector<vector<event>> new_motif;    //used to store the new motifs
    for (auto it = pre.begin(); it != pre.end();) {   //for each current prefix
        vector<event> key = *it;
        set<vertex> nodes = getNodes(key);
        if (e.first - key.front().first <= d_w && e.first - key.back().first <= d_c ) {  //check delta C and delta W
            if (nodes.find(u)!=nodes.end() || nodes.find(v)!=nodes.end()) {
                nodes.insert(u);
                nodes.insert(v);
                if (nodes.size() <= N_vtx) {    //check the number of vertices
                    vector<event> motif = key;
                    motif.push_back(e);
                    if (motif.size()==N_event && nodes.size()==N_vtx) {
                        string code = encodeMotif(motif);
                        motif_count[code] += 1;
                    } else if(motif.size()<N_event) {
                        new_motif.push_back(motif);
                    }
                }
            }
            ++it;
        } else {
            it = pre.erase(it);    //remove prefix if it exceeds the delta constrain
        }
    }
    //add the new motifs to the current prefix list
    if (!new_motif.empty()) {
        for (vector<event> const &mt: new_motif) {
            pre.insert(mt);
        }
    }
    vector<event> E;
    E.push_back(e);
    pre.insert(E); // add the new event to the current prefix list
    return;
}

void countSpecificmotif (event e, set<key>& pre, int& motif_count, string code_given, int N_vtx, int N_event, int d_c, int d_w){
    vector<vector<event>> new_motif;    //used to store the new motifs
    for (auto it = pre.begin(); it != pre.end();) {   //for each current prefix
        vector<event> key = *it;
        if (e.first - key.front().first <= d_w && e.first - key.back().first <= d_c) {  //check delta C and delta W
            vector<event> motif = key;
            motif.push_back(e);
            set<vertex> nodes = getNodes(motif);
            if (motif.size()==N_event && nodes.size()==N_vtx) {
                string code = encodeMotif(motif);
                int l = code.length();
                if (code==code_given.substr(0,l)) {
                    motif_count += 1;
                }
            } else if(motif.size()<N_event) {
                string code = encodeMotif(motif);
                int l = code.length();
                if (code==code_given.substr(0,l)) {
                    new_motif.push_back(motif);
                }
            }
            ++it;
        } else {
            it = pre.erase(it);    //remove prefix if it exceeds the delta constrain
        }
    }
    //add the new motifs to the current prefix list
    if (!new_motif.empty()) {
        for (vector<event> const &mt: new_motif) {
            pre.insert(mt);
        }
    }
    vector<event> E;
    E.push_back(e);
    pre.insert(E); // add the new event to the current prefix list
    return;
}

void removeIsomorphic (map<string, int>& motif_count){
    
    return;
}
