//
//  tmc.cpp
//  tmc
//
//  Created by Penghang Liu on 11/25/19.
//  Copyright © 2019 Penghang Liu. All rights reserved.
//

#include "tmc.hpp"

void Graph2motif(TGraph graph, adj_edges AE, TGraph graph_s, SGraph g, adj_edges BE, int d_c, int d_w, int N_vtx, int N_event, map<string, int>&  motif_count, bool multi, string method){
    for (auto it=graph.begin(); it!=graph.end(); ++it) {
        edge e = it->first;
        vector<timestamp> Tm = it->second;
        vertex u = e.first;
        vertex v = e.second;
        set<edge> edges(AE[u]);
        edges.insert(AE[v].begin(), AE[v].end());
        set<edge> edges_s(BE[u]);
        edges_s.insert(BE[v].begin(), BE[v].end());
        
        for (auto ap=edges.begin(); ap!=edges.end(); ++ap) {
            edge a = *ap;
            vector<timestamp> Ta = graph[a];
            string S;
            if (method == "v1") {
                S = induceEncode(a, e, g);
            } else {
                S = easyEncode(a, e);
            }
            for (int j=0; j<Tm.size(); j++) {
                for (int i=0; i<Ta.size(); i++) {
                    if(abs(Ta[i]-Tm[j])>d_c || Ta[i] > Tm[j]) continue;
                    motif_count[S] += 1;
                }
            }
            // 3n2e
            if (method == "v2") {
                for (auto cp=edges_s.begin(); cp!=edges_s.end(); ++cp) {
                    edge c = *cp;
                    if(checkNodes(a, e)==2 && checkNodes(a, c, e)>2) continue;
                    if(checkNodes(a, c, e)>3) continue;
                    vector<timestamp> Tc = graph_s[c];
                    string s = complexEncode(a, e, c);
                    for (int j=0; j<Tm.size(); j++) {
                        for (int i=0; i<Ta.size(); i++) {
                            if(abs(Ta[i]-Tm[j])>d_c || Ta[i] > Tm[j]) continue;
                            for (int l=0; l<Tc.size(); l++) {
                                if (Tc[l]>=(Ta[i]-d_c) && Tc[l]<=(Tm[j]+d_c)) {
                                    string plus = occurrence(Ta[i], Tm[j], Tc[l]);
                                    string m = s + plus;
                                    motif_count[m] += 1;
                                }
                            }
                        }
                    }
                }
            }
            // 3n3e
            for (auto bp=ap; bp!=edges.end(); ++bp) {
                edge b = *bp;
                if(checkNodes(a, b, e)>3) continue;
                vector<timestamp> Tb = graph[b];
                if (method == "v2") {
                    for (auto cp=edges_s.begin(); cp!=edges_s.end(); ++cp) {
                        edge c = *cp;
                        if(checkNodes(a, b, e)==2 && checkNodes(a, b, e, c)>2) continue;
                        if(checkNodes(a, b, e, c)>3) continue;
                        vector<timestamp> Tc = graph_s[c];
                        string s1 = complexEncode(a, e, b, c);
                        string s2 = complexEncode(b, e, a, c);
                        for (int j=0; j<Tm.size(); j++) {
                            for (int i=0; i<Ta.size(); i++) {
                                if(abs(Ta[i]-Tm[j])>d_c) continue;
                                for (int k=0; k<Tb.size(); k++) {
                                    if(abs(Tb[k]-Tm[j])>d_c) continue;
                                    if(abs(Ta[i]-Tb[k])>d_w) continue;
                                    for (int l=0; l<Tc.size(); l++) {
                                        if(Ta[i] < Tm[j] && Tm[j] < Tb[k]){
                                            if (Tc[l]>=(Ta[i]-d_c) && Tc[l]<=(Tb[k]+d_c)) {
                                                string plus = occurrence(Ta[i], Tm[j], Tb[k], Tc[l]);
                                                string m = s1 + plus;
                                                motif_count[m] += 1;
                                            }
                                        }
                                        if(a.first==b.first && a.second==b.second) continue;
                                        if (Ta[i] > Tm[j] && Tm[j] > Tb[k]) {
                                            if (Tc[l]>=(Tb[k]-d_c) && Tc[l]<=(Ta[i]+d_c)) {
                                                string plus = occurrence(Tb[k], Tm[j], Ta[i], Tc[l]);
                                                string m = s2 + plus;
                                                motif_count[m] += 1;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                string S1, S2;
                if (method == "v1") {
                    S1 = induceEncode(a, e, b, g);
                    S2 = induceEncode(b, e, a, g);
                } else {
                    S1 = easyEncode(a, e, b);
                    S2 = easyEncode(b, e, a);
                }
//                if (s1=="010101" || s2=="010101") {
//                    cout << "YES" << endl;
//                }
                for (int j=0; j<Tm.size(); j++) {
                    for (int i=0; i<Ta.size(); i++) {
                        if(abs(Ta[i]-Tm[j])>d_c) continue;
                        for (int k=0; k<Tb.size(); k++) {
                            if(abs(Tb[k]-Tm[j])>d_c) continue;
                            if(abs(Ta[i]-Tb[k])>d_w) continue;
                            if(Ta[i] < Tm[j] && Tm[j] < Tb[k]){
                                motif_count[S1] += 1;
                            }
                            if(a.first==b.first && a.second==b.second) continue;
                            if (Ta[i] > Tm[j] && Tm[j] > Tb[k]) {
                                motif_count[S2] += 1;
                            }
                        }
                    }
                }
            }
        }
    }
    return;
}

string occurrence(timestamp a, timestamp b, timestamp s){
    if (s <= a) {
        return "0";
    } else if (a < s && s <= b){
        return "1";
    } else {
        return "2";
    }
}

string occurrence(timestamp a, timestamp b, timestamp c, timestamp s){
    if (s <= a) {
        return "0";
    } else if (a < s && s <= b){
        return "1";
    } else if (b < s && s <= c){
        return "2";
    } else {
        return "3";
    }
}

set<edge> getLayer2 (set<vertex> nodes, SGraph g){
    set<edge> output;
    for (auto it=nodes.begin(); it!=nodes.end(); ++it) {
        for (auto itt = next(it, 1); itt!=nodes.end(); ++itt) {
            if (g[*it].find(*itt)!=g[*it].end()) {
                output.insert(make_pair(*it, *itt));
            }
        }
    }
    return output;
}

string induceEncode(edge a, edge b, SGraph g){
    string motif;
    map<vertex, string> code;
    code[a.first] = "0";
    code[a.second] = "1";
    motif.append("01");
    vector<vertex> temp;
    temp.push_back(b.first);
    temp.push_back(b.second);
//    temp.push_back(c.first);
//    temp.push_back(c.second);
    for (int i=0; i<temp.size(); i++) {
        if (code.find(temp[i])==code.end()){
            code[temp[i]] = "2";
        }
        motif.append(code[temp[i]]);
    }
    set<vertex> nodes;
    for (auto it=code.begin(); it!=code.end(); ++it) {
        vertex x = it->first;
        nodes.insert(x);
    }
    set<edge> edges = getLayer2(nodes, g);
    if (edges.size()>0) {
        motif.append("_");
    }
    for (auto ep=edges.begin(); ep!=edges.end(); ++ep) {
        edge e = *ep;
        string u = code[e.first];
        string v = code[e.second];
        if (stoi(u)<stoi(v)) {
            motif.append(u);
            motif.append(v);
        } else {
            motif.append(v);
            motif.append(u);
        }
    }
    return motif;
}

string induceEncode(edge a, edge b, edge c, SGraph g){
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
    set<vertex> nodes;
    for (auto it=code.begin(); it!=code.end(); ++it) {
        vertex x = it->first;
        nodes.insert(x);
    }
    set<edge> edges = getLayer2(nodes, g);
    if (edges.size()>0) {
        motif.append("_");
    }
    for (auto ep=edges.begin(); ep!=edges.end(); ++ep) {
        edge e = *ep;
        string u = code[e.first];
        string v = code[e.second];
        if (stoi(u)<stoi(v)) {
            motif.append(u);
            motif.append(v);
        } else {
            motif.append(v);
            motif.append(u);
        }
    }
    return motif;
}

string complexEncode(edge a, edge b, edge s){
    string motif;
    map<vertex, string> code;
    code[a.first] = "0";
    code[a.second] = "1";
    motif.append("01");
    vector<vertex> temp;
    temp.push_back(b.first);
    temp.push_back(b.second);
//    temp.push_back(c.first);
//    temp.push_back(c.second);
    for (int i=0; i<temp.size(); i++) {
        if (code.find(temp[i])==code.end()){
            code[temp[i]] = "2";
        }
        motif.append(code[temp[i]]);
    }
    motif.append("_");
    string u = code[s.first];
    string v = code[s.second];
    if (stoi(u)<stoi(v)) {
        motif.append(u);
        motif.append(v);
    } else {
        motif.append(v);
        motif.append(u);
    }
    motif.append("_");
    return motif;
}

string complexEncode(edge a, edge b, edge c, edge s){
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
    motif.append("_");
    string u = code[s.first];
    string v = code[s.second];
    if (stoi(u)<stoi(v)) {
        motif.append(u);
        motif.append(v);
    } else {
        motif.append(v);
        motif.append(u);
    }
    motif.append("_");
    return motif;
}

string easyEncode(edge a, edge b){
    string motif;
    map<vertex, string> code;
    code[a.first] = "0";
    code[a.second] = "1";
    motif.append("01");
    vector<vertex> temp;
    temp.push_back(b.first);
    temp.push_back(b.second);
    for (int i=0; i<temp.size(); i++) {
        if (code.find(temp[i])==code.end()){
            code[temp[i]] = "2";
        }
        motif.append(code[temp[i]]);
    }
    return motif;
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

//bool checkConnect(edge a, edge b, edge e){
//    set<vertex> V;
//    V.insert(a.first);
//    V.insert(a.second);
//    V.insert(b.first);
//    V.insert(b.second);
//    V.insert(e.first);
//    V.insert(e.second);
//    if (V.size()>3) {
//        return false;
//    }
//    return true;
//}

int checkNodes(edge a, edge b){
    set<vertex> V;
    V.insert(a.first);
    V.insert(a.second);
    V.insert(b.first);
    V.insert(b.second);
    return V.size();
}

int checkNodes(edge a, edge b, edge e, edge c){
    set<vertex> V;
    V.insert(a.first);
    V.insert(a.second);
    V.insert(b.first);
    V.insert(b.second);
    V.insert(e.first);
    V.insert(e.second);
    V.insert(c.first);
    V.insert(c.second);
    return V.size();
}

int checkNodes(edge a, edge b, edge e){
    set<vertex> V;
    V.insert(a.first);
    V.insert(a.second);
    V.insert(b.first);
    V.insert(b.second);
    V.insert(e.first);
    V.insert(e.second);
    return V.size();
}

void createGraph (string filename, SGraph& graph){
    ifstream in(filename);
    string line;
    
    while (getline(in, line)) {
        if (line[0] != '%' && line[0] != '#'){
            stringstream ss (line);
            vertex u, v;
            ss >> u >> v;
            if (u != v) {
                graph[u].insert(v);
                graph[v].insert(u);
            }
        }
    }
    cout << "layer2 nodes:" << graph.size() << endl;
    return;
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
    cout << "layer2 nodes:" << AE.size() << endl;
    cout << "layer2 edges:" << graph.size() << endl;
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
