//
//  tmc.cpp
//  tmc
//
//  Created by Penghang Liu on 11/25/19.
//  Copyright © 2019 Penghang Liu. All rights reserved.
//

#include "tmc.hpp"

void Graph2motif(WGraph graph, adj_edges AE, int N_vtx, int N_event, map<string, int>&  motif_count, FILE* oofp){
    int iter = 0;
    for (auto it=graph.begin(); it!=graph.end(); ++it) {
        if (iter % 100 == 0) {
            cout << iter << " of " << graph.size() << endl;
        }
        iter++;
        edge e = it->first;
        long We = it->second;
        vertex u = e.first;
        vertex v = e.second;
        set<edge> edges(AE[u]);
        edges.insert(AE[v].begin(), AE[v].end());
        
        for (auto ap=edges.begin(); ap!=edges.end(); ++ap) {
            // 3n2e
            edge a = *ap;
            long Wa = graph[a];
            if (a.first==u && a.second==v){
                if (We<=1) continue;
                fprintf(oofp, "%s,%s,%s,%ld\n", u.c_str(), v.c_str(), "R", We);
//                if(N_event==3){
//                    fprintf(oofp, "%s,%s,%s,%ld\n", u.c_str(), v.c_str(), "R", We*(We-1)*(We-2)/6);
//                }
            } else {
                string S1 = staticEncode(a, e);
//                string S2 = staticEncode(e, a);
                fprintf(oofp, "%s,%s,%s,%ld\n", a.first.c_str(), a.second.c_str(), S1.c_str(), min(We,Wa));
//                fprintf(oofp, "%s,%s,%s,%ld\n", u.c_str(), v.c_str(), S2.c_str(), We*Wa);
                // 3n3e
                if (N_event==3) {
                    for (auto bp=next(ap, 1); bp!=edges.end(); ++bp) {
                        edge b = *bp;
                        long Wb = graph[b];
                        if(isTri(a, b, e)){
                            fprintf(oofp, "%s,%s,%s,%ld\n", a.first.c_str(), a.second.c_str(), "T", min(We,min(Wa,Wb)));
                            fprintf(oofp, "%s,%s,%s,%ld\n", b.first.c_str(), b.second.c_str(), "T", min(We,min(Wa,Wb)));
                            fprintf(oofp, "%s,%s,%s,%ld\n", e.first.c_str(), e.second.c_str(), "T", min(We,min(Wa,Wb)));
                        }
                    }
                }
            }
        }
    }
    return;
}

string occurrence(timestamp a, timestamp b, vector<timestamp> T1, vector<timestamp> T2, vector<timestamp> T3, int d_c){
    string output;
    output.append(",");
    if (T1.size()>0) {
        set<string> temp;
        for (int i=0; i<T1.size(); i++) {
            if (T1[i]<a-d_c || T1[i]>b+d_c) continue;
            temp.insert(occurrence(a, b, T1[i]));
        }
        for (auto it=temp.begin(); it!=temp.end(); ++it) {
            output.append(*it);
        }
    }
    output.append(",");
    if (T2.size()>0) {
        set<string> temp;
        for (int i=0; i<T2.size(); i++) {
            if (T2[i]<a-d_c || T2[i]>b+d_c) continue;
            temp.insert(occurrence(a, b, T2[i]));
        }
        for (auto it=temp.begin(); it!=temp.end(); ++it) {
            output.append(*it);
        }
    }
    output.append(",");
    if (T3.size()>0) {
        set<string> temp;
        for (int i=0; i<T3.size(); i++) {
            if (T3[i]<a-d_c || T3[i]>b+d_c) continue;
            temp.insert(occurrence(a, b, T3[i]));
        }
        for (auto it=temp.begin(); it!=temp.end(); ++it) {
            output.append(*it);
        }
    }
    return output;
}

string occurrence(timestamp a, timestamp b, timestamp c, vector<timestamp> T1, vector<timestamp> T2, vector<timestamp> T3, int d_c){
    string output;
    output.append(",");
    if (T1.size()>0) {
        set<string> temp;
        for (int i=0; i<T1.size(); i++) {
            if (T1[i]<a-d_c || T1[i]>c+d_c) continue;
            temp.insert(occurrence(a, b, c, T1[i]));
        }
        for (auto it=temp.begin(); it!=temp.end(); ++it) {
            output.append(*it);
        }
    }
    output.append(",");
    if (T2.size()>0) {
        set<string> temp;
        for (int i=0; i<T2.size(); i++) {
            if (T2[i]<a-d_c || T2[i]>c+d_c) continue;
            temp.insert(occurrence(a, b, c, T2[i]));
        }
        for (auto it=temp.begin(); it!=temp.end(); ++it) {
            output.append(*it);
        }
    }
    output.append(",");
    if (T3.size()>0) {
        set<string> temp;
        for (int i=0; i<T3.size(); i++) {
            if (T3[i]<a-d_c || T3[i]>c+d_c) continue;
            temp.insert(occurrence(a, b, c, T3[i]));
        }
        for (auto it=temp.begin(); it!=temp.end(); ++it) {
            output.append(*it);
        }
    }
    return output;
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

bool isTri(edge a, edge b, edge c){
    map<vertex, int> code;
    vector<vertex> temp;
    temp.push_back(a.first);
    temp.push_back(a.second);
    temp.push_back(b.first);
    temp.push_back(b.second);
    temp.push_back(c.first);
    temp.push_back(c.second);
    for (int i=0; i<temp.size(); i++) {
        if (code.find(temp[i])==code.end()){
            code[temp[i]] = 1;
        } else {
            code[temp[i]] += 1;
        }
    }
    for (auto it=code.begin(); it!=code.end(); ++it) {
        if (it->second!=2) return false;
    }
    return true;
}

string staticEncode(edge a, edge b){
    string motif;
    if (a.second==b.first) {
        if (a.first==b.second) {
            motif="P";
        } else {
            motif="C-1";
        }
    } else {
        if (a.first==b.first) {
            motif="O";
        } else if (a.second==b.second){
            motif="I";
        } else if (b.second==a.first){
            motif = "C-2";
        }
    }
    return motif;
}

string induceEncode(edge a, edge b, SGraph g){
    string motif;
    map<vertex, string> code;
    code[a.first] = "0";
    code[a.second] = "1";
    motif.append("01");
    vector<vertex> temp;
    vector<vertex> nodes;
    nodes.push_back(a.first);
    nodes.push_back(a.second);
    temp.push_back(b.first);
    temp.push_back(b.second);
//    temp.push_back(c.first);
//    temp.push_back(c.second);
    for (int i=0; i<temp.size(); i++) {
        if (code.find(temp[i])==code.end()){
            code[temp[i]] = "2";
            nodes.push_back(temp[i]);
        }
        motif.append(code[temp[i]]);
    }
    string l2;
    for (auto it=nodes.begin(); it!=nodes.end(); ++it) {
        for (auto itt=next(it,1); itt!=nodes.end(); ++itt) {
            vertex x = *it;
            vertex y = *itt;
            if (g[x].find(y)!=g[x].end()) {
                l2.append(code[x]);
                l2.append(code[y]);
            }
        }
    }
    if (l2.size()>0) {
        motif.append("_");
        motif.append(l2);
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
    vector<vertex> nodes;
    nodes.push_back(a.first);
    nodes.push_back(a.second);
    temp.push_back(b.first);
    temp.push_back(b.second);
    temp.push_back(c.first);
    temp.push_back(c.second);
    for (int i=0; i<temp.size(); i++) {
        if (code.find(temp[i])==code.end()){
            code[temp[i]] = "2";
            nodes.push_back(temp[i]);
        }
        motif.append(code[temp[i]]);
    }
    string l2;
    for (auto it=nodes.begin(); it!=nodes.end(); ++it) {
        for (auto itt=next(it,1); itt!=nodes.end(); ++itt) {
            vertex x = *it;
            vertex y = *itt;
            if (g[x].find(y)!=g[x].end()) {
                l2.append(code[x]);
                l2.append(code[y]);
            }
        }
    }
    if (l2.size()>0) {
        motif.append("_");
        motif.append(l2);
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

string easyEncode(edge a, edge b, vector<vertex>& vlist){
    string motif;
    map<vertex, string> code;
    vlist.push_back(a.first);
    vlist.push_back(a.second);
    code[a.first] = "0";
    code[a.second] = "1";
    motif.append("01");
    vector<vertex> temp;
    temp.push_back(b.first);
    temp.push_back(b.second);
    for (int i=0; i<temp.size(); i++) {
        if (code.find(temp[i])==code.end()){
            code[temp[i]] = "2";
            vlist.push_back(temp[i]);
        }
        motif.append(code[temp[i]]);
    }
    return motif;
}

string easyEncode(edge a, edge b, edge c, vector<vertex>& vlist){
    string motif;
    map<vertex, string> code;
    vlist.push_back(a.first);
    vlist.push_back(a.second);
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
            vlist.push_back(temp[i]);
        }
        motif.append(code[temp[i]]);
    }
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

void createGraph (string filename, WGraph& graph, adj_edges& AE){
    ifstream in(filename);
    string line;
    
    while (getline(in, line)) {
        if (line[0] != '%' && line[0] != '#'){
            stringstream ss (line);
            vertex u, v;
            long w;
            edge e;
            event e_t;
            ss >> u >> v >> w;
            if (u != v) {
                e = make_pair(u, v);
                graph[e] = w;
                AE[u].insert(e);
                AE[v].insert(e);
            }
        }
    }
//    for (auto it=graph.begin(); it!=graph.end(); ++it) {
//        sort(it->second.begin(), it->second.end());
//        it->second.erase(unique(it->second.begin(), it->second.end()),it->second.end());
//    }
    cout << "graph nodes:" << AE.size() << endl;
    cout << "graph edges:" << graph.size() << endl;
    return;
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

void createUndirectedGraph (string filename, TGraph& graph){
    ifstream in(filename);
    string line;
    
    while (getline(in, line)) {
        if (line[0] != '%' && line[0] != '#'){
            stringstream ss (line);
            vertex u, v;
            timestamp t;
            edge e, e1;
            ss >> u >> v >> t;
            if (u != v) {
                e = make_pair(u, v);
                e1 = make_pair(v, u);
                graph[e].push_back(t);
                graph[e1].push_back(t);
            }
        }
    }
    for (auto it=graph.begin(); it!=graph.end(); ++it) {
        sort(it->second.begin(), it->second.end());
        it->second.erase(unique(it->second.begin(), it->second.end()),it->second.end());
    }
//    cout << "layer2 nodes:" << AE.size() << endl;
    cout << "layer2 edges:" << graph.size()/2 << endl;
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
