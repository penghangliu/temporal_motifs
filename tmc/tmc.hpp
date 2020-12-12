//
//  tmc.hpp
//  tmc
//
//  Created by Penghang Liu on 11/25/19.
//  Copyright © 2019 Penghang Liu. All rights reserved.
//

#ifndef tmc_hpp
#define tmc_hpp

#include <stdio.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <queue>
#include <stack>
#include <vector>
#include <set>
#include <map>
#include <unordered_map>
#include <algorithm>
#include <utility>
#include <string>
#include <initializer_list>

#include <assert.h>
#include <errno.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <random>
#include <chrono>
#include <sys/stat.h>

using namespace std;

typedef chrono::duration<double> tms;
typedef string vertex;
typedef long timestamp;
typedef pair<vertex, vertex> edge;
typedef pair<timestamp, edge> event;
typedef vector<event> key;  //a prefix or a motif
typedef pair<int, set<vertex>> counts;  //the count of the (prefix/motif) and the vertices in the (prefix/motif)
//typedef unordered_map<vector<event>, set<vertex>> prefix;
typedef unordered_map<vector<event>, pair<int, set<vertex>>> instancemap; //a hashtable of key and counts
typedef unordered_map<edge, vector<timestamp>> TGraph; //temporal graph
typedef unordered_map<vertex, set<edge>> adj_edges;
typedef unordered_map<vertex, set<vertex>> SGraph;

template <class T>
inline void hash_combine(std::size_t & seed, const T & v)
{
    std::hash<T> hasher;
    seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

namespace std
{
    template<typename S, typename T> struct hash<pair<S, T>>
    {
        inline size_t operator()(const pair<S, T> & v) const
        {
            size_t seed = 0;
            ::hash_combine(seed, v.first);
            ::hash_combine(seed, v.second);
            return seed;
        }
    };
    
    template<typename S, typename T> struct hash<vector<pair<S, T>>>
    {
        inline size_t operator()(const vector<pair<S, T>> &k) const
        {
            size_t seed = 0;
            for (auto it=k.begin(); it != k.end(); ++it) {
                hash_combine(seed, it->first);
                hash_combine(seed, it->second);
            }
            return seed;
        }
    };
}

inline void print_time (FILE* fp, const string& str, tms t) {
    fprintf (fp, "%s %.6lf\n", str.c_str(), t.count());
    fflush(fp);
}

set<edge> getLayer2 (set<vertex> nodes, SGraph g);
string induceEncode(edge a, edge b, SGraph g);
string induceEncode(edge a, edge b, edge c, SGraph g);
string occurrence(timestamp a, timestamp b, timestamp s);
string occurrence(timestamp a, timestamp b, timestamp c, timestamp s);
string occurrence(timestamp a, timestamp b, vector<timestamp> T1, vector<timestamp> T2, vector<timestamp> T3, int d_c);
string occurrence(timestamp a, timestamp b, timestamp c, vector<timestamp> T1, vector<timestamp> T2, vector<timestamp> T3, int d_c);
string easyEncode(edge a, edge b);
string easyEncode(edge a, edge b, edge c);
string easyEncode(edge a, edge b, vector<vertex>& vlist);
string easyEncode(edge a, edge b, edge c, vector<vertex>& vlist);
string complexEncode(edge a, edge b, edge s);
string complexEncode(edge a, edge b, edge c, edge s);
//bool checkConnect(edge a, edge b, edge e);
int checkNodes(edge a, edge b);
int checkNodes(edge a, edge b, edge e);
int checkNodes(edge a, edge b, edge e, edge c);
void createGraph (string filename, SGraph& graph);
void createGraph (string filename, TGraph& graph, adj_edges& AE);
void createUndirectedGraph (string filename, TGraph& graph);
void createEvents (string filename, vector<event>& events); //Load and sort the event list
void countInstance (event e, instancemap& imap, set<vector<event>>& keys, int N_vtx, int N_event, int d_c, int d_w);    //Increment the instance count and update the prefix type
string encodeMotif(vector<event> instance); //identify the type of motif
void countMotif (event e, set<key>& pre, map<string, int>& motif_count, int N_vtx, int N_event, int d_c, int d_w);
set<vertex> getNodes(vector<event> key);
void countSpecificmotif (event e, set<key>& pre, int& motif_count, string code_given, int N_vtx, int N_event, int d_c, int d_w);
char sconvert (int i);
void removeIsomorphic (map<string, int>& motif_count);
void Graph2motif(TGraph graph, adj_edges AE, TGraph graph_s, SGraph g, adj_edges BE, int d_c, int d_w, int N_vtx, int N_event, map<string, int>&  motif_count, bool multi, string method);
#endif /* tmc_hpp */
