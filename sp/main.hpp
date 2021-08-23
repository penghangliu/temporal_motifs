//
//  tmc.hpp
//  tmc
//
//  Created by Penghang Liu on 11/25/19.
//  Copyright © 2019 Penghang Liu. All rights reserved.
//

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

class Graph
{
    unordered_map<vertex, vector<vertex>> adj;
    
public:
    Graph(string filename);
    void addEdge(vertex u, vertex v);
    void removeEdge(vertex u, vertex v);
    int distance(vertex u, vertex v);
    bool check(vertex u, vertex v);
};

Graph::Graph(string filename)
{
    ifstream in(filename);
    string line;
    
    while (getline(in, line)) {
        if (line[0] != '%' && line[0] != '#'){
            stringstream ss (line);
            vertex u, v;
            ss >> u >> v;
            if (u != v) {
                adj[u].push_back(v);
                adj[v].push_back(u);
            }
        }
    }
    for (auto it=adj.begin(); it!=adj.end(); ++it) {
        it->second.erase(unique(it->second.begin(), it->second.end()),it->second.end());
    }
}

void Graph::addEdge(vertex u, vertex v)
{
    adj[u].push_back(v);
    adj[v].push_back(u);
}

void Graph::removeEdge(vertex u, vertex v)
{
    adj[u].erase(find(adj[u].begin(), adj[u].end(), v));
    adj[v].erase(find(adj[v].begin(), adj[v].end(), u));
}

int Graph::distance(vertex u, vertex v)
{
    queue<vertex> open;
    set<vertex> visited;
    int lvl = 1;
    int current = 1;
    
    open.push(u);
    while (!open.empty()) {
        vertex k = open.front();
        visited.insert(k);
        open.pop();
        current--;
        
        for (int i=0; i<adj[k].size(); i++) {
            vertex x = adj[k][i];
            if(x==v) return lvl;
            if (visited.find(x)==visited.end()) {
                open.push(x);
            }
        }
        
        if (current == 0) {
            lvl++;
            current = open.size();
        }
    }
    return -1;
}

bool Graph::check(vertex u, vertex v)
{
    if (find(adj[u].begin(), adj[u].end(), v) != adj[u].end()) {
        return true;
    }
    else return false;
}
