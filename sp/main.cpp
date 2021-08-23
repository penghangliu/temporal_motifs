//
//  main.cpp
//  tmc
//
//  Created by Penghang Liu on 11/24/19.
//  Copyright © 2019 Penghang Liu. All rights reserved.
//

#include <iostream>
#include "main.hpp"

int main(int argc, char * argv[]) {
    const auto t1 = chrono::steady_clock::now();
    
    string graph (argv[1]);
    string gname = graph.substr (graph.find_last_of("/") + 1, graph.find_last_of(".") - graph.find_last_of("/") -1);
    string target (argv[2]);
    string tname = target.substr (target.find_last_of("/") + 1, target.find_last_of(".") - target.find_last_of("/") -1);
    string out_file = "out_" + gname + "_" + tname;
    FILE* fp = fopen (out_file.c_str(), "w");
    
//Read file
    Graph G(graph);
        
    const auto t2 = chrono::steady_clock::now();
//    print_time (fp, "Read data time: ", t2 - t1);
    
//Count
    ifstream in(target);
    string line;
    
    int i = 0;
    while (getline(in, line)) {
        if (line[0] != '%' && line[0] != '#'){
            stringstream ss (line);
            vertex u, v;
            ss >> u >> v;
            bool isin = G.check(u,v);
            if (isin) {
                G.removeEdge(u,v);
            }
            int d = G.distance(u,v);
            fprintf(fp, "%s,%s,%d \n", u.c_str(), v.c_str(), d);
            if(isin) G.addEdge(u,v);
        }
        i++;
        if (i % 1000 == 0) {
            cout << i << endl;
        }
    }
    
    const auto t3 = chrono::steady_clock::now();
//    print_time (fp, "End-to-end Time: ", t3 - t1);
    
    fclose (fp);
    return 0;
}
