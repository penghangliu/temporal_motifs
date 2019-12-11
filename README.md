# Counting Temporal Motifs

Compile the code

`g++ -std=c++11 main.cpp tmc.cpp`

Run

version 1: enumerate all instance for all types of motifs

`./a.out v1 inputfile delta_C delta_W |V| |E|`

version 2: speed up version 1

`./a.out v2 inputfile delta_C delta_W |V| |E|`

version 3: specified motif query

`./a.out v3 inputfile delta_C delta_W motif`
