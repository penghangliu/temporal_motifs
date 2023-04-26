'make' to compile the code

'./TMC inputfile delta-C delta-W |V| |E| Resolution Consecutive(YES/NO) Constrained-dynamic-graphlet-counting(YES/NO)' to execute the code

    - for example, './TMC example.txt 15 30 3 3 5 NO YES' will count the three-node three-event motifs in example.txt, where delta-C is 15s, delta-W is 30s, based on then constrained dynamic graphlet countign algorithm. The dataset resolution is degraded to 5s. There are two output files:
    - '*_15_30_3_3' shows the number of each motif;
    - '*_15_0_3_3_time' contains the detail information, including intermediate event occurence and motif timespan.
    
