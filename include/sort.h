#ifndef SORT_H
#define SORT_H

#include "edgelist.h"
#include "graph.h"
// Order edges by id of a source vertex, 
// using the Counting Sort
// Complexity: O(E + V)

// struct Graph* radixSortEdgesBySourceOpenMP (struct Graph* graph);
// struct Graph* radixSortEdgesBySourceMPI (struct Graph* graph);
// struct Graph* radixSortEdgesBySourceHybrid (struct Graph* graph);
struct Graph* countSortEdgesBySource (struct Graph* graph);

extern int numThreads;

#endif