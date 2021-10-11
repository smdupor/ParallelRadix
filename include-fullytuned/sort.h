#ifndef SORT_H
#define SORT_H

#include "edgelist.h"
#include "graph.h"
#include "timer.h"
// Order edges by id of a source vertex, 
// using the Counting Sort
// Complexity: O(E + V)

struct Graph* radix_serial (struct Graph* graph, struct Timer *timer);
struct Graph* radixSortEdgesBySourceOpenMP (struct Graph* graph, struct Timer *timer);
struct Graph* radixSortEdgesBySourceMPI (struct Graph* graph, struct Timer *timer);
struct Graph* radixSortEdgesBySourceHybrid (struct Graph* graph, struct Timer *timer);
struct Graph* countSortEdgesBySource (struct Graph* graph, struct Timer *timer);
struct Graph* serial_count_sort (struct Graph* graph, struct Timer *timer);

extern int numThreads;

#endif