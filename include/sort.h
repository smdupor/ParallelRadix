/**
 * Sort.h Headers for sorting functions (See Sort.c) which parallelize Radix sort of edge lists within graph structures. Simple
 * OpenMP-based multiprocessor parallel radix, MPI-based HPC parallel radix, and OpenMP/MPI tandem parallel radix sorts
 * are implemented and explored. Several iterations of design and parallelization were explored, and some specific
 * sections were chosen to remain in serial format. Please consult accompanying paper for details on these design exper-
 * iments and the resultant design choices included in this code.
 *
 * Author: Stevan M. Dupor
 * Created: 09/15/21
 * Copyright 2021 by Stevan Dupor and NC State University. No unauthorized duplication permitted.
 */

#ifndef SORT_H
#define SORT_H

#include "edgelist.h"
#include "graph.h"
// Order edges by id of a source vertex,
// using the Counting Sort
// Complexity: O(E + V)

#ifdef OPENMP_HARNESS
struct Graph *radixSortEdgesBySourceOpenMP (struct Graph *graph);
#endif

#ifdef MPI_HARNESS
struct Graph *radixSortEdgesBySourceMPI (struct Graph *graph);
#endif

#ifdef HYBRID_HARNESS
struct Graph *radixSortEdgesBySourceHybrid (struct Graph *graph);
#endif

struct Graph *countSortEdgesBySource (struct Graph *graph);
struct Graph* serial_count_sort (struct Graph* graph);

extern int numThreads;

#endif