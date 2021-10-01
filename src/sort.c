#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

#ifdef OPENMP_HARNESS
#include <omp.h>
#endif

#ifdef MPI_HARNESS
#include <mpi.h>
#endif

#ifdef HYBRID_HARNESS
#include <omp.h>
#include <mpi.h>
#endif


#include "sort.h"
#include "graph.h"



#ifdef OPENMP_HARNESS
struct Graph *radixSortEdgesBySourceOpenMP (struct Graph *graph)
{

    printf("*** START Radix Sort Edges By Source OpenMP *** \n");
    return graph;
}
#endif

#ifdef MPI_HARNESS
struct Graph *radixSortEdgesBySourceMPI (struct Graph *graph)
{

    printf("*** START Radix Sort Edges By Source MPI*** \n");
    return graph;
}
#endif

#ifdef HYBRID_HARNESS
struct Graph *radixSortEdgesBySourceHybrid (struct Graph *graph)
{

    printf("*** START Radix Sort Edges By Source Hybrid*** \n");
    return graph;
}
#endif

/**
 * PARALLEL COUNT SORT for TESTING PURPOSES ONLY
 * @param graph
 * @return
 */
struct Graph *countSortEdgesBySource (struct Graph *graph)
{

   int i;
   int key;
   int pos;
   struct Edge *sorted_edges_array = newEdgeArray(graph->num_edges);

   // auxiliary arrays, allocated at the start up of the program
   int *vertex_count = (int *)malloc(graph->num_vertices * sizeof(int)); // needed for Counting Sort

   for(i = 0; i < graph->num_vertices; ++i)
   {
      vertex_count[i] = 0;
   }

   // count occurrence of key: id of a source vertex
   for(i = 0; i < graph->num_edges; ++i)
   {
      key = graph->sorted_edges_array[i].src;
      vertex_count[key]++;
   }

   // transform to cumulative sum
   for(i = 1; i < graph->num_vertices; ++i)
   {
      vertex_count[i] += vertex_count[i - 1];
   }

   // fill-in the sorted array of edges
   for(i = graph->num_edges - 1; i >= 0; --i)
   {
      key = graph->sorted_edges_array[i].src;
      pos = vertex_count[key] - 1;
      sorted_edges_array[pos] = graph->sorted_edges_array[i];
      vertex_count[key]--;
   }



   free(vertex_count);
   free(graph->sorted_edges_array);

   graph->sorted_edges_array = sorted_edges_array;

   return graph;

}