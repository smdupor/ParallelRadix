#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <omp.h>
#include <string.h>
// #include <mpi.h>

#include "sort.h"
#include "graph.h"


struct Graph* serial_count_sort (struct Graph* graph){
   int i;
   int key;
   int pos;
   struct Edge *sorted_edges_array = newEdgeArray(graph->num_edges);
   // auxiliary arrays, allocated at the start up of the program
   int *vertex_count = (int*)malloc(graph->num_vertices*sizeof(int)); // needed for Counting Sort
   for (i = 0; i < graph->num_vertices; ++i) {
      vertex_count[i] = 0;
   }
   // count occurrence of key: id of a source vertex
   for (i = 0; i < graph->num_edges; ++i) {
      key = graph->sorted_edges_array[i].src;
      vertex_count[key]++;
   }
   // transform to cumulative sum
   for (i = 1; i < graph->num_vertices; ++i) {
      vertex_count[i] += vertex_count[i - 1];
   }
   // fill-in the sorted array of edges
   for (i = graph->num_edges - 1; i >= 0; --i) {
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

struct Graph* countSortEdgesBySource (struct Graph* graph){
   int i;
   int key;
   int pos;
   struct Edge *sorted_edges_array = newEdgeArray(graph->num_edges);
   // auxiliary arrays, allocated at the start up of the program
   int *vertex_count = (int*)malloc(graph->num_vertices*sizeof(int)); // needed for Counting Sort

   /**************************** HOLDING ZONE **************************/
   // printf("Thread num: %i of %i \n", omp_get_thread_num(), omp_get_num_threads());

  // memset(vertex_count, 0, graph->num_vertices*sizeof(int));


   /**************************** HOLDING ZONE **************************/
   const int vertices = graph->num_vertices;
   const int edges = graph->num_edges;
#pragma omp parallel default(none) shared(graph, vertex_count, vertices, edges, sorted_edges_array) private (key, i, pos)
   {
#pragma omp for schedule(static, 8)
      for (i = 0; i < vertices; ++i) {
         vertex_count[i] = 0;
      }

      // count occurrence of key: id of a source vertex
#pragma omp for schedule(static, 8)
      for (i = 0; i < edges; ++i) {
         key = graph->sorted_edges_array[i].src;
         #pragma omp atomic
            vertex_count[key]++;
      }
#pragma omp single
      // transform to cumulative sum
      for (i = 1; i < vertices; ++i) {
         vertex_count[i] += vertex_count[i - 1];
      }

      // fill-in the sorted array of edges
#pragma omp single
      for (i = edges - 1; i >= 0; --i) {

             pos = --vertex_count[graph->sorted_edges_array[i].src];

           //  #pragma omp critical
         sorted_edges_array[pos] = graph->sorted_edges_array[i];
         //  #pragma omp atomic
         // vertex_count[key]--;
      }
   }

   free(vertex_count);
   free(graph->sorted_edges_array);
   graph->sorted_edges_array = sorted_edges_array;

   return graph;
}

struct Graph* radix_serial(struct Graph* graph){
   int i, key, pos;
   int max=0, granularity = 8, bitmask = 0xff; // 8-bit buckets, as masked by 0xff

   for (i = 0; i < graph->num_edges; ++i) {
      if (max < graph->sorted_edges_array[i].src)
         max = graph->sorted_edges_array[i].src;
   }

   struct Edge *sorted_edges_array = newEdgeArray(graph->num_edges);
   struct Edge *temp;
   int *vertex_count = (int*)malloc((2^granularity)*sizeof(int));

   for(int digits = 0; max>>digits > 0; digits += granularity) {
      // zero Out count array
      for (i = 0; i < graph->num_vertices; ++i) {
         vertex_count[i] = 0;
      }
      // count occurrence of key: id of a source vertex
      for (i = 0; i < graph->num_edges; ++i) {
         key = graph->sorted_edges_array[i].src;
         vertex_count[(key >> digits) & (bitmask)]++;
      }
      // transform to cumulative sum
      for (i = 1; i < graph->num_vertices; ++i) {
         vertex_count[i] += vertex_count[i - 1];
      }
      // fill-in the sorted array of edges
      for (i = graph->num_edges - 1; i >= 0; --i) {
         key = graph->sorted_edges_array[i].src;
         pos = vertex_count[(key >> digits) & bitmask] - 1;
         sorted_edges_array[pos] = graph->sorted_edges_array[i];
         vertex_count[(key >> digits) & bitmask]--;
      }
      // Swap the dirty and clean arrays
      temp = graph->sorted_edges_array;
      graph->sorted_edges_array = sorted_edges_array;
      sorted_edges_array = temp;
   }
   //free the extra arrays
   free(vertex_count);
   free(sorted_edges_array);
   return graph;
}

 struct Graph* radixSortEdgesBySourceOpenMP (struct Graph* graph){
    int i, key, pos;
    int max=0, granularity = 8, bitmask = 0xff; // 8-bit buckets, as masked by 0xff

    for (i = 0; i < graph->num_edges; ++i) {
       if (max < graph->sorted_edges_array[i].src)
          max = graph->sorted_edges_array[i].src;
    }

    struct Edge *sorted_edges_array = newEdgeArray(graph->num_edges);
    struct Edge *temp;
    int *vertex_count = (int*)malloc((2^granularity)*sizeof(int));

    for(int digits = 0; max>>digits > 0; digits += granularity) {
       // zero Out count array
       for (i = 0; i < graph->num_vertices; ++i) {
          vertex_count[i] = 0;
       }
       // count occurrence of key: id of a source vertex
       for (i = 0; i < graph->num_edges; ++i) {
          key = graph->sorted_edges_array[i].src;
          vertex_count[(key >> digits) & (bitmask)]++;
       }
       // transform to cumulative sum
       for (i = 1; i < graph->num_vertices; ++i) {
          vertex_count[i] += vertex_count[i - 1];
       }
       // fill-in the sorted array of edges
       for (i = graph->num_edges - 1; i >= 0; --i) {
          key = graph->sorted_edges_array[i].src;
          pos = vertex_count[(key >> digits) & bitmask] - 1;
          sorted_edges_array[pos] = graph->sorted_edges_array[i];
          vertex_count[(key >> digits) & bitmask]--;
       }
       // Swap the dirty and clean arrays
       temp = graph->sorted_edges_array;
       graph->sorted_edges_array = sorted_edges_array;
       sorted_edges_array = temp;
    }
    //free the extra arrays
    free(vertex_count);
    free(sorted_edges_array);
    return graph;
 }
 struct Graph* radixSortEdgesBySourceMPI (struct Graph* graph){

         // printf("*** START Radix Sort Edges By Source MPI*** \n");
   return graph;
 }
 struct Graph* radixSortEdgesBySourceHybrid (struct Graph* graph){

         // printf("*** START Radix Sort Edges By Source Hybrid*** \n");
   return graph;
 }
