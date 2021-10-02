/**
 * Sort.c Source code for sorting functions which parallelize Radix sort of edge lists within graph structures. Simple
 * OpenMP-based multiprocessor parallel radix, MPI-based HPC parallel radix, and OpenMP/MPI tandem parallel radix sorts
 * are implemented and explored. Several iterations of design and parallelization were explored, and some specific
 * sections were chosen to remain in serial format. Please consult accompanying paper for details on these design exper-
 * iments and the resultant design choices included in this code.
 *
 * Author: Stevan M. Dupor
 * Created: 09/15/21
 * Copyright 2021 by Stevan Dupor and NC State University. No unauthorized duplication permitted.
 */

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



/**
 * OpenMP-based parallelization of a radix sort. This flexible sort takes advantage of the configured number of processors,
 * but is designed with nested parallelism such that powers-of-two quantities of threads shall be used. Loop splitting
 * is the key source of paralellism, as parallelization of the final step in several methods proved to only yield perfor-
 * mance loss at benchmarking. Please consult the accompanying paper on these design choices and motivation for this
 * submitted design.
 */
#ifdef OPENMP_HARNESS
struct Graph *radixSortEdgesBySourceOpenMP (struct Graph *graph)
{
   printf("*** START Radix Sort Edges By Source OpenMP *** \n");
   int i, key, pos, digits=0;
   const int granularity = 8, bitmask = 0xff; // 8-bit buckets, as masked by 0xff
   const int iters = 32/granularity;

   // To support the nested parallelism below, and to ensure correct performance on the ARC cluster, force thread
   // size configuration to a power of two
   if(numThreads < 4)
      numThreads = 4;
   else if(numThreads % 4 != 0 && numThreads % 8 != 0 && numThreads % 16 != 0 && numThreads % 32 != 0 )
      numThreads = 8;

   // Ensure OMP is properly configured
   omp_set_nested(1);
   omp_set_num_threads(numThreads);

   // Calculate thread quantity for inner nested threads
   const int thread_per_bitbucket = numThreads/iters;

   // Initialize constants and data structures
   const int edges = graph->num_edges;
   const int blocksize = bitmask + 1;
   const int all_blocks = blocksize + (thread_per_bitbucket*blocksize);

   struct Edge *sorted_edges_array = newEdgeArray(graph->num_edges);
   struct Edge *temp;

   int vertex_count_tmp[iters][all_blocks];
   int vertex_count[iters][blocksize];

   // Initialize counting structures. Use serial for higher speed due to limited size.
   for (int k = 0; k < iters; ++k) {
      for (int q = 0; q < all_blocks; ++q)
         vertex_count_tmp[k][q] = 0;
      for (int q = 0; q < blocksize; ++q)
         vertex_count[k][q] = 0;
   }

   // Initialize offsets
   int thr=0, thr_offset=0;
   int offset=0;
   int end=0, row=0;

#pragma omp parallel for default(none) shared(vertex_count_tmp, graph, digits) private(i, thr, thr_offset, offset, end, key, row) num_threads(iters)
   for (digits = 0; digits < 32; digits = digits + granularity) {
#pragma omp parallel default(none) shared(vertex_count_tmp, graph, digits) private(i, thr, thr_offset, offset, end, key, row) num_threads(thread_per_bitbucket)
      {
         // Initialize offsets for each thread. Note, these run unsynchronized.
         thr = omp_get_thread_num();
         offset = thr * (edges / omp_get_num_threads());
         end = offset + (edges / omp_get_num_threads());
         thr_offset = blocksize*thr;
         row = digits/granularity;

         // Make sure we get all the edges in case there is a remainder
         if ((edges - end) < (edges/omp_get_num_threads()) && (edges - end) > 0)
            end = edges;

         // count occurrence of key: id of a source vertex
         for (i = offset; i < end; ++i) {
            key = graph->sorted_edges_array[i].src;
            vertex_count_tmp[row][((key >> digits) & bitmask) + thr_offset]++;
         }
      }
   }

   const int stop = bitmask+1;
   int joffset=0;
   // Collapse vertex count temp structure (expanded for multiple threads) into a single structure, vertex_count
#pragma omp parallel for default(none) shared(vertex_count, vertex_count_tmp) private(joffset, i) schedule(static) num_threads(4)
   for(int f = 0; f<iters; ++f){
      for (int j = 0; j < thread_per_bitbucket; ++j){
         joffset = j*stop;
         for (i = 0; i < stop; ++i) {
            if(j==0)
               vertex_count[f][i] = vertex_count_tmp[f][i];
            else
               vertex_count[f][i] += vertex_count_tmp[f][i + joffset];
         }
      }
   }

   // Perform the linear transformation step of the radix sort for all bucket sizes. Perform in serial for speed.
   for(int j=0;j < iters;++j)
      for (i = 1; i < blocksize; ++i)
         vertex_count[j][i] += vertex_count[j][i - 1];

      // Perform the sort action once for each bucket. Note, this step is performed serially as a result of extensive testing with pipelined
      // parallel methodologies. Please see the paper accompanying this code for details on this implementation choice.
   for (digits = 0; digits < 32; digits += granularity) {
      // Pre-calculate the bucket to sort from
      row = digits/granularity;
      // Sort this bucket
      for (i = edges - 1; i >= 0; --i) {
         key = graph->sorted_edges_array[i].src;
         key = (key >> digits) & bitmask;
         pos = --vertex_count[row][key];
         sorted_edges_array[pos] = graph->sorted_edges_array[i];
      }
      // Swap the dirty and clean arrays
      temp = graph->sorted_edges_array;
      graph->sorted_edges_array = sorted_edges_array;
      sorted_edges_array = temp;
   }

   // Free the un-used array
   free(sorted_edges_array);

   return graph;
}
#endif

/**
 * MPI-paralellized Radix sort of a graph's edgelist. This sort splits the buckets into two halves, with one MPI
 * host handling the most-singificant 16 bits and one process handling the least-significant 16 bits. The MPI protocol
 * is specifically designed to be extremely lightweight and keep dependence on high-latency gigabit ethernet-based TCP
 * to a minimum; as such, All init, count, crush, and transform steps are performed in parallel -- then, the hosts
 * synchronize and the root host performs the rest of the sort and graph functionality.
 *
 * PREREQUISITES: The initialization steps must be allowed to occur on all hosts. This is to prevent double-transmitting
 * graph data over TCP, and to place the load on the BeeGFS SAN instead of the inter-host TCP connection. Also, MPI_Init()
 * should have already been handled in the int main() function.
 *
 * POSTREQUISITE: int main() should also be prepared to handle return values from closing worker processes and freeing
 * un-used memory in these worker processes, while the root process finished handling the program tasks.
 */
#ifdef MPI_HARNESS
struct Graph *radixSortEdgesBySourceMPI (struct Graph *graph)
{
    printf("*** START Radix Sort Edges By Source MPI*** \n");

   int i, key, pos, digits;
   int granularity = 8, bitmask = 0xff; // 8-bit buckets, as masked by 0xff


   const int edges = graph->num_edges;

      int size, my_rank;
   MPI_Comm_size(MPI_COMM_WORLD, &size);
   MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
   if(size != 2 ) {
      printf("**************************************************************\n"
             "**************************************************************\n"
             "**************************************************************\n"
             "\n          INVALID MPI SIZING DETECTED. PLEASE RE-RUN     \n"
             "                     WITH TWO MPI HOSTS, eg:                   \n"
             "mpirun -n 2 ./bin/run-graph-mpi -f <graph file> -r [root] -n [num threads]\n\n"
             "**************************************************************\n"
             "**************************************************************\n"
             "**************************************************************\n");
         abort();
   }


   if (my_rank !=0)
      my_rank = 1;
   struct Edge *sorted_edges_array = newEdgeArray(graph->num_edges);
   struct Edge *temp;

   int vertex_count[4][bitmask+1];

   for(int k = 0 ; k<4;++k)
      for(int q = 0; q<=bitmask+1;++q)
         vertex_count[k][q] = 0;

   for (digits = my_rank * (2 * granularity);
        digits <= (my_rank * 2 * granularity) + granularity; digits = digits + granularity) {

      // zero Out count array
      for (i = 0; i < bitmask + 1; ++i) {
         vertex_count[(digits / granularity)][i] = 0;
      }

      // count occurrence of key: id of a source vertex
      for (i = 0; i < edges; ++i) {
         key = graph->sorted_edges_array[i].src;
         vertex_count[(digits / granularity)][(key >> digits) & (bitmask)]++;
      }

      // transform to cumulative sum
      for (i = 1; i < bitmask + 1; ++i) {
         vertex_count[(digits / granularity)][i] += vertex_count[(digits / granularity)][i - 1];
      }
   }

   if( my_rank == 0) {
      MPI_Recv(vertex_count[2], bitmask+1, MPI_INT, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vertex_count[3], bitmask+1, MPI_INT, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Ssend(vertex_count[0], bitmask+1, MPI_INT, 1, 0, MPI_COMM_WORLD);
      MPI_Ssend(vertex_count[1], bitmask+1, MPI_INT, 1, 0, MPI_COMM_WORLD);
   }
   else {
      MPI_Ssend(vertex_count[2], bitmask+1, MPI_INT, 0, 0, MPI_COMM_WORLD);
      MPI_Ssend(vertex_count[3], bitmask+1, MPI_INT, 0, 0, MPI_COMM_WORLD);
      MPI_Recv(vertex_count[0], bitmask+1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vertex_count[1], bitmask+1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
   }


   for (digits = 0; digits < 32; digits += granularity) {
      // fill-in the sorted array of edges
      for (i = edges - 1; i >= 0; --i) {
         key = graph->sorted_edges_array[i].src;
         key = (key >> digits) & bitmask;
         pos = --vertex_count[(digits / granularity)][key];
         sorted_edges_array[pos] = graph->sorted_edges_array[i];
      }
      // Swap the dirty and clean arrays
      temp = graph->sorted_edges_array;
      graph->sorted_edges_array = sorted_edges_array;
      sorted_edges_array = temp;
   }

   free(sorted_edges_array);

   return graph;
}
#endif

/**
 * MPI-with-OpenMP parallelized Radix sort of a graph's edgelist. This sort splits the buckets into two halves, with one MPI
 * host handling the most-singificant 16 bits and one process handling the least-significant 16 bits. On each host,
 * the init, count, crush, and transform tasks are divided with omp across as many processors as reasonably possible
 * (Currently tuned for the NCSU ARC cluster). Similar to the MPI section above, the MPI portion of the protocol
 * is specifically intended to be lightweight, with a small synchronization section between init-count-crush-transform and
 * sort-across-buckets.
 *
 * PREREQUISITES: The initialization steps must be allowed to occur on all hosts. This is to prevent double-transmitting
 * graph data over TCP, and to place the load on the BeeGFS SAN instead of the inter-host TCP connection. Also, MPI_Init()
 * should have already been handled in the int main() function.
 *
 * POSTREQUISITE: int main() should also be prepared to handle return values from closing worker processes and freeing
 * un-used memory in these worker processes, while the root process finished handling the program tasks.
 */
#ifdef HYBRID_HARNESS
struct Graph *radixSortEdgesBySourceHybrid (struct Graph *graph)
{
    printf("*** START Radix Sort Edges By Source Hybrid*** \n");
          int i, key, pos, digits=0;
      const int granularity = 8, bitmask = 0xff; // 8-bit buckets, as masked by 0xff
      const int iters = 32/granularity;

      const int thread_per_bitbucket = 8;///(8/granularity);
      omp_set_num_threads(thread_per_bitbucket);

      const int edges = graph->num_edges;
      const int blocksize = bitmask + 1;
      const int all_blocks = blocksize + (thread_per_bitbucket*blocksize);

   int size, my_rank;
   MPI_Comm_size(MPI_COMM_WORLD, &size);
   MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
   if(size != 2 ) {
      printf("**************************************************************\n"
             "**************************************************************\n"
             "**************************************************************\n"
             "\n          INVALID MPI SIZING DETECTED. PLEASE RE-RUN     \n"
             "                     WITH TWO MPI HOSTS, eg:                   \n"
             "mpirun -n 2 ./bin/run-graph-hybrid -f <graph file> -r [root] -n [num threads]\n\n"
             "**************************************************************\n"
             "**************************************************************\n"
             "**************************************************************\n");
         abort();
   }

      if (my_rank !=0)
         my_rank = 1;

      struct Edge *sorted_edges_array = newEdgeArray(graph->num_edges);
      struct Edge *temp;

      int vertex_count_tmp[iters][all_blocks];
      int vertex_count[iters][blocksize];


      for (int k = 0; k < iters; ++k) {
            for (int q = 0; q < all_blocks; ++q)
               vertex_count_tmp[k][q] = 0;
            for (int q = 0; q < blocksize; ++q)
               vertex_count[k][q] = 0;
         }

      int thr=0, thr_offset=0;
      int offset=0;
      int end=0, row=0;


      for (digits = my_rank * ((iters/2) * granularity); digits <= (my_rank * (iters/2) * granularity) + granularity; digits = digits + granularity) {
#pragma omp parallel default(none) shared(vertex_count_tmp, graph, digits, my_rank) private(i, thr, thr_offset, offset, end, key, row)
         {
            thr = omp_get_thread_num();
            offset = thr * (edges / omp_get_num_threads());
            end = offset + (edges / omp_get_num_threads());
            thr_offset = blocksize*thr;
            row = digits/granularity;
            if ((edges - end) < (edges/omp_get_num_threads()) && (edges - end) > 0)
               end = edges;

            // count occurrence of key: id of a source vertex
            for (i = offset; i < end; ++i) {

               key = graph->sorted_edges_array[i].src;
               vertex_count_tmp[row][((key >> digits) & bitmask) + thr_offset]++;
            }
         }
      }
      const int stop = bitmask+1;
      int joffset=0;

      for(int f = 0; f<iters; ++f){
         for (int j = 0; j < thread_per_bitbucket; ++j){
            joffset = j*stop;
            for (i = 0; i < stop; ++i) {
            if(j==0)
               vertex_count[f][i] = vertex_count_tmp[f][i];
            else
               vertex_count[f][i] += vertex_count_tmp[f][i + joffset];
            }
         }
      }

      for(int j=0;j < iters;++j)
         for (i = 1; i < blocksize; ++i)
            vertex_count[j][i] += vertex_count[j][i - 1];


      if( my_rank == 0) {
         if(iters == 4) {
            MPI_Recv(vertex_count[2], bitmask + 1, MPI_INT, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(vertex_count[3], bitmask + 1, MPI_INT, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Ssend(vertex_count[0], bitmask + 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
            MPI_Ssend(vertex_count[1], bitmask + 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
         }
         else if(iters == 2) {

            MPI_Recv(vertex_count[1], bitmask + 1, MPI_INT, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Ssend(vertex_count[0], bitmask + 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
         }
      }
      else {
         if(iters == 4) {
            MPI_Ssend(vertex_count[2], bitmask + 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
            MPI_Ssend(vertex_count[3], bitmask + 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
            MPI_Recv(vertex_count[0], bitmask + 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(vertex_count[1], bitmask + 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
         }
         else if (iters == 2) {
            MPI_Ssend(vertex_count[1], bitmask + 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
            MPI_Recv(vertex_count[0], bitmask + 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
         }
      }
      MPI_Barrier(MPI_COMM_WORLD);


      for (digits = 0; digits < 32; digits += granularity) {
         // fill-in the sorted array of edges
         row = digits/granularity;
         for (i = edges - 1; i >= 0; --i) {
            key = graph->sorted_edges_array[i].src;
            key = (key >> digits) & bitmask;
            pos = --vertex_count[row][key];
            sorted_edges_array[pos] = graph->sorted_edges_array[i];
         }
         // Swap the dirty and clean arrays
         temp = graph->sorted_edges_array;
         graph->sorted_edges_array = sorted_edges_array;
         sorted_edges_array = temp;
      }

      free(sorted_edges_array);

      return graph;

}
#endif

/** The original Serial count sort code. This function is ONLY used for datapoint generation, and validation of correct-
 * ness of the other sorting methods. In this submission, this function is NEVER called.
 *
 * @param graph the unsorted input graph.
 * @return a pointer to the sorted output graph.
 */
struct Graph *serial_count_sort (struct Graph *graph)
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

/** A simplistic parallel implementation of countsort, using built-in synchronization (atomic directives). This function
 * is used ONLY for datapoint generation and comparison purposes; for this submitted code, this function is NEVER called.
 *
 * @param graph The unsorted input graph.
 * @return  The sorted output graph.
 */
struct Graph * countSortEdgesBySource(struct Graph * graph) {
   int i, key, pos;
   struct Edge *sorted_edges_array = newEdgeArray(graph->num_edges);
   // auxiliary arrays, allocated at the start up of the program
   int *vertex_count = (int *) malloc(graph->num_vertices * sizeof(int)); // needed for Counting Sort
   const int vertices = graph->num_vertices;
   const int edges = graph->num_edges;

#pragma omp parallel for default(none) shared(vertex_count) schedule(static, 256)
   for (i = 0; i < vertices; ++i) {
      vertex_count[i] = 0;
   }

   // count occurrence of key: id of a source vertex
#pragma omp parallel for default(none) shared(vertex_count, graph) private(key, i) schedule(static, 256)
   for (i = 0; i < edges; ++i) {
      key = graph->sorted_edges_array[i].src;
#pragma omp atomic
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