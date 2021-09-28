#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <omp.h>
#include <string.h>
#include <mpi.h>

#include "sort.h"
#include "graph.h"
#include "timer.h"


struct Graph *serial_count_sort(struct Graph *graph, struct Timer *timer) {
   int i;
   int key;
   int pos;
   struct Edge *sorted_edges_array = newEdgeArray(graph->num_edges);
   // auxiliary arrays, allocated at the start up of the program
   int *vertex_count = (int *) malloc(graph->num_vertices * sizeof(int)); // needed for Counting Sort
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

struct Graph *countSortEdgesBySource(struct Graph *graph, struct Timer *timer) {
   int i;
   int key;
   int pos;
   struct Edge *sorted_edges_array = newEdgeArray(graph->num_edges);
   // auxiliary arrays, allocated at the start up of the program
   int *vertex_count = (int *) malloc(graph->num_vertices * sizeof(int)); // needed for Counting Sort

   /**************************** HOLDING ZONE **************************/
   // printf("Thread num: %i of %i \n", omp_get_thread_num(), omp_get_num_threads());

   // memset(vertex_count, 0, graph->num_vertices*sizeof(int));


   /**************************** HOLDING ZONE **************************/
   const int vertices = graph->num_vertices;
   const int edges = graph->num_edges;
   /*************** remove edges, vertices**/
#pragma omp parallel default(none) shared(graph, vertex_count, sorted_edges_array) private (key, i, pos)
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

struct Graph *radix_serial(struct Graph *graph, struct Timer *timer) {
   int i, key, pos;
   int max = 0, granularity = 8, bitmask = 0xff; // 8-bit buckets, as masked by 0xff

   const int edges = graph->num_edges;

   for (i = 0; i < edges; ++i) {
      if (max < graph->sorted_edges_array[i].src)
         max = graph->sorted_edges_array[i].src;
   }

   struct Edge *sorted_edges_array = newEdgeArray(graph->num_edges);
   struct Edge *temp;
   int *vertex_count = (int *) malloc((2 * bitmask) * sizeof(int));

   for (int digits = 0; max >> digits > 0; digits += granularity) {
      // zero Out count array
      for (i = 0; i < (2 * bitmask); ++i) {
         vertex_count[i] = 0;
      }
      // count occurrence of key: id of a source vertex
      for (i = 0; i < edges; ++i) {
         key = graph->sorted_edges_array[i].src;
         vertex_count[(key >> digits) & (bitmask)]++;
      }
      // transform to cumulative sum
      for (i = 1; i < (2 * bitmask); ++i) {
         vertex_count[i] += vertex_count[i - 1];
      }

      // fill-in the sorted array of edges
      for (i = edges - 1; i >= 0; --i) {
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
   return graph;
}

struct Graph *radixSortEdgesBySourceOpenMP(struct Graph *graph, struct Timer *timer) {
   int i, key, pos, digits=0;
   const int granularity = 8, bitmask = 0xff; // 8-bit buckets, as masked by 0xff
   const int iters = 32/granularity;

   const int thread_per_bitbucket = 2;
   omp_set_nested(1);
   omp_set_num_threads(thread_per_bitbucket*4);
   printf("OMP Num THR: %i", thread_per_bitbucket*4);

   const int edges = graph->num_edges;
   const int blocksize = bitmask + 1;
   const int all_blocks = blocksize + (thread_per_bitbucket*blocksize);

   struct Edge *sorted_edges_array = newEdgeArray(graph->num_edges);
   struct Edge *temp;

   int vertex_count_tmp[iters][all_blocks];
   int vertex_count[iters][blocksize];
   Start(&timer[INIT]);

   for (int k = 0; k < iters; ++k) {
      for (int q = 0; q < all_blocks; ++q)
         vertex_count_tmp[k][q] = 0;
      for (int q = 0; q < blocksize; ++q)
         vertex_count[k][q] = 0;
   }

   Stop(&timer[INIT]);

   int thr=0, thr_offset=0;
   int offset=0;
   int end=0, row=0;

   Start(&timer[COUNT]);
#pragma omp parallel for default(none) shared(vertex_count_tmp, graph, digits) private(i, thr, thr_offset, offset, end, key, row)
   for (digits = 0; digits < 32; digits = digits + granularity) {
#pragma omp parallel default(none) shared(vertex_count_tmp, graph, digits) private(i, thr, thr_offset, offset, end, key, row) num_threads(4)
      {
         thr = omp_get_thread_num();
         offset = thr * (edges / omp_get_num_threads());
         end = offset + (edges / omp_get_num_threads());
         thr_offset = blocksize*thr;
         row = digits/granularity;
        // printf("thread %i Row %i thr_off %i start %i end %i \n ", thr, row, thr_offset, offset, end);
         if ((edges - end) < (edges/omp_get_num_threads()) && (edges - end) > 0)
            end = edges;

         // count occurrence of key: id of a source vertex
         for (i = offset; i < end; ++i) {
            //  printf("i at kill was: %i of %i of %i on %i rnk %i addr %08xi\n", i, end, edges, thr, my_rank, *&graph->sorted_edges_array[i-1].src);
            key = graph->sorted_edges_array[i].src;
            vertex_count_tmp[row][((key >> digits) & bitmask) + thr_offset]++;

         }

      }
   }
   const int stop = bitmask+1;
   int joffset=0;
   Stop(&timer[COUNT]);

   Start(&timer[CRUSH]);

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
   Stop(&timer[CRUSH]);
   Start(&timer[XFORM]);

   for(int j=0;j < iters;++j)
      for (i = 1; i < blocksize; ++i)
         vertex_count[j][i] += vertex_count[j][i - 1];

   Stop(&timer[XFORM]);
   Start(&timer[SORT]);

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
   Stop(&timer[SORT]);
   free(sorted_edges_array);

   return graph;

   /*int i, key, pos, digits;
   int granularity = 8, bitmask = 0xff; // 8-bit buckets, as masked by 0xff

   const int edges = graph->num_edges;

   struct Edge *sorted_edges_array = newEdgeArray(graph->num_edges);
   struct Edge *temp;
   int vertex_count[4][bitmask + 1];

   omp_set_num_threads(4);
#pragma omp parallel default(none) shared(graph, vertex_count, sorted_edges_array, bitmask, granularity, temp, digits) private (key, i, pos)
   {
#pragma omp for schedule(static)
      for (digits = 0; digits < 32; digits += granularity) {
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
   }

   for (digits = 0; digits < 32; digits += granularity) {
      // fill-in the sorted array of edges

      for (i = edges - 1; i >= 0; --i) {
         key = graph->sorted_edges_array[i].src;
         key = (key >> digits) & bitmask;
         pos = --vertex_count[(digits / granularity)][key];
         sorted_edges_array[pos] = graph->sorted_edges_array[i];

      }
      printf("done with sort for digits: %i\n", digits);
      // Swap the dirty and clean arrays
      temp = graph->sorted_edges_array;
      graph->sorted_edges_array = sorted_edges_array;
      sorted_edges_array = temp;
   }

   //free the extra arrays

   free(sorted_edges_array);
   return graph;*/
}

struct Graph *radixSortEdgesBySourceMPI(struct Graph *graph, struct Timer *timer) {
      int i, key, pos, digits;
      int granularity = 8, bitmask = 0xff; // 8-bit buckets, as masked by 0xff

      const int edges = graph->num_edges;

      int my_rank;
      MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

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

   struct Graph *radixSortEdgesBySourceHybrid(struct Graph *graph, struct Timer *timer) {

      int i, key, pos, digits=0;
      const int granularity = 8, bitmask = 0xff; // 8-bit buckets, as masked by 0xff
      const int iters = 32/granularity;

      const int thread_per_bitbucket = 8;///(8/granularity);
      omp_set_num_threads(thread_per_bitbucket);

      const int edges = graph->num_edges;
      const int blocksize = bitmask + 1;
      const int all_blocks = blocksize + (thread_per_bitbucket*blocksize);
      Start(&timer[INIT]);
      Start(&timer[COUNT]);
      Start(&timer[CRUSH]);
      Start(&timer[XFORM]);
      Start((&timer[MPI_MSG]));
      Start(&timer[SORT]);

      int my_rank=0;
      MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

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

      Stop(&timer[INIT]);

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
              //  printf("i at kill was: %i of %i of %i on %i rnk %i addr %08xi\n", i, end, edges, thr, my_rank, *&graph->sorted_edges_array[i-1].src);
               key = graph->sorted_edges_array[i].src;
               vertex_count_tmp[row][((key >> digits) & bitmask) + thr_offset]++;
            }
         }
      }
      const int stop = bitmask+1;
      int joffset=0;
      Stop(&timer[COUNT]);



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
      Stop(&timer[CRUSH]);


      for(int j=0;j < iters;++j)
         for (i = 1; i < blocksize; ++i)
            vertex_count[j][i] += vertex_count[j][i - 1];

         Stop(&timer[XFORM]);


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
      Stop(&timer[MPI_MSG]);


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
      Stop(&timer[SORT]);
      free(sorted_edges_array);

      return graph;

   }
