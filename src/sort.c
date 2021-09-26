#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <omp.h>
#include <string.h>
#include <mpi.h>

#include "sort.h"
#include "graph.h"


struct Graph *serial_count_sort(struct Graph *graph) {
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

/*
struct Graph* countSortEdgesBySource (struct Graph* graph){
   int i;
   int key;
   int pos;
   struct Edge *sorted_edges_array = newEdgeArray(graph->num_edges);
   // auxiliary arrays, allocated at the start up of the program
   int *vertex_count = (int*)malloc(graph->num_vertices*sizeof(int)); // needed for Counting Sort

  //  **************************** HOLDING ZONE **************************
   // printf("Thread num: %i of %i \n", omp_get_thread_num(), omp_get_num_threads());

   // memset(vertex_count, 0, graph->num_vertices*sizeof(int));


   // **************************** HOLDING ZONE **************************
   const int vertices = graph->num_vertices;
   const int edges = graph->num_edges;
    // *************** remove edges, vertices**
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
}*/

struct Graph *countSortEdgesBySource(struct Graph *graph) {
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

struct Graph *radix_serial(struct Graph *graph) {
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

struct Graph *radixSortEdgesBySourceOpenMP(struct Graph *graph) {
   int i, key, pos, digits;
   int max = 0, granularity = 8, bitmask = 0xff; // 8-bit buckets, as masked by 0xff

   const int edges = graph->num_edges;
   const int jsize = edges / numThreads;
   const int jqty = edges / jsize;
   const int jrem = edges - (jsize * jqty);
   int maxes[jqty];
   printf("start maximizer section\n");

   /*********** remove jqty, jsize, edges*/

   for (i = 0; i < edges; ++i) {
      if (max < graph->sorted_edges_array[i].src)
         max = graph->sorted_edges_array[i].src;
   }
   /*
#pragma omp parallel default(none) shared(graph, maxes) private(i)
   {
#pragma omp single
      for (i = 0; i < jqty; ++i) {
         maxes[i] = 0;
      }
#pragma omp for schedule(static, 256)
      for (int j = 0; j < jqty; ++j) {
         for (i = j * jsize; i < (j + 1) * jsize; ++i) {
            if (maxes[j] < graph->sorted_edges_array[i].src)
               maxes[j] = graph->sorted_edges_array[i].src;
         }
      }
#pragma omp single
      for (i = jqty * jsize; i < edges; ++i) {
         if (maxes[jqty - 1] < graph->sorted_edges_array[i].src)
            maxes[jqty - 1] = graph->sorted_edges_array[i].src;
      }
   }
   max = 0;
   for (i = 0; i < jqty; ++i) {
      if (maxes[i] > max)
         max = maxes[i];
   }*/

   struct Edge *sorted_edges_array = newEdgeArray(graph->num_edges);
   struct Edge *temp;
   int vertex_count[4][bitmask + 1];

   omp_set_num_threads(4);
#pragma omp parallel default(none) shared(graph, vertex_count, sorted_edges_array, bitmask, granularity, max, temp, digits) private (key, i, pos)
   {
#pragma omp for schedule(static)
      for (digits = 0; digits < 32; digits += granularity) {
         // zero Out count array
         for (i = 0; i < bitmask + 1; ++i) {
            vertex_count[(digits / granularity)][i] = 0;
         }
         //  printf("Done with init for digits: %i,\n", digits);
         // count occurrence of key: id of a source vertex
         for (i = 0; i < edges; ++i) {
            key = graph->sorted_edges_array[i].src;
            vertex_count[(digits / granularity)][(key >> digits) & (bitmask)]++;
         }
         // printf("Done with count for digits: %i,\n", digits);
         // transform to cumulative sum
         for (i = 1; i < bitmask + 1; ++i) {
            vertex_count[(digits / granularity)][i] += vertex_count[(digits / granularity)][i - 1];
         }
         //printf("Done with linear xform for digits: %i,\n", digits);
      }
   }

   for (digits = 0; digits < 32; digits += granularity) {
      // fill-in the sorted array of edges
      //if(vertex_count[(digits / granularity)][1] == edges)
      // break;
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
   return graph;
}

struct Graph *radixSortEdgesBySourceMPI(struct Graph *graph) {
      int i, key, pos, digits;
      int max = 0, granularity = 8, bitmask = 0xff; // 8-bit buckets, as masked by 0xff

      const int edges = graph->num_edges;
      const int jsize = edges / numThreads;
      const int jqty = edges / jsize;
      const int jrem = edges - (jsize * jqty);
      int maxes[jqty];

      int my_rank;
      MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

      if (my_rank !=0)
         my_rank = 1;
   /*
      for (i = 0; i < edges; ++i) {
         if (max < graph->sorted_edges_array[i].src)
            max = graph->sorted_edges_array[i].src;
      }*/

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

/* MPI_Barrier(MPI_COMM_WORLD);
if(my_rank == 0)
 for(int k = 0 ; k<4;++k) {
    printf("Rank %i: %i: ", my_rank, k);
    for (int q = 0; q < 10; ++q)
       printf("%i.", vertex_count[k][q]);
    printf("\n");

 }
MPI_Barrier(MPI_COMM_WORLD);
if(my_rank == 1)
 for(int k = 0 ; k<4;++k) {
    printf("Rank %i: %i: ", my_rank, k);
    for (int q = 0; q < 10; ++q)
       printf("%i.", vertex_count[k][q]);
    printf("\n");
 }
MPI_Barrier(MPI_COMM_WORLD);*/
// }


   struct Graph *radixSortEdgesBySourceHybrid(struct Graph *graph) {
      int i, key, pos, digits=0;
      const int max = 0, granularity = 8, bitmask = 0xff; // 8-bit buckets, as masked by 0xff
      omp_set_num_threads(4);
      const int edges = graph->num_edges;
      const int thread_per_bitbucket = 4;///(8/granularity);
      // 2 threads = 2/4 = 0
      // 4 threads = 4/4 = 1
      // 8th = 8/4 =2
      // 16th = 16/4 = 4
      const int blocksize = bitmask + 1;
      const int all_blocks = blocksize + (thread_per_bitbucket*blocksize);

      int my_rank=0;
      MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

      if (my_rank !=0)
         my_rank = 1;
      /*
         for (i = 0; i < edges; ++i) {
            if (max < graph->sorted_edges_array[i].src)
               max = graph->sorted_edges_array[i].src;
         }*/

      struct Edge *sorted_edges_array = newEdgeArray(graph->num_edges);
      struct Edge *temp;

      int vertex_count_tmp[4][all_blocks];
      int vertex_count[4][blocksize];

#pragma omp parallel for default(none) shared (vertex_count_tmp) schedule(static)
      for(int k = 0 ; k<4;++k)
         for(int q = 0; q<=all_blocks;++q)
            vertex_count_tmp[k][q] = 0;

#pragma omp parallel for default(none) shared (vertex_count) schedule(static)
      for(int k = 0 ; k<4;++k)
         for(int q = 0; q<=blocksize;++q)
            vertex_count[k][q] = 0;

      int thr=0, thr_offset=0;
      int offset=0;
      int end=0, row=0;

      for (digits = my_rank * (2 * granularity); digits <= (my_rank * 2 * granularity) + granularity; digits = digits + granularity) {
#pragma omp parallel default(none) shared(vertex_count_tmp, graph, digits, my_rank) private(i, thr, thr_offset, offset, end, key, row)
         {
            thr = omp_get_thread_num();
/*            offset = thr * blocksize;
            end = offset + blocksize;
            //printf("This thread is # %i of %i digits %i granularity %i rank %i\n", thr, omp_get_num_threads(), digits, granularity, my_rank);
            // zero Out count array
#pragma omp for schedule(static)
            for (i = offset; i < end; ++i) {
               vertex_count[(digits / granularity)][i] = 0;
            }
*/
            offset = thr * (edges / omp_get_num_threads());
            end = offset + (edges / omp_get_num_threads()) - 1;
            thr_offset = blocksize*thr;
            row = digits/granularity;
            if ((edges - end) < (edges/omp_get_num_threads()) && (edges - end) > 0)
               end = edges-1;
          //  printf("Thread %i offset %i to end %i of %i of digits %i of rank %i on row %i w threadoff %i of allblock %i\n", thr, offset, end,edges, digits, my_rank, row, thr_offset, all_blocks);

            // count occurrence of key: id of a source vertex
            for (i = offset; i < end; ++i) {
              //  printf("i at kill was: %i of %i of %i on %i rnk %i addr %08xi\n", i, end, edges, thr, my_rank, *&graph->sorted_edges_array[i-1].src);
               key = graph->sorted_edges_array[i].src;
               vertex_count_tmp[row][((key >> digits) & bitmask) + thr_offset]++;
            }
         }

#pragma omp barrier

/************************************* issue with compression of arrays**********************/
      }
#pragma omp barrier
      MPI_Barrier(MPI_COMM_WORLD);
      for(int f = 0;f<4;++f){
         for (i = 0; i < 255; ++i) {
            for (int j = 0; j < 4; ++j)
               vertex_count[f][i] += vertex_count_tmp[f][i + (j * blocksize)];
         }
      }

         printf("EXIT PRAGMA LOOP ((EDGES ATM: %i EDGES REMT: %i))\n", edges, graph->num_edges);
      MPI_Barrier(MPI_COMM_WORLD);
      if(my_rank == 0)
         for(int k = 0 ; k<4;++k) {
            printf("Rank prexform %i: %i: ", my_rank, k);
            for (int q = 0; q < 255; ++q)
               printf("%i.", vertex_count[k][q]);
            printf("\n");

         }
      MPI_Barrier(MPI_COMM_WORLD);
      if(my_rank == 1)
         for(int k = 0 ; k<4;++k) {
            printf("Rankprexform %i: %i: ", my_rank, k);
            for (int q = 0; q < 255; ++q)
               printf("%i.", vertex_count[k][q]);
            printf("\n");
         }
      MPI_Barrier(MPI_COMM_WORLD);
      return NULL;
      // transform to cumulative sum
      for(int j=0;j < 4;++j)
         for (i = 1; i <= blocksize; ++i)
            vertex_count[j][i] += vertex_count[j][i - 1];


      MPI_Barrier(MPI_COMM_WORLD);
      if(my_rank == 0)
         for(int k = 0 ; k<4;++k) {
            printf("Rank %i: %i: ", my_rank, k);
            for (int q = 0; q < 255; ++q)
               printf("%i.", vertex_count[k][q]);
            printf("\n");

         }
      MPI_Barrier(MPI_COMM_WORLD);
      if(my_rank == 1)
         for(int k = 0 ; k<4;++k) {
            printf("Rank %i: %i: ", my_rank, k);
            for (int q = 0; q < 255; ++q)
               printf("%i.", vertex_count[k][q]);
            printf("\n");
         }
      MPI_Barrier(MPI_COMM_WORLD);
return NULL;
      MPI_Barrier(MPI_COMM_WORLD);
      printf("EXIT Linear Xform LOOP");

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
      MPI_Barrier(MPI_COMM_WORLD);
      if(my_rank == 0)
         for(int k = 0 ; k<4;++k) {
            printf("Rank %i: %i: ", my_rank, k);
            for (int q = 0; q < 10; ++q)
               printf("%i.", vertex_count[k][q]);
            printf("\n");

         }
      MPI_Barrier(MPI_COMM_WORLD);
      if(my_rank == 1)
         for(int k = 0 ; k<4;++k) {
            printf("Rank %i: %i: ", my_rank, k);
            for (int q = 0; q < 10; ++q)
               printf("%i.", vertex_count[k][q]);
            printf("\n");
         }
      MPI_Barrier(MPI_COMM_WORLD);
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
