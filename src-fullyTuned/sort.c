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

   Start(&timer[INIT]);
   Start(&timer[COUNT]);
   Start(&timer[CRUSH]);
   Start(&timer[XFORM]);
   Start(&timer[SORT]);
   Start(&timer[MPI_MSG]);
   Stop(&timer[MPI_MSG]);

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
   Stop(&timer[INIT]);
   // count occurrence of key: id of a source vertex
#pragma omp parallel for default(none) shared(vertex_count, graph) private(key, i) schedule(static, 256)
   for (i = 0; i < edges; ++i) {
      key = graph->sorted_edges_array[i].src;
#pragma omp atomic
      vertex_count[key]++;
   }
   Stop(&timer[COUNT]);
   Stop(&timer[CRUSH]);

   // transform to cumulative sum

   for (i = 1; i < graph->num_vertices; ++i) {
      vertex_count[i] += vertex_count[i - 1];
   }
   Stop(&timer[XFORM]);
   // fill-in the sorted array of edges
   for (i = graph->num_edges - 1; i >= 0; --i) {
      key = graph->sorted_edges_array[i].src;
      pos = vertex_count[key] - 1;
      sorted_edges_array[pos] = graph->sorted_edges_array[i];
      vertex_count[key]--;
   }
   Stop(&timer[SORT]);

   free(vertex_count);
   free(graph->sorted_edges_array);
   graph->sorted_edges_array = sorted_edges_array;
   return graph;






   /*
   *int i, key, pos;
   Start(&timer[INIT]);
   Start(&timer[COUNT]);
   Start(&timer[CRUSH]);
   Start(&timer[XFORM]);
   Start(&timer[SORT]);
   Start(&timer[MPI_MSG]);
   Stop(&timer[MPI_MSG]);

   const int thread_per_bitbucket = 4;
   omp_set_nested(1);
   omp_set_num_threads(thread_per_bitbucket);
   printf("ParalleCountsort Num THR: %i on OMP\n", thread_per_bitbucket);

   const int edges = graph->num_edges;
   const int blocksize = graph->num_vertices;
   const int all_blocks = blocksize * thread_per_bitbucket;

   printf("ALLBLOCKS %i bs %i edges %i\n", all_blocks, blocksize, edges);

   struct Edge *sorted_edges_array = newEdgeArray(graph->num_edges);
   struct Edge *temp;


   //int vertex_count_tmp[all_blocks];
   int * vertex_count_tmp = (int *) malloc((int) all_blocks*sizeof(int));


   int * vertex_count = (int *) malloc((int) blocksize*sizeof(int ));

   printf("Decl DONE.\n");
  // sleep(1);
#pragma omp parallel default(none) shared(vertex_count, vertex_count_tmp) num_threads(2)
   {
      if(omp_get_thread_num() == 0)
          for (int q = 0; q < all_blocks; ++q)
            vertex_count_tmp[q] = 0;
      if(omp_get_thread_num()==1)
        for (int q = 0; q < blocksize; ++q)
            vertex_count[q] = 0;
   }
printf("INIT DONE.\n");
   Stop(&timer[INIT]);

   int thr=0, thr_offset=0;
   int offset=0;
   int end=0;


#pragma omp parallel default(none) shared(vertex_count_tmp, graph) private(i, thr, thr_offset, offset, end, key) num_threads(thread_per_bitbucket)
      {
         thr = omp_get_thread_num();
         offset = edges/(omp_get_num_threads());
         end = offset + edges/(omp_get_num_threads());
         thr_offset = blocksize*thr;
     //    printf("thread %i thr_off %i start %i end %i \n ", thr, thr_offset, offset, end);


         if (((edges*omp_get_num_threads()) - end) < edges && ((edges*omp_get_num_threads()) - end) > 0)
            end = edges*omp_get_num_threads();
        // printf("thread %i thr_off %i start %i end %i \n ", thr, thr_offset, offset, end);
        // sleep(3);
         // count occurrence of key: id of a source vertex
       //  printf("Midsleep %i\n", thr);
         //if(thr<3)
        //    sleep(2);
       //  printf("Endsleep %i\n", thr);
         for (i = offset; i < end; ++i) {
        //      printf("i at kill was: %i of %i of %i on %i \n", i, end, edges, thr);

            key = graph->sorted_edges_array[i].src;
            vertex_count_tmp[key + thr_offset]++;

         }

      }
      printf("COUNTCRAP DONE\n");
   Stop(&timer[COUNT]);
   const int stop = blocksize;
   int joffset=0;
      for (int j = 0; j < thread_per_bitbucket; ++j){
         joffset = j*stop;
         for (i = 0; i < stop; ++i) {
            if(j==0)
               vertex_count[i] = vertex_count_tmp[i];
            else
               vertex_count[i] += vertex_count_tmp[i + joffset];
         }
      }
   printf("CRUSH IS DONE\n");
   Stop(&timer[CRUSH]);
free(vertex_count_tmp);


      for (i = 1; i < blocksize; ++i)
         vertex_count[i] += vertex_count[i - 1];

   Stop(&timer[XFORM]);
   printf("XFORM IS DONE\n");

      for (i = edges - 1; i >= 0; --i) {
         key = graph->sorted_edges_array[i].src;
         pos = --vertex_count[key];
         sorted_edges_array[pos] = graph->sorted_edges_array[i];
      }
      // Swap the dirty and clean arrays
   //   temp = graph->sorted_edges_array;
   printf("SORT1 IS DONE\n");

   free(vertex_count);
   free(graph->sorted_edges_array);
   graph->sorted_edges_array = sorted_edges_array;


   //free(graph->sorted_edges_array);
   printf("SORT 2IS DONE\n");
      graph->sorted_edges_array = sorted_edges_array;
   printf("SORT3 IS DONE\n");

   Stop(&timer[SORT]);
   printf("SORT 4IS DONE\n");
  // free(vertex_count);
   printf("SORT 5IS DONE\n");


   return graph;*/
}

struct Graph *radix_serial(struct Graph *graph, struct Timer *timer) {
   int i, key, pos;
   int max = 0, granularity = 8, bitmask = 0xff; // 8-bit buckets, as masked by 0xff

   Start(&timer[INIT]);
   Start(&timer[COUNT]);
   Stop(&timer[COUNT]);
   Start(&timer[CRUSH]);
   Stop(&timer[CRUSH]);
   Start(&timer[XFORM]);
   Stop(&timer[XFORM]);
   Start(&timer[SORT]);
   Start(&timer[MPI_MSG]);
   Stop(&timer[MPI_MSG]);

   const int edges = graph->num_edges;

   for (i = 0; i < edges; ++i) {
      if (max < graph->sorted_edges_array[i].src)
         max = graph->sorted_edges_array[i].src;
   }

   struct Edge *sorted_edges_array = newEdgeArray(graph->num_edges);
   struct Edge *temp;
   int *vertex_count = (int *) malloc((2 * bitmask) * sizeof(int));
   Stop(&timer[INIT]);

   for (int digits = 0; digits < 32; digits += granularity) {
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

   Stop(&timer[SORT]);
   return graph;
}

struct Graph *radixSortEdgesBySourceOpenMP(struct Graph *graph, struct Timer *timer) {
   int i, key, pos, digits=0;
   const int granularity = 8, bitmask = 0xff; // 8-bit buckets, as masked by 0xff
   const int iters = 32/granularity;
   Start(&timer[INIT]);
   Start(&timer[COUNT]);
   Start(&timer[CRUSH]);
   Start(&timer[XFORM]);
   Start(&timer[SORT]);
   Start(&timer[MPI_MSG]);
   Stop(&timer[MPI_MSG]);

   const int thread_per_bitbucket = 4;
   omp_set_nested(1);
   omp_set_num_threads(thread_per_bitbucket*iters);
   printf("OMP Num THR: %i on OMP RESOLUTION %i\n", thread_per_bitbucket*iters, granularity);

   const int edges = graph->num_edges;
   const int blocksize = bitmask + 1;
   const int all_blocks = blocksize + (thread_per_bitbucket*blocksize);

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


#pragma omp parallel for default(none) shared(vertex_count_tmp, graph, digits) private(i, thr, thr_offset, offset, end, key, row) num_threads(iters)
   for (digits = 0; digits < 32; digits = digits + granularity) {
#pragma omp parallel default(none) shared(vertex_count_tmp, graph, digits) private(i, thr, thr_offset, offset, end, key, row) num_threads(thread_per_bitbucket)
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
   Stop(&timer[CRUSH]);


   for(int j=0;j < iters;++j)
      for (i = 1; i < blocksize; ++i)
         vertex_count[j][i] += vertex_count[j][i - 1];

   Stop(&timer[XFORM]);


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

struct Graph *radixSortEdgesBySourceMPI(struct Graph *graph, struct Timer *timer) {
   Start(&timer[INIT]);
   Start(&timer[COUNT]);
   Start(&timer[CRUSH]);
   Start(&timer[XFORM]);
   Start((&timer[MPI_MSG]));
   Start(&timer[SORT]);
   int i, key, pos, digits;
      int granularity = 8, bitmask = 0xff; // 8-bit buckets, as masked by 0xff

      const int edges = graph->num_edges;

      int my_rank;
      MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

   //   if (my_rank !=0)
     //    my_rank = 1;

      struct Edge *sorted_edges_array = newEdgeArray(graph->num_edges);
      struct Edge *temp;

      int vertex_count[4][bitmask+1];

      for(int k = 0 ; k<4;++k)
         for(int q = 0; q<=bitmask+1;++q)
            vertex_count[k][q] = 0;
   Stop(&timer[INIT]);
   //   for (digits = my_rank * (2 * granularity);
     //      digits <= (my_rank * 2 * granularity) + granularity; digits = digits + granularity) {
   digits = my_rank * granularity;
         // zero Out count array
         for (i = 0; i < bitmask + 1; ++i) {
            vertex_count[(digits / granularity)][i] = 0;
         }
   Stop(&timer[COUNT]);
         // count occurrence of key: id of a source vertex
         for (i = 0; i < edges; ++i) {
            key = graph->sorted_edges_array[i].src;
            vertex_count[(digits / granularity)][(key >> digits) & (bitmask)]++;
         }
   Stop(&timer[CRUSH]);
         // transform to cumulative sum
         for (i = 1; i < bitmask + 1; ++i) {
            vertex_count[(digits / granularity)][i] += vertex_count[(digits / granularity)][i - 1];
         }
     // }


   Stop(&timer[XFORM]);

   if( my_rank == 0) {
      MPI_Recv(vertex_count[1], bitmask+1, MPI_INT, 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vertex_count[2], bitmask+1, MPI_INT, 2, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(vertex_count[3], bitmask+1, MPI_INT, 3, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
     // MPI_Ssend(vertex_count[0], bitmask+1, MPI_INT, 1, 0, MPI_COMM_WORLD);
      //MPI_Ssend(vertex_count[1], bitmask+1, MPI_INT, 1, 0, MPI_COMM_WORLD);
   }
   else {
      MPI_Ssend(vertex_count[my_rank], bitmask+1, MPI_INT, 0, my_rank, MPI_COMM_WORLD);
      //MPI_Ssend(vertex_count[3], bitmask+1, MPI_INT, 0, 0, MPI_COMM_WORLD);
      //MPI_Recv(vertex_count[0], bitmask+1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      //MPI_Recv(vertex_count[1], bitmask+1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

   }
   Stop((&timer[MPI_MSG]));
   return graph;

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
   Stop(&timer[SORT]);
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
