#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <omp.h>
#include <mpi.h>
#include <memory.h>

#include "graph.h"
#include "bfs.h"
#include "sort.h"
#include "edgelist.h" 
#include "vertex.h"
#include "timer.h"

int numThreads;



void file_dump(const struct Graph *graph, const char *fnp);
int validation_run(const struct Graph *serial, const struct Graph *test, short control);

void printMessageWithtime(const char * msg, double time){
    /*printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", msg);
    printf(" -----------------------------------------------------\n");
    printf("| %-51f | \n", time);
    printf(" -----------------------------------------------------\n");*/
    printf("%s :\t %f\n-----------------------------------------------------\n", msg, time);
}

static void usage(void) {
  printf("\nUsage: ./main -f <graph file> -r [root] -n [num threads]\n");
  printf("\t-f <graph file.txt>\n");
  printf("\t-h [help]\n");
  printf("\t-r [root/source]: BFS \n");
  printf("\t-n [num threads] default:max number of threads the system has\n");
  // _exit(-1);
}


int main(int argc, char **argv) {
  char *fvalue = NULL;
  //char *rvalue = NULL;
  char *nvalue = NULL;
  //int root = 0;

   int provided;
   MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &provided);

  numThreads = omp_get_max_threads();

  char *fnameb = NULL;
  int c;
  opterr = 0;

  while ((c = getopt (argc, argv, "f:r:n:h")) != -1){
    switch (c)
      {
      case 'h':
        usage();
        break;
      case 'f':
        fvalue = optarg;
        fnameb = fvalue;
        break;
      case 'r':
        //rvalue = optarg;
      //  root = atoi(rvalue);
        break;
       break;
      case 'n':
        nvalue = optarg;
        numThreads = atoi(nvalue);
        break;
      case '?':
        if (optopt == 'f')
          fprintf (stderr, "Option -%c <graph file> requires an argument  .\n", optopt);
        else if (optopt == 'r')
          fprintf (stderr, "Option -%c [root] requires an argument.\n", optopt);
        else if (optopt == 'n')
          fprintf (stderr, "Option -%c [num threads] requires an argument.\n", optopt);
        else if (isprint (optopt))
          fprintf (stderr, "Unknown option `-%c'.\n", optopt);
        else
          fprintf (stderr,
                   "Unknown option character `\\x%x'.\n",
                   optopt);
        usage();
        return 1;
      default:
        abort ();
      }
    }
    //Set number of threads for the program
    omp_set_nested(1);
    omp_set_num_threads(numThreads);
   struct Timer* timer2 = (struct Timer*) malloc(sizeof(struct Timer));
   struct Timer *timers = (struct Timer*) malloc(7*sizeof(struct Timer));

  /* int myrank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   if(myrank==0) {
*/

      //printf("Number of Threads: %i\n", numThreads);
     // printf(" -----------------------------------------------------\n");
      Start(timer2);
      struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));
      printf("initializing graphs...\n");
      fflush(stdout);
      struct Graph *graph = newGraph(fnameb);
      printf("loading array...\n");
      fflush(stdout);
      loadEdgeArray(fnameb, graph);
      printf("Copying graphs...\n");
      fflush(stdout);
      struct Graph *graph_ser = copyGraph(graph);
      printf("<1>");
      fflush(stdout);
      struct Graph *graph_omp = graph;
      printf("<2>");
      fflush(stdout);
     // struct Graph *graph_hyb = graph;
      //struct Graph* graph_mpi = graph;

      printf("<2>");
      fflush(stdout);

      // populate the edge array from file
      //loadEdgeArray(fnameb, graph);
      // loadEdgeArray(fnameb, graph_ser);
      //loadEdgeArray(fnameb, graph_omp);
      printf("\nDone. Sorting. Edges: %i Vertices: %i\n-----------------------------------------------------\n", graph->num_edges, graph->num_vertices);

      //char prefix[80] = "../";
      //char prefix[80] =  "/mnt/beegfs/smdupor/";

      Start(timer);
      graph_ser = serial_count_sort(graph_ser, timers); // you need to parallelize this function
      Stop(timer);
      printMessageWithtime("Serial Count Sort (Seconds)", Seconds(timer));
      //file_dump(graph_ser, strcat(prefix,"out-csser"));
      //  validation_run(graph_ser, graph, 1);
      //  printf("validation complete\n");

   graph_omp = countSortEdgesBySource(graph_omp, timers); // you need to parallelize this function
   Stop(timer);
   printf("Countsort Sorting OMP on thrQry %i (Seconds): %f \n",  numThreads, Seconds(timer));
   printf("Init: %f Count: %f Crush: %f Xform: %f MPImsg: %f Sort: %f\n", Millisecs(&timers[INIT]),
          Millisecs(&timers[COUNT]),
          Millisecs(&timers[CRUSH]), Millisecs(&timers[XFORM]), Millisecs(&timers[MPI_MSG]),
          Millisecs(&timers[SORT]));
   //file_dump(graph_omp, strcat(prefix,"out-omp"));
   if (validation_run(graph_ser, graph_omp, 2) == 0)
      printf("validation PASS\n");
   else
      printf("validation FAIL\n");

   freeGraph(graph_omp);

      // OMP Radix
      /*Start(timer);
      graph_omp = radixSortEdgesBySourceOpenMP(graph_omp, timers); // you need to parallelize this function
      Stop(timer);
      printMessageWithtime("Radix Sorting OMP (Seconds)", Seconds(timer));
      printf("Init: %f Count: %f Crush: %f Xform: %f MPImsg: %f Sort: %f\n", Millisecs(&timers[INIT]),
             Millisecs(&timers[COUNT]),
             Millisecs(&timers[CRUSH]), Millisecs(&timers[XFORM]), Millisecs(&timers[MPI_MSG]),
             Millisecs(&timers[SORT]));
      //file_dump(graph_omp, strcat(prefix,"out-omp"));
      if (validation_run(graph_ser, graph_omp, 1) == 0)
         printf("validation PASS\n");
      freeGraph(graph_omp);*/

      // Serial Radix
/*   Start(timer);
   graph_omp = radix_serial(graph_omp, timers); // you need to parallelize this function
   Stop(timer);
   printMessageWithtime("Radix Sorting OMP (Seconds)", Seconds(timer));
   printf("Init: %f Count: %f Crush: %f Xform: %f MPImsg: %f Sort: %f\n", Millisecs(&timers[INIT]),
          Millisecs(&timers[COUNT]),
          Millisecs(&timers[CRUSH]), Millisecs(&timers[XFORM]), Millisecs(&timers[MPI_MSG]),
          Millisecs(&timers[SORT]));
   //file_dump(graph_omp, strcat(prefix,"out-omp"));
   if (validation_run(graph_ser, graph_omp, 1) == 0)
      printf("validation PASS\n");
   else
      printf("validation FAIL\n");

   freeGraph(graph_omp);
*/

      //  freeGraph(graph);
/*
        Start(timer);
       graph_mpi = radixSortEdgesBySourceMPI(graph_mpi, timers); // you need to parallelize this function
        Stop(timer);
   int my_rank;
   MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
          if (graph_mpi != NULL && my_rank==0) {
            printMessageWithtime("Radix Sorting MPI (Seconds)", Seconds(timer));
            printf("InitCountCrushXformMPImsgSort: %f, %f, %f, %f, %f, %f \n", Millisecs(&timers[INIT]),
                   Millisecs(&timers[COUNT]),
                   Millisecs(&timers[CRUSH]), Millisecs(&timers[XFORM]), Millisecs(&timers[MPI_MSG]),
                   Millisecs(&timers[SORT]));
               if (validation_run(graph_ser, graph_mpi, 1) == 0)
                  printf("VALIDATION PASS\n");
               else
                  printf("VALIDATION FAIL\n");

         } else {
            printf("Slave process, exit, Radix MPI.\n");

             printf("InitCountCrushXformMPImsgSort: %f, %f, %f, %f, %f \n", Millisecs(&timers[INIT]),
                    Millisecs(&timers[COUNT]),
                    Millisecs(&timers[CRUSH]), Millisecs(&timers[XFORM]), Millisecs(&timers[MPI_MSG]));
         }*/
      //if (myrank == 1) {
     /*    MPI_Barrier(MPI_COMM_WORLD);

         Start(timer);
         graph_hyb = radixSortEdgesBySourceHybrid(graph_hyb, timers); // you need to parallelize this function
         Stop(timer);
         if (graph_hyb != NULL) {
            printMessageWithtime("Radix Sorting Hybrid (Seconds)", Seconds(timer));
            printf("InitCountCrushXformMPImsgSort: %f, %f, %f, %f, %f, %f \n", Millisecs(&timers[INIT]),
                   Millisecs(&timers[COUNT]),
                   Millisecs(&timers[CRUSH]), Millisecs(&timers[XFORM]), Millisecs(&timers[MPI_MSG]),
                   Millisecs(&timers[SORT]));
            int myrank;
            MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
            if (myrank == 0) {
               if (validation_run(graph_ser, graph_hyb, 1) == 0)
                  printf("VALIDATION PASS\n");
               else
                  printf("VALIDATION FAIL\n");
            }
         } else {
            printf("Slave process, exit.");
         }
         MPI_Barrier(MPI_COMM_WORLD);
     // }*/
      //freeGraph(graph_mpi);
      if (graph_ser)
        freeGraph(graph_ser);
      //if (graph_hyb)
      //   freeGraph(graph_hyb);
     // if (graph)
    //      freeGraph(graph);

         Stop(timer2);
      printMessageWithtime("Total Sim Time: ", Seconds(timer2));

   //}///////////// TODO MPI DISABLATION
   MPI_Finalize();

    return 0;
}

void file_dump(const struct Graph *graph, const char *fnp) {
   FILE *fpp;
   fpp = fopen(fnp, "w");
   for(int i = 0; i < graph->num_edges - 1 ; ++i) {
      //printf("Edge: %i  -  %i\n", graph->sorted_edges_array[i].src, graph->sorted_edges_array[i].dest);
      fprintf(fpp,"Edge: %i  -  %i\n", graph->sorted_edges_array[i].src, graph->sorted_edges_array[i].dest);
   }
   fclose(fpp);
}

int validation_run(const struct Graph *serial, const struct Graph *test, short control) {
   // Control = 0 NO printing
   // Control = 1 Print SUMS
   // Control = 2 Verbose
   int success = -1;
   switch(control){
      case 0:
         if(serial->num_edges != test->num_edges)
            return success;
         success = 0;
         for(int i=0; i< serial->num_edges;++i){
            if(serial->sorted_edges_array[i].src != test->sorted_edges_array[i].src &&
                     serial->sorted_edges_array[i].dest != test->sorted_edges_array[i].dest) {
               ++success;
            }
         }
         return success;
      case 1:
         if(serial->num_edges != test->num_edges) {
          printf("Edge Quantity Failure: Input: %i  Output: %i\n", serial->num_edges, test->num_edges);
            return success;
         }
         success = 0;
         for(int i=0; i< serial->num_edges;++i){
            if(serial->sorted_edges_array[i].src != test->sorted_edges_array[i].src &&
               serial->sorted_edges_array[i].dest != test->sorted_edges_array[i].dest) {
               ++success;
            }
         }
         if(success>0)
         printf("Interal Graph Failures Qty: %i",success);
         return success;
      case 2:
         if(serial->num_edges != test->num_edges) {
            printf("Edge Quantity Failure: Input: %i  Output: %i\n", serial->num_edges, test->num_edges);
            return success;
         }
         success = 0;
         for(int i=0; i< serial->num_edges;++i){
            if(serial->sorted_edges_array[i].src != test->sorted_edges_array[i].src &&
               serial->sorted_edges_array[i].dest != test->sorted_edges_array[i].dest) {
               ++success;
               printf("Internal Graph Failure on Row: %i, Contents %i - %i  VS %i - %i\n", i,serial->sorted_edges_array[i].src ,
                      test->sorted_edges_array[i].src , serial->sorted_edges_array[i].dest, test->sorted_edges_array[i].dest);
            }
         }
         if(success>0)
         printf("Interal Graph Failures Qty: %i\n",success);
         return success;
   }
   return success;
}

