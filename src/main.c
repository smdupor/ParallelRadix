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

    printf("Number of Threads: %i\n", numThreads);
    printf(" -----------------------------------------------------\n");
   Start(timer2);
    struct Timer* timer = (struct Timer*) malloc(sizeof(struct Timer));
   printf("initializing graphs...");
   fflush(stdout);
    struct Graph* graph = newGraph(fnameb);
   printf("\nloading array...");
   fflush(stdout);
    loadEdgeArray(fnameb, graph);
   printf("\nCopying graphs...");
   fflush(stdout);
    struct Graph* graph_ser = copyGraph(graph);
   printf("<1>");
   fflush(stdout);
 //   struct Graph* graph_omp = copyGraph(graph);
   printf("<2>");
   fflush(stdout);
   struct Graph* graph_hyb = copyGraph(graph);
   struct Graph* graph_mpi = graph;

   printf("<2>");
   fflush(stdout);

    // populate the edge array from file
    //loadEdgeArray(fnameb, graph);
  // loadEdgeArray(fnameb, graph_ser);
   //loadEdgeArray(fnameb, graph_omp);
   printf("\nDone. Sorting.\n-----------------------------------------------------\n");

   //char prefix[80] = "../";
   char prefix[80] =  "/mnt/beegfs/smdupor/";

    // Serial
  //  Start(timer);
 //   graph = countSortEdgesBySource(graph); // you need to parallelize this function
   //   Stop(timer);
   //file_dump(graph, strcat(prefix, "out-cspar"));
  // printMessageWithtime("Parallel Count Sort (Seconds)",Seconds(timer));

   Start(timer);
   graph_ser = serial_count_sort(graph_ser); // you need to parallelize this function
   Stop(timer);
  printMessageWithtime("Serial Count Sort (Seconds)",Seconds(timer));
  //file_dump(graph_ser, strcat(prefix,"out-csser"));
 //  validation_run(graph_ser, graph, 1);
 //  printf("validation complete\n");

   // OMP Radix
 //  Start(timer);
 //   graph_omp = radixSortEdgesBySourceOpenMP(graph_omp); // you need to parallelize this function
//   Stop(timer);
 //  printMessageWithtime("Radix Sorting OMP (Seconds)",Seconds(timer));
   //file_dump(graph_omp, strcat(prefix,"out-omp"));
 //  validation_run(graph_ser, graph_omp, 1);
 //  printf("validation complete\n");
//  freeGraph(graph_omp);
 //  freeGraph(graph);

   Start(timer);
  graph_mpi = radixSortEdgesBySourceMPI(graph_mpi); // you need to parallelize this function
   Stop(timer);
   if(graph_mpi != NULL) {
      printMessageWithtime("Radix Sorting MPI (Seconds)", Seconds(timer));
     // freeGraph(graph_mpi);
      int myrank;
      MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
      if(myrank == 0)
      if(validation_run(graph_ser, graph_mpi, 1) == 0)
         printf("VALIDATION PASS\n");
      else
         printf("VALIDATION FAIL\n");
   }
   else {
      printf("Slave process, exit.");
      //freeGraph(graph_mpi);
   }
   MPI_Barrier(MPI_COMM_WORLD);

   Start(timer);
   graph_hyb = radixSortEdgesBySourceHybrid(graph_hyb); // you need to parallelize this function
   Stop(timer);
   if(graph_hyb != NULL) {
      printMessageWithtime("Radix Sorting Hybrid (Seconds)", Seconds(timer));
      // freeGraph(graph_mpi);
      int myrank;
      MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
      if(myrank == 0)
         if(validation_run(graph_ser, graph_hyb, 1) == 0)
            printf("VALIDATION PASS\n");
         else
            printf("VALIDATION FAIL\n");
   }
   else {
      printf("Slave process, exit.");
      //freeGraph(graph_mpi);
   }
   MPI_Barrier(MPI_COMM_WORLD);

 freeGraph(graph_mpi);
   freeGraph(graph_ser);
 //Stop(timer2);
 //  printMessageWithtime("Total Sim Time: ", Seconds(timer2) );


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

