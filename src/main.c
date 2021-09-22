#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <omp.h>
#include <memory.h>

#include "graph.h"
#include "bfs.h"
#include "sort.h"
#include "edgelist.h" 
#include "vertex.h"
#include "timer.h"

int numThreads;

void file_dump(const struct Graph *graph, const char *fnp);

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
  char *rvalue = NULL;
  char *nvalue = NULL;
  int root = 0;

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
        rvalue = optarg;
        root = atoi(rvalue);
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

    printf("Number of Threads: %i\n", numThreads);
    printf(" -----------------------------------------------------\n");

    struct Timer* timer = (struct Timer*) malloc(sizeof(struct Timer));

    struct Graph* graph = newGraph(fnameb);
    struct Graph* graph_ser = newGraph(fnameb);
    struct Graph* graph_omp = newGraph(fnameb);

    printf("loading arrays...\n-----------------------------------------------------\n");
    // populate the edge array from file
    loadEdgeArray(fnameb, graph);
   loadEdgeArray(fnameb, graph_ser);
   loadEdgeArray(fnameb, graph_omp);

    // Serial
    Start(timer);
    graph = countSortEdgesBySource(graph); // you need to parallelize this function
      Stop(timer);
   file_dump(graph, "../out-cspar");
   freeGraph(graph);
   printMessageWithtime("Parallel Count Sort (Seconds)",Seconds(timer));

   Start(timer);
   graph_ser = serial_count_sort(graph_ser); // you need to parallelize this function
   Stop(timer);
   file_dump(graph_ser, "../out-csser");
   freeGraph(graph_ser);
   printMessageWithtime("Serial Count Sort (Seconds)",Seconds(timer));

   // OMP Radix
   Start(timer);
    graph_omp = radixSortEdgesBySourceOpenMP(graph_omp); // you need to parallelize this function
   Stop(timer);
   file_dump(graph_omp, "../out-omp");
   freeGraph(graph_omp);
   printMessageWithtime("Radix Sorting OMP (Seconds)",Seconds(timer));
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


