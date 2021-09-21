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
    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", msg);
    printf(" -----------------------------------------------------\n");
    printf("| %-51f | \n", time);
    printf(" -----------------------------------------------------\n");
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
/*
    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "File Name");
    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", fnameb);
    printf(" -----------------------------------------------------\n");
*/    printf("| %-51s | \n", "Number of Threads");
    printf(" -----------------------------------------------------\n");
    printf("| %-51u | \n", numThreads);
    printf(" -----------------------------------------------------\n");

    struct Timer* timer = (struct Timer*) malloc(sizeof(struct Timer));
  /*
    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "New graph calculating size");
    printf(" -----------------------------------------------------\n");
    Start(timer);*/
    struct Graph* graph = newGraph(fnameb);
    struct Graph* graph_ser = newGraph(fnameb);
 /*   Stop(timer);
    printMessageWithtime("New Graph Created",Seconds(timer));

    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Populate Graph with edges");
    printf(" -----------------------------------------------------\n");
   * Start(timer);
   */ // populate the edge array from file
    loadEdgeArray(fnameb, graph);
   loadEdgeArray(fnameb, graph_ser);
   // Stop(timer);
    //printMessageWithtime("Time load edges to graph (Seconds)",Seconds(timer));

// you need to parallelize this function
    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "COUNT Sort Graph");
    printf(" -----------------------------------------------------\n");


   char *fns = "out-csser";
   char *fnp = "out-cspar";

    // Serial
    Start(timer);
    graph = countSortEdgesBySource(graph); // you need to parallelize this function
    // graph = radixSortEdgesBySource(graph); // you need to parallelize this function
    Stop(timer);
   file_dump(graph, "../out-cspar");

   printMessageWithtime("Parallel Sorting (Seconds)",Seconds(timer));

   Start(timer);
   graph_ser = serial_count_sort(graph_ser); // you need to parallelize this function
   // graph = radixSortEdgesBySource(graph); // you need to parallelize this function
   Stop(timer);
   file_dump(graph_ser, "../out-csser");

   printMessageWithtime("Serial Sorting (Seconds)",Seconds(timer));

    /*
    // For testing purpose.
    
    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Map vertices to Edges");
    printf(" -----------------------------------------------------\n");
    Start(timer);
    mapVertices(graph);
    Stop(timer);
    printMessageWithtime("Time Mapping (Seconds)",Seconds(timer));

    printf(" *****************************************************\n");
    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "BFS Algorithm (PUSH/PULL)");
    printf(" -----------------------------------------------------\n");

    printf("| %-51s | \n", "PUSH");
  
    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "ROOT/SOURCE");
    printf(" -----------------------------------------------------\n");
    printf("| %-51u | \n", root);
    printf(" -----------------------------------------------------\n");
    Start(timer);
    
    breadthFirstSearchGraphPush(root, graph);

    // for (int p = 0; p < 10; ++p)
    // {
    //   breadthFirstSearchGraphPush(p, graph); 
    // }

    Stop(timer);
    printMessageWithtime("Time BFS (Seconds)",Seconds(timer));

    printf("CONTENTS OF GRAPH FOR DEBUG:\n");



  //  Start(timer);
    freeGraph(graph);
   // Stop(timer);
    //printMessageWithtime("Free Graph (Seconds)",Seconds(timer));*/






    return 0;
}

void file_dump(const struct Graph *graph, const char *fnp) {
   FILE *fpp;
   fpp = fopen(fnp, "w");
   for(int i = 0; i <graph->num_edges - 1 ; ++i) {
      //printf("Edge: %i  -  %i\n", graph->sorted_edges_array[i].src, graph->sorted_edges_array[i].dest);
      fprintf(fpp, "Edge: %i  -  %i\n", graph->sorted_edges_array[i].src, graph->sorted_edges_array[i].dest);
   }
   close(fpp);
}


