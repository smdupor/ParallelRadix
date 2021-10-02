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

int numThreads;

void printMessageWithtime(const char *msg, double time)
{

    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", msg);
    printf(" -----------------------------------------------------\n");
    printf("| %-51f | \n", time);
    printf(" -----------------------------------------------------\n");

}

static void usage(void)
{
    printf("\nUsage: ./main -f <graph file> -r [root] -n [num threads]\n");
    printf("\t-f <graph file.txt>\n");
    printf("\t-h [help]\n");
    printf("\t-r [root/source]: BFS \n");
    printf("\t-n [num threads] default:max number of threads the system has\n");
    // _exit(-1);
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


int main(int argc, char **argv)
{

#ifdef MPI_HARNESS
   int provided;
   MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &provided);
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
#endif

#ifdef HYBRID_HARNESS
   int provided;
   MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &provided);
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
#endif

    char *fvalue = NULL;
    char *rvalue = NULL;
    char *nvalue = NULL;

    int root = 0;


    numThreads = omp_get_max_threads();
    char *fnameb = NULL;

    int c;
    opterr = 0;

    while ((c = getopt (argc, argv, "f:r:n:h")) != -1)
    {
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
           // usage();
            return 1;
        default:
            abort ();
        }
    }


    //Set number of threads for the program
    omp_set_nested(1);
    omp_set_num_threads(numThreads);
   struct Timer *timer = (struct Timer *) malloc(sizeof(struct Timer));
   struct Graph *graph;
#ifdef OPENMP_HARNESS
    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "OPENMP Implementation");
#endif

#ifdef MPI_HARNESS

   if(my_rank == 0) {
      printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "MPI Implementation");

#endif

#ifdef HYBRID_HARNESS
   if(my_rank == 0) {
    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Hybrid (OMP+MPI) Implementation");

#endif

   printf(" -----------------------------------------------------\n");
   printf("| %-51s | \n", "File Name");
   printf(" -----------------------------------------------------\n");
   printf("| %-51s | \n", fnameb);
   printf(" -----------------------------------------------------\n");
   printf("| %-51s | \n", "Number of Threads");
   printf(" -----------------------------------------------------\n");
   printf("| %-51u | \n", numThreads);
   printf(" -----------------------------------------------------\n");


   printf(" -----------------------------------------------------\n");
   printf("| %-51s | \n", "New graph calculating size");
   printf(" -----------------------------------------------------\n");
   Start(timer);
   graph = newGraph(fnameb);
   Stop(timer);
   printMessageWithtime("New Graph Created", Seconds(timer));


   printf(" -----------------------------------------------------\n");
   printf("| %-51s | \n", "Populate Graph with edges");
   printf(" -----------------------------------------------------\n");
   Start(timer);
   // populate the edge array from file
   loadEdgeArray(fnameb, graph);
   Stop(timer);
   printMessageWithtime("Time load edges to graph (Seconds)", Seconds(timer));
#ifdef MPI_HARNESS
   }
    else {
       graph = newGraph(fnameb);
          loadEdgeArray(fnameb, graph);
    }
#endif

#ifdef HYBRID_HARNESS
   }
    else {
       graph = newGraph(fnameb);
          loadEdgeArray(fnameb, graph);
    }
#endif


   /////////////////////////////////////////////////////////////////////////
    struct Graph *copy = copyGraph(graph);
      copy = countSortEdgesBySource(copy);
///////////////////////////////////////////////////////////////////////
#ifdef OPENMP_HARNESS
   // you need to parallelize this function
   printf(" -----------------------------------------------------\n");
   printf("| %-51s | \n", "OPENMP RADIX Sort Graph");
   printf(" -----------------------------------------------------\n");
   Start(timer);
   graph = radixSortEdgesBySourceOpenMP(graph); // you need to parallelize this function
   Stop(timer);
   printMessageWithtime("Time Sorting (Seconds)", Seconds(timer));
#endif

#ifdef MPI_HARNESS
   // you need to parallelize this function
if(my_rank == 0) {
   printf(" -----------------------------------------------------\n");
   printf("| %-51s | \n", "MPI RADIX Sort Graph");
   printf(" -----------------------------------------------------\n");
}
   Start(timer);
   graph = radixSortEdgesBySourceMPI(graph); // you need to parallelize this function
   Stop(timer);
   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Finalize();
   if(my_rank > 0) {
      printf("*************** WORKER PROCESS %i EXITING\n", my_rank);
      if(graph)
         freeGraph(graph);
      return 0;
   }

   printMessageWithtime("Time Sorting (Seconds)", Seconds(timer));
#endif

#ifdef HYBRID_HARNESS
   // you need to parallelize this function
   if(my_rank == 0) {
   printf(" -----------------------------------------------------\n");
   printf("| %-51s | \n", "HYBRID OMP/MPI RADIX Sort Graph");
   printf(" -----------------------------------------------------\n");
   }
   Start(timer);
   graph = radixSortEdgesBySourceHybrid(graph); // you need to parallelize this function
   Stop(timer);
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Finalize();
      if(my_rank > 0) {
      printf("*************** WORKER PROCESS %i EXITING\n", my_rank);
      freeGraph(graph);
      return 0;
   }
   printMessageWithtime("Time Sorting (Seconds)", Seconds(timer));
#endif



    // For testing purpose.

    printf(" -----------------------------------------------------\n");
    printf("| %-51s | \n", "Map vertices to Edges");
    printf(" -----------------------------------------------------\n");
    Start(timer);
    mapVertices(graph);
    Stop(timer);
    printMessageWithtime("Time Mapping (Seconds)", Seconds(timer));

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
    printMessageWithtime("Time BFS (Seconds)", Seconds(timer));

    if(validation_run(copy, graph,1) ==0)
       printf("VALIDATION PASS");
    else
       printf("VALIDATION FAIL");

   //freeGraph(copy);

    Start(timer);
    freeGraph(graph);
    Stop(timer);
    printMessageWithtime("Free Graph (Seconds)", Seconds(timer));

    return 0;
}


