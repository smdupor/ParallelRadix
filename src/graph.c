#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>

#include "graph.h"
#include "edgelist.h"
#include "vertex.h"


void   loadGraphInfo(const char * fname, struct Graph *graph){
		 loadEdgeArrayInfo(fname, &graph->num_vertices, &graph->num_edges);
}

// initialize a new graph from file
struct Graph * newGraph(const char * fname){
	int i;
	struct Graph* graph = (struct Graph*) malloc(sizeof(struct Graph));

	loadGraphInfo(fname, graph);

	graph->parents  = (int*) malloc( graph->num_vertices *sizeof(int));
    graph->vertices = newVertexArray(graph->num_vertices);
    graph->sorted_edges_array = newEdgeArray(graph->num_edges);

    for(i = 0; i < graph->num_vertices; i++){
        graph->parents[i] = -1;
    }

    graph->iteration = 0;
    graph->processed_nodes = 0;
	return graph;
}

struct Graph * copyGraph(const struct Graph * orig){
   //int i;

   struct Graph* graph = (struct Graph*) malloc(sizeof(*orig));

   graph->num_vertices = orig->num_vertices;
   graph->num_edges = orig->num_edges;

   graph->parents = orig->parents;
   graph->vertices = newVertexArray(graph->num_vertices);
   graph->sorted_edges_array = newEdgeArray(graph->num_edges);

   for(int i=0;i<graph->num_vertices;++i){
      graph->vertices[i].edges_idx = orig->vertices[i].edges_idx;
      graph->vertices[i].in_degree = orig->vertices[i].in_degree;
      graph->vertices[i].out_degree = orig->vertices[i].out_degree;
   }

   for(int i=0; i<graph->num_edges;++i){
      graph->sorted_edges_array[i].dest = orig->sorted_edges_array[i].dest;
      graph->sorted_edges_array[i].src = orig->sorted_edges_array[i].src;
   }

   graph->iteration = orig->iteration;
   graph->processed_nodes = orig->processed_nodes;

   return graph;
}


void freeGraph(struct Graph *graph){
    free(graph->vertices);
    free(graph->sorted_edges_array);
	free(graph);
}