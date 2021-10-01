#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

#include "edgelist.h"
#include "graph.h"

inline int maxTwoIntegers(int num1, int num2) {
   if (num1 >= num2)
      return num1;
   else
      return num2;
}

void printEdgeArray(struct Edge *edgeArray, int numOfEdges) {
   int i;
   printf("\n");
   for (i = 0; i < numOfEdges; i++) {
      printf("%d->%d\n", edgeArray[i].src, edgeArray[i].dest);
   }
}


// initialize a new edge array
struct Edge *newEdgeArray(int numOfEdges) {
   struct Edge *edgeArray = (struct Edge *) malloc(numOfEdges * sizeof(struct Edge));
   int i;
   for (i = 0; i < numOfEdges; i++) {
      edgeArray[i].dest = 0;
      edgeArray[i].src = 0;
   }
   return edgeArray;
}

struct Graph *loadEdgeArray(const char *fname, struct Graph *graph) {
   FILE *pText;
   int size = 0, i;
   int src = 0, dest = 0;

   pText = fopen(fname, "r");

   if (pText == NULL) {
      return NULL;
   }

   while (1) {
      i = fscanf(pText, "%d\t%d\n", &src, &dest);
      if (i == EOF)
         break;

      graph->sorted_edges_array[size].src = src;
      graph->sorted_edges_array[size].dest = dest;
     /* if(size % 0xfffff == 0) {
         printf(".");
         fflush(stdout);
      }*/

      size++;
   }
   fclose(pText);
   return graph;
}

void loadEdgeArrayInfo(const char *fname, int *numOfVertices, int *numOfEdges) {
   FILE *pText;
   int size = 0, i;
   int src = 0, dest = 0, num_vertices = 0;
   pText = fopen(fname, "r");

   if (pText == NULL) {
      // err(1, "open: %s", fname);
      return;
   }
   while (1) {
      size++;
     /* if(size % 0xfffff == 0) {
         printf(".");
         fflush(stdout);
      }*/
      i = fscanf(pText, "%d\t%d\n", &src, &dest);
      if (i == EOF)
         break;
     /* if(src>dest)
         tmp = src;
      else
         tmp = dest;
      if(tmp>num_vertices)
         num_vertices = tmp;
*/
      num_vertices = maxTwoIntegers(num_vertices, maxTwoIntegers(src, dest));
   }
   fclose(pText);

   *numOfVertices = num_vertices + 1;
   *numOfEdges = size - 1;
  // printf("Number Vertices: %i\n", num_vertices+1);
   //printf("Number Edges: %i\n", size);
}

/*
void loadEdgeArrayInfo(const char *fname, int *numOfVertices, int *numOfEdges) {
   FILE *pText;
   int size = 0, i;
   int src = 0, dest = 0, num_vertices = 0, tmp=0;
   pText = fopen(fname, "r");

   if (pText == NULL) {
      // err(1, "open: %s", fname);
      return;
   }
   while (1) {
      size++;
      if(size % 0xfffff == 0) {
         printf(".");
         fflush(stdout);
      }
      i = fscanf(pText, "%d\t%d\n", &src, &dest);
      if (i == EOF)
         break;
      if(src>dest)
         tmp = src;
      else
         tmp = dest;
      if(tmp>num_vertices)
         num_vertices = tmp;

      //num_vertices = maxTwoIntegers(num_vertices, maxTwoIntegers(src, dest));
   }
   fclose(pText);

   *numOfVertices = num_vertices + 1;
   *numOfEdges = size - 1;
}*/


