#ifndef BFS_H
#define BFS_H

#include "edgelist.h" 
#include "vertex.h"
#include "graph.h"
#include "bitmap.h"
#include "arrayQueue.h"

#define TRUE 1
#define FALSE 0

//breasth first search algorithm
void breadthFirstSearchGraphPush(int source, struct Graph* graph);
int topDownStepGraphCSR(struct Graph* graph, struct ArrayQueue* sharedFrontierQueue, struct ArrayQueue** localFrontierQueues);


extern int numThreads;

#endif