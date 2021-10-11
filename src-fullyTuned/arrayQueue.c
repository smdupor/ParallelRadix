#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <stdint.h>
#include <omp.h>

#include "arrayQueue.h"
#include "bitmap.h"

struct ArrayQueue *newArrayQueue(uint32_t size){

    struct ArrayQueue* arrayQueue = (struct ArrayQueue*) malloc( sizeof(struct ArrayQueue));
	

    arrayQueue->head = 0;
    arrayQueue->tail = 0;
    arrayQueue->tail_next = 0;
    arrayQueue->size = size;
  


    
    arrayQueue->queue = (uint32_t*) malloc(size*sizeof(uint32_t));
	

    arrayQueue->q_bitmap = newBitmap(size);

    arrayQueue->q_bitmap_next = newBitmap(size);

    return arrayQueue;

}


void resetArrayQueue(struct ArrayQueue *q){

	q->head = 0;
    q->tail = 0;
    q->tail_next = 0;
    clearBitmap(q->q_bitmap);

}

void freeArrayQueue(struct ArrayQueue *q){

	freeBitmap(q->q_bitmap_next);
	freeBitmap(q->q_bitmap);

	free(q->queue);
	free(q);

}



void enArrayQueueWithBitmap (struct ArrayQueue *q, uint32_t k){

	q->queue[q->tail] = k;
	setBit(q->q_bitmap, k);
	q->tail = q->tail_next;
	q->tail++;
	q->tail_next++;

}


void enArrayQueueDelayed (struct ArrayQueue *q, uint32_t k){

	q->queue[q->tail_next] = k;
	q->tail_next++;

}

void enArrayQueueDelayedWithBitmap (struct ArrayQueue *q, uint32_t k){

	q->queue[q->tail_next] = k;
	setBit(q->q_bitmap_next, k);
	q->tail_next++;

}


void slideWindowArrayQueue (struct ArrayQueue *q){

		q->head = q->tail;
		q->tail = q->tail_next;
	
}

void slideWindowArrayQueueBitmap (struct ArrayQueue *q){

	q->head = q->tail;
	q->tail = q->tail_next;
	swapBitmaps(&q->q_bitmap, &q->q_bitmap_next);
	clearBitmap(q->q_bitmap_next);

}

uint32_t deArrayQueue(struct ArrayQueue *q){

	uint32_t k = q->queue[q->head];
	clearBit(q->q_bitmap,k);
	q->head = (q->head+1)%q->size;

	return k;

}


uint32_t frontArrayQueue (struct ArrayQueue *q){

	uint32_t k = q->queue[q->head];

	return k;

}

uint8_t isEmptyArrayQueueCurr (struct ArrayQueue *q){

  if((q->tail > q->head))
  	return 0;
  else
  	return 1;

}

uint8_t isEmptyArrayQueue (struct ArrayQueue *q){

  if(!isEmptyArrayQueueCurr(q) || !isEmptyArrayQueueNext(q))
  	return 0;
  else
  	return 1;

}

uint8_t isEmptyArrayQueueNext (struct ArrayQueue *q){

  if((q->tail_next > q->head))
  	return 0;
  else
  	return 1;

}

uint8_t  isEnArrayQueued 	(struct ArrayQueue *q, uint32_t k){


	return getBit(q->q_bitmap, k);

}

uint8_t  isEnArrayQueuedNext 	(struct ArrayQueue *q, uint32_t k){


	return getBit(q->q_bitmap_next, k);

}

uint32_t sizeArrayQueueCurr(struct ArrayQueue *q){

	return q->tail - q->head;

}

uint32_t sizeArrayQueueNext(struct ArrayQueue *q){

	return q->tail_next - q->tail;
}


uint32_t sizeArrayQueue(struct ArrayQueue *q){

	return q->tail_next - q->head;

}


/*

	Flush only works for single thread make sure is multi threaded safe

*/

void flushArrayQueueToShared(struct ArrayQueue *local_q, struct ArrayQueue *shared_q){

uint32_t shared_q_tail_next = shared_q->tail_next;
	shared_q->tail_next += local_q->tail;
uint32_t local_q_size = local_q->tail - local_q->head;

//reset the local queue state the ruse purpose is alloc/dealoc local queus create overhead.

memcpy(&shared_q->queue[shared_q_tail_next],&local_q->queue[local_q->head],local_q_size*(sizeof(uint32_t)));

	local_q->head = 0;
    local_q->tail = 0;
    local_q->tail_next = 0;

}


/*

you need to implement this if needed

*/
void enArrayQueueAtomic (struct ArrayQueue *q, uint32_t k){

/*

	atomic operation

*/
	
}



void enArrayQueue (struct ArrayQueue *q, uint32_t k){

	q->queue[q->tail] = k;
	q->tail = (q->tail+1)%q->size;
	q->tail_next = q->tail;

}
