#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <stdint.h>

#include "bitmap.h"

struct Bitmap* newBitmap( uint32_t size){


	
        struct Bitmap* bitmap = (struct Bitmap*) malloc( sizeof(struct Bitmap));
        bitmap->bitarray = (uint32_t*) malloc(sizeof(uint32_t)*((size+kBitsPerWord - 1)/kBitsPerWord));
	

      
    memset(bitmap->bitarray, 0, (sizeof(uint32_t)*((size+kBitsPerWord - 1)/kBitsPerWord)));
	bitmap->size =  size;
	bitmap->numSetBits =  0;

	return bitmap;
}


void freeBitmap( struct Bitmap* bitmap){

        free(bitmap->bitarray);
        free(bitmap);
	
}

void clearBitmap(struct Bitmap* bitmap){

	memset(bitmap->bitarray, 0, (sizeof(uint32_t)*((bitmap->size+kBitsPerWord - 1)/kBitsPerWord)));
	bitmap->numSetBits =  0;

}

void setBit(struct Bitmap* bitmap, uint32_t pos){

	bitmap->bitarray[word_offset(pos)] |= (uint32_t) (1 << bit_offset(pos));

}

void setBitRange(struct Bitmap* bitmap, uint32_t start,uint32_t end){

 uint32_t pos;

 for (pos = start; pos < end; ++pos)
 {
 	setBit(bitmap, pos);
 }

}



uint32_t getBit(struct Bitmap* bitmap, uint32_t pos){

	return (bitmap->bitarray[word_offset(pos)] >> bit_offset(pos)) & 1l;;

}


void clearBit(struct Bitmap* bitmap, uint32_t pos){

	bitmap->bitarray[word_offset(pos)] &= ((uint32_t) (~(1l << bit_offset(pos))));

}




void swapBitmaps (struct Bitmap** bitmap1, struct Bitmap** bitmap2){

	
	struct Bitmap* temp_bitmap = *bitmap1;
	*bitmap1 = *bitmap2;
	*bitmap2 = temp_bitmap;

}



uint32_t getNumOfSetBits (struct Bitmap* bitmap){

	uint32_t i;
	uint32_t numSetBits = 0;

	#pragma omp parallel for reduction(+:numSetBits) schedule(dynamic,256)
	for(i= 0 ; i < (bitmap->size); i++){
		if(getBit(bitmap, i))
			numSetBits++;
	}

	return numSetBits;
}

void printSetBits (struct Bitmap* bitmap){

	uint32_t i;

	for(i= 0 ; i < (bitmap->size); i++){
		if(getBit(bitmap, i)){
			printf("**%u \n", i);
		}
	}

}


/***********************************
you need to implement this function
***********************************/
void setBitAtomic(struct Bitmap* bitmap, uint32_t pos){

/*
	atomic operatoin here
*/
	
}