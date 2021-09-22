# globals
APP                = main

# dirs
APP_DIR           	= .
SRC_DIR           	= src
OBJ_DIR			  	= obj
INC_DIR			  	= include
BIN_DIR				= bin

# compilers
CC				  = gcc
# CC				  = mpicc

INC = 	-I$(APP_DIR)/$(INC_DIR)/
		
# flags
CFLAGS            = -O3 -Wall -m64 -fopenmp -s

.PHONY: all
all: test

$(APP_DIR)/$(OBJ_DIR)/$(APP).o: $(APP_DIR)/$(SRC_DIR)/$(APP).c directories
	@echo 'making $(APP) <- $(APP).o'
	@$(CC) $(CFLAGS) $(INC) -c -o $(APP_DIR)/$(OBJ_DIR)/$(APP).o $(APP_DIR)/$(SRC_DIR)/$(APP).c

$(APP_DIR)/$(OBJ_DIR)/edgelist.o: $(APP_DIR)/$(SRC_DIR)/edgelist.c $(APP_DIR)/$(INC_DIR)/edgelist.h directories
	@echo 'making $(APP) <- edgelist.o'
	@$(CC) $(CFLAGS) $(INC) -c -o $(APP_DIR)/$(OBJ_DIR)/edgelist.o $(APP_DIR)/$(SRC_DIR)/edgelist.c

$(APP_DIR)/$(OBJ_DIR)/sort.o: $(APP_DIR)/$(SRC_DIR)/sort.c $(APP_DIR)/$(INC_DIR)/sort.h directories
	@echo 'making $(APP) <- sort.o'
	@$(CC) $(CFLAGS) $(INC) -c -o $(APP_DIR)/$(OBJ_DIR)/sort.o $(APP_DIR)/$(SRC_DIR)/sort.c

$(APP_DIR)/$(OBJ_DIR)/vertex.o: $(APP_DIR)/$(SRC_DIR)/vertex.c $(APP_DIR)/$(INC_DIR)/vertex.h directories
	@echo 'making $(APP) <- vertex.o'
	@$(CC) $(CFLAGS) $(INC) -c -o $(APP_DIR)/$(OBJ_DIR)/vertex.o $(APP_DIR)/$(SRC_DIR)/vertex.c

$(APP_DIR)/$(OBJ_DIR)/timer.o: $(APP_DIR)/$(SRC_DIR)/timer.c $(APP_DIR)/$(INC_DIR)/timer.h directories
	@echo 'making $(APP) <- timer.o'
	@$(CC) $(CFLAGS) $(INC) -c -o $(APP_DIR)/$(OBJ_DIR)/timer.o $(APP_DIR)/$(SRC_DIR)/timer.c

$(APP_DIR)/$(OBJ_DIR)/graph.o: $(APP_DIR)/$(SRC_DIR)/graph.c $(APP_DIR)/$(INC_DIR)/graph.h directories
	@echo 'making $(APP) <- graph.o'
	@$(CC) $(CFLAGS) $(INC) -c -o $(APP_DIR)/$(OBJ_DIR)/graph.o $(APP_DIR)/$(SRC_DIR)/graph.c

$(APP_DIR)/$(OBJ_DIR)/bfs.o: $(APP_DIR)/$(SRC_DIR)/bfs.c $(APP_DIR)/$(INC_DIR)/bfs.h directories
	@echo 'making $(APP) <- bfs.o'
	@$(CC) $(CFLAGS) $(INC) -c -o $(APP_DIR)/$(OBJ_DIR)/bfs.o $(APP_DIR)/$(SRC_DIR)/bfs.c

$(APP_DIR)/$(OBJ_DIR)/arrayQueue.o: $(APP_DIR)/$(SRC_DIR)/arrayQueue.c $(APP_DIR)/$(INC_DIR)/arrayQueue.h directories
	@echo 'making $(APP) <- arrayQueue.o'
	@$(CC) $(CFLAGS) $(INC) -c -o $(APP_DIR)/$(OBJ_DIR)/arrayQueue.o $(APP_DIR)/$(SRC_DIR)/arrayQueue.c

$(APP_DIR)/$(OBJ_DIR)/bitmap.o: $(APP_DIR)/$(SRC_DIR)/bitmap.c $(APP_DIR)/$(INC_DIR)/bitmap.h directories
	@echo 'making $(APP) <- bitmap.o'
	@$(CC) $(CFLAGS) $(INC) -c -o $(APP_DIR)/$(OBJ_DIR)/bitmap.o $(APP_DIR)/$(SRC_DIR)/bitmap.c

.PHONY: app
app: $(APP_DIR)/$(OBJ_DIR)/$(APP).o

.PHONY: edgelist
edgelist: $(APP_DIR)/$(OBJ_DIR)/edgelist.o

.PHONY: sort
sort: $(APP_DIR)/$(OBJ_DIR)/sort.o

.PHONY: vertex
vertex: $(APP_DIR)/$(OBJ_DIR)/vertex.o

.PHONY: timer
timer: $(APP_DIR)/$(OBJ_DIR)/timer.o

.PHONY: bfs
bfs: $(APP_DIR)/$(OBJ_DIR)/bfs.o

.PHONY: graph
graph: $(APP_DIR)/$(OBJ_DIR)/graph.o

.PHONY: arrayQueue
arrayQueue: $(APP_DIR)/$(OBJ_DIR)/arrayQueue.o

.PHONY: bitmap
bitmap: $(APP_DIR)/$(OBJ_DIR)/bitmap.o

.PHONY: directories
directories :
	@mkdir -p $(APP_DIR)/$(BIN_DIR)
	@mkdir -p $(APP_DIR)/$(OBJ_DIR)


.PHONY: test
test: app edgelist sort vertex timer bfs graph bitmap arrayQueue
	@echo 'linking $(APP) <- $(APP).o edgelist.o sort.o vertex.o timer.o bfs.o graph.o'
	@$(CC) 	$(APP_DIR)/$(OBJ_DIR)/$(APP).o 			\
			$(APP_DIR)/$(OBJ_DIR)/edgelist.o 		\
			$(APP_DIR)/$(OBJ_DIR)/sort.o 			\
			$(APP_DIR)/$(OBJ_DIR)/vertex.o 			\
			$(APP_DIR)/$(OBJ_DIR)/timer.o  			\
			$(APP_DIR)/$(OBJ_DIR)/bfs.o 			\
			$(APP_DIR)/$(OBJ_DIR)/graph.o  			\
			$(APP_DIR)/$(OBJ_DIR)/arrayQueue.o  	\
			$(APP_DIR)/$(OBJ_DIR)/bitmap.o  		\
		 -o $(APP_DIR)/$(BIN_DIR)/$(APP)			\
			$(CFLAGS)

n=1
# r=3009230
# f=./datasets/RMAT/RMAT22
# r=6
# f=./datasets/test/test.txt
r=3120
f=./datasets/facebook/facebook_combined.txt

.PHONY: run
run: test
	#$(APP_DIR)/$(BIN_DIR)/$(APP) -f $(f) -n $(n) -r $(r) -h

.PHONY: debug
debug: test	
	gdb $(APP_DIR)/$(BIN_DIR)/$(APP)

.PHONY: clean
clean:
	@rm -rf  $(APP_DIR)/$(BIN_DIR)
	@rm -rf  $(APP_DIR)/$(OBJ_DIR)
	
