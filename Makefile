
#########################################################
#       		 GENERAL DIRECTOIRES   	    			#
#########################################################
# globals binaary /bin
APP                 = run-graph

# dirs Root app
APP_DIR             = .

#dir root/managed_folders
SRC_DIR           	= src
OBJ_DIR_OMP			= obj-openmp
OBJ_DIR_MPI			= obj-mpi
OBJ_DIR_HYB			= obj-hybrid
INC_DIR			  	= include
BIN_DIR			  	= bin

##################################################
##################################################

#########################################################
#       		         GRAPH ARGUMENTS    			#
#########################################################

BENCHMARKS_DIR    	= /mnt/beegfs/smdupor/

#Tests
# GRAPH_NAME = wiki-vote
GRAPH_NAME = asdf
# GRAPH_NAME = test

# GONG # https://gonglab.pratt.duke.edu/google-dataset
# GRAPH_NAME = gplus

# GAP # https://sparse.tamu.edu/MM/GAP/
# GRAPH_NAME = road

# SNAP # https://snap.stanford.edu/data/
# GRAPH_NAME = SNAP-cit-Patents
# GRAPH_NAME = SNAP-com-Orkut

# KONECT # http://konect.cc/networks/wikipedia_link_en/
# export GRAPH_NAME = KONECT-wikipedia_link_en

FILE_BIN_TYPE = 64m.txt
FILE_BIN = $(BENCHMARKS_DIR)/$(FILE_BIN_TYPE)


NUM_THREADS   = 4
ROOT 		  = 1

ARGS = -n $(NUM_THREADS) -r $(ROOT) -h
##################################################
##################################################

CC		    = gcc
CC_MPI		= mpicc
APP_INC     = -I$(APP_DIR)/$(INC_DIR)

# flags
CFLAGS   =  -O3 -Wall -m64 -fopenmp -g
LFLAGS 	 = -lm 

##################################################
##################################################

##############################################
#                  COMPILATION VARIABLES     #
##############################################

SRC_FILES_ALGO			=  $(wildcard $(APP_DIR)/$(SRC_DIR)/*.c)
OBJ_FILES_ALGO_OMP      =  $(patsubst $(APP_DIR)/$(SRC_DIR)/%.c,$(APP_DIR)/$(OBJ_DIR_OMP)/%.o,$(SRC_FILES_ALGO))
OBJ_FILES_ALGO_MPI      =  $(patsubst $(APP_DIR)/$(SRC_DIR)/%.c,$(APP_DIR)/$(OBJ_DIR_MPI)/%.o,$(SRC_FILES_ALGO))
OBJ_FILES_ALGO_HYB      =  $(patsubst $(APP_DIR)/$(SRC_DIR)/%.c,$(APP_DIR)/$(OBJ_DIR_HYB)/%.o,$(SRC_FILES_ALGO))
INC_FILES_ALGO       	=  $(wildcard $(APP_DIR)/$(INC_DIR)/*.h)
ALL_HEADER_FILES        =   $(INC_FILES_ALGO) 

##################################################
##################################################


#########################################################
#       		 OPEN GRAPH LIBRARY     				#
#########################################################

.PHONY: all
all: app-openmp app-mpi app-hybrid

.PHONY: directories
directories :
	@mkdir -p $(APP_DIR)/$(BIN_DIR)
	@mkdir -p $(APP_DIR)/$(OBJ_DIR_OMP)
	@mkdir -p $(APP_DIR)/$(OBJ_DIR_MPI)
	@mkdir -p $(APP_DIR)/$(OBJ_DIR_HYB)

##################################################
##################################################

.PHONY: app-openmp
app-openmp: MODE = -DOPENMP_HARNESS
app-openmp : directories $(APP_DIR)/$(BIN_DIR)/$(APP)-openmp
	@echo "\n ******************************************************************************  "
	@echo " * DONE!! NOTHING ELSE TO COMPILE ---> openmp: ./$(word 2,$^)"
	@echo " ******************************************************************************  \n"

.PHONY: app-mpi
app-mpi: MODE = -DMPI_HARNESS
app-mpi : directories $(APP_DIR)/$(BIN_DIR)/$(APP)-mpi
	@echo "\n ******************************************************************************  "
	@echo " * DONE!! NOTHING ELSE TO COMPILE ---> mpi: ./$(word 2,$^)"
	@echo " ******************************************************************************  \n"

.PHONY: app-hybrid
app-hybrid: MODE = -DHYBRID_HARNESS
app-hybrid : directories $(APP_DIR)/$(BIN_DIR)/$(APP)-hybrid
	@echo "\n ******************************************************************************  "
	@echo " * DONE!! NOTHING ELSE TO COMPILE ---> hybrid: ./$(word 2,$^)"
	@echo " ******************************************************************************  \n"


$(APP_DIR)/$(BIN_DIR)/$(APP)-openmp :  $(OBJ_FILES_ALGO_OMP)
	@$(CC) $(CFLAGS) -o $@ $^  $(LFLAGS)

$(APP_DIR)/$(BIN_DIR)/$(APP)-mpi :  $(OBJ_FILES_ALGO_MPI)
	@$(CC_MPI) $(CFLAGS) -o $@ $^  $(LFLAGS)

$(APP_DIR)/$(BIN_DIR)/$(APP)-hybrid :  $(OBJ_FILES_ALGO_HYB)
	@$(CC_MPI) $(CFLAGS) -o $@ $^  $(LFLAGS)

$(APP_DIR)/$(OBJ_DIR_OMP)/%.o : $(APP_DIR)/$(SRC_DIR)/%.c $(INC_FILES_ALGO)
	$(CC) $(CFLAGS) $(APP_INC) $(MODE) -c -o $@ $<

$(APP_DIR)/$(OBJ_DIR_MPI)/%.o : $(APP_DIR)/$(SRC_DIR)/%.c $(INC_FILES_ALGO)
	$(CC_MPI) $(CFLAGS) $(APP_INC) $(MODE) -c -o $@ $<

$(APP_DIR)/$(OBJ_DIR_HYB)/%.o : $(APP_DIR)/$(SRC_DIR)/%.c $(INC_FILES_ALGO)
	$(CC_MPI) $(CFLAGS) $(APP_INC) $(MODE) -c -o $@ $<

##############################################
#              TOP LEVEL RULES               #
##############################################

.PHONY: run
run: run-openmp

.PHONY: clean-openmp
clean-openmp:
	@rm -fr $(APP_DIR)/$(OBJ_DIR_OMP)
	@rm -fr $(APP_DIR)/$(BIN_DIR)/$(APP)-openmp

.PHONY: clean-mpi
clean-mpi:
	@rm -fr $(APP_DIR)/$(OBJ_DIR_OMP)
	@rm -fr $(APP_DIR)/$(BIN_DIR)/$(APP)-mpi

.PHONY: clean-hybrid
clean-hybrid:
	@rm -fr $(APP_DIR)/$(OBJ_DIR_OMP)
	@rm -fr $(APP_DIR)/$(BIN_DIR)/$(APP)-hybrid

.PHONY: clean
clean: clean-all

.PHONY: clean-all
clean-all:
	@rm -fr $(APP_DIR)/$(OBJ_DIR_OMP)
	@rm -fr $(APP_DIR)/$(OBJ_DIR_MPI)
	@rm -fr $(APP_DIR)/$(OBJ_DIR_HYB)
	@rm -fr $(APP_DIR)/$(BIN_DIR)


.PHONY: run-openmp
run-openmp: app-openmp
	./$(APP_DIR)/$(BIN_DIR)/$(APP)-openmp -f $(FILE_BIN) $(ARGS)

.PHONY: run-mpi
run-mpi: app-mpi
	./$(APP_DIR)/$(BIN_DIR)/$(APP)-mpi -f $(FILE_BIN) $(ARGS)

.PHONY: run-hybrid
run-hybrid: app-hybrid
	./$(APP_DIR)/$(BIN_DIR)/$(APP)-hybrid -f $(FILE_BIN) $(ARGS)

.PHONY: debug-openmp
debug-openmp: app-openmp
	gdb -ex=r --args ./$(APP_DIR)/$(BIN_DIR)/$(APP)-openmp -f $(FILE_BIN) $(ARGS)

.PHONY: debug-mpi
debug-mpi: app-mpi
	gdb -ex=r --args ./$(APP_DIR)/$(BIN_DIR)/$(APP)-mpi -f $(FILE_BIN) $(ARGS)

.PHONY: debug-hybrid
debug-hybrid: app-hybrid
	gdb -ex=r --args ./$(APP_DIR)/$(BIN_DIR)/$(APP)-hybrid -f $(FILE_BIN) $(ARGS)
