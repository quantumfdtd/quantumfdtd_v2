TARGET = mpisolve
OBJ = intde2.o random.o paramreader.o grid.o outputroutines.o initialconditions.o specialfunctions.o potential.o mpisolve.o

PROFILE = #-pg
DEBUG = -Wno-deprecated #-W -Wall -g
INCLUDES = -I include
OPTIMIZATION = -O2 -funroll-loops -finline-functions
CXXFLAGS = $(DEBUG) $(PROFILE) $(OPTIMIZATION) $(FLOWTRACE) $(INCLUDES)
MPIFLAGS = $(CXXFLAGS)
CC = mpicxx

all: $(TARGET)

mpisolve: $(OBJ)

run:	
	mpirun -np 3 mpisolve

run1:   
	mpirun -np 9 mpisolve

run2:   
	mpirun -np 17 mpisolve

run4:   
	mpirun -np 33 mpisolve

run8:   
	mpirun -np 65 mpisolve

cleandata:
	./scripts/cleandatafiles.sh

clean: cleandata
	rm -f *\.o *~

mrproper: clean
	rm -f $(TARGET)

.c.o:
	$(CC) $(MPIFLAGS) -c $*.c

help:
	@echo 'Build targets:'
	@echo '  all          - Builds mpisolve standalone'
	@echo '  solve        - Builds mpisolve standalone'
	@echo 'Cleanup targets:'
	@echo '  clean        - Remove generated files + Emacs leftovers'
	@echo '  mrproper     - Removes generated targets + all aboves'
	@echo 'Exec targets:'
	@echo '  run          - Run mpisolve (with first cleaning up data)'
