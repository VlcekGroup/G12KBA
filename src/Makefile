

.PHONY = all clean


APP := G1_G2_main
All: $(APP)

CC=CC
#DEFS=-DNOTHREADS=1
#DEBUG=-g -O0 
OPT=-O3 -fast
CFLAGS= $(DEBUG) $(OPT) $(DEFS) -fopenmp 
LDFLAGS= -lm  $(DEBUG) $(OPT) $(DEFS) -fopenmp -lpthread -lm -ldl 

DEPS:=collision_integral.h  commutators.h  G1_G2_t_step.h  HF.h  print_matrix.h  psi.h  write_to_file.h

%.o: %.cpp $(DEPS)
	$(CC) -c $(CFLAGS)  $<

G1_G2_main: Main.o
	$(CC) $(LDFLAGS) $< -o $@


clean:
	rm -rf *.o $(APP)
