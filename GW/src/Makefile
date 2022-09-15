.PHONY = all clean


APP := G1_G2_main
All: $(APP)

CC=g++
DEFS=-DNOTHREADS=1
OPT=-O3
CFLAGS=$(OPT) $(DEFS) -fopenmp
LDFLAGS= $(DEFS) -fopenmp -lmkl_intel_lp64 -lmkl_gnu_thread -lpthread -lm -ldl -lmkl_core -O3

DEPS:=collision_integral.h  commutators.h  G1_G2_t_step.h  HF.h  print_matrix.h  psi.h  write_to_file.h  RK_solver.h pi.h

G1_G2_main:
        $(CC) Main.cpp $(LDFLAGS) $< -o G1_G2_main


clean:
        rm -rf *.o $(APP)
