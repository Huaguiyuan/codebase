CC      = g++
CFLAGS  = -Wall -O3 -std=c++11 -fopenmp -frounding-math
LDFLAGS = -lgsl -lgslcblas -lCGAL

\OBJ = hamiltonian.o fermi_surface.o brillouin_zone.o symmetry.o particle.o coupling_grid.o loop.o test_prog.o

OBJ = hamiltonian.o fermi_surface.o brillouin_zone.o symmetry.o particle.o coupling_grid.o loop.o prog.o

prog: $(OBJ) 
	$(CC) $(CFLAGS) $(OBJ) $(LDFLAGS) -o wcrg_prog 

%.o: %.cpp
	$(CC) $(CFLAGS) -c $< 
