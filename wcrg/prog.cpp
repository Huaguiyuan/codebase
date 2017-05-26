//g++ -O3 -lgsl -o prog hamiltonian.o fermi_surface.o couplings.o loop.o prog.cpp
//g++ -std=c++11 -O3 -lgsl -lgslcblas -fopenmp -o prog hamiltonian.o fermi_surface.o couplings.o loop.o prog.cpp
#include "hamiltonian.hpp"
#include "brillouin_zone.hpp"
#include "symmetry.hpp"
#include "fermi_surface.hpp"
#include "particle.hpp"
#include "coupling_grid.hpp"
#include "loop.hpp"


int main () 
{
   brillouin_zone bz("square");
   hamiltonian hh(-1.47,0.0);
   fermi_surface fs(bz,hh,24);

   coupling_grid WW(800,800,bz,hh,fs);

   std::map<int, std::string> jobs = { {200,"loops_soc0_grid200_100915.dat"}, {300,"loops_soc0_grid300_100915.dat"}, 
                              				 {400,"loops_soc0_grid400_100915.dat"}, {500,"loops_soc0_grid500_100915.dat"}};

   std::map<int, std::string>::iterator it;

   int grid_res;
   std::string file;
   for (it = jobs.begin(); it != jobs.end(); it++)
	 {  
     grid_res = it->first;
		 file = it->second;
		 loop LL(1.0e-4,grid_res,grid_res,0,bz,hh,fs,WW);
   	 LL.set_grid_flag(true);
   	 LL.calc_loops();
   	 LL.print_loop_integrals(file.c_str());
	 }
}
