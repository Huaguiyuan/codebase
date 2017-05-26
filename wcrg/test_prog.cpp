//g++ -O3 -lgsl -o prog hamiltonian.o fermi_surface.o couplings.o loop.o prog.cpp
//g++ -std=c++11 -O3 -lgsl -lgslcblas -fopenmp -o prog hamiltonian.o fermi_surface.o couplings.o loop.o prog.cpp
#include "brillouin_zone.hpp"
#include "symmetry.hpp"
#include "hamiltonian.hpp"
#include "fermi_surface.hpp"
#include "particle.hpp"
#include "coupling_grid.hpp"
#include "loop.hpp"


int main () 
{
   hamiltonian pb(-1.47,0.1);

   cout << "test hamiltonian class" << endl;
   Vector3d kk(1.0,2.0,3.0); 
   cout << pb.get_ham_mat(kk) << endl;
   cout << "---------------" << endl;
   cout << pb.get_ham_eigmat(kk) << endl;
   cout << "---------------" << endl;
	 cout << pb.get_ham_eigvec(kk, 0) << endl;
   cout << "---------------" << endl;
   cout << pb.get_ham_eigval(kk, 0) << endl;
   cout << "---------------" << endl;
   cout << pb.get_ham_eigval(kk) << endl;
   cout << "---------------" << endl;
   cout << pb.get_dim() << endl;
   cout << "---------------" << endl;
   //double tp[7] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
   //cout << "---------------" << endl;
   //pb.set_model_pars(0.0, tp, 0);

   brillouin_zone bz("square");

   cout << "test brillouin_zone class" << endl;
   double ang = 4.057890511;
   Vector3d pt(0.0,0.0,0.0);
   cout << bz.get_U_pt(ang) << endl;
   cout << "---------------" << endl;
   cout << bz.get_K_pt(ang) << endl;
   cout << "---------------" << endl;
   cout << bz.get_patch_len(ang) << endl;
   cout << "---------------" << endl;
   cout << bz.r2xy(0.4,ang,pt) << endl;
   cout << "---------------" << endl;
   cout << bz.dr2dx(0.4,ang) << endl;
   cout << "---------------" << endl;
   cout << bz.dr2dx(2.4,ang) << endl;
   cout << "---------------" << endl;
   cout << bz.get_bz() << endl;
   cout << "---------------" << endl;
   cout << bz.get_um() << endl;
   bz.print_discrt("bz_data.dat",24,40);
   
   cout << "test symmetry class" << endl;
   kk << 1.0, 2.0, 3.0; 
   sym_rep sym;
   const Matrix3d_vec &k_ops = sym.get_k_ops();
   cout << k_ops[1];
   cout << "---------------" << endl;
   sym.sym_test(kk,pb);
   //---------------------------------------------------
   cout << "test fermi_surface class" << endl;
   fermi_surface fs(bz,pb,24);
   cout << "---------------" << endl;
   fs.print_fs("test_fs.dat");
   cout << "---------------" << endl;
   fs.print_states("test_states.dat");
   cout << "---------------" << endl;
   fs.print_sym_ops("test_sym_ops.dat");
   //---------------------------------------------------
   cout << "test particle class" << endl;
   Vector3d k1, k2, k3;
   k1 << -2.607318021, 2.460566224, 0.0;
   k2 <<  3.089218476, -0.52559591, 0.0;
   k3 <<  1.324518476, -1.52559591, 0.0;
   cout << "before particle creation" << endl;
   particle p1(0,0,0,fs,pb); 
   particle p2(10,0,0,fs,pb); 
   particle p3(k3,0,0,fs,pb); 
   cout << p1 << endl;
   cout << Tinv(p1) << endl;
   cout << Pinv(p1) << endl;
   particle p5 = p3;
   p1 = p3;
   cout << p5 << endl;
   cout << "---------------" << endl;
   cout << p1.get_patch() << endl;
   cout << "---------------" << endl;
   cout << p1.get_kk() << endl;
   cout << "---------------" << endl;
   cout << p1.get_bix() << endl;
   cout << "---------------" << endl;
   cout << p1.get_sp() << endl;
   cout << "---------------" << endl;
   cout << p1.get_state() << endl;;
   cout << "---------------" << endl;
 	 cout << p1.get_en() << endl;
   cout << "---------------" << endl;
 	 cout << p1 << endl;
   cout << "---------------" << endl;
   for (int op = 0; op < sym.get_sym_num(); op++)
   {
 	 		cout << Uop(p1,op) << endl;
 	 		cout << charc(p1,op) << endl;
	 }
   //---------------------------------------------------
   cout << "test coupling grid class" << endl;
   coupling_grid WW(100,100,bz,pb,fs);
   cout << "---------------" << endl;
   WW.VV_check_sym("test_VV_sym.dat");
   cout << "---------------" << endl;
   WW.print_VV("test_VV_exct.dat",0);
   cout << "---------------" << endl;
   WW.print_VV("test_VV_grid.dat",1);
   //---------------------------------------------------
   cout << "test loop class" << endl;
   //loop LL(1.0e-4,100,100,0,bz,pb,fs,WW);
   loop LL(1.0e-4,100,100,0,bz,pb,fs,WW);
   LL.set_grid_flag(false);
   LL.print_loops("test_loops_2b_exact_func2.dat",0);
   LL.print_loops("test_loops_2e_exact_func2.dat",3);
   LL.calc_loops();
   LL.print_loop_integrals("loop_result.dat");
}
