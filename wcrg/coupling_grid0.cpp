#include "coupling_grid.hpp"

coupling_grid::coupling_grid(int nx, int ny, brillouin_zone bz, hamiltonian H, fermi_surface fs) : bz_(bz), fs_(fs), H_(H)
{
	grid_nx = nx;
	grid_ny = ny;
	b_pars_ = fs_.get_band_params();
	calc_grid(grid_nx, grid_ny);
  sym_ = sym_rep();
}

void coupling_grid::calc_grid(int n1, int n2)
{
	int band_lim = b_pars_.band_num;
	int spin_lim = b_pars_.spin_num;

  // calculates enclosing rectangular, kx_lo ..etc can also be hard coded 
  const Polygon &bz_zone = bz_.get_bz();
  double kx_lo = bz_zone.left_vertex()->x();
  double kx_up = bz_zone.right_vertex()->x();
  double ky_up = bz_zone.top_vertex()->y();
  double ky_lo = bz_zone.bottom_vertex()->y();
  //-------

	Vector3d kk;
	double kx, ky;
	MatrixXcd mat;
	en_grid = new VectorXd*[n1];
	states_grid = new VectorXcd***[n1];
   	for ( int i = 0; i < n1; i++ )
   	{
     	kx = kx_lo + (i + 0.5)*(kx_up - kx_lo)/n1;
   		en_grid[i]     = new VectorXd[n2];
  		states_grid[i] = new VectorXcd**[n2];
      for ( int j = 0; j < n2; j++ )
      {
      	ky = ky_lo + (j + 0.5)*(ky_up - ky_lo)/n2;
				kk << kx, ky, 0.0;
        en_grid[i][j] = H_.get_ham_eigval(kk);
				mat = H_.get_ham_eigmat(kk);
				states_grid[i][j] = new VectorXcd*[band_lim];
   			for ( int b = 0; b < band_lim; b++ )
   			{
					states_grid[i][j][b] = new VectorXcd[spin_lim];
      		for ( int s = 0; s < spin_lim; s++ )
					{
						states_grid[i][j][b][s] = mat.col(fs_.ind(b,s));	
					}
				}
      }
   	}
	
}

cmplx coupling_grid::onsite_interaction(VectorXcd u1, VectorXcd u2, VectorXcd u3, VectorXcd u4) const
{
  int dim = H_.get_dim();

	double U_intra = 0.3;
	double U_inter = 0.0;
  double J_h     = 0.0;
  double J_p     = 0.0;

	cmplx ww = 0.0;
  int orb_num = dim/2;

  // intra orbital
  for (int i = 0; i < orb_num; i++)
	{
	 	ww += U_intra*( conj(u1(2*i))*conj(u2(2*i + 1))*u3(2*i)*u4(2*i + 1)
	                 -conj(u1(2*i+1))*conj(u2(2*i))*u3(2*i)*u4(2*i + 1)
                   -conj(u1(2*i))*conj(u2(2*i + 1))*u3(2*i+1)*u4(2*i)
                   +conj(u1(2*i+1))*conj(u2(2*i))*u3(2*i+1)*u4(2*i) );
	}
	return ww;
}

cmplx coupling_grid::VV_exct(const particle &p1, const particle &p2, const particle &p3, const particle &p4) const
{
	cmplx res;
	
	VectorXcd u1 = get_state_exct(p1);
	VectorXcd u2 = get_state_exct(p2);
	VectorXcd u3 = get_state_exct(p3);
	VectorXcd u4 = get_state_exct(p4);

	res = onsite_interaction(u1, u2, u3, u4);
	return res;
}

cmplx coupling_grid::VV_grid(const particle &p1, const particle &p2, const particle &p3, const particle &p4) const
{
  cmplx res;

	VectorXcd u1 = get_state_grid(p1);
	VectorXcd u2 = get_state_grid(p2);
	VectorXcd u3 = get_state_grid(p3);
	VectorXcd u4 = get_state_grid(p4);

	res = onsite_interaction(u1, u2, u3, u4);
	return res;
}

VectorXcd coupling_grid::get_state_exct(const particle &pp) const
{	
	int patch   = pp.get_patch();
  Vector3d kk = pp.get_kk();
	int bix     = pp.get_bix();
	int sp      = pp.get_sp();
	if (patch == -1)
	{
		Vector3d fold_vec = bz_.mapin(kk);
    return H_.get_ham_eigvec(fold_vec, fs_.ind(bix,sp));
	}
	else
	{
		return fs_.get_state(patch,bix,sp);
	}
}

Vector2i coupling_grid::get_grid_point(Vector3d kk) const
{
  Vector2i ik;
	ik(0) = (int) rint((kk(0) + 1.0*M_PI)*grid_nx/(2.0*M_PI) - 0.5);       
  ik(1) = (int) rint((kk(1) + 1.0*M_PI)*grid_ny/(2.0*M_PI) - 0.5);       
  if (ik(0) >= grid_nx) ik(0) = grid_nx - 1;
  if (ik(0) < 0) ik(0) = 0;
  if (ik(1) >= grid_ny) ik(1) = grid_ny - 1;
  if (ik(1) < 0) ik(1) = 0;
  return ik;
}

VectorXcd coupling_grid::get_state_grid(const particle &pp) const
{	
	int patch   = pp.get_patch();
  Vector3d kk = pp.get_kk();
	int bix     = pp.get_bix();
	int sp      = pp.get_sp();

	if (patch == -1)
	{
		Vector3d fold_vec = bz_.mapin(kk);
    Vector2i ki = get_grid_point(fold_vec);
   	return states_grid[ki(0)][ki(1)][bix][sp];
	}
	else
	{
		return fs_.get_state(patch,bix,sp);
	}
}

void coupling_grid::print_VV(const char* filename, int flag) const
{
   ofstream file_output(filename);

   int band_lim = b_pars_.band_num;
   int spin_lim = b_pars_.spin_num;

   int npatch = fs_.get_npatch(); 
   int num = 5;
   int width1 = 6;
   int width2 = 18;
   file_output.precision(10); 
   cout.precision(10);
   
   particle p1(fs_,H_), p2(fs_,H_), p3(fs_,H_), p4(fs_,H_);
   Vector3d p4_fold;

   int i1, i2, i3;
   for (int bix1 = 0; bix1 < band_lim; bix1++)
   {
    for (int sp1 = 0; sp1 < spin_lim; sp1++)
    { 
     for (int m1 = 0; m1 < num; m1++)
     { 
      i1 = rand() % npatch;
	    p1.set_particle(i1,bix1,sp1);
      for (int bix2 = 0; bix2 < band_lim; bix2++)
      {
       for (int sp2 = 0; sp2 < spin_lim; sp2++)
       { 
        for (int m2 = 0; m2 < num; m2++)
        { 
         i2 = rand() % npatch;
	       p2.set_particle(i2,bix2,sp2);
         for (int bix3 = 0; bix3 < band_lim; bix3++)
         {
          for (int sp3 = 0; sp3 < spin_lim; sp3++)
          { 
           for (int m3 = 0; m3 < num; m3++)
           { 
            i3 = rand() % npatch;
	          p3.set_particle(i3,bix3,sp3);
       	    for (int bix4 = 0; bix4 < band_lim; bix4++)
            {
             for (int sp4 = 0; sp4 < spin_lim; sp4++)
             { 
						 p4_fold = bz_.mapin(p1.get_kk() + p2.get_kk() - p3.get_kk());
	           p4.set_particle(p4_fold, bix4, sp4);
	           file_output << setw(width1) << right << bix1;
   	         file_output << setw(width1) << right << sp1;
   	         file_output << setw(width1) << right << i1;
	           file_output << setw(width1) << right << bix2;
   	         file_output << setw(width1) << right << sp2;
   	         file_output << setw(width1) << right << i2;
	           file_output << setw(width1) << right << bix3;
   	         file_output << setw(width1) << right << sp3;
   	         file_output << setw(width1) << right << i3;
	           file_output << setw(width1) << right << bix4;
   	         file_output << setw(width1) << right << sp4;
	           if (flag == 0)
	           {
	            file_output << setw(width2) << right << VV_exct(p1,p2,p3,p4).real();
	            file_output << setw(width2) << right << VV_exct(p1,p2,p3,p4).imag();
	           }
	           else
	           {
	            file_output << setw(width2) << right << VV_grid(p1,p2,p3,p4).real();
	            file_output << setw(width2) << right << VV_grid(p1,p2,p3,p4).imag();
	           }
	            file_output << "\n";
	          }
           }
	        }
         }
	      }
	     }
      }
     }
    }
   }
  }
  file_output.close();
}

void coupling_grid::VV_check_sym(const char* filename) const
{
   ofstream file_output(filename);

   int band_lim = b_pars_.band_num;
   int spin_lim = b_pars_.spin_num;

   int npatch = fs_.get_npatch(); 
   int num = 5;
   int width1 = 6;
   int width2 = 18;
   file_output.precision(6); 
   cout.precision(10);
   
   particle p1(fs_,H_),  p2(fs_,H_),  p3(fs_,H_),  p4(fs_,H_);
   particle p1s(fs_,H_), p2s(fs_,H_), p3s(fs_,H_), p4s(fs_,H_);
   Vector3d p4_fold, p4s_fold;
   
   const Matrix3d_vec &k_ops = sym_.get_k_ops();

   cmplx err, ch1, ch2, ch3, ch4;

   int i1, i2, i3;
   int i1s, i2s, i3s;

	 for (int op = 0; op < sym_.get_sym_num(); op++)
   {
   	for (int bix1 = 0; bix1 < band_lim; bix1++)
   	{
     for (int sp1 = 0; sp1 < spin_lim; sp1++)
     { 
      for (int m1 = 0; m1 < num; m1++)
      { 
       i1  = rand() % npatch;
	     p1.set_particle(i1,bix1,sp1);
			 i1s = fs_.sym_patch(i1, bix1, sp1, op);
	     p1s.set_particle(i1s,bix1,sp1);
       for (int bix2 = 0; bix2 < band_lim; bix2++)
       {
        for (int sp2 = 0; sp2 < spin_lim; sp2++)
        { 
         for (int m2 = 0; m2 < num; m2++)
         { 
          i2 = rand() % npatch;
	        p2.set_particle(i2,bix2,sp2);
			    i2s = fs_.sym_patch(i2, bix2, sp2, op);
	        p2s.set_particle(i2s,bix2,sp2);
          for (int bix3 = 0; bix3 < band_lim; bix3++)
          {
           for (int sp3 = 0; sp3 < spin_lim; sp3++)
           { 
            for (int m3 = 0; m3 < num; m3++)
            { 
             i3 = rand() % npatch;
	           p3.set_particle(i3,bix3,sp3);
			       i3s = fs_.sym_patch(i3, bix3, sp3, op);
	           p3s.set_particle(i3s,bix3,sp3);
       	     for (int bix4 = 0; bix4 < band_lim; bix4++)
             {
              for (int sp4 = 0; sp4 < spin_lim; sp4++)
              { 
	  					 p4_fold = bz_.mapin(p1.get_kk() + p2.get_kk() - p3.get_kk());
	             p4.set_particle(p4_fold, bix4, sp4);
						   //p4s_fold = bz_.mapin(p1s.get_kk() + p2s.get_kk() - p3s.get_kk());
						   p4s_fold = k_ops[op]*p4_fold;
	             p4s.set_particle(p4s_fold, bix4, sp4);

               //characters not here in purple bronze model
							 //transformation behavior is put in by hand
							 //ch1 = fs_.get_char(i1, bix1, sp1, op);
							 //ch2 = fs_.get_char(i2, bix2, sp2, op);
							 //ch3 = fs_.get_char(i3, bix3, sp3, op);
							 //ch4 = fs_.get_char(p4_fold, bix4, sp4, op);
							 //err = VV_exct(p1s,p2s,p3s,p4s) - conj(ch1)*conj(ch2)*ch3*ch4*VV_exct(p1,p2,p3,p4);

							 switch (op)
							 {
							 	case 0:
							 		err = VV_exct(p1s,p2s,p3s,p4s) - VV_exct(p1,p2,p3,p4);
									break;
							 	case 1:
							 		err = VV_exct(p1s,p2s,p3s,p4s) - conj(VV_exct(Sp(p1),Sp(p2),Sp(p3),Sp(p4)));
									break;
							 	case 2:
							 		err = VV_exct(p1s,p2s,p3s,p4s) - conj(VV_exct(p1,p2,p3,p4));
									break;
							 	case 3:
							 		err = VV_exct(p1s,p2s,p3s,p4s) - VV_exct(Sp(p1),Sp(p2),Sp(p3),Sp(p4));
									break;
               }

							 if (abs(err) > 1e-8)
							 {
  	            file_output << setw(width1) << right << op;
	              file_output << setw(width2) << right << err;
                //file_output << setw(width2) << right << p4_fold;
                file_output << setw(width2) << right << VV_exct(p1s,p2s,p3s,p4s);
                file_output << setw(width2) << right << VV_exct(p1,p2,p3,p4);
                file_output << setw(width2) << right << ch1;
                file_output << setw(width2) << right << ch2;
                file_output << setw(width2) << right << ch3;
                file_output << setw(width2) << right << ch4;
	              //file_output << setw(width1) << right << bix1;
   	            //file_output << setw(width1) << right << sp1;
   	            //file_output << setw(width1) << right << i1;
	              //file_output << setw(width1) << right << bix2;
   	            //file_output << setw(width1) << right << sp2;
   	            //file_output << setw(width1) << right << i2;
	              //file_output << setw(width1) << right << bix3;
   	            //file_output << setw(width1) << right << sp3;
   	            //file_output << setw(width1) << right << i3;
	              //file_output << setw(width1) << right << bix4;
   	            //file_output << setw(width1) << right << sp4;
	              file_output << "\n";
							 }
	            }
             }
	          }
           }
	        }
	       }
        }
       }
      }
     }
    }
	 }
   file_output.close();
}

