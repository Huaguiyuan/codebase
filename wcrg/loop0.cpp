#include "loop.hpp"
#include "bz_funcs.hpp"

loop::loop(double T, int nx, int ny, int int_flag, const brillouin_zone &bz, const hamiltonian &H, const fermi_surface &fs, const couplings &cp): bz_(bz), fs_(fs), H_(H), cp_(cp)
{
	int_grid_nx_ = nx;
	int_grid_ny_ = ny;
  integration_flag_ = int_flag; // 0 for kartesian, 1 for radial along r2xy paths
	temperature_ = T;
	tflag_       = false;
	grid_flag_   = false;
  b_pars_      = fs_.get_band_params();
  npatch_      = fs_.get_npatch();
  sb_lim_      = bpars_.spin_lim;
  //sb_lim_      = bpars_.band_lim; number of extra degrees of freedom on a outer diagram line (depends on model)
	initialize_loops();
}

void loop::set_tflag(bool tflag)
{
	tflag_ = tflag;
}

void loop::set_grid_flag(bool grid_flag)
{
	grid_flag_ = grid_flag;
}

double loop::n_F(const particle &p1)
{
	double occ;
	double en;
	en = p1.get_en();
	if (tflag_) 
	{
		if (en >= 0.0)
		{
			occ = 0.0;
		}
		else
		{
			occ = 1.0;
		}	
	}
	else
	{
		occ = 1.0/(1.0 + exp(en/temperature_));
	}
	return occ;
}

double loop::frequency_integral(const particle &p1, const particle &p2)
{
	double  res, numer, denom;
	numer = n_F(p1) - n_F(p2);
	if (tflag_)
	{
		if (numer == 0.0)
		{
			res = 0.0;
		}
		else
		{
			denom = p1.get_en() - p2.get_en();
			res = numer/denom;
		}
	}
	else
	{
		denom = p1.get_en() - p2.get_en();
		if (abs(denom) > 1.0e-10)
		{
			res = numer/denom;
		}
		else 
		{
			res = (-1.0/temperature_)/(2.0 + 2.0*cosh(p1.get_en()/temperature_));
		}

	}
	return res;
}

cmplx loop::loop_1a_func(const particle &p1, const particle &p2)
{
	return cp_.VV_exct(p1,Tinv(p1),p2,Tinv(p2));
}

cmplx loop::loop_2b_func(const particle &p1, const particle &p2, const particle &p3, int bix4, int sp4)
{
	cmplx V1, V2;

  Vector3d p4_fold = bz::mapin(p3.get_kk() - p1.get_kk() - p2.get_kk());
  particle p4(p4_fold,bix4,sp4,fs_,H_);

	double loop = frequency_integral(p3, p4);

	if (abs(loop) > 1.0e-10)
	{
		if (grid_flag_) 
		{
			V1 = cp_.VV_grid(p1,p4,p3,Tinv(p2));
			V2 = cp_.VV_grid(p3,Tinv(p1),p2,p4);
		}
		else
		{
			V1 = cp_.VV_exct(p1,p4,p3,Tinv(p2));
			V2 = cp_.VV_exct(p3,Tinv(p1),p2,p4);
		}
		return loop*V1*V2;
	}
	else
	{
		return 0.0;
	}
}

cmplx loop::loop_2c_func(const particle &p1, const particle &p2, const particle &p3, int bix4, int sp4)
{
	cmplx V1, V2;

  Vector3d p4_fold = bz::mapin(p3.get_kk() - p1.get_kk() + p2.get_kk());
  particle p4(p4_fold,bix4,sp4,fs_,H_);

	double loop = frequency_integral(p3, p4);

	if (abs(loop) > 1.0e-10)
	{
		if (grid_flag_) 
		{
			V1 = cp_.VV_grid(p1,p4,p3,p2);
			V2 = cp_.VV_grid(p3,Tinv(p1),p4,Tinv(p2));
		}
		else
		{
			V1 = cp_.VV_exct(p1,p4,p3,p2);
			V2 = cp_.VV_exct(p3,Tinv(p1),p4,Tinv(p2));
		}
		return loop*V1*V2;
	}
	else
	{
		return 0.0;
	}
}

cmplx loop::loop_2d_func(const particle &p1, const particle &p2, const particle &p3, int bix4, int sp4)
{
	cmplx V1, V2;

  Vector3d p4_fold = bz::mapin(p3.get_kk() + p1.get_kk() - p2.get_kk());
  particle p4(p4_fold,bix4,sp4,fs_,H_);

	double loop = frequency_integral(p3, p4);

	if (abs(loop) > 1.0e-10)
	{
		if (grid_flag_) 
		{
			V1 = cp_.VV_grid(p1,p3,p2,p4);
			V2 = cp_.VV_grid(Tinv(p1),p4,p3,Tinv(p2));
		}
		else
		{
			V1 = cp_.VV_exct(p1,p3,p2,p4);
			V2 = cp_.VV_exct(Tinv(p1),p4,p3,Tinv(p2));
		}
		return loop*V1*V2;
	}
	else
	{
		return 0.0;
	}
}

cmplx loop::loop_2e_func(const particle &p1, const particle &p2, const particle &p3, int bix4, int sp4)
{
	cmplx V1, V2;

  Vector3d p4_fold = bz::mapin(p3.get_kk() - p1.get_kk() + p2.get_kk());
  particle p4(p4_fold,bix4,sp4,fs_,H_);

	double loop = frequency_integral(p3, p4);

	if (abs(loop) > 1.0e-10)
	{
		if (grid_flag_) 
		{
			V1 = cp_.VV_grid(p1,p4,p2,p3);
			V2 = cp_.VV_grid(Tinv(p1),p3,Tinv(p2),p4);
		}
		else
		{
			V1 = cp_.VV_exct(p1,p4,p2,p3);
			V2 = cp_.VV_exct(Tinv(p1),p3,Tinv(p2),p4);
		}
		return loop*V1*V2;
	}
	else
	{
		return 0.0;
	}
}

cmplx loop::loop_int(const particle &p1, const particle &p2, cmplx (loop::*func) (const particle &, const particle &, const particle &, int, int))
{
   int band_lim = b_pars_.band_num;
   int spin_lim = b_pars_.spin_num;

   cmplx res;
   double q3x,q3y;
   Vector3d q3;
   double ang, r, rlim;
   Vector3d pt;
   particle p3(fs_,H_);

   double dinf = (2*M_PI/int_grid_n1_)*(2*M_PI/int_grid_n2_)/(4*M_PI*M_PI); 
   pt << 0.0, 0.0, 0.0;

   for (int bix3 = 0; bix3 < band_lim; bix3++)
   {
    for (int sp3 = 0; sp3 < spin_lim; sp3++)
    {
     for (int bix4 = 0; bix4 < band_lim; bix4++)
     {
      for (int sp4 = 0; sp4 < spin_lim; sp4++)
      {
       for (int ii = 0; ii < int_grid_n1_; ii++)
       {
				switch (integration_flag_)
				{
					case 0:
        		q3x  = -M_PI + (ii - 0.5)*2*M_PI/int_grid_n1_;
						break;
					case 1;
		    		ang  = (2*ii + 1)*M_PI/int_grid_n1_; 
		    		rlim = bz_.get_patch_len(ang);
						break;
				}
        for (int jj = 0; jj < int_grid_n2_; jj++)
        {
					switch (integration_flag_)
					{
						case 0:
         			q3y  = -M_PI + (jj - 0.5)*2*M_PI/int_grid_n2_;
         			q3 << q3x, q3y, 0.0;
							break;
						case 1:
			   			r    = jj*rlim/int_grid_n2_;
			   			q3   = r2xy(r,ang,pt);
							dinf = bz_.dr2dx(r,ang);
							break;
					}
         p3.set_particle(q3,bix3,sp3);
         res = res + dinf*(this->*func)(p1, p2, p3, bix4, sp4);
        }
       }
      }
     }
    }
   }
   return res;
}	

cmplx loop::loop_int(const particle &p1, const particle &p2, cmplx (loop::*func) (const particle &, const particle &, const particle &, int, int), int sp3, int sp4)
{
   int band_lim = b_pars_.band_num;

   cmplx res;
   double q3x,q3y;
   Vector3d q3;
   double ang, r, rlim;
   Vector3d pt;
   particle p3(fs_,H_);

   double dinf = (2*M_PI/int_grid_nx_)*(2*M_PI/int_grid_ny_)/(4*M_PI*M_PI); 
   pt << 0.0, 0.0, 0.0;

   for (int bix3 = 0; bix3 < band_lim; bix3++)
   {
    for (int bix4 = 0; bix4 < band_lim; bix4++)
    {
     for (int ii = 0; ii < int_grid_n1_; ii++)
     {
			switch (integration_flag_)
			{
				case 0:
        	q3x  = -M_PI + (ii - 0.5)*2*M_PI/int_grid_n1_;
					break;
				case 1;
		   		ang  = (2*ii + 1)*M_PI/int_grid_n1_; 
		   		rlim = bz_.get_patch_len(ang);
					break;
			}
      for (int jj = 0; jj < int_grid_n2_; jj++)
      {
				switch (integration_flag_)
				{
					case 0:
       			q3y  = -M_PI + (jj - 0.5)*2*M_PI/int_grid_n2_;
       			q3 << q3x, q3y, 0.0;
						break;
					case 1:
			   		r    = jj*rlim/int_grid_n2_;
			 			q3   = r2xy(r,ang,pt);
						dinf = bz_.dr2dx(r,ang);
						break;
				}
        p3.set_particle(q3,bix3,sp3);
        res = res + dinf*(this->*func)(p1, p2, p3, bix4, sp4);
      }
     }
    }
   }
   return res;
}	

void loop::loop_main()
{
   int i1, i1s, i2s;

   int_vec::const_iterator ic;
   int_vec red_ind;

   particle p1(fs_,H_), p2(fs_,H_);

   // calculate symmetry reduced diagrams
   for (int sb1 = 0; sb1 < sb_lim_; sb1++)
   {
    red_ind = fs_.sym_red_ind(sb1,0);
    for (ic = red_ind.begin(); ic!=red_ind.end(); ic++)
    {
     i1 = *ic;
     cout << (npatch_*sb1 + i1)/(sb_lim_*npatch_)*100 << "% calculated" << endl;
     for (int sb2 = 0; sb2 < sb_lim_; sb2++)
     {
      for (int i2  = 0; i2 < npatch_; i2++)
      {

       // purple bronze case
       //p1.set_particle(i1,sb1,0);
       //p2.set_particle(i2,sb2,0);
       //loop_1a_[sb1][i1][sb2][i2] = loop_1a_func(p1, p2);
       //------------------------------	       
       //p1.set_particle(i1,sb1,0);
       //p2.set_particle(i2,sb2,0);
       //loop_2b_[sb1][i1][sb2][i2] = loop_int(p1, p2, &loop::loop_2b_func, 0, 1);
       //------------------------------
       //p1.set_particle(i1,sb1,0);
       //p2.set_particle(i2,sb2,0);
       //loop_2e_[sb1][i1][sb2][i2] = loop_int(p1, p2, &loop::loop_2e_func, 1, 1);

       //LaOStO sum encloses all internal degrees of freedom (pseudospin here)
       p1.set_particle(i1,0,sb1);
       p2.set_particle(i2,0,sb2);
       loop_1a_[sb1][i1][sb2][i2] = loop_1a_func(p1, p2);
       loop_2b_[sb1][i1][sb2][i2] = loop_int(p1, p2, &loop::loop_2b_func);
       loop_2c_[sb1][i1][sb2][i2] = loop_int(p1, p2, &loop::loop_2c_func);
       loop_2d_[sb1][i1][sb2][i2] = loop_int(p1, p2, &loop::loop_2d_func);
       loop_2e_[sb1][i1][sb2][i2] = loop_int(p1, p2, &loop::loop_2e_func);
      }
     }
    }
   }
   // reecover remaining entries by symmetry
   for (int sb1 = 0; sb1 < sb_lim_; sb1++)
   {
    red_ind = fs_.sym_red_ind(sb1,0);
    for (ic = red_ind.begin(); ic!=red_ind.end(); ic++)
    {
     i1 = *ic;
     for (int sb2 = 0; sb2 < sb_lim_; sb2++)
     {
      for (int i2  = 0; i2 < npatch_; i2++)
      {
			 for (int op = 0; op < 4; op++)
			 {
				i1s = fs_.sym_patch(op,i1,0,sb1);
				i2s = fs_.sym_patch(op,i2,0,sb2);
				if (op > 1)
				{
					loop_1a_[sb1][i1s][sb2][i2s] = conj(loop_1a_[sb1][i1][sb2][i2]);
					loop_2b_[sb1][i1s][sb2][i2s] = conj(loop_2b_[sb1][i1][sb2][i2]);
					loop_2e_[sb1][i1s][sb2][i2s] = conj(loop_2e_[sb1][i1][sb2][i2]);
				}
				else
				{
					loop_1a_[sb1][i1s][sb2][i2s] = loop_1a_[sb1][i1][sb2][i2];
					loop_2b_[sb1][i1s][sb2][i2s] = loop_2b_[sb1][i1][sb2][i2];
					loop_2e_[sb1][i1s][sb2][i2s] = loop_2e_[sb1][i1][sb2][i2];
				}
			 }
			}
		 }
		}
	 }
}

void loop::initialize_loops()
{
   loop_1a_ = new cmplx***[sb_lim_];
   loop_2b_ = new cmplx***[sb_lim_];
   loop_2c_ = new cmplx***[sb_lim_];
   loop_2d_ = new cmplx***[sb_lim_];
   loop_2e_ = new cmplx***[sb_lim_];
   for (int i=0; i < sb_lim_; i++)
   {
	   loop_1a_[i] = new cmplx**[npatch_];
	   loop_2b_[i] = new cmplx**[npatch_];
	   loop_2c_[i] = new cmplx**[npatch_];
	   loop_2d_[i] = new cmplx**[npatch_];
	   loop_2e_[i] = new cmplx**[npatch_];
     for (int j=0; j < npatch_; j++)
     {
	   	loop_1a_[i][j] = new cmplx*[sb_lim_];
	   	loop_2b_[i][j] = new cmplx*[sb_lim_];
	   	loop_2c_[i][j] = new cmplx*[sb_lim_];
	   	loop_2d_[i][j] = new cmplx*[sb_lim_];
	   	loop_2e_[i][j] = new cmplx*[sb_lim_];
     	for (int k=0; k < sb_lim_; k++)
     	{
	   		loop_1a_[i][j][k] = new cmplx[npatch_];
	   		loop_2b_[i][j][k] = new cmplx[npatch_];
	   		loop_2c_[i][j][k] = new cmplx[npatch_];
	   		loop_2d_[i][j][k] = new cmplx[npatch_];
	   		loop_2e_[i][j][k] = new cmplx[npatch_];
			}
   	 }
    }
}

void loop::print_loop_integrals(const char* file_prefix)
{

   int width1 = 6;
   int width2 = 18;

   string pre_str(file_prefix);
   string loop_1a_str("loop_1a.dat");
   string loop_2b_str("loop_2b.dat");
   string loop_2b_str("loop_2c.dat");
   string loop_2b_str("loop_2d.dat");
   string loop_2e_str("loop_2e.dat");

   ofstream out_1a((pre_str + loop_1a_str).c_str());
   ofstream out_2b((pre_str + loop_2b_str).c_str());
   ofstream out_2c((pre_str + loop_2c_str).c_str());
   ofstream out_2d((pre_str + loop_2d_str).c_str());
   ofstream out_2e((pre_str + loop_2e_str).c_str());

   out_1a.precision(10); 
   out_2b.precision(10); 
   out_2c.precision(10); 
   out_2d.precision(10); 
   out_2e.precision(10); 
   out_1a << setw(width1) << right << "sb1" << setw(width1) << right << "p1" << setw(width1) << right << "sb2" << setw(width1) << right << "p2";
   out_1a << setw(width2+width2) << right << "loop_int" << "\n";
   out_2b << setw(width1) << right << "sb1" << setw(width1) << right << "p1" << setw(width1) << right << "sb2" << setw(width1) << right << "p2";
   out_2b << setw(width2+width2) << right << "loop_int" << "\n";
   out_2c << setw(width1) << right << "sb1" << setw(width1) << right << "p1" << setw(width1) << right << "sb2" << setw(width1) << right << "p2";
   out_2c << setw(width2+width2) << right << "loop_int" << "\n";
   out_2d << setw(width1) << right << "sb1" << setw(width1) << right << "p1" << setw(width1) << right << "sb2" << setw(width1) << right << "p2";
   out_2d << setw(width2+width2) << right << "loop_int" << "\n";
   out_2e << setw(width1) << right << "sb1" << setw(width1) << right << "p1" << setw(width1) << right << "sb2" << setw(width1) << right << "p2";
   out_2e << setw(width2+width2) << right << "loop_int" << "\n";
   for (int i=0; i < sb_lim_; i++)
   {
    for (int j=0; j < npatch_; j++)
    {
 	   for (int k=0; k < sb_lim_; k++)
 	   {
      for (int l=0; l < npatch_; l++)
      {
       out_1a << setw(width1) << right << i ;
       out_1a << setw(width1) << right << j ;
       out_1a << setw(width1) << right << k ;
       out_1a << setw(width1) << right << l ;
       out_1a << setw(width2) << right << loop_1a_[i][j][k][l].real();
       out_1a << setw(width2) << right << loop_1a_[i][j][k][l].imag();
       out_1a << "\n";
       out_2b << setw(width1) << right << i ;
       out_2b << setw(width1) << right << j ;
       out_2b << setw(width1) << right << k ;
       out_2b << setw(width1) << right << l ;
       out_2b << setw(width2) << right << loop_2b_[i][j][k][l].real();
       out_2b << setw(width2) << right << loop_2b_[i][j][k][l].imag();
       out_2b << "\n";
       out_2c << setw(width1) << right << i ;
       out_2c << setw(width1) << right << j ;
       out_2c << setw(width1) << right << k ;
       out_2c << setw(width1) << right << l ;
       out_2c << setw(width2) << right << loop_2c_[i][j][k][l].real();
       out_2c << setw(width2) << right << loop_2c_[i][j][k][l].imag();
       out_2c << "\n";
       out_2d << setw(width1) << right << i ;
       out_2d << setw(width1) << right << j ;
       out_2d << setw(width1) << right << k ;
       out_2d << setw(width1) << right << l ;
       out_2d << setw(width2) << right << loop_2d_[i][j][k][l].real();
       out_2d << setw(width2) << right << loop_2d_[i][j][k][l].imag();
       out_2d << "\n";
       out_2e << setw(width1) << right << i ;
       out_2e << setw(width1) << right << j ;
       out_2e << setw(width1) << right << k ;
       out_2e << setw(width1) << right << l ;
       out_2e << setw(width2) << right << loop_2e_[i][j][k][l].real();
       out_2e << setw(width2) << right << loop_2e_[i][j][k][l].imag();
       out_2e << "\n";
      }
     }
    }
   }
}


void loop::print_loops(const char* filename, int flag)
{
   ofstream file_output(filename);

   int band_lim = b_pars_.band_num;
   int spin_lim = b_pars_.spin_num;

   int num = 20;
   int width1 = 6;
   int width2 = 18;
   file_output.precision(10); 
   cout.precision(10);
   
   particle p1(fs_,H_), p2(fs_,H_), p3(fs_,H_);
   double q3x, q3y;
   Vector3d q3;
   cmplx ll;

   default_random_engine generator;
   uniform_real_distribution<double> distribution(-M_PI,M_PI);

   for (int bix1 = 0; bix1 < band_lim; bix1++)
   {
    for (int sp1 = 0; sp1 < spin_lim; sp1++)
    { 
     for (int m1 = 0; m1 < npatch_; m1++)
     { 
      cout << bix1 << " " << sp1 << " " << m1 << endl;
      p1.set_particle(m1,bix1,sp1);
      for (int bix2 = 0; bix2 < band_lim; bix2++)
      {
       for (int sp2 = 0; sp2 < spin_lim; sp2++)
       { 
        for (int m2 = 0; m2 < npatch_; m2++)
        {
         p2.set_particle(m2,bix2,sp2);
         for (int bix3 = 0; bix3 < band_lim; bix3++)
         {
          for (int sp3 = 0; sp3 < spin_lim; sp3++)
          { 
           for (int m3 = 0; m3 < num; m3++)
	   			 {
	          q3x = distribution(generator);  
	          q3y = distribution(generator);  
	          q3 << q3x, q3y, 0.0;
	          p3.set_particle(q3,bix3,sp3);
       	    for (int bix4 = 0; bix4 < band_lim; bix4++)
            {
             for (int sp4 = 0; sp4 < spin_lim; sp4++)
	     			 {
	       			switch (flag)
	       			{
	       				case 0:
             			ll = loop_2b_func(p1,p2,p3,bix4,sp4);
                  break;
	       				case 1:
             			ll = loop_2c_func(p1,p2,p3,bix4,sp4);
                  break;
	       				case 2:
             			ll = loop_2d_func(p1,p2,p3,bix4,sp4);
                  break;
	       				case 3:
             			ll = loop_2e_func(p1,p2,p3,bix4,sp4);
                  break;
	            }
	       			if (abs(ll) > 1.0e-7)
               { 
	       				file_output << setw(width1) << right << bix1 ;
   	       			file_output << setw(width1) << right << sp1 ;
   	       			file_output << setw(width1) << right << m1 ;
	       				file_output << setw(width1) << right << bix2 ;
   	       			file_output << setw(width1) << right << sp2 ;
   	       			file_output << setw(width1) << right << m2 ;
	       				file_output << setw(width1) << right << bix3 ;
   	       			file_output << setw(width1) << right << sp3 ;
   	       			file_output << setw(width2) << right << q3x ;
   	       			file_output << setw(width2) << right << q3y ;
	       				file_output << setw(width1) << right << bix4 ;
   	       			file_output << setw(width1) << right << sp4 ;
   	       			file_output << setw(width2) << right << ll.real() ;
   	       			file_output << setw(width2) << right << ll.imag() ;
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
   file_output.close();
}
