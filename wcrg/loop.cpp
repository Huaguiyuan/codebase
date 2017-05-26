#include "loop.hpp"
#include "brillouin_zone.hpp"

loop::loop(double T, int nx, int ny, int int_flag, const brillouin_zone &bz, const hamiltonian &H, const fermi_surface &fs, const coupling_grid &cp): bz_(bz), fs_(fs), H_(H), cp_(cp)
{
	int_grid_n1_ = nx;
	int_grid_n2_ = ny;
  integration_flag_ = int_flag; // 0 for kartesian, 1 for radial along r2xy paths
	temperature_ = T;
	tflag_       = false; //true for zero temperature
	grid_flag_   = true; // uses states on precalculated grid
  b_pars_      = fs_.get_band_params();
  npatch_      = fs_.get_npatch();
  spin_lim_    = b_pars_.spin_num;
  band_lim_    = b_pars_.band_num;
  sym_         = sym_rep();
  set_in_particle_list();
  set_out_particle_list();
}

void loop::set_in_particle_list()
{
  particle pp(fs_,H_);
  int_vec red_ind;
  int_vec::const_iterator ic;
  for (int bix = 0; bix < band_lim_; bix++)
  {
 		red_ind = fs_.sym_red_ind(bix,0);
 		for (ic = red_ind.begin(); ic!=red_ind.end(); ic++)
 		{
  		pp.set_particle(*ic,bix,0);
			in_particles_.push_back(pp);
		}
	}
}

void loop::set_out_particle_list()
{
  particle pp(fs_,H_);
  for (int bix = 0; bix < band_lim_; bix++)
  {
		for (int p = 0; p < npatch_; p++)
 		{
   		pp.set_particle(p,bix,0);
			out_particles_.push_back(pp);
		}
	}
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

  Vector3d p4_fold = bz_.mapin(p3.get_kk() - p1.get_kk() - p2.get_kk());
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

  Vector3d p4_fold = bz_.mapin(p3.get_kk() - p1.get_kk() + p2.get_kk());
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

  Vector3d p4_fold = bz_.mapin(p3.get_kk() + p1.get_kk() - p2.get_kk());
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

  Vector3d p4_fold = bz_.mapin(p3.get_kk() - p1.get_kk() + p2.get_kk());
  particle p4(p4_fold,bix4,sp4,fs_,H_);

	double loop = frequency_integral(p3, p4);

	if (abs(loop) > 1.0e-10)
	{
		if (grid_flag_) 
		{
			V1 = cp_.VV_grid(p1,p4,p2,p3);
			V2 = cp_.VV_grid(Pinv(p1),p3,Pinv(p2),p4);
		}
		else
		{
			V1 = cp_.VV_exct(p1,p4,p2,p3);
			V2 = cp_.VV_exct(Pinv(p1),p3,Pinv(p2),p4);
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

   double dinf ; 
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
					case 1:
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
              dinf = (2*M_PI/int_grid_n1_)*(2*M_PI/int_grid_n2_); //compared to weejee there is another 1/(2pi)^2 factor
							break;
						case 1:
			   			r    = jj*rlim/int_grid_n2_;
			   			q3   = bz_.r2xy(r,ang,pt);
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

   double dinf;
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
				case 1:
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
            dinf = (2*M_PI/int_grid_n1_)*(2*M_PI/int_grid_n2_)/(4*M_PI*M_PI); 
						break;
					case 1:
			   		r    = jj*rlim/int_grid_n2_;
			 			q3   = bz_.r2xy(r,ang,pt);
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

void loop::calc_loops()
{
   particle p1(fs_,H_), p2(fs_,H_);


   result rr;
   std::list<result> loop_list_;

   // calculate symmetry reduced diagrams
   for (int i = 0; i < in_particles_.size(); i++)
   {
    //cout << "p1: " << std::distance(in_particles_.begin(), it1) << endl;
# pragma omp parallel for firstprivate(p1,p2,rr)
   	for (int j = 0; j < out_particles_.size(); j++)
   	{
     //cout << "p2: " << std::distance(out_particles_.begin(), it2) << endl;
		 p1           = in_particles_[i];
	   p2           = out_particles_[j];
     rr.in_patch  = p1.get_patch();
     rr.in_bix    = p1.get_bix();
     rr.in_sp     = p1.get_sp();
     rr.out_patch = p2.get_patch();
     rr.out_bix   = p2.get_bix();
     rr.out_sp    = p2.get_sp();
     rr.loop_1a   = loop_1a_func(p1, p2);
     rr.loop_2b   = loop_int(p1, p2, &loop::loop_2b_func, 0, 1);
     rr.loop_2c   = 0.0;
     rr.loop_2d   = 0.0;
     rr.loop_2e   = loop_int(p1, p2, &loop::loop_2e_func, 1, 1);
# pragma omp critical
     loop_list_.push_back(rr);
    }
   }

   // symmetry application
   for (int i = 0; i < in_particles_.size(); i++)
   {
   	for (int j = 0; j < out_particles_.size(); j++)
   	{
		 p1           = in_particles_[i];
	   p2           = out_particles_[j];
		 for (int op = 0; op < sym_.get_sym_num(); op++)
		 {
			p1 = Uop(p1, op);
			p2 = Uop(p2, op);
     	rr.in_patch  = p1.get_patch();
     	rr.in_bix    = p1.get_bix();
     	rr.in_sp     = p1.get_sp();
     	rr.out_patch = p2.get_patch();
     	rr.out_bix   = p2.get_bix();
     	rr.out_sp    = p2.get_sp();
      //characters not needed here, symmetry properties are applied by hand
     	//rr.loop_1a   = rr.loop_1a*charc(p1,op)*charc(Tinv(p1),op)*conj(charc(p2,op))*conj(charc(Tinv(p2),op));
     	//rr.loop_2b   = rr.loop_2b*charc(p1,op)*charc(Tinv(p1),op)*conj(charc(p2,op))*conj(charc(Tinv(p2),op));
     	//rr.loop_2c   = rr.loop_2c*charc(p1,op)*charc(Tinv(p1),op)*conj(charc(p2,op))*conj(charc(Tinv(p2),op));
     	//rr.loop_2d   = rr.loop_2d*charc(p1,op)*charc(Tinv(p1),op)*conj(charc(p2,op))*conj(charc(Tinv(p2),op));
     	//rr.loop_2e   = rr.loop_2e*charc(p1,op)*charc(Tinv(p1),op)*conj(charc(p2,op))*conj(charc(Tinv(p2),op));
      switch (op)
			{
			case 0:
     		rr.loop_1a   = rr.loop_1a;
     		rr.loop_2b   = rr.loop_2b;
     		rr.loop_2c   = rr.loop_2c;
     		rr.loop_2d   = rr.loop_2d;
     		rr.loop_2e   = rr.loop_2e;
			case 1:
     		rr.loop_1a   = rr.loop_1a;
     		rr.loop_2b   = rr.loop_2b;
     		rr.loop_2c   = rr.loop_2c;
     		rr.loop_2d   = rr.loop_2d;
     		rr.loop_2e   = rr.loop_2e;
			case 2:
     		rr.loop_1a   = conj(rr.loop_1a);
     		rr.loop_2b   = conj(rr.loop_2b);
     		rr.loop_2c   = conj(rr.loop_2c);
     		rr.loop_2d   = conj(rr.loop_2d);
     		rr.loop_2e   = conj(rr.loop_2e);
			case 3:
     		rr.loop_1a   = conj(rr.loop_1a);
     		rr.loop_2b   = conj(rr.loop_2b);
     		rr.loop_2c   = conj(rr.loop_2c);
     		rr.loop_2d   = conj(rr.loop_2d);
     		rr.loop_2e   = conj(rr.loop_2e);
			}
     	loop_list_.push_back(rr);
		 }
    }
   }
}

void loop::print_loop_integrals(const char* filename)
{
   ofstream file_output(filename);
   int width = 10;
   file_output.precision(8); 
   file_output << setw(width) << right << "p1_band" << setw(width) << right << "p1_spin" << setw(width) << right << "p1_patch";
   file_output << setw(width) << right << "p2_band" << setw(width) << right << "p2_spin" << setw(width) << right << "p2_patch";
   file_output << setw(width) << right << "Re_1a"   << setw(width) << right << "Im_1a";  
   file_output << setw(width) << right << "Re_2b"   << setw(width) << right << "Im_2b";  
   file_output << setw(width) << right << "Re_2c"   << setw(width) << right << "Im_2c";  
   file_output << setw(width) << right << "Re_2d"   << setw(width) << right << "Im_2d";  
   file_output << setw(width) << right << "Re_2e"   << setw(width) << right << "Im_2e";  
   file_output << "\n";

   std::list<result>::iterator it;
   for (it = loop_list_.begin(); it != loop_list_.end(); it++)
   {
       file_output << setw(width) << right << it->in_bix;
       file_output << setw(width) << right << it->in_sp;
       file_output << setw(width) << right << it->in_patch;
       file_output << setw(width) << right << it->out_bix;
       file_output << setw(width) << right << it->out_sp;
       file_output << setw(width) << right << it->out_patch;

       file_output << setw(width) << right << it->loop_1a.real();
       file_output << setw(width) << right << it->loop_1a.imag();
       file_output << setw(width) << right << it->loop_2b.real();
       file_output << setw(width) << right << it->loop_2b.imag();
       file_output << setw(width) << right << it->loop_2c.real();
       file_output << setw(width) << right << it->loop_2c.imag();
       file_output << setw(width) << right << it->loop_2d.real();
       file_output << setw(width) << right << it->loop_2d.imag();
       file_output << setw(width) << right << it->loop_2e.real();
       file_output << setw(width) << right << it->loop_2e.imag();
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
