#include "fermi_surface.hpp"

fermi_surface::fermi_surface(brillouin_zone BZ, hamiltonian H, int npatch) : BZ_(BZ), H_(H)
{
	npatch_         = npatch;
	bands_.lo_band  = 0; 
	bands_.band_num = 1;
	bands_.spin_num = 2;
  sym_ = c4v_dxz_dyz();
  calc_FS();
  calc_minus();
  calc_chars();
}

fermi_surface::fermi_surface(const fermi_surface &other) : BZ_(other.BZ_), H_(other.H_), sym_(other.sym_)
{
	npatch_          = other.npatch_;
	bands_.lo_band   = other.bands_.lo_band;
	bands_.band_num  = other.bands_.band_num;
	bands_.spin_num  = other.bands_.spin_num;

	// allocate memory for FS and FS_states
	int nn = npatch_;
	int band_lim = bands_.band_num;
	int spin_lim = bands_.spin_num;
 	FS = new Vector3d**[band_lim];
 	FS_states = new VectorXcd**[band_lim];
 	for ( int i = 0; i < band_lim; i++ )
 	{
 		FS[i] = new Vector3d*[spin_lim];
   	FS_states[i] = new VectorXcd*[band_lim];
    for ( int j = 0; j < spin_lim; j++ )
    	{
      	FS[i][j] = new Vector3d[nn];
        FS_states[i][j] = new VectorXcd[nn];
      }
  }
	// copy FS and FS_states
  for ( int i = 0; i < band_lim; i++ )
  {
  	for ( int j = 0; j < spin_lim; j++ )
    {
    	for ( int k = 0; k < nn; k++ )
      {
				FS[i][j][k]        = other.FS[i][j][k];
				FS_states[i][j][k] = other.FS_states[i][j][k];
      }
    }
  }

 minus = new int**[npatch_];
 for ( int i = 0; i < npatch_; i++ )
 {
 		minus[i] = new int*[band_lim];
 		for ( int j = 0; j < band_lim; j++ )
 		{
     	 minus[i][j] = new int[spin_lim];
			 for (int k = 0; k < spin_lim; k++)
			 {
					 minus[i][j][k] = other.minus[i][j][k];
			 }
	 	}
	}

 chars = new cmplx***[npatch_];
 for ( int i = 0; i < npatch_; i++ )
 {
 		chars[i] = new cmplx**[band_lim];
 		for ( int j = 0; j < band_lim; j++ )
 		{
     	 chars[i][j] = new cmplx*[spin_lim];
			 for (int k = 0; k < spin_lim; k++)
			 {
     	    chars[i][j][k] = new cmplx[sym_.get_sym_num()];
			 		for (int l = 0; l < sym_.get_sym_num(); l++)
			 		{
						chars[i][j][k][l] = other.chars[i][j][k][l];
			 		}
			 }
	 	}
	}
}

fermi_surface& fermi_surface :: operator= (const fermi_surface &other)
{
	// avoid self-assignment
	if ( this != &other )
  {       
    BZ_  = other.BZ_;
		H_   = other.H_;
    sym_ = other.sym_;
		npatch_          = other.npatch_;
		bands_.lo_band   = other.bands_.lo_band;
		bands_.band_num  = other.bands_.band_num;
		bands_.spin_num  = other.bands_.spin_num;

		// allocate memory for FS and FS_states
		int nn = npatch_;
		int band_lim = bands_.band_num;
		int spin_lim = bands_.spin_num;
   	FS = new Vector3d**[band_lim];
   	FS_states = new VectorXcd**[band_lim];
   	for ( int i = 0; i < band_lim; i++ )
   	{
      FS[i] = new Vector3d*[spin_lim];
   		FS_states[i] = new VectorXcd*[band_lim];
      for ( int j = 0; j < spin_lim; j++ )
      	{
        	FS[i][j] = new Vector3d[nn];
         	FS_states[i][j] = new VectorXcd[nn];
      	}
   	}
		// copy FS and FS_states
   	for ( int i = 0; i < band_lim; i++ )
   	{
    	for ( int j = 0; j < spin_lim; j++ )
      {
      	for ( int k = 0; k < nn; k++ )
      	{
					FS[i][j][k]        = other.FS[i][j][k];
					FS_states[i][j][k] = other.FS_states[i][j][k];
      	}
      }
   	}
    // copy minus array
 		minus = new int**[npatch_];
 		for ( int i = 0; i < npatch_; i++ )
 		{
 			minus[i] = new int*[band_lim];
 			for ( int j = 0; j < band_lim; j++ )
 			{
     	 	minus[i][j] = new int[spin_lim];
			 	for (int k = 0; k < spin_lim; k++)
			 	{
					 minus[i][j][k] = other.minus[i][j][k];
			 	}
	 		}
		}

	 // copy characters
 	chars = new cmplx***[npatch_];
 	for ( int i = 0; i < npatch_; i++ )
 	{
 		 chars[i] = new cmplx**[band_lim];
 		 for ( int j = 0; j < band_lim; j++ )
 		 {
     	 chars[i][j] = new cmplx*[spin_lim];
			 for (int k = 0; k < spin_lim; k++)
			 {
     	    chars[i][j][k] = new cmplx[sym_.get_sym_num()];
			 		for (int l = 0; l < sym_.get_sym_num(); l++)
			 		{
						chars[i][j][k][l] = other.chars[i][j][k][l];
			 		}
		   }
	 	 }
	 }
	}
	return *this;
}

Vector3d fermi_surface::get_kpt(int patch, int bix, int sp) const
{
	if ( 0 <= patch && patch < npatch_ && 0 <= bix && bix < bands_.band_num && 0 <= sp && sp < bands_.spin_num)
  {
		return FS[bix][sp][patch];
	}
	else
  {
		std::invalid_argument("invalid indices in get_kpt");  
		return FS[0][0][0];
	}
}

VectorXcd fermi_surface::get_state(int patch, int bix, int sp) const
{
	if ( 0 <= patch && patch < npatch_ && 0 <= bix && bix < bands_.band_num && 0 <= sp && sp < bands_.spin_num)
	{
		return FS_states[bix][sp][patch];
	}
	else
  {
		std::invalid_argument("invalid indices in get_state");  
		return FS_states[0][0][0];
	}
}

cmplx fermi_surface::get_char(int patch, int bix, int sp, int op) const
{
	if ( 0 <= op && op < sym_.get_sym_num() && 0 <= patch && patch < npatch_ && 0 <= bix && bix < bands_.band_num && 0 <= sp && sp < bands_.spin_num)
	{
		return chars[patch][bix][sp][op];
	}
	else
  {
		std::invalid_argument("invalid indices in get_chars");  
    cout << " invalid indices in get_chars" << endl;
		return 1.0;
	}
}

void fermi_surface::set_band_params(int lo_band, int band_num, int spin_num)
{
	bands_.lo_band = lo_band;
	bands_.band_num = band_num;
	bands_.spin_num = spin_num;
}

int fermi_surface::get_npatch() const
{
	return npatch_;
}

band_params fermi_surface::get_band_params() const
{
	struct band_params b_pars;
	b_pars.lo_band  = bands_.lo_band;
	b_pars.band_num = bands_.band_num;
	b_pars.spin_num = bands_.spin_num;
	return b_pars;
}

int_vec fermi_surface::sym_red_ind(int bix, int sp) const
{
	int_vec red_ind;
 	Vector3d fs_vec;
 	double ang;
 	Vector2d ang_lim(sym_.get_sym_range());

	if (0 <= bix && bix < bands_.band_num && 0 <= sp && sp < bands_.spin_num)
	{
  	for (int p = 0; p < npatch_; p++)
  	{ 
	   	fs_vec = FS[bix][sp][p];
     	ang    = atan2(fs_vec(1), fs_vec(0));
     	if (ang_lim(0) < ang && ang < ang_lim(1) > 0.0)
     	{
       	red_ind.push_back(p);
     	}
  	}
  	return red_ind;
	}
	else
	{
		std::invalid_argument("invalid indices in sym_red_ind");  
		return red_ind;
	}
}

int fermi_surface::ind(int bix,int sp) const
{
	if (0 <= bix && bix < bands_.band_num && 0 <= sp && sp < bands_.spin_num)
	{
		return (bands_.lo_band + bix)*bands_.spin_num + sp;
	}
	else
	{
		std::invalid_argument("invalid indices in ind");  
		return 0;
	}
}

int fermi_surface::find_patch(Vector3d kk, int bix, int sp) const
{
	Vector3d ki, ktemp; 
	int mpatch;
	double knorm;
	double dist = 10000;

	if (0 <= bix && bix < bands_.band_num && 0 <= sp && sp < bands_.spin_num)
  {
		for (int i = 0; i < npatch_; i++)
		{
			ki    = FS[bix][sp][i];
			ktemp = ki - kk;
			knorm = ktemp.norm();
			if (knorm < dist)
			{
				mpatch = i;
				dist = knorm; 
			}		
		}
		return mpatch;
	}
	else
	{
		std::invalid_argument("invalid radius in find_patch");  
		return 0;
	}
}

int fermi_surface::minus_patch(int patch, int bix, int sp) const
{
	if ( 0 <= patch && patch < npatch_ && 0 <= bix && bix < bands_.band_num && 0 <= sp && sp < bands_.spin_num)
	{
		Vector3d kf = FS[bix][sp][patch];
		return find_patch(-kf,bix,sp);
	}
	else
	{
		std::invalid_argument("invalid indices in minus_patch");  
    return 0;
	}
}

int fermi_surface::sym_patch(int op, int patch, int bix, int sp) const
{
	if (0 <= op && op < sym_.get_sym_num() && 0 <= patch && patch < npatch_ && 0 <= bix && bix < bands_.band_num && 0 <= sp && sp < bands_.spin_num)
	{
  	const Matrix3d_vec &operation = sym_.get_k_ops();
		Vector3d kf = FS[bix][sp][patch];
  	return find_patch(operation[op]*kf, bix, sp);
	}
	else
	{
		std::invalid_argument("invalid indices or operation in sym_patch");  
		return 0;
	}
}

double get_enrad_wrapper(double rr, void *params)
{
	struct add_params * p = (struct add_params *) params;
	Vector3d pt = p->pt;
	double ang = p->ang;
	int bdx = p->bdx;
	int sp = p->sp;
	int m = p->fs->ind(bdx,sp);
	Vector3d kk = p->bz->r2xy(rr,ang,pt);
	return p->clptr->get_ham_eigval(kk,m);
}

cmplx fermi_surface::get_char(Vector3d kk, int bix, int sp, int op) const
{
 int dim = H_.get_dim();

 const Matrix4cd_vec &state_ops = sym_.get_state_ops();
 const Matrix3d_vec  &k_ops     = sym_.get_k_ops();

 VectorXcd vec0, vec1;
 VectorXcd diff; 
 ArrayXd vec_abs(dim);

 int index;
 cmplx res;
 cmplx max_val;

 vec0 = state_ops[op]*H_.get_ham_eigvec(kk, ind(bix,sp));
 vec1 = H_.get_ham_eigvec(k_ops[op]*kk, ind(bix,sp));

 vec_abs = vec0.array().abs();

 max_val = vec_abs.maxCoeff(&index);
 res  = vec1(index)/vec0(index);
 diff = vec1 - res*vec0;

 if (diff.norm() > 1e-8)
 {
  cout << "error in character calculation:" << endl;
  cout << diff.norm() << endl;
  std::invalid_argument("error in character calculation");  
  return res;
 }
 else
 {
  return res;
 }
}

void fermi_surface::calc_chars()
{
 int sym_num = sym_.get_sym_num();
 int band_lim = bands_.band_num;
 int spin_lim = bands_.spin_num;

 chars = new cmplx***[npatch_];
 for ( int i = 0; i < npatch_; i++ )
 {
 		chars[i] = new cmplx**[band_lim];
 		for ( int j = 0; j < band_lim; j++ )
 		{
     	 chars[i][j] = new cmplx*[spin_lim];
			 for (int k = 0; k < spin_lim; k++)
			 {
     	 		chars[i][j][k] = new cmplx[sym_num];
			 		for (int l = 0; l < sym_num; l++)
			 		{
							chars[i][j][k][l] = get_char(FS[j][k][i],j,k,l);
							//cout << i << " " << j << " " << k << " " << l << " " << abs(res) << " " << arg(res) << endl;
			 		}
			 }
	 	}
	}
}

void fermi_surface::calc_minus()
{
 int band_lim = bands_.band_num;
 int spin_lim = bands_.spin_num;

 minus = new int**[npatch_];
 for ( int i = 0; i < npatch_; i++ )
 {
 		minus[i] = new int*[band_lim];
 		for ( int j = 0; j < band_lim; j++ )
 		{
     	 minus[i][j] = new int[spin_lim];
			 for (int k = 0; k < spin_lim; k++)
			 {
					 minus[i][j][k] = minus_patch(i,j,k);
			 }
	 	}
	}
}

void fermi_surface::calc_FS()
{
	int band_lim = bands_.band_num;
	int spin_lim = bands_.spin_num;
  Vector3d kk;

 	FS = new Vector3d**[band_lim];
 	FS_states = new VectorXcd**[band_lim];
 	for ( int i = 0; i < band_lim; i++ )
 	{
   		FS[i] = new Vector3d*[spin_lim];
   		FS_states[i] = new VectorXcd*[band_lim];
   		for ( int j = 0; j < spin_lim; j++ )
   		{
       		FS[i][j] = new Vector3d[npatch_];
       		FS_states[i][j] = new VectorXcd[npatch_];
			 		calc_single_FS(i,j,FS[i][j]);
					for (int k = 0; k < npatch_; k++)
					{
					  	kk = FS[i][j][k];
							FS_states[i][j][k] = H_.get_ham_eigvec(kk, ind(i,j));
					}
      }
   }
}

void fermi_surface::calc_single_FS(int bdx, int sp, Vector3d *fs_pts)
{
 	int max_iter = 100;
	double abs_err = 1e-10;
	double rel_err = 1e-10;
  Vector3d kpt(0.0,0.0,0.0);

	struct add_params pars = {kpt, 0.0, bdx, sp, &BZ_, this, &H_};

	int count = 0;
	int status, iter;
	double r0, r1, rsol;
	double sign;
	const gsl_root_fsolver_type *T;
	gsl_root_fsolver *s;
	gsl_function F;
  F.function = &get_enrad_wrapper;
  T = gsl_root_fsolver_brent;
  s = gsl_root_fsolver_alloc (T);

	for ( int i = 0; i < npatch_; i++)
	{
			pars.ang = (2*i + 1)*M_PI/npatch_; 
     	F.params = &pars;
			iter = 0;
			r0 = 0;
			r1 = BZ_.get_patch_len(pars.ang);
			sign = GSL_FN_EVAL(&F,r0)*GSL_FN_EVAL(&F,r1);
			assert(sign < 0.0 && "no root bracketing");
     	gsl_root_fsolver_set (s, &F, r0, r1);
      do
      {
				iter++;
        status = gsl_root_fsolver_iterate (s);
     		rsol   = gsl_root_fsolver_root (s);
				r0     = gsl_root_fsolver_x_lower (s);
				r1     = gsl_root_fsolver_x_upper (s);
				status = gsl_root_test_interval (r0, r1, abs_err, rel_err);
			}
			while (status == GSL_CONTINUE && iter < max_iter);
			fs_pts[count] << BZ_.r2xy(rsol,pars.ang,pars.pt);
			count = count + 1;
	}
  gsl_root_fsolver_free (s);
}

void fermi_surface::print_fs(const char* filename) const
{
   ofstream file_output(filename);

   int band_lim = bands_.band_num;
   int spin_lim = bands_.spin_num;
   int width1 = 6;
   int width2 = 15;
   file_output.precision(10); 
   file_output << setw(width1) << right << "band" << setw(width1) << right << "spin" << setw(width1) << right << "patch";
   file_output << setw(width2) << right << "px" << setw(width2) << right << "py" << setw(width2) << right << "pz";  
   file_output << "\n";

   for (int bix = 0; bix < band_lim; bix++)
   {
	  for (int sp = 0; sp < spin_lim; sp++)
	  { 
     for ( int i = 0; i < npatch_; i++)
   	 {
			file_output << setw(width1) << right << bix ;
   		file_output << setw(width1) << right << sp ;
   		file_output << setw(width1) << right << i ;
			file_output << setw(width2) << right << FS[bix][sp][i](0);
			file_output << setw(width2) << right << FS[bix][sp][i](1);
			file_output << setw(width2) << right << FS[bix][sp][i](2);
			file_output << "\n";
		 }
	  }
   }
   file_output.close();
}	

void fermi_surface::print_states(const char* filename) const
{
   ofstream file_output(filename);
   int band_lim = bands_.band_num;
   int spin_lim = bands_.spin_num;
   int dim = H_.get_dim();
   int width1 = 6;
   int width2 = 18;
   file_output.precision(10); 
   file_output << setw(width1) << right << "band" << setw(width1) << right << "spin" << setw(width1) << right << "patch";
   file_output << setw(width2) << left << "FS_pt";  
   file_output << setw(3*width2) << left << "state";  
   file_output << "\n";
   cout.precision(10);
   
   for (int bix = 0; bix < band_lim; bix++)
   {
	  for (int sp = 0; sp < spin_lim; sp++)
	  { 
   	 	for ( int i = 0; i < npatch_; i++)
   		{
			file_output << setw(width1) << right << bix ;
   		file_output << setw(width1) << right << sp ;
   		file_output << setw(width1) << right << i ;
			file_output << setw(width2) << right << FS[bix][sp][i](0);
			file_output << setw(width2) << right << FS[bix][sp][i](1);
			file_output << setw(width2) << right << FS[bix][sp][i](2);
			for (int kk = 0; kk < dim; kk++)
			{
				file_output << setw(width2) << right << (FS_states[bix][sp][i](kk)).real();
				file_output << setw(width2) << right << (FS_states[bix][sp][i](kk)).imag();
			}	
			file_output << "\n";
		}
	  }
   }
   file_output.close();
}	

void fermi_surface::print_sym_ops(const char* filename) const
{
   ofstream file_output(filename);
   int band_lim = bands_.band_num;
   int spin_lim = bands_.spin_num;
   int width = 6;
   file_output.precision(10); 
   file_output << setw(width) << right << "band" << setw(width) << right << "spin" << setw(width) << right << "patch";
   file_output << setw(width) << right << "op" << setw(2*width) << right << "op*patch";  
   file_output << "\n";
   cout.precision(10);

   for (int bix = 0; bix < band_lim; bix++)
   {
	  for (int sp = 0; sp < spin_lim; sp++)
	  { 
   	 for ( int i = 0; i < npatch_; i++)
   	 {
   	 	for ( int op = 0; op < sym_.get_sym_num(); op++)
   		{
			file_output << setw(width) << right << bix ;
   		file_output << setw(width) << right << sp ;
   		file_output << setw(width) << right << i ;
			file_output << setw(width) << right << op;
			file_output << setw(2*width) << right << sym_patch(op, i, bix, sp);
      file_output << "\n";
      }
     }
    }
   }
}
