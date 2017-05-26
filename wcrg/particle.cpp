#include "particle.hpp"

particle::particle(const fermi_surface &fs, const hamiltonian &H): H_(H), fs_(fs)
{
	patch_ = -1;
	kk_    << 0.0, 0.0, 0.0;
	bix_   = 0;
	sp_    = 0;
}

particle::particle(int patch, int bix, int sp, const fermi_surface &fs, const hamiltonian &H): H_(H), fs_(fs)
{
	patch_ = patch;
	kk_    = fs_.get_kpt(patch,bix,sp);
	bix_   = bix;
	sp_    = sp;
}

particle::particle(Vector3d kk, int bix, int sp, const fermi_surface &fs, const hamiltonian &H): H_(H), fs_(fs)
{
	patch_ = -1;
	kk_    = kk;
	bix_   = bix;
	sp_    = sp;
}

particle& particle::operator= (const particle &other)
{
	// avoid self-assignment
	if ( this != &other )
  {       
		patch_ = other.patch_;
		kk_    = other.kk_;
		bix_   = other.bix_;
		sp_    = other.sp_;
  }
  return *this;
}

particle Pinv(const particle &other)
{
	particle pp(other);
  pp.Pinv();
  return pp;
}

particle Tinv(const particle &other)
{
	particle pp(other);
  pp.Tinv();
  return pp;
}

particle Sp(const particle &other)
{
	particle pp(other);
  pp.Sp();
  return pp;
}

particle Uop(const particle &other, int op)
{
	particle pp(other);
  pp.Uop(op);
  return pp;
}

cmplx charc(const particle &other, int op)
{
	particle pp(other);
  return pp.get_char(op);
}

particle::particle(const particle &other) : H_(other.H_), fs_(other.fs_)
{
	patch_  = other.patch_;
	kk_     = other.kk_;
	bix_    = other.bix_;
	sp_     = other.sp_;
}

ostream& operator<<(ostream& out, const particle& pp)
{
  out.precision(8);
	out << "patch = " << pp.get_patch() << " kk = " << setw(10) << right << (pp.get_kk()).transpose(); 
  out << " bix = " << pp.get_bix() << " sp = " << pp.get_sp();
  return out; 
}

void particle::set_particle(int patch, int bix, int sp)
{
	patch_ = patch;
	kk_    = fs_.get_kpt(patch, bix, sp);
	bix_   = bix;
	sp_    = sp;
}

void particle::set_particle(Vector3d kk, int bix, int sp)
{
	patch_ = -1;
	kk_    = kk;
	bix_   = bix;
	sp_    = sp;
}

int particle::get_patch() const
{
	return patch_;
}

Vector3d particle::get_kk() const
{
	return kk_;
}

int particle::get_bix() const
{
	return bix_;
}

int particle::get_sp() const
{
	return sp_;
}

VectorXcd particle::get_state() const
{	
	MatrixXcd mat  = H_.get_ham_eigmat(kk_);
	return mat.col(fs_.ind(bix_,sp_));
}

double particle::get_en() const
{	
	return H_.get_ham_eigval(kk_, fs_.ind(bix_,sp_));
}

void particle::Tinv()
{
	int k1m;
	sp_ = (sp_ + 1) % 2;
	if (patch_ == -1)
	{
		set_particle(-kk_, bix_, sp_);
	}
	else
	{
		k1m = fs_.get_minus(patch_, bix_, sp_);
		set_particle(k1m, bix_, sp_);
	}	
}

void particle::Pinv()
{
	int k1m;
	if (patch_ == -1)
	{
		set_particle(-kk_, bix_, sp_);
	}
	else
	{
		k1m = fs_.get_minus(patch_, bix_, sp_);
		set_particle(k1m, bix_, sp_);
	}	
}

void particle::Sp()
{
	sp_ = (sp_ + 1) % 2;
	set_particle(kk_, bix_, sp_);
}

void particle::Uop(int op)
{
	int k1s;
	if (patch_ == -1)
	{
  	sym_rep sym;
		set_particle(sym.get_k_op(op)*kk_, bix_, sp_);
	}
	else
	{
		k1s = fs_.sym_patch(patch_, bix_, sp_, op);
		set_particle(k1s, bix_, sp_);
	}	
}
 
cmplx particle::get_char(int op)
{
	if (patch_ == -1)
	{
	  return fs_.get_char(kk_, bix_, sp_, op);
	}
	else
	{
	  return fs_.get_char(patch_, bix_, sp_, op);
	}	
}
