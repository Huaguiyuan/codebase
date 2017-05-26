#ifndef PARTICLE
#define PARTICLE

#include "fermi_surface.hpp"
#include "symmetry.hpp"

class particle
{

  private:
		int patch_; //-1 if particle does not correspond to patch
		Vector3d kk_;  // teilchenimpuls
		int bix_; // bandindex (starts at 0)
		int sp_ ; //spinindex (0 oder 1)
  
  private:
		const hamiltonian &H_;
    const fermi_surface &fs_;

  private:
    int ind();

	public:
		particle(const fermi_surface &fs, const hamiltonian &H);
		particle(int patch, int bix, int sp, const fermi_surface &fs, const hamiltonian &H);
		particle(Vector3d kk, int bix, int sp, const fermi_surface &fs, const hamiltonian &H);
    particle& operator= (const particle &other);
    particle(const particle &other);
    
	  void set_particle(int patch, int bix, int sp);
	  void set_particle(Vector3d kk, int bix, int sp);

    int get_patch() const;
    Vector3d get_kk() const;
    int get_bix() const;
    int get_sp() const;

  	VectorXcd get_state() const;
  	double get_en() const;

	  void Pinv(); 
	  void Tinv();
	  void Sp();
    void Uop(int op);
 
		cmplx get_char(int op);
};

particle Pinv(const particle &other);
particle Tinv(const particle &other);
particle Sp(const particle &other);
particle Uop(const particle &other, int op);
cmplx charc(const particle &other, int op);
ostream& operator<<(ostream& out, const particle& pp);

#endif
