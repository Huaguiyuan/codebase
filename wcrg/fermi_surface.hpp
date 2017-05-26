#ifndef FERMI_SURFACE
#define FERMI_SURFACE

#include "hamiltonian.hpp"
#include "brillouin_zone.hpp"
#include "symmetry.hpp"

struct band_params
{
	int lo_band; 
	int band_num;
	int spin_num;
};

class fermi_surface 
{
  private:
  int npatch_;
  brillouin_zone BZ_;
	hamiltonian H_;
  sym_rep sym_;
	struct band_params bands_;
	Vector3d ***FS_; 
	VectorXcd ***FS_states_; 
  int ***minus_;
  cmplx ****chars_;

  public:
	fermi_surface(brillouin_zone BZ, hamiltonian H, int npatch);
  fermi_surface(const fermi_surface &other);
	fermi_surface& operator=(const fermi_surface& other);

  void calc_single_FS(int bdx, int p, Vector3d *fs_pts);
  void calc_FS();
  void calc_minus();
  void calc_chars();

  Vector3d  get_kpt(int patch, int bix, int sp) const;
  VectorXcd get_state(int patch, int bix, int sp) const;
	cmplx get_char(int patch, int bix, int sp, int op) const;
	cmplx get_char(Vector3d kk, int bix, int sp, int op) const;
  int get_minus(int patch, int bix, int sp) const;

	void print_fs(const char* filename) const;
	void print_states(const char* filename) const;
  void print_sym_ops(const char* filename) const;

  public:
	void set_band_params(int lo_band, int band_num, int spin_num);
	band_params get_band_params() const;
  int get_npatch() const;

  int_vec sym_red_ind(int bix, int sp) const;
  int ind(int bix,int sp) const; //maps band and spin index to eigenvalue index

  int minus_patch(int patch, int bix, int sp) const; //calculates minus patch
  int find_patch(Vector3d kk, int bix, int sp) const;
  int sym_patch(int patch, int bix, int sp, int op) const;
};

struct add_params
{
	Vector3d pt;
	double ang;
	int bdx;
	int sp;
  brillouin_zone* bz;
	fermi_surface* fs;
	hamiltonian* clptr;
};
#endif
