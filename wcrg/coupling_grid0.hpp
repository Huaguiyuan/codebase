#ifndef COUPLING_GRID
#define COUPLING_GRID

#include "particle.hpp"
#include "symmetry.hpp"

class coupling_grid 
{
  private:
  brillouin_zone bz_;
	fermi_surface fs_;
	hamiltonian H_;
	struct band_params b_pars_;
  sym_rep sym_;

  public:
	VectorXd **en_grid;
	VectorXcd ****states_grid;
	int grid_nx;
	int grid_ny;

  public:
	coupling_grid(int nx, int ny, brillouin_zone bz, hamiltonian H, fermi_surface fs);

  void calc_grid(int n1, int n2); //calc rectangular grid around BZ
  cmplx onsite_interaction(VectorXcd u1, VectorXcd u2, VectorXcd u3, VectorXcd u4) const;

	cmplx VV_exct(const particle &p1, const particle &p2, const particle &p3, const particle &p4) const;
	cmplx VV_grid(const particle &p1, const particle &p2, const particle &p3, const particle &p4) const;
	void print_VV(const char* filename, int flag) const;
	void VV_check_sym(const char* filename) const;

	VectorXcd get_state_grid(const particle &pp) const;
	VectorXcd get_state_exct(const particle &pp) const;
  Vector2i get_grid_point(const Vector3d kk) const;
};
#endif
