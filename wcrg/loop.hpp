#include "coupling_grid.hpp"

struct result
{
	int in_patch;
  int in_bix;
  int in_sp;
	int out_patch;
  int out_bix;
  int out_sp;
  cmplx  loop_1a;
  cmplx  loop_2b;
  cmplx  loop_2c;
  cmplx  loop_2d;
  cmplx  loop_2e;
};

class loop 
{
  private:
    const brillouin_zone &bz_;
		const fermi_surface &fs_;
		const hamiltonian &H_;
		const coupling_grid &cp_;
		bool tflag_;
		bool grid_flag_;
  	struct band_params b_pars_;
    sym_rep sym_;

  private:
		double temperature_; // 0 for kartesian, 1 for radial along r2xy paths
    int integration_flag_;
		int int_grid_n1_;
		int int_grid_n2_;
    int npatch_;
    int spin_lim_;   
    int band_lim_;   

  private:
		std::vector<particle> out_particles_;
		std::vector<particle> in_particles_;
    std::list<result> loop_list_;

  public:
		loop(double T, int nx, int ny, int int_flag, const brillouin_zone &bz, const hamiltonian &H, const fermi_surface &fs, const coupling_grid &cp);

    void set_in_particle_list();
    void set_out_particle_list();
		void calc_loops();

  	//sp3, sp4 spin of inner propagators
  	cmplx loop_int(const particle &p1, const particle &p2, cmplx (loop::*func) (const particle &, const particle &, const particle &, int, int), int sp3, int sp4); 
  	//integration over all spin indices sp3 and sp4 in inner propagators
  	cmplx loop_int(const particle &p1, const particle &p2, cmplx (loop::*func) (const particle &, const particle &, const particle &, int, int)); 
  	cmplx loop_1a_func(const particle &p1, const particle &p2);
  	cmplx loop_2b_func(const particle &p1, const particle &p2, const particle &p3, int bix4, int sp4);
  	cmplx loop_2c_func(const particle &p1, const particle &p2, const particle &p3, int bix4, int sp4);
  	cmplx loop_2d_func(const particle &p1, const particle &p2, const particle &p3, int bix4, int sp4);
  	cmplx loop_2e_func(const particle &p1, const particle &p2, const particle &p3, int bix4, int sp4);

  	double frequency_integral(const particle &p1, const particle &p2);
  	double n_F(const particle &p1);
		void set_tflag(bool tflag);
  	void set_grid_flag(bool grid_flag);
		void print_loop_integrals(const char* file_prefix);
  	void print_loops(const char* filename, int flag);
};

