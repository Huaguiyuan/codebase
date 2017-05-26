#ifndef HAMILTONIAN 
#define HAMILTONIAN

#include "./general_header.hpp"
#include "./symmetry.hpp"

struct model_pars
{
	double mu;
	double tp[7] ;
	int gauge;
};

class hamiltonian 
{
  private:
	int dim_ham_;
  MatrixXcd Hmat_help(Vector3d kk, int flag) const; // flag = 0/1 for zero/nonzero spin orbit 
	struct model_pars mp_;
  sym_rep sym_;

  public:
	hamiltonian(double mu);
  hamiltonian(const hamiltonian &other);
	hamiltonian& operator=(const hamiltonian& other);

  int get_dim() const;
  struct model_pars get_model_pars() const;
  void set_model_pars(double mu, double tp[7], int gauge);

  MatrixXcd get_ham_mat(Vector3d kk) const;
  MatrixXcd get_ham_mat(double kx, double ky, double kz) const;
	MatrixXcd get_ham_eigmat(Vector3d kk) const;
	VectorXcd get_ham_eigvec(Vector3d kk, int m) const;

  double    get_ham_eigval(Vector3d kk, int m) const;
  VectorXd  get_ham_eigval(Vector3d kk) const;

  void ham_test(const char* filename, int ngrid, int bix) const;
  void sym_test(Vector3d kk) const; 
};

#endif
