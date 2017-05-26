#ifndef HAMILTONIAN 
#define HAMILTONIAN

#include "./general_header.hpp"

struct model_pars
{
	double mu;
	double soc;
	double t;
	double eta;
	double ts;
	double tp;
	int gauge;
};

class hamiltonian 
{
  private:
	int dim_ham_;
  MatrixXcd Hmat_help(Vector3d kk, int flag) const; // flag = 0/1 for zero/nonzero spin orbit 

  protected:
	struct model_pars mp_;

  public:
	hamiltonian(double mu, double soc);
  hamiltonian(const hamiltonian &other);
	hamiltonian& operator=(const hamiltonian& other);

  int get_dim() const;
  void set_model_pars(double mu, double soc, double t, double eta, double ts, double tp, int gauge);

  MatrixXcd get_ham_mat(Vector3d kk) const;
	MatrixXcd get_ham_eigmat(Vector3d kk) const;
	VectorXcd get_ham_eigvec(Vector3d kk, int m) const;

  double    get_ham_eigval(Vector3d kk, int m) const;
  VectorXd  get_ham_eigval(Vector3d kk) const;

};

#endif
