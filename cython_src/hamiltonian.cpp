#include "hamiltonian.hpp"
using namespace constants;


hamiltonian::hamiltonian(double mu)
{
	mp_.mu    = mu;
	mp_.tp[0] = -0.610;
	mp_.tp[1] = 0.300;
	mp_.tp[2] = 0.010;
	mp_.tp[3] = 0.040;
	mp_.tp[4] = 0.140;
	mp_.tp[5] = 0.020;
	mp_.tp[6] = 0.250;
	mp_.gauge = 0;
	dim_ham_ = 4;
  sym_ = sym_rep();
}

hamiltonian::hamiltonian(const hamiltonian &other) : sym_(other.sym_)
{
	mp_      = other.mp_;
	dim_ham_ = other.dim_ham_;
}

hamiltonian& hamiltonian :: operator= (const hamiltonian &other)
{
	// avoid self-assignment
	if ( this != &other )
  {       
		mp_      = other.mp_;
		dim_ham_ = other.dim_ham_;
    sym_     = other.sym_;
	}
	return *this;
}

int hamiltonian::get_dim() const
{
	return dim_ham_;
}

void hamiltonian::set_model_pars(double mu, double tp[7], int gauge)
{
	mp_.mu    = mu;
	mp_.tp[0] = tp[0];
	mp_.tp[1] = tp[1];
	mp_.tp[2] = tp[2];
	mp_.tp[3] = tp[3];
	mp_.tp[4] = tp[4];
	mp_.tp[5] = tp[5];
	mp_.tp[6] = tp[6];
  mp_.gauge = gauge;
}

struct model_pars hamiltonian::get_model_pars() const
{
	return mp_;
}

MatrixXcd hamiltonian::get_ham_mat(Vector3d kk) const
{	
	return Hmat_help(kk, 1);
}

MatrixXcd hamiltonian::get_ham_mat(double k1, double k2, double k3) const
{	
	Vector3d kk(k1,k2,k3);
	return Hmat_help(kk, 1);
}

VectorXcd hamiltonian::get_ham_eigvec(Vector3d kk, int m) const
{
	MatrixXcd mat = get_ham_eigmat(kk);
  assert (0 <= m && m < dim_ham_ && "invalid dimension in get_ham_eigvec");
	return mat.col(m);
}

VectorXd hamiltonian::get_ham_eigval(Vector3d kk) const
{	
	MatrixXcd mat  = get_ham_mat(kk);
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> hermitian_matrix(mat);
	assert( hermitian_matrix.info() == 0 ); 
	return hermitian_matrix.eigenvalues();
}

double hamiltonian::get_ham_eigval(Vector3d kk, int m) const
{	
	MatrixXcd mat  = get_ham_mat(kk);
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> hermitian_matrix(mat);
	assert( hermitian_matrix.info() == 0 ); 
  assert (0 <= m && m < dim_ham_ && "invalid dimension in get_ham_eigval");
	return hermitian_matrix.eigenvalues()(m);
}


MatrixXcd hamiltonian::get_ham_eigmat(Vector3d kk) const
{	
  cmplx ii(0.0,1.0);
	MatrixXcd eigmat(dim_ham_,dim_ham_);

	double phi;
  Vector3d kks; 
  kks << copysign(1.0, kk(1))*kk(0), abs(kk(1)), 0.0;
   	
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> herm_mat(Hmat_help(kks,1));
	assert( herm_mat.info() == 0 ); 
  MatrixXcd eigvecs = herm_mat.eigenvectors();

  MatrixXcd T = sym_.get_Tinv();

	assert (mp_.gauge == 0 && "invalid gauge parameter in get_ham_eigmat");
	for ( int j = 0; j < eigvecs.cols(); j++ )
	{	
   	phi = arg(eigvecs(j,j));
   	eigvecs.col(j) = eigvecs.col(j)*exp(-ii*phi);
 	}

  if (kk(1) < 0.0)
  {
   	eigvecs = T*(eigvecs.conjugate());
 	}

	return eigvecs;
}


MatrixXcd hamiltonian::Hmat_help(Vector3d kk, int flag) const
{	
	double kx = kk(0);
	double ky = kk(1);
	cmplx ii(0.0,1.0);

  double tp[7];

	tp[0] = mp_.tp[0] - mp_.mu;
	tp[1] = mp_.tp[1];
	tp[2] = mp_.tp[2];
	tp[3] = mp_.tp[3];
	tp[4] = mp_.tp[4];
	tp[5] = mp_.tp[5];
	tp[6] = mp_.tp[6];

	MatrixXcd mat(dim_ham_,dim_ham_);
	MatrixXcd mat0(dim_ham_/2,dim_ham_/2);
	mat.setZero();
	mat0.setZero();

	double H0_s  = pow(tp[5],2.0)/(4*tp[6]) - (tp[1]+tp[2])*(cos(kx)+cos(ky))
	             + pow(tp[3],2.0)/(2*tp[6])*(pow(sin(kx),2.0)+pow(sin(ky),2.0)) - tp[0];

	double H0_xy = (tp[4]+pow(tp[3],2)/tp[6])*sin(kx)*sin(ky);
	double H0_x2_y2 = (-tp[1]+tp[2] + pow(tp[3],2.0)/(2*tp[6])*(cos(kx)+cos(ky)))*(cos(kx)-cos(ky));
	
	MatrixXcd H0, HH;
	H0.resize(2,2);
	HH.resize(4,4);
	H0 << cmplx(+H0_s+H0_x2_y2,0.0), cmplx(+H0_xy,0.0),
	      cmplx(+H0_xy,0.0),         cmplx(+H0_s-H0_x2_y2,0.0);

  assert ( 0 <= flag && flag <= 1 && "invalid flag in Hmat_help"); 

	if (flag == 0)
	{
		return H0;
	}
	else if (flag == 1)
	{
		HH = kroneckerProduct(H0,constants::pauli_0);
		
		double alpha = tp[3]*tp[5]/(2*tp[6]);
		
		MatrixXcd H_soc = MatrixXcd::Zero(4,4);
		
		double Hso0 = ( tp[5]/2.0 + pow(tp[5],2.0)/(4*tp[6]) );
		
		H_soc = Hso0*kroneckerProduct(constants::pauli_y,constants::pauli_z)
		      + alpha*sin(kx)*kroneckerProduct(constants::pauli_0,constants::pauli_y)
		      - alpha*sin(ky)*kroneckerProduct(constants::pauli_0,constants::pauli_x)
		      - alpha*sin(kx)*kroneckerProduct(constants::pauli_x,constants::pauli_x)
		      + alpha*sin(ky)*kroneckerProduct(constants::pauli_x,constants::pauli_y)
		      - alpha*sin(kx)*kroneckerProduct(constants::pauli_z,constants::pauli_y)
		      - alpha*sin(ky)*kroneckerProduct(constants::pauli_z,constants::pauli_x);
		
		HH = HH + H_soc;
		return HH;
	}
	else
	{
		return HH;
	}
}

void hamiltonian::ham_test(const char* filename, int ngrid, int bix) const
{
   ofstream file_output(filename);
   int dim = dim_ham_;
   int width1 = 13;
   int width2 = 18;
   double kx, ky, res;
   Vector3d kk;
   MatrixXcd HH, eigvecs;
   VectorXcd vv;
   cmplx lam;
   file_output.precision(8); 
   file_output << setw(width1) << right << "kx" << setw(width1) << right << "ky" << setw(4) << right << "bix";
   file_output << setw(width2) << left << " state:";  
   file_output << "\n";
   
   for (int i = 0; i < ngrid; i++)
   {
    kx = -M_PI + M_PI*(2*i + 1)/ngrid;
    for (int j = 0; j < ngrid; j++)
	  { 
     ky = -M_PI + M_PI*(2*j + 1)/ngrid;
     kk << kx, ky, 0.0;
	   HH  = get_ham_mat(kk);
	   vv  = get_ham_eigvec(kk, bix);
     lam = get_ham_eigval(kk, bix);
     res = (HH*vv - lam*vv).norm();
     //cout << kx << " " << ky << " " << res << endl;
     assert (res < 1.0e-8 && "eigvec test failed in ham_test");
   	 file_output << setw(width1) << right << kx;
		 file_output << setw(width1) << right << ky;
   	 file_output << setw(4) << right << bix ;
		 for (int kk = 0; kk < dim; kk++)
		 {
		  file_output << setw(width2) << right << (vv(kk)).real();
			file_output << setw(width2) << right << (vv(kk)).imag();
		 }	
		 file_output << "\n";
		}
   }
   file_output.close();
}

void hamiltonian::sym_test(Vector3d kk) const
{
  double err;
	Matrix3d kmat;
	Matrix4cd UU, res_mat;
	for (int i = 0; i < sym_.get_sym_num(); i++)
	{
		kmat    = sym_.get_k_op(i);
    UU      = sym_.get_state_op(i);
		res_mat = UU*get_ham_mat(kk)*UU.inverse() - get_ham_mat(kmat*kk);
    err     = res_mat.norm();
    if (err > 1.0E-8) 
		{
			cout << "in operation: " << i << " symmetry check gives an error diff_norm = " << err << endl;
      cout << "kmat: " << kmat << endl;
			cout << "state_mat: " << UU << endl;
		}
	}
}
