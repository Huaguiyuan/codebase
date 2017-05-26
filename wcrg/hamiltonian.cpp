#include "hamiltonian.hpp"
using namespace constants;

hamiltonian::hamiltonian(double mu, double soc)
{
	mp_.mu  = mu;
	mp_.soc = soc;
	mp_.t   = 1.0;
	mp_.eta = 1.2;
	mp_.ts  =  0.072*mp_.eta*mp_.t;
	mp_.tp  = -0.048*mp_.eta*mp_.t;
	mp_.gauge = 2;
	dim_ham_ = 8;
}

hamiltonian::hamiltonian(const hamiltonian &other)
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
	}
	return *this;
}

int hamiltonian::get_dim() const
{
	return dim_ham_;
}

void hamiltonian::set_model_pars(double mu, double soc, double t, double eta, double ts, double tp, int gauge)
{
	mp_.mu    = mu;
	mp_.soc   = soc;
	mp_.t     = t;
  mp_.eta   = eta;
	mp_.ts    = ts;
  mp_.tp    = tp;
  mp_.gauge = gauge;
}

MatrixXcd hamiltonian::get_ham_mat(Vector3d kk) const
{	
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
	MatrixXcd eigmat(dim_ham_,dim_ham_);
	cmplx ii(0.0,1.0);
	int real_ind = 0;
	double phi;
	Vector3d kk_abs;
  kk_abs << abs(kk(0)), abs(kk(1)), abs(kk(2));
	MatrixXcd mat0 = Hmat_help(kk_abs,0);
	MatrixXcd mat  = Hmat_help(kk_abs,1);
	VectorXd  eigvals0, eigvals;
	MatrixXcd eigvecs0, eigvecs;
	MatrixXcd VV(dim_ham_,mat0.cols());

	Matrix4cd Mo_0, Mo_1;
	Mo_0.setIdentity();
	Mo_1.setZero();
	Mo_1(0,2) = 1.0;
	Mo_1(1,3) = 1.0;
	Mo_1(2,0) = 1.0;
	Mo_1(3,1) = 1.0;

	MatrixXcd MM, M;
	MM = kroneckerProduct(Mo_1,pauli_0);
	M  = kroneckerProduct(Mo_0,pauli_x);

	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> herm_mat0(mat0);
	assert( herm_mat0.info() == 0 ); 
	eigvals0 = herm_mat0.eigenvalues();
  eigvecs0 = herm_mat0.eigenvectors();
   	
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> herm_mat(mat);
	assert( herm_mat.info() == 0 ); 
	eigvals = herm_mat.eigenvalues();
  eigvecs = herm_mat.eigenvectors();

	VV.setZero();
	for ( int j = 0; j < mat0.cols(); j++ )
	{	
    phi = arg(eigvecs0(real_ind,j));
    for ( int i = 0; i < mat0.rows(); i++ )
    {
 	 		 VV(2*i,j) = exp(-ii*phi)*eigvecs0(i,j);
    }
  } 
        
	MatrixXcd lstmat(dim_ham_,2);
	VectorXcd res(2);
	VectorXcd v1(dim_ham_), v2(dim_ham_);

	assert (mp_.gauge == 0 || mp_.gauge == 1 || mp_.gauge == 2 || "invalid gauge parameter in get_ham_eigmat");

	for ( int i = 0; i < dim_ham_/2; i++)
	{
	    if (mp_.gauge == 0)
	    {
	    	lstmat.col(0) = eigvecs.col(2*i);
	    	lstmat.col(1) = eigvecs.col(2*i+1);
	    	res = lstmat.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(VV.col(i));
	    	v1  = lstmat*res;
	    	v1  = v1/v1.norm();
	    	phi = arg(v1(0));
            	v1  = exp(-ii*phi)*v1;
	    	v2  = M*v1;
				if ( kk(0) >= 0.0 && kk(1) >= 0.0)
				{
					eigmat.col(2*i)   = v1;
					eigmat.col(2*i+1) = MM*v2.conjugate();
				}
				else if ( kk(0) < 0.0 && kk(1) >= 0.0)
				{
					eigmat.col(2*i)   = v1.conjugate();
					eigmat.col(2*i+1) = MM*v2;
				}
				else if ( kk(0) < 0.0 && kk(1) < 0.0)
				{
					eigmat.col(2*i)   = MM*v1;
					eigmat.col(2*i+1) = v2.conjugate();
				}
				else if ( kk(0) >= 0.0 && kk(1) < 0.0)
				{
					eigmat.col(2*i)   = MM*v1.conjugate();
					eigmat.col(2*i+1) = v2;
				}
	    }
	    if (mp_.gauge == 1)
	    {
	    	lstmat.col(0) = eigvecs.col(2*i);
	    	lstmat.col(1) = eigvecs.col(2*i+1);

	    	res = lstmat.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(VV.col(i));
	    	v1  = lstmat*res;
	    	v1  = v1/v1.norm();
	    	phi = arg(v1(0));
       	v1  = exp(-ii*phi)*v1;

	    	res = lstmat.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(M*(VV.col(i)));
	    	v2  = lstmat*res;
	    	v2  = v2/v2.norm();
	    	phi = arg(v2(1));
        v2  = exp(-ii*phi)*v2;

				if ( kk(0) >= 0.0 && kk(1) >= 0.0)
				{
					eigmat.col(2*i)   = v1;
					eigmat.col(2*i+1) = v2;
				}
				else if ( kk(0) < 0.0 && kk(1) >= 0.0)
				{
					eigmat.col(2*i)   = MM*M*v2;
					eigmat.col(2*i+1) = MM*M*v1;
				}
				else if ( kk(0) < 0.0 && kk(1) < 0.0)
				{
					eigmat.col(2*i)   = M*v2.conjugate();
					eigmat.col(2*i+1) = M*v1.conjugate();
				}
				else if ( kk(0) >= 0.0 && kk(1) < 0.0)
				{
					eigmat.col(2*i)   = MM*v1.conjugate();
					eigmat.col(2*i+1) = MM*v2.conjugate();
				}
	    }
	    if (mp_.gauge == 2)
	    {
	    	lstmat.col(0) = eigvecs.col(2*i);
	    	lstmat.col(1) = eigvecs.col(2*i+1);

	    	res = lstmat.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(VV.col(i));
	    	v1  = lstmat*res;
	    	v1  = v1/v1.norm();
	    	phi = arg(v1(0));
       	v1  = exp(-ii*phi)*v1;

	    	res = lstmat.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(M*(VV.col(i)));
	    	v2  = lstmat*res;
	    	v2  = v2/v2.norm();
	    	phi = arg(v2(1));
       	v2  = exp(-ii*phi)*v2;

				if ( kk(0) >= 0.0 && kk(1) >= 0.0)
				{
					eigmat.col(2*i)   = v1;
					eigmat.col(2*i+1) = v2;
				}
				else if ( kk(0) < 0.0 && kk(1) >= 0.0)
				{
					eigmat.col(2*i)   = v1.conjugate();
					eigmat.col(2*i+1) = v2.conjugate();
				}
				else if ( kk(0) < 0.0 && kk(1) < 0.0)
				{
					eigmat.col(2*i)   = (M*v2).conjugate();
					eigmat.col(2*i+1) = (M*v1).conjugate();
				}
				else if ( kk(0) >= 0.0 && kk(1) < 0.0)
				{
					eigmat.col(2*i)   = M*v2;
					eigmat.col(2*i+1) = M*v1;
				}
	    }

	} 
	return eigmat;
}

MatrixXcd hamiltonian::Hmat_help(Vector3d kk, int flag) const
{	
	double kx = kk(0);
	double ky = kk(1);
	cmplx ii(0.0,1.0);

	double mu, soc, t, eta, ts, tp;

	mu  = mp_.mu;
	soc = mp_.soc;
	t   = mp_.t;
	eta = mp_.eta;
	ts  = mp_.ts;
	tp  = mp_.tp;

	MatrixXcd id0(dim_ham_/2,dim_ham_/2);
	MatrixXcd id(dim_ham_,dim_ham_);
	id0.setIdentity();
	id.setIdentity();

	MatrixXcd mat(dim_ham_,dim_ham_);
	MatrixXcd mat0(dim_ham_/2,dim_ham_/2);
	mat.setZero();
	mat0.setZero();

  mat0(0,1) = tp*exp(ii*kx/2.0);
  mat0(1,0) = conj(mat0(0,1)); 
  mat0(0,2) = 2.0*ts*exp(-ii*kx/2.0)*cos(ky/2.0);
  mat0(2,0) = conj(mat0(0,2));
  mat0(0,3) = 2.0*t*cos(ky/2.0);
  mat0(3,0) = conj(mat0(0,3));
  mat0(1,2) = 2.0*t*cos(ky/2.0);
  mat0(2,1) = conj(mat0(1,2));
  mat0(2,3) = tp*exp(-ii*kx/2.0);
  mat0(3,2) = conj(mat0(2,3));

  assert ( 0 <= flag && flag <= 1 && "invalid flag in Hmat_help"); 

	if (flag == 0)
	{
	
		mat0 = -1.0*mat0;
		mat0 = mat0 - mu*id0;
		return mat0;
	}
	else if (flag == 1)
	{
		mat.block(0,2,2,2) = mat0(0,1)*pauli_0;
		mat.block(2,0,2,2) = mat.block(0,2,2,2).conjugate();

		mat.block(0,4,2,2) = mat0(0,2)*pauli_0;
		mat.block(4,0,2,2) = mat.block(0,4,2,2).conjugate();
    
		mat.block(0,6,2,2) = mat0(0,3)*pauli_0 - soc*sin(ky/2.0)*pauli_z;
		mat.block(6,0,2,2) = mat.block(0,6,2,2).conjugate();

		mat.block(2,4,2,2) = mat0(1,2)*pauli_0 + soc*sin(ky/2.0)*pauli_z;
		mat.block(4,2,2,2) = mat.block(2,4,2,2).conjugate();

		mat.block(4,6,2,2) = mat0(2,3)*pauli_0;
		mat.block(6,4,2,2) = mat.block(4,6,2,2).conjugate();

		mat = -1.0*mat;
		mat = mat - mu*id;
		return mat;
	}
	else
	{
		return mat;
	}
}
