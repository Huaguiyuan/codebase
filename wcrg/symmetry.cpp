#include "symmetry.hpp"

sym_rep::sym_rep():sym()
{
	irred_angle_range_ << 0.0, M_PI/1;
  cmplx I(0.0,1.0);
  Matrix3d kmat;
  Matrix2cd smat;
  OrbMat omat;
  Vector3d kk;
  kk << 0.0, 0.0, 1.0;
  //identity
  kmat = sym::rot_mat(0.0,kk);
  smat = sym::spin_rot_mat(0.0,kk);
  omat << 1, 0, 0, 0,
          0, 1, 0, 0,
          0, 0, 1, 0,
          0, 0, 0, 1;
  k_ops_.push_back(kmat);
  spin_ops_.push_back(smat);
  orb_ops_.push_back(omat);
  state_ops_.push_back(kroneckerProduct(omat,smat));
  // pi rotation
  kmat = sym::rot_mat(M_PI,kk);
  smat = sym::spin_rot_mat(M_PI,kk);
  omat << 0, 0, 1, 0,
          0, 0, 0, 1,
          1, 0, 0, 0,
          0, 1, 0, 0;
  k_ops_.push_back(kmat);
  spin_ops_.push_back(smat);
  orb_ops_.push_back(omat);
  state_ops_.push_back(kroneckerProduct(omat,smat));
  // x reflection
  kk << 1.0, 0.0, 0.0;
  kmat = sym::ref_mat(kk);
  smat = sym::spin_ref_mat(kk);
  omat << 0, 0, 1, 0,
          0, 0, 0, 1,
          1, 0, 0, 0,
          0, 1, 0, 0;
  k_ops_.push_back(kmat);
  spin_ops_.push_back(smat);
  orb_ops_.push_back(omat);
  state_ops_.push_back(kroneckerProduct(omat,smat));
  // y reflection
  kk << 0.0, 1.0, 0.0;
  kmat = sym::ref_mat(kk);
  smat = sym::spin_ref_mat(kk);
  omat << 1, 0, 0, 0,
          0, 1, 0, 0,
          0, 0, 1, 0,
          0, 0, 0, 1;
  k_ops_.push_back(kmat);
  spin_ops_.push_back(smat);
  orb_ops_.push_back(omat);
  state_ops_.push_back(kroneckerProduct(omat,smat));

  sym_num_ = k_ops_.size();
}

int sym_rep::get_sym_num() const
{
	return sym_num_;
}

Vector2d sym_rep::get_irred_angle_range() const
{
	return irred_angle_range_;
}

Matrix3d_vec sym_rep::get_k_ops() const
{
	return k_ops_;
}

Matrix2cd_vec sym_rep::get_spin_ops() const
{
	return spin_ops_;
}

OrbMat_vec sym_rep::get_orb_ops() const
{
	return orb_ops_;
}

StateMat_vec sym_rep::get_state_ops() const
{
	return state_ops_;
}

Matrix3d sym_rep::get_k_op(int op) const
{
	return k_ops_[op];
}

Matrix2cd sym_rep::get_spin_op(int op) const
{
	return spin_ops_[op];
}

OrbMat sym_rep::get_orb_op(int op) const
{
	return orb_ops_[op];
}

StateMat sym_rep::get_state_op(int op) const
{
	return state_ops_[op];
}

void sym_rep::sym_test(Vector3d kk, const hamiltonian &H) const
{
  double err;
	Matrix3d kmat;
	StateMat UU, res_mat;
	for (int i = 0; i < sym_num_; i++)
	{
		kmat    = k_ops_[i];
    UU      = state_ops_[i];
		res_mat = UU*H.get_ham_mat(kk)*UU.inverse() - H.get_ham_mat(kmat*kk);
    err     = res_mat.norm();
    if (err > 1.0E-8) 
		{
			cout << "in operation: " << i << " symmetry check gives an error diff_norm = " << err << endl;
      cout << "kmat: " << kmat << endl;
			cout << "state_mat: " << UU << endl;
		}
	}
}

sym::sym()
{
}

Matrix3d sym::rot_mat(double ang, Vector3d nvec)
{
	nvec = nvec.normalized();
  Matrix3d Id, Ux, UU;
  Id.setIdentity();
  
  double &x = nvec(0);
  double &y = nvec(1);
  double &z = nvec(2);
 
	Ux << 0.0,  -z,  y, 
          z, 0.0, -x,
         -y,   x, 0.0;

  UU << x*x, x*y, x*z, 
        x*y, y*y, y*z,
        x*z, y*z, z*z;

  return cos(ang)*Id + sin(ang)*Ux + (1.0 - cos(ang))*UU;
}

Matrix2cd sym::spin_rot_mat(double ang, Vector3d nvec)
{
  cmplx I(0.0,1.0); 
	nvec = nvec.normalized();
  Matrix2cd rot_mat;

  rot_mat << cos(ang/2) - I*nvec(2)*sin(ang/2), (-I*nvec(0) - nvec(1))*sin(ang/2),
                   (-I*nvec(0) + nvec(1))*sin(ang/2), cos(ang/2) + I*nvec(2)*sin(ang/2); 

  return rot_mat;
}

Matrix3d sym::ref_mat(Vector3d nvec)
{
	nvec = nvec.normalized();
  Matrix3d refmat;

  double &x = nvec(0);
  double &y = nvec(1);
  double &z = nvec(2);

	refmat  << 1.0-2.0*x*x,    -2.0*x*y,    -2.0*x*z, 
			          -2.0*x*y, 1.0-2.0*y*y,    -2.0*y*z,
				        -2.0*x*z,    -2.0*y*z, 1.0-2.0*z*z;
  return refmat;
}

Matrix2cd sym::spin_ref_mat(Vector3d nvec)
{
	nvec = nvec.normalized();
  Vector3d ez(0.0,0.0,1.0);
  double ang = acos(nvec.dot(ez));
  Vector3d rvec = nvec.cross(ez);

  return spin_rot_mat(-ang,rvec)*constants::pauli_z*spin_rot_mat(ang,rvec);
}
