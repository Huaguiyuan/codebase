#ifndef SYM
#define SYM

#include "general_header.hpp"

class sym
{
	protected:
		Matrix3d rot_mat(double ang, Vector3d nvec);
		Matrix3d ref_mat(Vector3d nvec);
		Matrix2cd spin_rot_mat(double ang, Vector3d nvec);
		Matrix2cd spin_ref_mat(Vector3d nvec);

  public:
		sym();
};
#endif

#ifndef SYM_REP
#define SYM_REP

typedef Matrix<cmplx,2,2> OrbMat;
typedef std::vector<OrbMat> OrbMat_vec;
typedef Matrix<cmplx,4,4> StateMat;
typedef std::vector<StateMat> StateMat_vec;

class sym_rep : sym
{
	private:
    int sym_num_;
    Matrix3d_vec  k_ops_;
    Matrix2cd_vec spin_ops_;
    Matrix2cd_vec orb_ops_;
    Matrix4cd_vec state_ops_;

	public:
    sym_rep();
		int get_sym_num() const;

    Matrix3d_vec  get_k_ops() const;
    Matrix2cd_vec get_spin_ops() const;
    Matrix2cd_vec get_orb_ops() const;
    Matrix4cd_vec get_state_ops() const;

    Matrix3d  get_k_op(int op) const;
    Matrix2cd get_spin_op(int op) const;
    OrbMat    get_orb_op(int op) const;
    StateMat  get_state_op(int op) const;
 
    StateMat  get_Tinv() const;
};
#endif
