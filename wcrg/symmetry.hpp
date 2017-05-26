#ifndef SYM
#define SYM

#include "general_header.hpp"
#include "hamiltonian.hpp"

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

typedef Matrix<cmplx,4,4> OrbMat;
typedef std::vector<OrbMat> OrbMat_vec;
typedef Matrix<cmplx,8,8> StateMat;
typedef std::vector<StateMat> StateMat_vec;

class sym_rep : sym
{
	private:
		Vector2d irred_angle_range_;
    int sym_num_;
    Matrix3d_vec  k_ops_;
    Matrix2cd_vec spin_ops_;
    OrbMat_vec orb_ops_;
    StateMat_vec state_ops_;

	public:
    sym_rep();
		Vector2d get_irred_angle_range() const;
		int get_sym_num() const;

    Matrix3d_vec  get_k_ops() const;
    Matrix2cd_vec get_spin_ops() const;
    OrbMat_vec get_orb_ops() const;
    StateMat_vec get_state_ops() const;

    Matrix3d  get_k_op(int op) const;
    Matrix2cd get_spin_op(int op) const;
    OrbMat get_orb_op(int op) const;
    StateMat get_state_op(int op) const;

    void sym_test(Vector3d kk, const hamiltonian &H) const; 
};
#endif
