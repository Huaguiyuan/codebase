#include "./eigen/Eigen/Dense"
using namespace Eigen;

#ifndef MATRIXXCD_PY_H
#define MATRIXXCD_PY_H
 
class MatrixXcd_py : public MatrixXcd { 
    public: 
        MatrixXcd_py() : MatrixXcd() { }
        MatrixXcd_py(int rows,int cols) : MatrixXcd(rows,cols) { }
        MatrixXcd_py(const MatrixXcd other) : MatrixXcd(other) { } 
};
#endif

#ifndef VECTOR3D_PY_H
#define VECTOR3D_PY_H
 
class Vector3d_py : public Vector3d { 
    public: 
        Vector3d_py() : Vector3d() { }
        Vector3d_py(double kx,double ky, double kz) : Vector3d(kx,ky,kz) { }
        Vector3d_py(const Vector3d other) : Vector3d(other) { } 
};
#endif

#ifndef VECTORXCD_PY_H
#define VECTORXCD_PY_H
 
class VectorXcd_py : public VectorXcd { 
    public: 
        VectorXcd_py() : VectorXcd() { }
        VectorXcd_py(const VectorXcd other) : VectorXcd(other) { } 
        VectorXcd_py(int ind) : VectorXcd(ind) { }
};
#endif

#ifndef VECTORXD_PY_H
#define VECTORXD_PY_H
 
class VectorXd_py : public VectorXd { 
    public: 
        VectorXd_py() : VectorXd() { }
        VectorXd_py(const VectorXd other) : VectorXd(other) { } 
        VectorXd_py(int ind) : VectorXd(ind) { }
};
#endif