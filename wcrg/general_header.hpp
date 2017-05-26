#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <cstdlib>
#include <stdexcept>
#include <complex>
#include <cmath>
#include <random>
#include <vector>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/Ray_2.h>
#include <CGAL/Direction_2.h>
#include <CGAL/Segment_2.h>
#include <CGAL/intersections.h>
#include "./eigen/Eigen/Eigenvalues"
#include "./eigen/Eigen/Dense"
#include "./eigen/unsupported/Eigen/KroneckerProduct"

using namespace std;
using namespace Eigen;
typedef std::complex<double> cmplx;
typedef std::vector<int> int_vec;
typedef std::vector<Matrix3d> Matrix3d_vec;
typedef std::vector<Matrix2cd> Matrix2cd_vec;

typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_2 Point;
typedef CGAL::Polygon_2<K> Polygon;
typedef CGAL::Direction_2<K> Direct;
typedef CGAL::Ray_2<K> Ray;
typedef Polygon::Edge_const_iterator EdgeIterator;
typedef Polygon::Vertex_iterator VertexIterator;

#ifndef CONSTANTS
#define CONSTANTS
namespace constants
{
	const Matrix2cd pauli_0 = (Matrix2cd() << cmplx(+1.0,+0.0), cmplx(+0.0,+0.0),
	                                            cmplx(+0.0,+0.0), cmplx(+1.0,+0.0)).finished();
	const Matrix2cd pauli_x = (Matrix2cd() << cmplx(+0.0,+0.0), cmplx(+1.0,+0.0),
	                                            cmplx(+1.0,+0.0), cmplx(+0.0,+0.0)).finished();
	const Matrix2cd pauli_y = (Matrix2cd() << cmplx(+0.0,+0.0), cmplx(+0.0,-1.0),
	                                            cmplx(+0.0,+1.0), cmplx(+0.0,+0.0)).finished();
	const Matrix2cd pauli_z = (Matrix2cd() << cmplx(+1.0,+0.0), cmplx(+0.0,+0.0),
	                                            cmplx(+0.0,+0.0), cmplx(-1.0,+0.0)).finished();

}
#endif
