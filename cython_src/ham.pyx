import numpy as np
cdef extern from "cpp_MatrixXcd_py.h":
			cdef cppclass VectorXd_py:
				VectorXd_py()
				VectorXd_py(VectorXd_py other)
				int size()
				complex coeff(int)

cdef extern from "cpp_MatrixXcd_py.h":
			cdef cppclass VectorXcd_py:
				VectorXcd_py()
				VectorXcd_py(VectorXcd_py other)
				int size()
				complex coeff(int)

cdef extern from "cpp_MatrixXcd_py.h":
			cdef cppclass Vector3d_py:
				Vector3d_py()
				Vector3d_py(float kx, float ky, float kz)
				Vector3d_py(Vector3d_py other)

cdef extern from "cpp_MatrixXcd_py.h":
			cdef cppclass MatrixXcd_py:
				MatrixXcd_py()
				MatrixXcd_py(int d1, int s2)
				MatrixXcd_py(MatrixXcd_py other)
				int rows()
				int cols()
				complex coeff(int,int)

cdef extern from "hamiltonian.hpp": 
		cdef cppclass hamiltonian:
			hamiltonian(float mu)
			int get_dim() const
			MatrixXcd_py get_ham_mat(Vector3d_py kk)
			MatrixXcd_py get_ham_eigmat(Vector3d_py kk)
			VectorXd_py  get_ham_eigval(Vector3d_py kk)

cdef class ham_py:
		cdef hamiltonian* thisptr # hold a C++ instance

		def __cinit__(self, mu):
			self.thisptr = new hamiltonian( mu)

		def __dealloc__(self):
			del self.thisptr

		def get_ham_mat(self,kx,ky,kz):
			cdef Vector3d_py kk   = Vector3d_py(kx,ky,kz)
			cdef MatrixXcd_py mat = self.thisptr.get_ham_mat(kk)
			res_mat = np.zeros((mat.rows(),mat.cols()), dtype = complex)
			for row in range(mat.rows()):
				for col in range(mat.cols()):
					res_mat[row, col] = mat.coeff(row, col)
			return res_mat 

		def get_ham_eigmat(self,kx,ky,kz):
			cdef Vector3d_py kk   = Vector3d_py(kx,ky,kz)
			cdef MatrixXcd_py mat = self.thisptr.get_ham_eigmat(kk)
			res_mat = np.zeros((mat.rows(),mat.cols()), dtype = complex)
			for row in range(mat.rows()):
				for col in range(mat.cols()):
					res_mat[row, col] = mat.coeff(row, col)
			return res_mat 

		def get_ham_eigval(self,kx,ky,kz):
			cdef Vector3d_py kk   = Vector3d_py(kx,ky,kz)
			cdef VectorXd_py vec  = self.thisptr.get_ham_eigval(kk)
			res_vec = np.zeros(vec.size(), dtype = float)
			for jj in range(vec.size()):
				dum = vec.coeff(jj)
				res_vec[jj] = np.real(dum)
			return res_vec 