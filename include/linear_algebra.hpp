#ifndef LINEAR_ALGEBRA_H_
#define LINEAR_ALGEBRA_H_

#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/Eigen>
#include <Eigen/IterativeLinearSolvers>

typedef float Scalar;

typedef Eigen::VectorX<Scalar> Vector;

typedef Eigen::MatrixX<Scalar> Mat;

typedef Eigen::Vector3<Scalar> Vec3;

typedef Eigen::Matrix3<Scalar> Mat3;

typedef Eigen::Matrix4<Scalar> Mat4;

typedef Eigen::SparseMatrix<Scalar> SparseMatrix;

typedef Eigen::Triplet<Scalar> Triplet;

#endif // LINEAR_ALGEBRA_H_
