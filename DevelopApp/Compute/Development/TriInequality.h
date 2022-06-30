#pragma once
#include <iostream>
#include <Eigen/Eigen>
using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::MatrixXi;
using SPMatrixXd = Eigen::SparseMatrix<double, Eigen::ColMajor>;

class TriInequality
{
public:
	TriInequality();
	TriInequality(const VectorXd& l0_, const MatrixXi& E_, double similar, const VectorXd& eps_);
	void init(const VectorXd& l0_, const MatrixXi& E_, double similar_, const VectorXd& eps_);
	VectorXd opt(const VectorXd& lt_, int maxIters = 1e3);

private:
	VectorXd m_l, m_lt, m_l0;
	MatrixXi m_E;

	MatrixXd m_SL, m_lambda;
	Eigen::SimplicialLDLT<SPMatrixXd> m_solver;

public:
	VectorXd m_eps;
	double m_s;

private:
	void updateL();
	void updateSL();
	void updateLambda();
	bool istri();
	bool converge();
};



