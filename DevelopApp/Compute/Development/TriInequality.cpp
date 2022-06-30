#include "TriInequality.h"

TriInequality::TriInequality()
{

}

TriInequality::TriInequality(const VectorXd& l0_, const MatrixXi& E_, double similar_, const VectorXd& eps_)
{
	init(l0_, E_, similar_, eps_);
	m_SL = MatrixXd::Zero(9, m_E.cols());
	m_lambda = MatrixXd::Zero(9, m_E.cols());
}

void TriInequality::init(const VectorXd& l0_, const MatrixXi& E_, double similar_, const VectorXd& eps_)
{
	m_l0 = l0_;
	//m_l0 = l0_.mean() * VectorXd::Ones(l0_.size());
	m_E = E_;
	m_s = similar_;
	m_eps = eps_;
	m_SL = MatrixXd::Zero(9, m_E.cols());
	m_lambda = MatrixXd::Zero(9, m_E.cols());


	std::vector<Eigen::Triplet<double>> trips;
	trips.reserve(9 * E_.cols());
	for (int i = 0; i < E_.cols(); ++i)
	{
		int e0 = E_(0, i);
		int e1 = E_(1, i);
		int e2 = E_(2, i);

		trips.emplace_back(e0, e0, 3 + m_s * (m_l0(e1) + m_l0(e2)) / m_l0(e0));
		trips.emplace_back(e1, e1, 3 + m_s * (m_l0(e0) + m_l0(e2)) / m_l0(e1));
		trips.emplace_back(e2, e2, 3 + m_s * (m_l0(e0) + m_l0(e1)) / m_l0(e2));

		trips.emplace_back(e0, e1, -m_s);
		trips.emplace_back(e1, e0, -m_s);
		trips.emplace_back(e0, e2, -m_s);
		trips.emplace_back(e2, e0, -m_s);
		trips.emplace_back(e1, e2, -m_s);
		trips.emplace_back(e2, e1, -m_s);
	}

	for (int i = 0; i < m_l0.size(); ++i)
	{
		trips.emplace_back(i, i, 1.0);
	}

	SPMatrixXd spm(m_l0.size(), m_l0.size());
	spm.setFromTriplets(trips.begin(), trips.end());
	m_solver.compute(spm);
}

VectorXd TriInequality::opt(const VectorXd& lt_, int maxIters)
{
	m_l = m_lt = lt_;
	m_SL.setZero();
	m_lambda.setZero();

	for (int i = 0; i < maxIters && !converge(); i++)
	{
		updateSL();
		updateL();
		updateLambda();
	}

	return m_l;
}

bool TriInequality::istri()
{
	for (int i = 0; i < m_E.cols(); i++)
	{
		if (m_l(m_E(0, i)) + m_l(m_E(1, i)) - m_l(m_E(2, i)) < m_eps(i) ||
			m_l(m_E(1, i)) + m_l(m_E(2, i)) - m_l(m_E(0, i)) < m_eps(i) ||
			m_l(m_E(2, i)) + m_l(m_E(0, i)) - m_l(m_E(1, i)) < m_eps(i))
			return false;
	}

	return false;
}

bool TriInequality::converge()
{
	//double res = 0;
	//for (int i = 0; i < m_E.cols(); i++)
	//{
	//	res += pow(m_l(m_E(0, i)) - m_SL(0, i), 2.0);
	//	res += pow(m_l(m_E(0, i)) - m_SL(1, i), 2.0);
	//	res += pow(m_l(m_E(0, i)) - m_SL(2, i), 2.0);

	//	res += pow(m_l(m_E(1, i)) - m_SL(3, i), 2.0);
	//	res += pow(m_l(m_E(1, i)) - m_SL(4, i), 2.0);
	//	res += pow(m_l(m_E(1, i)) - m_SL(5, i), 2.0);

	//	res += pow(m_l(m_E(2, i)) - m_SL(6, i), 2.0);
	//	res += pow(m_l(m_E(2, i)) - m_SL(7, i), 2.0);
	//	res += pow(m_l(m_E(2, i)) - m_SL(8, i), 2.0);
	//}

	// ===================================
	double res = 0;
	for (int i = 0; i < m_E.cols(); i++)
	{
		const double& m_l_0 = m_l(m_E(0, i));
		const double& m_l_1 = m_l(m_E(1, i));
		const double& m_l_2 = m_l(m_E(2, i));

		const Eigen::VectorXd m_SL_i = m_SL.col(i);

		res += (m_l_0 - m_SL_i(0)) * (m_l_0 - m_SL_i(0))
			+ (m_l_0 - m_SL_i(1)) * (m_l_0 - m_SL_i(1))
			+ (m_l_0 - m_SL_i(2)) * (m_l_0 - m_SL_i(2))
			+ (m_l_1 - m_SL_i(3)) * (m_l_1 - m_SL_i(3))
			+ (m_l_1 - m_SL_i(4)) * (m_l_1 - m_SL_i(4))
			+ (m_l_1 - m_SL_i(5)) * (m_l_1 - m_SL_i(5))
			+ (m_l_2 - m_SL_i(6)) * (m_l_2 - m_SL_i(6))
			+ (m_l_2 - m_SL_i(7)) * (m_l_2 - m_SL_i(7))
			+ (m_l_2 - m_SL_i(8)) * (m_l_2 - m_SL_i(8));
	}
	// ==================================

	return (res < m_eps.squaredNorm());
}

void TriInequality::updateSL()
{
	VectorXd r = VectorXd::Zero(m_l0.size());
	r = m_l.cwiseQuotient(m_l0);

	//for (int i = 0; i < m_E.cols(); i++)
	//{
	//	m_SL(0, i) = m_l(m_E(0, i)) + m_lambda(0, i);
	//	m_SL(1, i) = m_l(m_E(0, i)) + m_lambda(1, i);
	//	m_SL(2, i) = m_l(m_E(0, i)) + m_lambda(2, i);
	//	m_SL(0, i) = fmax(m_SL(0, i), m_l(m_E(2, i)) - m_l(m_E(1, i)) + m_eps(i));
	//	m_SL(1, i) = fmin(m_SL(1, i), m_l(m_E(1, i)) + m_l(m_E(2, i)) - m_eps(i));
	//	m_SL(2, i) = fmax(m_SL(2, i), m_l(m_E(1, i)) - m_l(m_E(2, i)) + m_eps(i));

	//	m_SL(3, i) = m_l(m_E(1, i)) + m_lambda(3, i);
	//	m_SL(4, i) = m_l(m_E(1, i)) + m_lambda(4, i);
	//	m_SL(5, i) = m_l(m_E(1, i)) + m_lambda(5, i);
	//	m_SL(3, i) = fmax(m_SL(3, i), m_l(m_E(2, i)) - m_l(m_E(0, i)) + m_eps(i));
	//	m_SL(4, i) = fmax(m_SL(4, i), m_l(m_E(0, i)) - m_l(m_E(2, i)) + m_eps(i));
	//	m_SL(5, i) = fmin(m_SL(5, i), m_l(m_E(0, i)) + m_l(m_E(2, i)) - m_eps(i));

	//	m_SL(6, i) = m_l(m_E(2, i)) + m_lambda(6, i);
	//	m_SL(7, i) = m_l(m_E(2, i)) + m_lambda(7, i);
	//	m_SL(8, i) = m_l(m_E(2, i)) + m_lambda(8, i);
	//	m_SL(6, i) = fmin(m_SL(6, i), m_l(m_E(0, i)) + m_l(m_E(1, i)) - m_eps(i));
	//	m_SL(7, i) = fmax(m_SL(7, i), m_l(m_E(0, i)) - m_l(m_E(1, i)) + m_eps(i));
	//	m_SL(8, i) = fmax(m_SL(8, i), m_l(m_E(1, i)) - m_l(m_E(0, i)) + m_eps(i));
	//}

	// ================================================
	m_SL = m_lambda;
	for (int i = 0; i < m_E.cols(); i++)
	{
		double& m_SL_0 = m_SL(0, i);
		double& m_SL_1 = m_SL(1, i);
		double& m_SL_2 = m_SL(2, i);
		double& m_SL_3 = m_SL(3, i);
		double& m_SL_4 = m_SL(4, i);
		double& m_SL_5 = m_SL(5, i);
		double& m_SL_6 = m_SL(6, i);
		double& m_SL_7 = m_SL(7, i);
		double& m_SL_8 = m_SL(8, i);

		const double& m_l_0 = m_l(m_E(0, i));
		const double& m_l_1 = m_l(m_E(1, i));
		const double& m_l_2 = m_l(m_E(2, i));

		const double& m_eps_i = m_eps(i);

		m_SL_0 += m_l_0;
		m_SL_1 += m_l_0;
		m_SL_2 += m_l_0;
		m_SL_3 += m_l_1;
		m_SL_4 += m_l_1;
		m_SL_5 += m_l_1;
		m_SL_6 += m_l_2;
		m_SL_7 += m_l_2;
		m_SL_8 += m_l_2;

		m_SL_0 = fmax(m_SL_0, m_l_2 - m_l_1 + m_eps_i);
		m_SL_1 = fmin(m_SL_1, m_l_1 + m_l_2 - m_eps_i);
		m_SL_2 = fmax(m_SL_2, m_l_1 - m_l_2 + m_eps_i);

		m_SL_3 = fmax(m_SL_3, m_l_2 - m_l_0 + m_eps_i);
		m_SL_4 = fmax(m_SL_4, m_l_0 - m_l_2 + m_eps_i);
		m_SL_5 = fmin(m_SL_5, m_l_2 + m_l_0 - m_eps_i);
		
		m_SL_6 = fmin(m_SL_6, m_l_0 + m_l_1 - m_eps_i);
		m_SL_7 = fmax(m_SL_7, m_l_0 - m_l_1 + m_eps_i);
		m_SL_8 = fmax(m_SL_8, m_l_1 - m_l_0 + m_eps_i);
	}
	//=============================================
}

void TriInequality::updateL()
{
	m_l = m_lt;
	VectorXd ln = VectorXd::Ones(m_l.size());

	//for (int i = 0; i < m_E.cols(); i++)
	//{
	//	m_l(m_E(0, i)) += m_SL(0, i) - m_lambda(0, i);
	//	m_l(m_E(0, i)) += m_SL(1, i) - m_lambda(1, i);
	//	m_l(m_E(0, i)) += m_SL(2, i) - m_lambda(2, i);
	//	ln(m_E(0, i)) += 3;

	//	m_l(m_E(1, i)) += m_SL(3, i) - m_lambda(3, i);
	//	m_l(m_E(1, i)) += m_SL(4, i) - m_lambda(4, i);
	//	m_l(m_E(1, i)) += m_SL(5, i) - m_lambda(5, i);
	//	ln(m_E(1, i)) += 3;

	//	m_l(m_E(2, i)) += m_SL(6, i) - m_lambda(6, i);
	//	m_l(m_E(2, i)) += m_SL(7, i) - m_lambda(7, i);
	//	m_l(m_E(2, i)) += m_SL(8, i) - m_lambda(8, i);
	//	ln(m_E(2, i)) += 3;
	//}


	// ===================================
	Eigen::MatrixXd temp_matrix = m_SL - m_lambda;
	for (int i = 0; i < m_E.cols(); i++)
	{
		m_l(m_E(0, i)) += temp_matrix(0, i) + temp_matrix(1, i) + temp_matrix(2, i);
		ln(m_E(0, i)) += 3;

		m_l(m_E(1, i)) += temp_matrix(3, i) + temp_matrix(4, i) + temp_matrix(5, i);
		ln(m_E(1, i)) += 3;

		m_l(m_E(2, i)) += temp_matrix(6, i) + temp_matrix(7, i) + temp_matrix(8, i);
		ln(m_E(2, i)) += 3;
	}
	// ====================================

	//m_l = m_l.cwiseQuotient(ln);
	m_l = m_solver.solve(m_l);
}

void TriInequality::updateLambda()
{
	//for (int i = 0; i < m_E.cols(); i++)
	//{
	//	m_lambda(0, i) += m_l(m_E(0, i)) - m_SL(0, i);
	//	m_lambda(1, i) += m_l(m_E(0, i)) - m_SL(1, i);
	//	m_lambda(2, i) += m_l(m_E(0, i)) - m_SL(2, i);

	//	m_lambda(3, i) += m_l(m_E(1, i)) - m_SL(3, i);
	//	m_lambda(4, i) += m_l(m_E(1, i)) - m_SL(4, i);
	//	m_lambda(5, i) += m_l(m_E(1, i)) - m_SL(5, i);

	//	m_lambda(6, i) += m_l(m_E(2, i)) - m_SL(6, i);
	//	m_lambda(7, i) += m_l(m_E(2, i)) - m_SL(7, i);
	//	m_lambda(8, i) += m_l(m_E(2, i)) - m_SL(8, i);
	//}

	m_lambda = -m_SL;
	for (int i = 0; i < m_E.cols(); i++)
	{
		const double& m_l_0 = m_l(m_E(0, i));
		const double& m_l_1 = m_l(m_E(1, i));
		const double& m_l_2 = m_l(m_E(2, i));

		m_lambda(0, i) += m_l_0;
		m_lambda(1, i) += m_l_0;
		m_lambda(2, i) += m_l_0;

		m_lambda(3, i) += m_l_1;
		m_lambda(4, i) += m_l_1;
		m_lambda(5, i) += m_l_1;

		m_lambda(6, i) += m_l_2;
		m_lambda(7, i) += m_l_2;
		m_lambda(8, i) += m_l_2;
	}
}