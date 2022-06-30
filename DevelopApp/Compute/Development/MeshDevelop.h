#pragma once
#include "MeshLaplace.h"
#include "AndersonAcceleration.h"
#include "TriInequality.h"
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Vector3d;
using Eigen::Matrix3d;
using std::vector;
using SPMatrixXd = Eigen::SparseMatrix<double, Eigen::ColMajor>;

#define USE_FAST_SVD

class MeshDevelop
{
public:
	MeshDevelop();
	void oriMesh(const Mesh& mesh);
	void tarMesh(const Mesh& mesh, double dist);
	void outMesh(Mesh& mesh);

	double develop(int maxIter = 2e4);
	
	void GetDevelopEnergy(std::vector<double>& energy) {
		energy = dev_energy;
	}


private:
	Mesh tmp;
	MeshTopo m_topo;
	MatrixXd m_vpfnen;
	MatrixXd m_vpfnen0;
	MatrixXd m_eb0, m_vt;
	VectorXd m_el, m_ew, m_ed, m_es, m_eA, m_fA, m_vA;
	MatrixXd m_fDhDeReU; // 9 * size

	std::vector<bool> m_fixV;
	std::vector<Eigen::JacobiSVD<MatrixXd>> m_svd;
	std::vector<Matrix3d> m_svdM, m_svdU, m_svdV;
	std::vector<Vector3d> m_svdA;

	std::vector<Eigen::SelfAdjointEigenSolver<Matrix3d>> m_eig;

	double m_wD = 1.0;
	//double m_wD = 0.1;
	double m_wP = 1e3;
	double m_similar = 0;
	double m_dist;
	double m_scale;

	SPMatrixXd m_lap, m_elap;
	Eigen::SimplicialLDLT<SPMatrixXd> m_ldlt;
	std::unique_ptr<AndersonAcceleration> AA;
	TriInequality triNeq;

	// Develop energy
	std::vector<double> dev_energy;


private:
	void localStep();
	void globalStep();
	void updateN();
	void updateV();
	void updatePara(int iter);

	void calcEN();
	void calcFN();

	void updateD();
	void updateEU();
	void updateER();
	void updateVP();
	void updateVT();
	void updateES();

	double energy();
	double energyD();
	double energyFN();
	double energyDP();
	double energyRP();

	bool inBound();


private:
	std::vector<double> EN0_i_00;
	std::vector<double> EN0_i_01;
	std::vector<double> EN0_i_02;
	std::vector<double> EN0_i_03;
	std::vector<double> EN0_i_10;
	std::vector<double> EN0_i_11;
	std::vector<double> EN0_i_12;
	std::vector<double> EN0_i_13;
	std::vector<double> EN0_i_20;
	std::vector<double> EN0_i_21;
	std::vector<double> EN0_i_22;
	std::vector<double> EN0_i_23;

	std::vector<double> EN_i_00;
	std::vector<double> EN_i_01;
	std::vector<double> EN_i_02;
	std::vector<double> EN_i_03;
	std::vector<double> EN_i_10;
	std::vector<double> EN_i_11;
	std::vector<double> EN_i_12;
	std::vector<double> EN_i_13;
	std::vector<double> EN_i_20;
	std::vector<double> EN_i_21;
	std::vector<double> EN_i_22;
	std::vector<double> EN_i_23;


	std::vector<double> e0fn0_i_0;
	std::vector<double> e0fn0_i_1;
	std::vector<double> e0fn0_i_2;
	std::vector<double> e0fn_i_0;
	std::vector<double> e0fn_i_1;
	std::vector<double> e0fn_i_2;
	std::vector<double> e1fn0_i_0;
	std::vector<double> e1fn0_i_1;
	std::vector<double> e1fn0_i_2;
	std::vector<double> e1fn_i_0;
	std::vector<double> e1fn_i_1;
	std::vector<double> e1fn_i_2;

	std::vector<double> diag_i_0;
	std::vector<double> diag_i_1;
	std::vector<double> diag_i_2;
	std::vector<double> diag_i_3;
};