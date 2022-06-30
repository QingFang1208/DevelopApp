#pragma once
#include "../../MeshDefinition.h"
#include "matlab_utils.h"

typedef Eigen::SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of double
typedef Eigen::Triplet<double> T;

namespace FiedlerPara
{
	const double min_val_ = 0.1;
	const double max_val_ = 0.99;
}

class FiedlerVector
{
public:
	FiedlerVector(Mesh& mesh, std::vector<int>& glo2loc, std::vector<int>& loc2glo, std::vector<double>& sub_fiedler);

	void Run();

private:
	void ComputeGlobalDistance();

	void ComputeLocalMeanDistance();

	void ConstructAffinity();

	void ComputeEigenFuntion();

private:
	Mesh& mesh_;
	std::vector<int>& glo2loc_;
	std::vector<int>& loc2glo_;
	std::vector<double>& sub_fiedler_;	// fiedler of sub seg

	double epsilon_1_ = 0.0001, epsilon_2_ = 0.1;

	std::vector<double> angle_dist_;
	double loc_dis_mean_;

	SpMat aff_matrix_;
};

