#pragma once
#include <Eigen/Sparse>
#include "../SegMesh.h"

typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Triplet<double> T;

class GaussMin
{
public:
	GaussMin(SegMesh& seg_mesh);

	void Reset();

	void SetOriMesh(Mesh& mesh);

	void UpdateWeight()
	{
		lambda_dev_ *= update_lambda_precent_ * update_lambda_precent_;
		lambda_seam_ *= update_lambda_precent_;
		lambda_smooth_ *= update_lambda_precent_;
	}

	void Run(double max_cur);

private:
	void Init();

	void SetPosition();
	void SetStepH();

	void SetGloToLoc();
	void SetIdentityMatrix();
	void SetInnerSmoothMatrix();
	void SetSeamSmoothMatrix();

	//void SetInitJTJ();
	void SetInitTriplet();

	void Opt();

	void ComputeNormK();
	void ComputeR();
	void UpdateEnergy();
	void ComputeJ();

	bool LinearSearch(double& alpha_k, double c2 = 0.9);

	//void Projection()
	//{
	//	for (int i = 0; i < new_pos_.rows(); ++i)
	//	{
	//		new_pos_.row(i) = ori_pos_.row(i) + (new_pos_.row(i) - ori_pos_.row(i))
	//		/ std::max(1.0, (new_pos_.row(i) - ori_pos_.row(i)).norm() / m_dist_);
	//	}
	//}

	std::vector<T> ComputeTempTriplet(const int& start_idx,
		const int& thread_id);

	template<typename Scalar>
	inline Scalar AngleDefect(const VH& v_h, const Eigen::Matrix<Scalar, -1, 3>& temp_pos);

private:
	SegMesh& seg_mesh_;
	Mesh& mesh_;

	Eigen::MatrixX3d ori_pos_;

	double m_dist_;

	std::vector<double> step_h_;

	//int n_inner_, n_seam_;
	std::vector<int> glo2inner_;
	//std::vector<int> glo2seam_;
	std::vector<int> inner2glo_;
	std::vector<int> seam2glo_;

	//SpMat m_init_JTJ_;
	std::vector<T> init_triplet_;

	Eigen::VectorXd r_;
	SpMat J_;
	SpMat m_JTJ_;

	SpMat m_sm_;
	SpMat m_se_;
	SpMat m_id_;
	SpMat m_de_;

	Eigen::MatrixX3d new_pos_;
	Eigen::VectorXd norm_K_;
	
	//SpMat partial_E_;
	//SpMat m_JTJ_;
	//Eigen::VectorXd v_JTd_;
	Eigen::VectorXd v_direction_;
	Eigen::MatrixX3d m_direction_;
	Eigen::VectorXd v_gradient_;

	double total_energy_;
	double smooth_energy_;
	double seam_energy_;
	double close_energy_;
	double dev_energy_;

	Eigen::SimplicialLDLT<SpMat> m_ldlt_;
	//Eigen::BiCGSTAB<SpMat_Row, Eigen::IncompleteLUT<double>> m_bicg_;

	double cur_bound_ = DBL_MAX;

	double lambda_exp_close_ = 1.0;
	double lambda_exp_smooth_ = 1.0;
	double lambda_exp_seam_ = 1.0;
	Eigen::VectorXd r_exp_close_;
	Eigen::VectorXd r_exp_smooth_;
	Eigen::VectorXd r_exp_seam_;

	Eigen::VectorXd m_sm_sq_;
	Eigen::VectorXd m_se_sq_;

	int iter_num_ = 500;
	int output_num_ = 100;
	int update_lambda_num_ = 100;
	double small_gradient_precent_ = 0.001;
	double target_dist_ = 0.05;//0.05
	double step_precent_ = 0.01;

	double lambda_smooth_ = 2.0;//2.0
	double lambda_seam_ = 2.0;//2.0
	double lambda_close_ = 1.0;//1.0
	double lambda_dev_ = 1.0;//1.0

	double update_lambda_precent_ = (1 + 1e-1);

	double break_diff_;

	int p_ = 1;
};

