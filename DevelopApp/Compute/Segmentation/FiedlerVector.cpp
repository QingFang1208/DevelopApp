#include "FiedlerVector.h"

FiedlerVector::FiedlerVector(Mesh& mesh, std::vector<int>& glo2loc,
	std::vector<int>& loc2glo, std::vector<double>& sub_fiedler)
	:mesh_(mesh), glo2loc_(glo2loc), loc2glo_(loc2glo), sub_fiedler_(sub_fiedler)
{
}

void FiedlerVector::Run()
{
	mesh_.update_face_normals();
	ComputeGlobalDistance();
	ComputeLocalMeanDistance();
	ConstructAffinity();
	ComputeEigenFuntion();
}

void FiedlerVector::ComputeGlobalDistance()
{
	angle_dist_.clear();
	angle_dist_.resize(mesh_.n_edges(), 0);

	for (const EH& e_h : mesh_.edges())
	{
		HEH he_h_0 = mesh_.halfedge_handle(e_h, 0);
		HEH he_h_1 = mesh_.halfedge_handle(e_h, 1);

		FH f_0 = mesh_.face_handle(he_h_0);
		FH f_1 = mesh_.face_handle(he_h_1);

		if (!mesh_.is_valid_handle(f_0) || !mesh_.is_valid_handle(f_1)) continue;

		OpenMesh::Vec3d normal_0 = mesh_.normal(f_0);
		OpenMesh::Vec3d normal_1 = mesh_.normal(f_1);

		double diff_normal = dot(normal_0 - normal_1, normal_0 - normal_1);
		//if (diff_normal > 0.05)
		//{
		//	diff_normal *= 10;
		//}

		angle_dist_[e_h.idx()] = diff_normal;
	}
}

void FiedlerVector::ComputeLocalMeanDistance()
{
	double sum_dist = 0;
	int sum_num = 0;
	for (const EH& e_h : mesh_.edges())
	{
		HEH he_h_0 = mesh_.halfedge_handle(e_h, 0);
		HEH he_h_1 = mesh_.halfedge_handle(e_h, 1);

		FH f_0 = mesh_.face_handle(he_h_0);
		FH f_1 = mesh_.face_handle(he_h_1);

		if (!mesh_.is_valid_handle(f_0) || !mesh_.is_valid_handle(f_1)) continue;

		int loc_0_idx = glo2loc_[f_0.idx()];
		int loc_1_idx = glo2loc_[f_1.idx()];

		if (loc_0_idx != -1 && loc_1_idx != -1)
		{
			sum_dist += angle_dist_[e_h.idx()];
			sum_num++;
		}
	}

	loc_dis_mean_ = sum_dist / sum_num;
}

void FiedlerVector::ConstructAffinity()
{
	// matrix
	std::vector<T> coeff_triplet;
	for (const EH& e_h : mesh_.edges())
	{
		HEH he_h_0 = mesh_.halfedge_handle(e_h, 0);
		HEH he_h_1 = mesh_.halfedge_handle(e_h, 1);

		FH f_0 = mesh_.face_handle(he_h_0);
		FH f_1 = mesh_.face_handle(he_h_1);

		if (!mesh_.is_valid_handle(f_0) || !mesh_.is_valid_handle(f_1)) continue;

		int loc_0_idx = glo2loc_[f_0.idx()];
		int loc_1_idx = glo2loc_[f_1.idx()];

		double e_len = mesh_.calc_edge_length(e_h);

		if (loc_0_idx != -1 && loc_1_idx != -1)
		{
			double factor = e_len * exp(-angle_dist_[e_h.idx()] / loc_dis_mean_);

			if (factor > epsilon_1_ && factor < epsilon_2_) {
				coeff_triplet.emplace_back(T(loc_0_idx, loc_1_idx, -factor));
				coeff_triplet.emplace_back(T(loc_1_idx, loc_0_idx, -factor));
				coeff_triplet.emplace_back(T(loc_0_idx, loc_0_idx, factor));
				coeff_triplet.emplace_back(T(loc_1_idx, loc_1_idx, factor));
			}
			else if (factor < epsilon_1_) {
				coeff_triplet.emplace_back(T(loc_0_idx, loc_1_idx, -epsilon_1_));
				coeff_triplet.emplace_back(T(loc_1_idx, loc_0_idx, -epsilon_1_));
				coeff_triplet.emplace_back(T(loc_0_idx, loc_0_idx, epsilon_1_));
				coeff_triplet.emplace_back(T(loc_1_idx, loc_1_idx, epsilon_1_));
			}
			else {
				coeff_triplet.emplace_back(T(loc_0_idx, loc_1_idx, -epsilon_2_));
				coeff_triplet.emplace_back(T(loc_1_idx, loc_0_idx, -epsilon_2_));
				coeff_triplet.emplace_back(T(loc_0_idx, loc_0_idx, epsilon_2_));
				coeff_triplet.emplace_back(T(loc_1_idx, loc_1_idx, epsilon_2_));
			}
		}
	}
	aff_matrix_.resize(loc2glo_.size(), loc2glo_.size());
	aff_matrix_.setFromTriplets(coeff_triplet.begin(), coeff_triplet.end());
}

void FiedlerVector::ComputeEigenFuntion()
{
	getMatEngine().connect("");
	eigen2matlab("L", aff_matrix_);
	matlabEval("[V,D] = eigs(L, 2 ,'smallestabs','Tolerance',1e-7);");

	Eigen::VectorXd	fiedler_vector;
	matlabEval("V2 = V(:, 2);");

	// [a, b] -> [min, max]
	matlabEval("V2 = " 
		+ std::to_string(FiedlerPara::max_val_ - FiedlerPara::min_val_)
		+ " * (V2 - min(V2)) / (max(V2) - min(V2)) + "
		+ std::to_string(FiedlerPara::min_val_) + ";");
	matlab2eigen("V2", fiedler_vector);
	matlabEval("[B,I] = sort(V2);");
	//matlabEval("figure(1); plot(1:size(B, 1), B);");
	//matlabEval("figure(2); plot(1:(size(B, 1) - 1), diff(B));");

	sub_fiedler_.clear();
	sub_fiedler_.resize(mesh_.n_faces(), 0.1 * FiedlerPara::min_val_);
	for (size_t i = 0; i < loc2glo_.size(); i++)
	{
		sub_fiedler_[loc2glo_[i]] = fiedler_vector(i);
	}
}
