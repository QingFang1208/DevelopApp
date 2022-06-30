#include "SegCut.h"

SegCut::SegCut(SegMesh& seg_mesh,
	std::vector<int>& glo2loc,
	std::vector<double>& sub_fiedler)
	:seg_mesh_(seg_mesh), 
	glo2loc_(glo2loc), 
	sub_fiedler_(sub_fiedler)
{
}

//bool SegCut::DiffRun()
//{
//	int active_f_n = 0;
//	for (size_t i = 0; i < glo2loc_.size(); i++)
//	{
//		if (glo2loc_[i] == -1) continue;
//		active_f_n++;
//	}
//
//	min_cut_diff_ = (0.99 - 0.1) * 10 / active_f_n;
//
//	double cut_bound;
//	if (ComputeMaxDiff(cut_bound))
//	{
//		SegBound(cut_bound);
//
//		//std::cout << "Successful segmentation!" << std::endl;
//
//		return true;
//	}
//	else
//	{
//		std::cout << "Small diff or seg size!" << std::endl;
//
//		return false;
//	}
//}
//
//bool SegCut::CurRun()
//{
//	min_cut_diff_ = 0;
//
//
//}

bool SegCut::Run()
{
	// Segment the current patch
	ComputeMinCutDiff();

	if (ComputeCutBound())
	{
		SegBound();

		std::cout << "Successful segmentation!" << std::endl;

		return true;
	}
	else
	{
		std::cout << "Small diff or seg size!" << std::endl;

		return false;
	}

	//int active_f_n = 0;
	//for (size_t i = 0; i < glo2loc_.size(); i++)
	//{
	//	if (glo2loc_[i] == -1) continue;
	//	active_f_n++;
	//}

	//min_cut_diff_ = (0.99 - 0.1) * 10 / active_f_n;

	//double cut_bound;
	//if (ComputeMaxDiff(cut_bound))
	//{
	//	SegBound(cut_bound);

	//	//std::cout << "Successful segmentation!" << std::endl;
	//	
	//	return true;
	//}
	//else
	//{
	//	std::cout << "Small diff or seg size!" << std::endl;

	//	return false;
	//}
}

void SegCut::ComputeDiff()
{
	// Fiedler for active faces
	sort_fiedler_.clear();
	for (size_t i = 0; i < sub_fiedler_.size(); i++)
	{
		if (glo2loc_[i] == -1) continue;

		sort_fiedler_.push_back(sub_fiedler_[i]);
	}
	std::sort(sort_fiedler_.begin(), sort_fiedler_.end(), std::greater<double>());

	// Compute the difference of the sorted Fiedler vector
	diff_bound_.clear();
	for (size_t i = 0; i < sort_fiedler_.size() - 1; i++)
	{
		double temp_diff = sort_fiedler_[i] - sort_fiedler_[i + 1];

		if (temp_diff < min_cut_diff_) continue;

		diff_bound_.push_back({ temp_diff , sort_fiedler_[i + 1] + DBL_EPSILON });
	}

	std::sort(diff_bound_.begin(), diff_bound_.end());
}

//bool SegCut::ComputeMaxDiff(double& cut_bound)
//{
//	ComputeDiff();
//
//	// cur bound
//	while (diff_bound_.size() > 0)
//	{
//		double temp_diff = diff_bound_.back().first;
//		double temp_bound = diff_bound_.back().second;
//		diff_bound_.pop_back();
//
//		int max_num = 0, min_num = 0;
//		for (size_t i = 0; i < sort_fiedler_.size(); i++)
//		{
//			if (sort_fiedler_[i] > temp_bound)
//			{
//				max_num++;
//			}
//			else
//			{
//				min_num++;
//			}
//		}
//
//		if (max_num > size_precent_ * sub_fiedler_.size() && min_num > size_precent_ * sub_fiedler_.size())
//		{
//			cut_bound = temp_bound;
//			return true;
//		}
//	}
//	return false; 
//}

void SegCut::SegBound()
{
	Mesh& mesh = seg_mesh_.GetMesh();
	std::vector<bool>& seam_status = seg_mesh_.GetSeam();

	for (EH e_h : mesh.edges())
	{
		if (mesh.is_boundary(e_h)) continue;

		HEH he_h = mesh.halfedge_handle(e_h, 0);
		FH f_0 = mesh.face_handle(he_h);
		FH f_1 = mesh.opposite_face_handle(he_h);

		if (glo2loc_[f_0.idx()] == -1 || glo2loc_[f_1.idx()] == -1) continue;

		double sign = (sub_fiedler_[f_0.idx()] - cut_bound_)
			* (sub_fiedler_[f_1.idx()] - cut_bound_);

		if (sign > 0) continue;

		seam_status[e_h.idx()] = true;
	}

	seg_mesh_.BoundToIdx();
}
