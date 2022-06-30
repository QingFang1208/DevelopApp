#include "CurSegCut.h"

CurSegCut::CurSegCut(SegMesh& seg_mesh,
	std::vector<int>& glo2loc,
	std::vector<double>& sub_fiedler)
	:SegCut(seg_mesh, glo2loc, sub_fiedler)
{
}

void CurSegCut::ComputeMinCutDiff()
{
	min_cut_diff_ = 0;
}

bool CurSegCut::ComputeCutBound()
{
	Mesh& mesh = seg_mesh_.GetMesh();
	const auto& seam_status = seg_mesh_.GetSeam();
	auto abs_v_gauss = seg_mesh_.GetModifiedGauss();

	ComputeDiff();

	// cur bound

	while (diff_bound_.size() > 0)
	{
		double temp_diff = diff_bound_.back().first;
		double temp_bound = diff_bound_.back().second;
		diff_bound_.pop_back();

		if (IsAcrossHighCur(mesh, seam_status, abs_v_gauss, temp_bound)
			&& IsProperSize(mesh, temp_bound))
		{
			cut_bound_ = temp_bound;
			return true;
		}
	}
	return false;
}

bool CurSegCut::IsAcrossHighCur(const Mesh& mesh,
	const std::vector<bool>& seam_status,
	const std::vector<double>& abs_v_gauss,
	double temp_bound)
{
	for (const EH& e_h : mesh.edges())
	{
		if (mesh.is_boundary(e_h)) continue;

		HEH he_h = mesh.halfedge_handle(e_h, 0);
		FH f_0 = mesh.face_handle(he_h);
		FH f_1 = mesh.opposite_face_handle(he_h);

		if (glo2loc_[f_0.idx()] == -1 || glo2loc_[f_1.idx()] == -1) continue;

		double sign = (sub_fiedler_[f_0.idx()] - temp_bound)
			* (sub_fiedler_[f_1.idx()] - temp_bound);

		if (sign > 0) continue;

		VH to_v = mesh.to_vertex_handle(he_h);
		VH from_v = mesh.from_vertex_handle(he_h);

		if (abs_v_gauss[to_v.idx()] > DBL_EPSILON 
			|| abs_v_gauss[from_v.idx()] > DBL_EPSILON)
		{
			return true;
		}
	}

	return false;
}

bool CurSegCut::IsProperSize(const Mesh& mesh, const double& temp_bound)
{
	int less_cout = 0, larger_cout = 0;
	size_t n_f = mesh.n_faces();
	for (size_t f_id = 0; f_id < n_f; ++f_id)
	{
		if (glo2loc_[f_id] == -1) continue;

		if (sub_fiedler_[f_id] > temp_bound)
		{
			++larger_cout;
		}
		else
		{
			++less_cout;
		}
	}

	if (larger_cout > size_precent_ * glo2loc_.size() 
		&& less_cout > size_precent_ * glo2loc_.size())
	{
		return true;
	}
	else
	{
		return false;
	}
}