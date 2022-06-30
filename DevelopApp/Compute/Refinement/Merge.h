#pragma once
#include "../SegMesh.h"

class Merge
{
public:
	Merge(SegMesh& seg_mesh);

	void RemoveOneTriangle(const double& merge_angle_bound = M_PI / 180 * 5);

	void SegMerge(
		const double& merge_angle_bound = M_PI / 180 * 5, 
		const double& len_rate = 0.05);

private:
	void ComputeStartandEndHalfedge(const std::vector<int>& v_cout,
		std::vector<std::vector<std::vector<HEH>>>& start_he_array,
		std::vector<std::vector<std::vector<HEH>>>& end_he_array);

	void ComputeBoundNumAndFeatureStatus(
		const double& merge_angle_bound,
		const std::vector<int>& v_cout,
		std::vector<std::vector<int>>& seg_bound_num,
		std::vector<std::vector<bool>>& seg_feature_status);

	void ComputeLargeIdx(
		const double& len_rate,
		const int& min_idx,
		const std::vector<std::vector<int>>& seg_bound_num,
		const std::vector<std::vector<bool>>& seg_feature_status,
		const std::vector<std::vector<std::vector<HEH>>>& start_he_array,
		const std::vector<std::vector<std::vector<HEH>>>& end_he_array,
		int& large_seg_idx);

private:
	SegMesh& seg_mesh_;
};

