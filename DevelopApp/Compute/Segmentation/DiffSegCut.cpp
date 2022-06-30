#include "DiffSegCut.h"

DiffSegCut::DiffSegCut(SegMesh& seg_mesh, std::vector<int>& glo2loc, std::vector<double>& sub_fiedler)
	:SegCut(seg_mesh, glo2loc, sub_fiedler)
{
}

void DiffSegCut::ComputeMinCutDiff()
{
	int active_f_n = 0;
	for (size_t i = 0; i < glo2loc_.size(); i++)
	{
		if (glo2loc_[i] == -1) continue;
		active_f_n++;
	}

	// Threshold of the Fiedler vector difference 
	min_cut_diff_ = (FiedlerPara::max_val_ - FiedlerPara::min_val_) * max_rate_ / active_f_n;
}

bool DiffSegCut::ComputeCutBound()
{
	ComputeDiff();

	// Determine whether to segment
	while (diff_bound_.size() > 0)
	{
		double temp_diff = diff_bound_.back().first;
		double temp_bound = diff_bound_.back().second;
		diff_bound_.pop_back();

		int max_num = 0, min_num = 0;
		for (size_t i = 0; i < sort_fiedler_.size(); i++)
		{
			if (sort_fiedler_[i] > temp_bound)
			{
				max_num++;
			}
			else
			{
				min_num++;
			}
		}

		if (max_num > 0 * sub_fiedler_.size() && min_num > 0 * sub_fiedler_.size())
		{
			cut_bound_ = temp_bound;
			return true;
		}
	}

	return false;
}