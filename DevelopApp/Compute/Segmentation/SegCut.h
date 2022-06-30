#pragma once
#include <iostream>
#include "../SegMesh.h"

class SegCut
{
public:
	SegCut(SegMesh& seg_mesh,
		std::vector<int>& glo2loc,
		std::vector<double>& sub_fiedler);

	//bool DiffRun();

	//bool CurRun();

	bool Run();

protected:
	virtual void ComputeMinCutDiff() = 0;

	virtual bool ComputeCutBound() = 0;

	void SegBound();

	void ComputeDiff();

protected:
	SegMesh& seg_mesh_;
	std::vector<int>& glo2loc_;
	std::vector<double>& sub_fiedler_;

	std::vector<double> sort_fiedler_;
	std::vector<std::pair<double, double>> diff_bound_;
	double cut_bound_;

	double min_cut_diff_ = 0.0001;
	const double size_precent_ = 0.01;
};

