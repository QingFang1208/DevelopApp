#pragma once

#include "SegCut.h"
#include "FiedlerVector.h"

class DiffSegCut : public SegCut
{
public:
	DiffSegCut(SegMesh& seg_mesh,
		std::vector<int>& glo2loc,
		std::vector<double>& sub_fiedler);

private:
	void ComputeMinCutDiff();

	bool ComputeCutBound();

private:
	// =========== 2021-9-27 ============== //
	double max_rate_ = 2;
	// =========== 2021-9-27 ============== //
};

