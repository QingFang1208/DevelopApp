#pragma once
#include "SegCut.h"

class CurSegCut : public SegCut
{
public:
	CurSegCut(SegMesh& seg_mesh,
		std::vector<int>& glo2loc,
		std::vector<double>& sub_fiedler);


private:
	void ComputeMinCutDiff();

	bool ComputeCutBound();

	bool IsAcrossHighCur(const Mesh& mesh, 
		const std::vector<bool>& seam_status,
		const std::vector<double>& abs_v_gauss,
		double temp_bound);

	bool IsProperSize(const Mesh& mesh, const double& temp_bound);
};

