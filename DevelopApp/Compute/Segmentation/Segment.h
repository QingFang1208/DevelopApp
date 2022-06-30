#pragma once
#include <iostream>
#include "FiedlerVector.h"

#include "DiffSegCut.h"
#include "CurSegCut.h"

class Segment
{
public:
	enum SegMode { DIFF, CUR };

public:
	Segment(SegMesh& seg_mesh, SegMode seg_mode);

	~Segment();

	void Init(const std::vector<bool>& f_status);

	bool Run();

	//void AddNewSeam(std::vector<int>& seam_e);

	//void AddNewSeam(std::vector<bool>& seam_status);

	//void GetSegBound(std::vector<int>& seg_bound);

	void GetFiedler(std::vector<double>& sub_fiedler);

private:
	SegMesh& seg_mesh_;

	std::vector<int> glo2loc_, loc2glo_;

	FiedlerVector* fiedler_;
	std::vector<double> sub_fiedler_;

	SegCut* seg_cut_;
};

