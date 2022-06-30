#include "Segment.h"

Segment::Segment(SegMesh& seg_mesh, SegMode seg_mode)
	:seg_mesh_(seg_mesh)
{
	Mesh& mesh = seg_mesh.GetMesh();

	fiedler_ = new FiedlerVector(mesh, glo2loc_, loc2glo_, sub_fiedler_);

	if (seg_mode == DIFF)
	{
		seg_cut_ = new DiffSegCut(seg_mesh_, glo2loc_, sub_fiedler_);
	}
	if (seg_mode == CUR)
	{
		seg_cut_ = new CurSegCut(seg_mesh_, glo2loc_, sub_fiedler_);
	}

	glo2loc_.resize(mesh.n_faces());
	loc2glo_.resize(mesh.n_faces());
	for (size_t i = 0; i < mesh.n_faces(); i++)
	{
		glo2loc_[i] = i;
		loc2glo_[i] = i;
	}
}

Segment::~Segment()
{
	if (fiedler_)
	{
		delete fiedler_;
		fiedler_ = NULL;
	}

	if (seg_cut_)
	{
		delete seg_cut_;
		seg_cut_ = NULL;
	}
}

bool Segment::Run()
{
	std::cout << "------ Fiedler ------" << std::endl;
	fiedler_->Run();

	std::cout << "------ Seg Cut ------" << std::endl;
	return seg_cut_->Run();
}

void Segment::Init(const std::vector<bool>& f_status)
{
	std::cout << "------ Glo2Loc Init ------" << std::endl;

	Mesh& mesh = seg_mesh_.GetMesh();

	// Find the current patch
	// Local to global
	int count = 0;
	glo2loc_ = std::vector<int>(mesh.n_faces(), -1);
	loc2glo_.resize(0);
	for (size_t i = 0; i < mesh.n_faces(); i++)
	{
		if (f_status[i])
		{
			glo2loc_[i] = count;
			loc2glo_.push_back(i);
			count++;
		}
	}
}

//void Segment::AddNewSeam(std::vector<int>& seam_e)
//{
//	seam_e = seg_bound_;
//}
//
//void Segment::AddNewSeam(std::vector<bool>& seam_status)
//{
//	if (seg_bound_.size() < 20)
//	{
//		seam_status.assign(mesh_.n_edges(), false);
//		return;
//	}
//
//	for (int i : seg_bound_)
//	{
//		seam_status[i] = true;
//	}
//}
//
//void Segment::GetSegBound(std::vector<int>& seg_bound)
//{
//	seg_bound = seg_bound_;
//}
//
void Segment::GetFiedler(std::vector<double>& sub_fiedler)
{
	sub_fiedler = sub_fiedler_;
}