#include "SegMesh.h"

SegMesh::SegMesh(Mesh& mesh)
	:mesh_(mesh)
{
	Init();
}

void SegMesh::ReadSeg(const std::string& file_name)
{
	std::ifstream seg_file(file_name);
	for (size_t i = 0; i < mesh_.n_faces(); i++)
	{
		FH f_h = mesh_.face_handle(i);
		int temp_id;
		seg_file >> temp_id;

		seg_id_[i] = temp_id;
	}
	seg_file.close();

	IdxToBound();
}

void SegMesh::Init()
{
	seg_num_ = 1;

	seam_status_.assign(mesh_.n_edges(), false);

	seg_id_.assign(mesh_.n_faces(), 0);
}

Mesh& SegMesh::GetMesh()
{
	return mesh_;
}

std::vector<bool>& SegMesh::GetSeam()
{
	return seam_status_;
}

void SegMesh::WriteSeg(const std::string& file_name)
{
	std::ofstream seg_cout(file_name);
	for (size_t i = 0; i < mesh_.n_faces(); i++)
	{
		if (i != mesh_.n_faces() - 1)
		{
			seg_cout << seg_id_[i] << '\n';
		}
		else
		{
			seg_cout << seg_id_[i];
		}

	}
	seg_cout.close();
}

void SegMesh::WriteMesh(const std::string& file_name)
{
	MeshTools::WriteMesh(mesh_, file_name);
}

std::vector<double> SegMesh::GetAbsGauss(bool remove_seam)
{
	MeshTools::ComputeGaussianCurvature(mesh_, abs_v_gauss_);

	for (size_t i = 0; i < abs_v_gauss_.size(); i++)
	{
		abs_v_gauss_[i] = abs(abs_v_gauss_[i]);
	}

	if (remove_seam)
	{
		for (const EH& e_h : mesh_.edges())
		{
			if (!seam_status_[e_h.idx()]) continue;

			HEH he_0 = mesh_.halfedge_handle(e_h, 0);

			abs_v_gauss_[mesh_.to_vertex_handle(he_0).idx()] = 0;
			abs_v_gauss_[mesh_.from_vertex_handle(he_0).idx()] = 0;
		}
	}

	modified_gauss_ = abs_v_gauss_;

	return abs_v_gauss_;
}

void SegMesh::SetModifiedGauss(std::vector<double> v_gauss)
{
	modified_gauss_ = v_gauss;
}

std::vector<double>& SegMesh::GetModifiedGauss()
{
	return modified_gauss_;
}

void SegMesh::ComputeVertexCount(std::vector<int>& v_cout)
{
	v_cout.assign(mesh_.n_vertices(), 0);
	HEH he_0, he_1;
	for (const EH& e_h : mesh_.edges())
	{
		if (!seam_status_[e_h.idx()]) continue;

		he_0 = mesh_.halfedge_handle(e_h, 0);
		he_1 = mesh_.halfedge_handle(e_h, 1);

		v_cout[mesh_.to_vertex_handle(he_0).idx()]++;
		v_cout[mesh_.to_vertex_handle(he_1).idx()]++;
	}
}

void SegMesh::ComputeVertexSeg(std::vector<int>& v_seg)
{
	v_seg.assign(mesh_.n_vertices(), -1);
	for (const FH& f_h : mesh_.faces())
	{
		int f_seg = seg_id_[f_h.idx()];
		for (const VH& fv_h : mesh_.fv_range(f_h))
		{
			v_seg[fv_h.idx()] = f_seg;
		}
	}

	for (const EH& e_h : mesh_.edges())
	{
		if (!mesh_.is_boundary(e_h) && !seam_status_[e_h.idx()]) continue;

		VH v_0 = mesh_.to_vertex_handle(mesh_.halfedge_handle(e_h, 0));
		VH v_1 = mesh_.to_vertex_handle(mesh_.halfedge_handle(e_h, 1));

		v_seg[v_0.idx()] = -1;
		v_seg[v_1.idx()] = -1;
	}
}


void SegMesh::BoundToIdx()
{
	seg_id_.assign(mesh_.n_faces(), -1);

	seg_num_ = 0;

	do
	{
		int start_idx = -1;
		for (int i = 0; i < mesh_.n_faces(); i++)
		{
			if (seg_id_[i] != -1) continue;
			
			start_idx = i; break;
		}

		if (start_idx == -1) return;

		std::vector<int> f_stack;
		f_stack.push_back(start_idx);
		do
		{
			int top_idx = f_stack.back(); f_stack.pop_back();
			seg_id_[top_idx] = seg_num_;
			
			for (HEH fh_h : mesh_.fh_range(mesh_.face_handle(top_idx)))
			{
				FH oppo_f_h = mesh_.opposite_face_handle(fh_h);

				if (!mesh_.is_valid_handle(oppo_f_h)) continue;
				if (seg_id_[oppo_f_h.idx()] != -1) continue;

				int adj_e_idx = mesh_.edge_handle(fh_h).idx();

				if (!seam_status_[adj_e_idx])
				{
					f_stack.push_back(oppo_f_h.idx());
				}
			}
		} while (f_stack.size() != 0);

		++seg_num_;

	} while (true);
}

std::vector<int>& SegMesh::GetSegId()
{
	return seg_id_;
}

void SegMesh::IdxToBound()
{
	seam_status_.assign(mesh_.n_edges(), false);

	int idx_0, idx_1;
	for (const EH& e_h : mesh_.edges())
	{
		if (mesh_.is_boundary(e_h)) continue;

		idx_0 = mesh_.face_handle(mesh_.halfedge_handle(e_h, 0)).idx();
		idx_1 = mesh_.face_handle(mesh_.halfedge_handle(e_h, 1)).idx();

		if (seg_id_[idx_0] == seg_id_[idx_1]) continue;
		
		seam_status_[e_h.idx()] = true;
	}

	BoundToIdx();
}

void SegMesh::ComputeSegFaceNum(std::vector<int>& f_num)
{
	f_num.assign(seg_num_, 0);
	for (int temp_id:seg_id_)
	{
		++f_num[temp_id];
	}
}

int& SegMesh::GetSegNum()
{
	return seg_num_;
}

void SegMesh::ComputeVertexStatus(std::vector<bool>& seam_v_status)
{
	seam_v_status.assign(mesh_.n_vertices(), false);
	for (size_t i = 0; i < seam_status_.size(); i++)
	{
		if (!seam_status_[i]) continue;
		
		HEH he_h = mesh_.halfedge_handle(EH(i), 0);

		seam_v_status[mesh_.to_vertex_handle(he_h).idx()] = true;
		seam_v_status[mesh_.from_vertex_handle(he_h).idx()] = true;
	}
}
