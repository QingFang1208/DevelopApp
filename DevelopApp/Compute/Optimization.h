#pragma once
#include <queue>
#include <fstream>
#include "SegMesh.h"
#include "Development\MeshDevelop.h"
#include "Segmentation\Segment.h"
#include "Refinement\Merge.h"
#include "Refinement\GaussMin.h"

const std::string FILE_PATH = "mesh_data/";

class DeforInfo
{
public:
	void UpdateCurvatureEnergy(SegMesh& seg_mesh, bool remove_neibor = false)
	{
		const Mesh& mesh  = seg_mesh.GetMesh();
		const auto& seam_status = seg_mesh.GetSeam();
		const auto& seg_id = seg_mesh.GetSegId();
		const int& seg_num = seg_mesh.GetSegNum();

		high_cur_.assign(seg_num, 0);

		std::vector<double> abs_v_gauss = seg_mesh.GetAbsGauss(true);

		if (remove_neibor)
		{
			for (EH e_h : mesh.edges())
			{
				if (!mesh.is_boundary(e_h) && !seam_status[e_h.idx()]) continue;

				HEH he_0 = mesh.halfedge_handle(e_h, 0);

				VH to_v = mesh.to_vertex_handle(he_0);
				VH from_v = mesh.from_vertex_handle(he_0);

				for (const VH& adj_v:mesh.vv_range(to_v))
				{
					abs_v_gauss[adj_v.idx()] = 0;
				}

				for (const VH& adj_v : mesh.vv_range(from_v))
				{
					abs_v_gauss[adj_v.idx()] = 0;
				}
			}
		}

		for (EH e_h : mesh.edges())
		{
			if (mesh.is_boundary(e_h)) continue;

			if (seam_status[e_h.idx()]) continue;

			HEH he_0 = mesh.halfedge_handle(e_h, 0);

			int to_idx = mesh.to_vertex_handle(he_0).idx();
			int from_idx = mesh.from_vertex_handle(he_0).idx();

			int temp_seg_id = seg_id[mesh.face_handle(he_0).idx()];

			high_cur_[temp_seg_id] = std::max(high_cur_[temp_seg_id], abs_v_gauss[to_idx]);
			high_cur_[temp_seg_id] = std::max(high_cur_[temp_seg_id], abs_v_gauss[from_idx]);
		}

		max_cur_ = *std::max_element(high_cur_.begin(), high_cur_.end());
	}

	bool IsSmallCurvature()
	{
		return max_cur_ < M_PI / 180 / 2;
	}

	void UpdateSegEnergy(SegMesh& seg_mesh)
	{
		const Mesh& mesh = seg_mesh.GetMesh();
		const auto& seam_status = seg_mesh.GetSeam();
		const auto& seg_id = seg_mesh.GetSegId();
		const int& seg_num = seg_mesh.GetSegNum();

		seg_energy_.assign(seg_num, 0);

		// seg total energy
		for (EH e_h : mesh.edges())
		{
			if (mesh.is_boundary(e_h)) continue;

			if (seam_status[e_h.idx()]) continue;

			HEH he_0 = mesh.halfedge_handle(e_h, 0);
			HEH he_1 = mesh.halfedge_handle(e_h, 1);

			EH prev_0 = mesh.edge_handle(mesh.prev_halfedge_handle(he_0));
			EH next_0 = mesh.edge_handle(mesh.next_halfedge_handle(he_0));
			EH prev_1 = mesh.edge_handle(mesh.prev_halfedge_handle(he_1));
			EH next_1 = mesh.edge_handle(mesh.next_halfedge_handle(he_1));

			if (seam_status[prev_0.idx()] ||
				seam_status[next_0.idx()] ||
				seam_status[prev_1.idx()] ||
				seam_status[next_1.idx()]) continue;

			FH f_0 = mesh.face_handle(he_0);
			FH f_1 = mesh.face_handle(he_1);

			if (seg_id[f_0.idx()] != seg_id[f_1.idx()]) continue;

			seg_energy_[seg_id[f_0.idx()]] = std::max(dev_energy_[e_h.idx()], seg_energy_[seg_id[f_0.idx()]]);
		}
	}

	void ResetDiffBound(int n_e)
	{
		energy_diff_bound_ *= 4 * n_e;
	}

	bool IsSmallDiff()
	{
		return abs(prev_energy_ - cur_energy_) < energy_diff_bound_;
	}

	void UpdateCurEnergy(double cur_energy)
	{
		prev_energy_ = cur_energy_;
		cur_energy_ = cur_energy;
	}

	void OutputInfo()
	{
		std::ofstream file_infor(FILE_PATH + "deform_and_seg_information.txt");

		file_infor << "Iter Num: " << cout_ << std::endl;
		file_infor << "Prev energy: " << prev_energy_ << std::endl;
		file_infor << "Cur energy: " << cur_energy_ << std::endl;
		file_infor << "Energy diff: " << prev_energy_ - cur_energy_ << std::endl;
		//file_infor << "Max Energy: " << max_energy_ << std::endl;

		file_infor << "Max Curvature: " << max_cur_ << std::endl;

		file_infor << "Seg Num: " << high_cur_.size() << std::endl;

		file_infor << "Seg Energy and Curvature: " << std::endl;
		for (size_t i = 0; i < high_cur_.size(); i++)
		{
			file_infor << i << "th energy: " << seg_energy_[i] << "  high curvature: " << high_cur_[i] << std::endl;
		}

		file_infor.close();

		// Output developability energy
		std::ofstream file_deform(FILE_PATH + "deformation_energy.txt");
		for (size_t i = 0; i < dev_energy_.size(); i++)
		{
			file_deform << dev_energy_[i] << std::endl;
		}
		file_deform.close();
	}

public:
	int cout_ = 0;
	double prev_energy_ = DBL_MAX;
	double cur_energy_ = DBL_MAX / 2;
	double max_energy_ = DBL_MAX;
	double max_cur_ = DBL_MAX;

	std::vector<double> dev_energy_;

	std::vector<double> high_cur_;
	std::vector<double> seg_energy_;

	double energy_diff_bound_ = 1e-8;
};

class Optimization
{
public:
	Optimization(SegMesh& seg_mesh);

	~Optimization();

	void SetOriMesh(const Mesh& mesh);

	void Deformation();

	void Segmentation();

	void Run();

	void Refinement();

private:
	int GetLargeSegFaces(
		const std::vector<int>& false_seg_idx,
		std::vector<bool>& f_status,
		int& max_idx);

	void Smoothing();

	void RemoveSmallSeg();


private:
	SegMesh& seg_mesh_;
	Mesh ori_mesh_;

	DeforInfo dev_info_;
	//Mesh mesh_;

	// Deformation
	const int sub_iter_ = 2000;
	const double target_dist_ = 0.005;
	const int max_iter_ = 10;
	//const int max_iter_ = 50; // More iterations if necessary
	const double min_seg_energy_ = 1e-5;
	const double min_cur_bound_ = M_PI / 180;
	
	// Segmentation
	const double delta_len = 0.1;
	const double epsilon_merge = M_PI / 180 * 5;
	const int iter_step_ = 10;

	// Refinement
	const double cur_bound_ = 5 * 1e-4;

	bool seg_status_ = false;
};


