#include "Optimization.h"

Optimization::Optimization(SegMesh& seg_mesh)
	:seg_mesh_(seg_mesh)
{
	Mesh& mesh = seg_mesh.GetMesh();

	dev_info_.ResetDiffBound(mesh.n_edges());

	SetOriMesh(mesh);
}

Optimization::~Optimization()
{
}

void Optimization::SetOriMesh(const Mesh& mesh)
{
	ori_mesh_ = mesh;
}

void Optimization::Deformation()
{
	// Developability-encouraged deformation
	std::cout << "======== Deformation ========" << std::endl;
	Mesh& mesh = seg_mesh_.GetMesh();

	auto& dev_energy = dev_info_.dev_energy_;
	int& cout = dev_info_.cout_;

	MeshDevelop deve;
	deve.oriMesh(ori_mesh_);
	deve.tarMesh(mesh, target_dist_);

	cout = 0;
	do
	{
		std::cout << cout << " th developability deformation:" << std::endl;
		double d_energy = deve.develop(sub_iter_);
		dev_info_.UpdateCurEnergy(d_energy);

		deve.outMesh(mesh);
		deve.GetDevelopEnergy(dev_energy);

		seg_mesh_.WriteMesh(FILE_PATH + "deform_" + std::to_string(cout) + ".obj");

		if (dev_info_.IsSmallDiff())
		{
			std::cout << "------ Small Diff ------" << std::endl;
			break;
		}

		cout++;

		if (cout >= max_iter_)
		{
			std::cout << "------ Max Iter ------" << std::endl;
			break;
		}

	} while (true);

	mesh.update_face_normals();

	std::cout << "------ Output Dev Info ------ " << std::endl;
	dev_info_.UpdateCurvatureEnergy(seg_mesh_);
	dev_info_.UpdateSegEnergy(seg_mesh_);
	dev_info_.OutputInfo();

	std::cout << "------ Output Mesh ------ " << std::endl;
	seg_mesh_.WriteMesh(FILE_PATH + "first_deform.obj");
}

void Optimization::Segmentation()
{
	Mesh& mesh = seg_mesh_.GetMesh();

	// Load Deformation Energy
	auto& dev_energy = dev_info_.dev_energy_;
	dev_energy.resize(mesh.n_edges());

	std::ifstream file_deform(FILE_PATH + "deformation_energy.txt");
	for (size_t i = 0; i < dev_energy.size(); i++)
	{
		file_deform >> dev_energy[i];
	}
	file_deform.close();


	// Step-by-step segmentation
	std::cout << "======== Segmentation ========" << std::endl;
	Segment seg(seg_mesh_, Segment::SegMode::DIFF);

	std::vector<int> false_seg_idx(0);
	do
	{
		std::cout << "------ Get Large ------" << std::endl;

		std::vector<bool> f_status;
		int max_idx;
		if (!GetLargeSegFaces(false_seg_idx, f_status, max_idx)) {
			seg_status_ = false;
			false_seg_idx.clear();
			//continue;
			break;
		}

		std::cout << "------ Segment ------" << std::endl;
		seg.Init(f_status);
		seg_status_ = seg.Run();

		if (!seg_status_)
			break;

	} while (true);

	// Merging
	std::cout << "------ Merge ------" << std::endl;
	Merge merge(seg_mesh_);
	for (int i = 0; i < iter_step_; i++)
	{
		merge.RemoveOneTriangle(M_PI);
		merge.SegMerge(epsilon_merge / iter_step_ * (i + 1), delta_len);
		merge.RemoveOneTriangle(M_PI);
	}

	RemoveSmallSeg();

	std::cout << "------ Output Dev Info ------ " << std::endl;
	seg_mesh_.BoundToIdx();
	dev_info_.UpdateCurvatureEnergy(seg_mesh_);
	dev_info_.UpdateSegEnergy(seg_mesh_);
	dev_info_.OutputInfo();

	std::cout << "------ Output Seg ------ " << std::endl;
	seg_mesh_.WriteSeg(FILE_PATH + "first_deform_seg.txt");
}

void Optimization::Run()
{
	Deformation();

	Segmentation();

	Refinement();
}

void Optimization::Refinement()
{
	Smoothing();

	seg_mesh_.WriteMesh(FILE_PATH + "smooth.obj");
	seg_mesh_.WriteSeg(FILE_PATH + "smooth_seg.txt");

	//Merge merge(seg_mesh_);

	//merge.RemoveOneTriangle();
	//merge.SegMerge(epsilon_merge, delta_len);
	//merge.RemoveOneTriangle();

	//seg_mesh_.WriteMesh(FILE_PATH + "start.obj");
	//seg_mesh_.WriteSeg(FILE_PATH + "start_seg.txt");

	int cout = 0;
	GaussMin cur_min(seg_mesh_);
	cur_min.SetOriMesh(ori_mesh_);

	do
	{
		dev_info_.UpdateCurvatureEnergy(seg_mesh_);
		double prev_cur = dev_info_.max_cur_;

		cur_min.Run(prev_cur);

		seg_mesh_.WriteMesh(FILE_PATH + "min_" + std::to_string(cout) + ".obj");
		seg_mesh_.WriteSeg(FILE_PATH + "min_seg_" + std::to_string(cout) + ".txt");

		dev_info_.UpdateCurvatureEnergy(seg_mesh_);
		double max_cur = dev_info_.max_cur_;
		std::cout << "refinement: " << cout << " max_cur: " << max_cur << std::endl;

		if (max_cur < cur_bound_)
		{
			seg_mesh_.WriteMesh(FILE_PATH + "final.obj");
			seg_mesh_.WriteSeg(FILE_PATH + "final_seg.txt");

			return;
		}

		if (max_cur > prev_cur || max_cur < prev_cur * (1 - .1))
		{
			std::cout << "large change!" << std::endl;
		}
		else
		{
			std::cout << "not high!" << std::endl;
			cur_min.UpdateWeight();
		}
		cout++;

	} while (cout < 100);

}

int Optimization::GetLargeSegFaces(
	const std::vector<int>& false_seg_idx,
	std::vector<bool>& f_status,
	int& max_idx)
{
	// Seg idx
	std::vector<int>& seg_id = seg_mesh_.GetSegId();
	Mesh& mesh = seg_mesh_.GetMesh();

	// Def info
	dev_info_.UpdateCurvatureEnergy(seg_mesh_, true);
	auto high_cur = dev_info_.high_cur_;

	dev_info_.UpdateSegEnergy(seg_mesh_);
	auto seg_energy = dev_info_.seg_energy_;

	for (int i:false_seg_idx)
	{
		seg_energy[i] = 0;
	}

	double& max_energy = dev_info_.max_energy_;

	// Seg with largest energy
	max_idx = -1;
	do
	{
		auto max_iter = std::max_element(seg_energy.begin(), seg_energy.end());

		if ((*max_iter) < min_seg_energy_) break;

		int temp_idx = std::distance(seg_energy.begin(), max_iter);
		if (high_cur[temp_idx] > min_cur_bound_)
		{
			//std::cout << "max energy: " << *max_iter << std::endl;
			max_energy = *max_iter;
			max_idx = temp_idx;
			break;
		}

		*max_iter = 0;

	} while (true);

	if (max_idx == -1)
	{
		return false;
	}
	else
	{
		f_status.assign(mesh.n_faces(), false);
		for (const FH& f_h : mesh.faces())
		{
			if (seg_id[f_h.idx()] == max_idx)
			{
				f_status[f_h.idx()] = true;
			}
		}

		return true;
	}
}

void Optimization::Smoothing()
{
	Mesh& mesh = seg_mesh_.GetMesh();
	std::vector<double> abs_v_gauss = seg_mesh_.GetAbsGauss(true);

	std::vector<std::pair<double, int>> gauss_idx;
	for (size_t i = 0; i < abs_v_gauss.size(); i++)
	{
		gauss_idx.push_back({ abs_v_gauss[i], i });
	}

	std::sort(gauss_idx.begin(), gauss_idx.end(), std::greater<std::pair<double, int>>());

	for (size_t i = 0; i < gauss_idx.size(); i++)
	{
		VH v_h = VH(i);

		if (mesh.is_boundary(v_h)) continue;

		OpenMesh::Vec3d temp_direction(0, 0, 0);
		int cont = 0;
		for (const VH& adj_v : mesh.vv_range(v_h))
		{
			temp_direction += mesh.point(adj_v) - mesh.point(v_h);
			cont++;
		}
		temp_direction /= cont;

		mesh.set_point(v_h, mesh.point(v_h) + temp_direction);
	}
}

void Optimization::RemoveSmallSeg()
{
	const Mesh& mesh = seg_mesh_.GetMesh();
	std::vector<bool>& seg_status = seg_mesh_.GetSeam();
	std::vector<int>& seg_id = seg_mesh_.GetSegId();
	int& seg_num = seg_mesh_.GetSegNum();

	std::vector<double> v_gauss_abs = seg_mesh_.GetAbsGauss();

	double max_cur = *std::max_element(v_gauss_abs.begin(), v_gauss_abs.end()) * mesh.n_faces();

	do
	{
		std::vector<int> f_num;
		seg_mesh_.ComputeSegFaceNum(f_num);

		int small_max_id = -1;
		int small_max_num = -1;
		for (size_t i = 0; i < f_num.size(); i++)
		{
			if (f_num[i] > 0.005 * mesh.n_faces()) continue;

			if (f_num[i] > small_max_num)
			{
				small_max_id = i;
				small_max_num = f_num[i];
			}
		}

		if (small_max_id == -1) break;

		double total_gauss, gauss_0, gauss_1;
		std::vector<double> adj_seg_cur(seg_num, max_cur * 2);
		for (const FH& f_h : mesh.faces())
		{
			if (seg_id[f_h.idx()] != small_max_id) continue;

			for (const HEH& adj_he : mesh.fh_range(f_h))
			{
				EH e_h = mesh.edge_handle(adj_he);

				if (!seg_status[e_h.idx()] || mesh.is_boundary(e_h)) continue;

				int seg_1 = seg_id[mesh.opposite_face_handle(adj_he).idx()];

				gauss_0 = v_gauss_abs[mesh.to_vertex_handle(adj_he).idx()];
				gauss_1 = v_gauss_abs[mesh.from_vertex_handle(adj_he).idx()];

				total_gauss = gauss_0 + gauss_1;

				adj_seg_cur[seg_1] = adj_seg_cur[seg_1] > max_cur + DBL_EPSILON ?
					total_gauss : total_gauss + adj_seg_cur[seg_1];
			}
		}

		int merge_id = std::min_element(adj_seg_cur.begin(), adj_seg_cur.end()) - adj_seg_cur.begin();

		for (const FH& f_h : mesh.faces())
		{
			if (seg_id[f_h.idx()] != small_max_id) continue;

			seg_id[f_h.idx()] = merge_id;
		}

		seg_mesh_.IdxToBound();

	} while (true);
}
