#include "Merge.h"

Merge::Merge(SegMesh& seg_mesh)
	:seg_mesh_(seg_mesh)
{
}

void Merge::RemoveOneTriangle(const double& merge_angle_bound)
{
	//  ___/\___ ---> ________ 

	const Mesh& mesh = seg_mesh_.GetMesh();
	int& seg_num = seg_mesh_.GetSegNum();
	std::vector<int>& seg_id = seg_mesh_.GetSegId();

	std::vector<double> abs_v_gauss = seg_mesh_.GetAbsGauss();

	bool ok;
	int i;
	HEH he_0, he_1, he_2;
	FH f_0, f_1, f_2;
	VH to_0, to_1, to_2;
	for (ok = false, i = 0; !ok && i < 100; ++i) {
		ok = true;

		for (const FH& f_h : mesh.faces())
		{
			he_0 = mesh.halfedge_handle(f_h);
			he_1 = mesh.next_halfedge_handle(he_0);
			he_2 = mesh.prev_halfedge_handle(he_0);

			f_0 = mesh.opposite_face_handle(he_0);
			f_1 = mesh.opposite_face_handle(he_1);
			f_2 = mesh.opposite_face_handle(he_2);

			to_0 = mesh.to_vertex_handle(he_0);
			to_1 = mesh.to_vertex_handle(he_1);
			to_2 = mesh.to_vertex_handle(he_2);

			if (f_0.is_valid() && f_1.is_valid()
				&& seg_id[f_0.idx()] == seg_id[f_1.idx()] 
				&& abs_v_gauss[to_0.idx()] < merge_angle_bound)
			{
				seg_id[f_h.idx()] = seg_id[f_0.idx()];
				ok = false;
				continue;
			}

			if (f_1.is_valid() && f_2.is_valid()
				&& seg_id[f_1.idx()] == seg_id[f_2.idx()]
				&& abs_v_gauss[to_1.idx()] < merge_angle_bound)
			{
				seg_id[f_h.idx()] = seg_id[f_1.idx()];
				ok = false;
				continue;
			}

			if (f_2.is_valid() && f_0.is_valid()
				&& seg_id[f_2.idx()] == seg_id[f_0.idx()]
				&& abs_v_gauss[to_2.idx()] < merge_angle_bound)
			{
				seg_id[f_h.idx()] = seg_id[f_2.idx()];
				ok = false;
				continue;
			}
		}
	}

	seg_mesh_.IdxToBound();
}

void Merge::SegMerge(const double& merge_angle_bound, const double& len_rate)
{
	const Mesh& mesh = seg_mesh_.GetMesh();
	int& seg_num = seg_mesh_.GetSegNum();
	std::vector<int>& seg_id = seg_mesh_.GetSegId();
	std::vector<bool>& seam_status = seg_mesh_.GetSeam();

	std::vector<int> f_num;
	std::vector<int> v_cout;
	std::vector<std::vector<std::vector<HEH>>> start_he_array;
	std::vector<std::vector<std::vector<HEH>>> end_he_array;
	std::vector<std::vector<int>> seg_bound_num;
	std::vector<std::vector<bool>> seg_feature_status;
	do
	{
		// F_num of each seg
		seg_mesh_.ComputeSegFaceNum(f_num);

		// Seam and v_cout
		seg_mesh_.ComputeVertexCount(v_cout);

		// Bound start and end
		ComputeStartandEndHalfedge(v_cout, start_he_array, end_he_array);

		// Seg_bound_num and seg_feature status
		ComputeBoundNumAndFeatureStatus(merge_angle_bound, v_cout, seg_bound_num, seg_feature_status);

		int large_seg_idx = -1;
		do
		{
			int min_idx = std::distance(f_num.begin(), std::min_element(f_num.begin(), f_num.end()));

			if (f_num[min_idx] == mesh.n_faces()) break;
			f_num[min_idx] = mesh.n_faces();

			ComputeLargeIdx(len_rate, min_idx, seg_bound_num, seg_feature_status,
				start_he_array, end_he_array, large_seg_idx);

			if (large_seg_idx != -1)
			{
				for (const EH& e_h : mesh.edges())
				{
					if (mesh.is_boundary(e_h)) continue;

					if (!seam_status[e_h.idx()]) continue;

					FH f_0 = mesh.face_handle(mesh.halfedge_handle(e_h, 0));
					FH f_1 = mesh.face_handle(mesh.halfedge_handle(e_h, 1));

					if ((seg_id[f_0.idx()] == min_idx && seg_id[f_1.idx()] == large_seg_idx)
						|| (seg_id[f_1.idx()] == min_idx && seg_id[f_0.idx()] == large_seg_idx))
					{
						seam_status[e_h.idx()] = false;
					}
				}

				seg_mesh_.BoundToIdx();

				break;
			}

		} while (true);


		// Update seg idx
		if (large_seg_idx == -1)
		{
			break;
		}
	} while (true);
}

void Merge::ComputeStartandEndHalfedge(const std::vector<int>& v_cout, 
	std::vector<std::vector<std::vector<HEH>>>& start_he_array,
	std::vector<std::vector<std::vector<HEH>>>& end_he_array)
{
	const Mesh& mesh = seg_mesh_.GetMesh();
	int& seg_num = seg_mesh_.GetSegNum();
	std::vector<int>& seg_id = seg_mesh_.GetSegId();
	std::vector<bool>& seam_status = seg_mesh_.GetSeam();

	start_he_array.assign(seg_num, std::vector<std::vector<HEH>>(seg_num));
	end_he_array.assign(seg_num, std::vector<std::vector<HEH>>(seg_num));

	for (const HEH& he_h : mesh.halfedges())
	{
		if (mesh.is_boundary(mesh.edge_handle(he_h))) continue;

		if (!seam_status[mesh.edge_handle(he_h).idx()]) continue;

		int seg_0_idx = seg_id[mesh.face_handle(he_h).idx()];
		int seg_1_idx = seg_id[mesh.opposite_face_handle(he_h).idx()];

		if (mesh.is_boundary(mesh.to_vertex_handle(he_h)) 
			|| v_cout[mesh.to_vertex_handle(he_h).idx()] != 2)
		{
			end_he_array[seg_0_idx][seg_1_idx].push_back(he_h);
		}

		if (mesh.is_boundary(mesh.from_vertex_handle(he_h))
			|| v_cout[mesh.from_vertex_handle(he_h).idx()] != 2)
		{
			start_he_array[seg_0_idx][seg_1_idx].push_back(he_h);
		}
	}
}

void Merge::ComputeBoundNumAndFeatureStatus(
	const double& merge_angle_bound,
	const std::vector<int>& v_cout, 
	std::vector<std::vector<int>>& seg_bound_num, 
	std::vector<std::vector<bool>>& seg_feature_status)
{
	const Mesh& mesh = seg_mesh_.GetMesh();
	int& seg_num = seg_mesh_.GetSegNum();
	std::vector<int>& seg_id = seg_mesh_.GetSegId();
	std::vector<bool>& seam_status = seg_mesh_.GetSeam();

	std::vector<double> abs_v_gauss = seg_mesh_.GetAbsGauss();

	seg_bound_num.assign(seg_num, std::vector<int>(seg_num, 0));
	seg_feature_status.assign(seg_num, std::vector<bool>(seg_num, true));
	for (const EH& e_h : mesh.edges())
	{
		if (mesh.is_boundary(e_h)) continue;
		
		if (!seam_status[e_h.idx()]) continue;

		HEH he_0 = mesh.halfedge_handle(e_h, 0);
		HEH he_1 = mesh.halfedge_handle(e_h, 1);

		FH f_0 = mesh.face_handle(he_0);
		FH f_1 = mesh.face_handle(he_1);

		int min_idx = std::min(seg_id[f_0.idx()], seg_id[f_1.idx()]);
		int max_idx = std::max(seg_id[f_0.idx()], seg_id[f_1.idx()]);

		// Seg_feature_status
		bool all_small_cur = true;
		if (v_cout[mesh.to_vertex_handle(he_0).idx()] == 2
			&& abs_v_gauss[mesh.to_vertex_handle(he_0).idx()] > merge_angle_bound)
		{
			all_small_cur = false;
		}
		if (v_cout[mesh.to_vertex_handle(he_1).idx()] == 2
			&& abs_v_gauss[mesh.to_vertex_handle(he_1).idx()] > merge_angle_bound)
		{
			all_small_cur = false;
		}
		seg_feature_status[min_idx][max_idx] = !((!seg_feature_status[min_idx][max_idx]) || (!all_small_cur));
		seg_feature_status[max_idx][min_idx] = !((!seg_feature_status[max_idx][min_idx]) || (!all_small_cur));

		// Seg_bound_num
		seg_bound_num[min_idx][max_idx]++;
		seg_bound_num[max_idx][min_idx]++;
	}
}

void Merge::ComputeLargeIdx(
	const double& len_rate,
	const int& min_idx,
	const std::vector<std::vector<int>>& seg_bound_num, 
	const std::vector<std::vector<bool>>& seg_feature_status,
	const std::vector<std::vector<std::vector<HEH>>>& start_he_array,
	const std::vector<std::vector<std::vector<HEH>>>& end_he_array,
	int& large_seg_idx)
{
	const Mesh& mesh = seg_mesh_.GetMesh();
	int& seg_num = seg_mesh_.GetSegNum();
	std::vector<int>& seg_id = seg_mesh_.GetSegId();
	std::vector<bool>& seam_status = seg_mesh_.GetSeam();

	const std::vector<int>& temp_total_num = seg_bound_num[min_idx];
	const int total_num = std::accumulate(temp_total_num.begin(), temp_total_num.end(), 0);
	double min_cor = DBL_MAX, temp_cor;
	for (int i = 0; i < seg_num; i++)
	{
		if (seg_feature_status[min_idx][i] && temp_total_num[i] == total_num)
		{
			large_seg_idx = i;
			return;
		}
		else if (seg_feature_status[min_idx][i] && start_he_array[min_idx][i].size() == 1 && temp_total_num[i] > len_rate * total_num)
		{
			//std::cout << "temp: " << temp_total_num[i] << std::endl;

			HEH start_he = start_he_array[min_idx][i][0];


			//std::cout << mesh.face_handle(start_he).idx() << std::endl;


			HEH end_he = end_he_array[min_idx][i][0];

			double min_start = DBL_MAX;
			for (const HEH& adj_0 : mesh.voh_range(mesh.from_vertex_handle(start_he)))
			{
				if (adj_0 == start_he) continue;

				if (!seam_status[mesh.edge_handle(adj_0).idx()] && mesh.is_boundary(mesh.edge_handle(adj_0))) continue;

				for (const HEH& adj_1 : mesh.voh_range(mesh.from_vertex_handle(start_he)))
				{
					if (adj_1 == start_he) continue;

					if (!seam_status[mesh.edge_handle(adj_1).idx()] && mesh.is_boundary(mesh.edge_handle(adj_1))) continue;

					temp_cor = dot(mesh.calc_edge_vector(adj_0).normalized(), mesh.calc_edge_vector(adj_1).normalized());

					if (temp_cor < min_start)
					{
						min_start = temp_cor;
					}
				}
			}

			double min_end = DBL_MAX;
			for (const HEH& adj_0 : mesh.vih_range(mesh.to_vertex_handle(end_he)))
			{
				if (adj_0 == end_he) continue;

				if (!seam_status[mesh.edge_handle(adj_0).idx()] && mesh.is_boundary(mesh.edge_handle(adj_0))) continue;

				for (const HEH& adj_1 : mesh.vih_range(mesh.to_vertex_handle(end_he)))
				{
					if (adj_1 == end_he) continue;

					if (!seam_status[mesh.edge_handle(adj_1).idx()] && mesh.is_boundary(mesh.edge_handle(adj_1))) continue;

					temp_cor = dot(mesh.calc_edge_vector(adj_0).normalized(), mesh.calc_edge_vector(adj_1).normalized());

					if (temp_cor < min_end)
					{
						min_end = temp_cor;
					}
				}
			}

			//std::cout << "min: " << min_start << "   " << min_end << "   " << min_cor << std::endl;

			if (min_cor > min_end + min_start)
			{
				large_seg_idx = i;
				min_cor = min_end + min_start;
			}
		}
	}
}
