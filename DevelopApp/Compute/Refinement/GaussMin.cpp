#include "GaussMin.h"

GaussMin::GaussMin(SegMesh& seg_mesh)
	:seg_mesh_(seg_mesh), mesh_(seg_mesh.GetMesh())
{
	SetPosition();
	SetStepH();

	//lambda_close_ /= m_dist_;

	lambda_dev_ *= m_dist_;
}

void GaussMin::Reset()
{
	//lambda_dev_ *= 0.9;
	//lambda_smooth_ *= 0.9;
	//lambda_seam_ *= 0.9;
}

void GaussMin::SetOriMesh(Mesh& mesh)
{
	int n_v = mesh.n_vertices();

	// position
	ori_pos_.resize(n_v, 3);
	for (const VH& v_h : mesh.vertices())
	{
		ori_pos_.row(v_h.idx()) = Eigen::Map<Eigen::Vector3d>(mesh.point(v_h).data(), 3);
	}
}

void GaussMin::Run(double max_cur)
{
	//lambda_dev_ = m_dist_ / 0.01;

	if (cur_bound_ == DBL_MAX)
	{
		cur_bound_ = max_cur;
	}

	//lambda_seam_ *= pow(cur_bound_ / max_cur, p_);
	//lambda_smooth_ *= pow(cur_bound_ / max_cur, p_);

	cur_bound_ = max_cur;

	//break_diff_ = 1e-3 * pow(0.001 / cur_bound_, p_);
	break_diff_ = 0;

	for (const VH& v_h : mesh_.vertices())
	{
		new_pos_.row(v_h.idx()) = Eigen::Map<Eigen::Vector3d>(mesh_.point(v_h).data(), 3);
	}

	Init();

	Opt();

	for (size_t i = 0; i < new_pos_.rows(); i++)
	{
		mesh_.set_point(VH(i),
			OpenMesh::Vec3d(new_pos_(i, 0), new_pos_(i, 1), new_pos_(i, 2)));
	}
}

void GaussMin::Init()
{
	std::cout << "Coeff close: " << lambda_close_
		<< " coeff dev: " << lambda_dev_ << std::endl;

	SetGloToLoc();
	
	SetIdentityMatrix();
	SetInnerSmoothMatrix();
	SetSeamSmoothMatrix();

	if (m_sm_.rows() > 0)
	{
		m_sm_sq_ = (m_sm_ * ori_pos_).rowwise().squaredNorm();
	}

	if (m_se_.rows() > 0)
	{
		m_se_sq_ = (m_se_ * ori_pos_).rowwise().squaredNorm();
	}

	SetInitTriplet();
}

void GaussMin::SetPosition()
{
	int n_v = mesh_.n_vertices();

	// position
	ori_pos_.resize(n_v, 3);
	for (const VH& v_h : mesh_.vertices())
	{
		ori_pos_.row(v_h.idx()) = Eigen::Map<Eigen::Vector3d>(mesh_.point(v_h).data(), 3);
	}

	new_pos_ = ori_pos_;

	m_direction_.resize(ori_pos_.rows(), ori_pos_.cols());

	OpenMesh::Vec3d boxMin(DBL_MAX), boxMax(-DBL_MAX);
	for (const VH& v_h : mesh_.vertices())
	{
		boxMin.minimize(mesh_.point(v_h));
		boxMax.maximize(mesh_.point(v_h));
	}
	m_dist_ = (boxMax - boxMin).norm() * target_dist_;
}

void GaussMin::SetStepH()
{
	step_h_.resize(mesh_.n_vertices(), 0);

	double min_len = DBL_MAX;
	for (const EH& e_h : mesh_.edges())
	{
		min_len = std::min(min_len, mesh_.calc_edge_length(e_h));
	}

	for (size_t i = 0; i < mesh_.n_vertices(); i++)
	{
		step_h_[i] = step_precent_ * min_len;
	}
}


void GaussMin::ComputeNormK()
{
#pragma omp parallel for
	for (int i = 0; i < inner2glo_.size(); i++)
	{
		norm_K_(i) = AngleDefect(VH(inner2glo_[i]), new_pos_) / cur_bound_;
	}

}

void GaussMin::ComputeR()
{
	int n_v = mesh_.n_vertices();
	int n_inner = inner2glo_.size();
	int n_seam = seam2glo_.size();
	Eigen::MatrixX3d temp_matrix;

	int start_idx = 0;

	// close
	temp_matrix = lambda_close_ * (new_pos_ - ori_pos_);
	r_.segment(start_idx, n_v) = temp_matrix.col(0); start_idx += n_v;
	r_.segment(start_idx, n_v) = temp_matrix.col(1); start_idx += n_v;
	r_.segment(start_idx, n_v) = temp_matrix.col(2); start_idx += n_v;

	// smooth
	temp_matrix = lambda_smooth_ * (m_sm_ * new_pos_);
	r_.segment(start_idx, n_inner) = temp_matrix.col(0); start_idx += n_inner;
	r_.segment(start_idx, n_inner) = temp_matrix.col(1); start_idx += n_inner;
	r_.segment(start_idx, n_inner) = temp_matrix.col(2); start_idx += n_inner;

	// seam
	temp_matrix = lambda_seam_ * (m_se_ * new_pos_);
	r_.segment(start_idx, n_seam) = temp_matrix.col(0); start_idx += n_seam;
	r_.segment(start_idx, n_seam) = temp_matrix.col(1); start_idx += n_seam;
	r_.segment(start_idx, n_seam) = temp_matrix.col(2); start_idx += n_seam;

	// curvature
	ComputeNormK();
	r_.segment(start_idx, n_inner) = lambda_dev_ * Eigen::pow(norm_K_.array(), p_).matrix();
	start_idx += n_inner;

	// exp close
	r_exp_close_ = Eigen::exp(lambda_exp_close_ * (
		((new_pos_ - ori_pos_) / m_dist_).rowwise().squaredNorm().array() - 1.0)).matrix();
	r_.segment(start_idx, n_v) = r_exp_close_;
	start_idx += n_v;

	// exp smooth
	r_exp_smooth_ = Eigen::exp(lambda_exp_smooth_ * (
		(m_sm_ * new_pos_).rowwise().squaredNorm() - m_sm_sq_).array()).matrix();
	r_.segment(start_idx, n_inner) = r_exp_smooth_;
	start_idx += n_inner;

	// exp seam
	r_exp_seam_ = Eigen::exp(lambda_exp_seam_ * (
		(m_se_ * new_pos_).rowwise().squaredNorm() - m_se_sq_).array()).matrix();
	r_.segment(start_idx, n_seam) = r_exp_seam_;
	start_idx += n_seam;
}

void GaussMin::UpdateEnergy()
{
	ComputeR();

	// energy
	close_energy_ = (m_id_ * (new_pos_ - ori_pos_)).squaredNorm();
	smooth_energy_ = (m_sm_ * new_pos_).squaredNorm();
	seam_energy_ = (m_se_ * new_pos_).squaredNorm();
	dev_energy_ = Eigen::pow(norm_K_.array(), p_).matrix().squaredNorm();

	total_energy_ = 0.5 * r_.squaredNorm();
}

void GaussMin::ComputeJ()
{
	int n_v = mesh_.n_vertices();
	int n_inner = inner2glo_.size();
	int n_seam = seam2glo_.size();

	// close, smooth, seam
	std::vector<T> temp_triplet = init_triplet_;
	int start_idx = 6 * n_v;

	// K
	std::vector<std::vector<T>> K_triplet(6);

#pragma omp parallel for num_threads(6)
	for (int i = 0; i < 6; i++)
	{
		//printf("i = % d, ThreadId = % d\n", i, omp_get_thread_num());
		//i_triplet[i] = ComputeTempTriplet(i);
		K_triplet[i] = ComputeTempTriplet(start_idx, i);
	}


	for (int i = 0; i < 6; i++)
	{
		temp_triplet.insert(temp_triplet.begin(), K_triplet[i].begin(), K_triplet[i].end());
	}

	start_idx += n_inner;

	// exp close
	Eigen::MatrixX3d temp_matrix = lambda_exp_close_ * r_exp_close_.asDiagonal()
		* (new_pos_ - ori_pos_) / m_dist_ / m_dist_;
	std::vector<T> exp_close_triplet;
	for (size_t i = 0; i < n_v; i++)
	{
		exp_close_triplet.push_back(T(start_idx + i, i, temp_matrix(i, 0)));
		exp_close_triplet.push_back(T(start_idx + i, i + n_v, temp_matrix(i, 1)));
		exp_close_triplet.push_back(T(start_idx + i, i + 2 * n_v, temp_matrix(i, 2)));
	}
	start_idx += n_v;

	// exp smooth
	temp_matrix = lambda_exp_smooth_ * r_exp_smooth_.asDiagonal() * m_sm_ * new_pos_;
	std::vector<T> exp_smooth_triplet;
	for (size_t i = 0; i < n_inner; i++)
	{
		exp_smooth_triplet.push_back(T(start_idx + i, i, temp_matrix(i, 0)));
		exp_smooth_triplet.push_back(T(start_idx + i, i + n_v, temp_matrix(i, 1)));
		exp_smooth_triplet.push_back(T(start_idx + i, i + 2 * n_v, temp_matrix(i, 2)));
	}
	start_idx += n_inner;

	// exp seam
	temp_matrix = lambda_exp_seam_ * m_se_.transpose() * r_exp_seam_.asDiagonal() * m_se_ * new_pos_;
	std::vector<T> exp_seam_triplet;
	for (size_t i = 0; i < n_seam; i++)
	{
		exp_seam_triplet.push_back(T(start_idx + i, i, 0));
		exp_seam_triplet.push_back(T(start_idx + i, i + n_v, 0));
		exp_seam_triplet.push_back(T(start_idx + i, i + 2 * n_v, 0));
	}
	start_idx += n_seam;

	temp_triplet.insert(temp_triplet.begin(), exp_close_triplet.begin(), exp_close_triplet.end());
	temp_triplet.insert(temp_triplet.begin(), exp_smooth_triplet.begin(), exp_smooth_triplet.end());
	temp_triplet.insert(temp_triplet.begin(), exp_seam_triplet.begin(), exp_seam_triplet.end());

	J_.resize(start_idx, 3 * n_v);
	J_.setFromTriplets(temp_triplet.begin(), temp_triplet.end());
}

void GaussMin::SetGloToLoc()
{
	int n_v = mesh_.n_vertices();
	const auto& seam_status = seg_mesh_.GetSeam();
	
	std::vector<bool> seam_v_status(n_v, false);

	HEH temp_he;
	int to_idx, from_idx;
	for (const EH& e_h : mesh_.edges())
	{
		if (!seam_status[e_h.idx()]) continue;

		temp_he = mesh_.halfedge_handle(e_h, 0);

		to_idx = mesh_.to_vertex_handle(temp_he).idx();
		from_idx = mesh_.from_vertex_handle(temp_he).idx();

		seam_v_status[to_idx] = true;
		seam_v_status[from_idx] = true;
	}

	int n_inner = 0;
	glo2inner_.assign(n_v, -1);

	seam2glo_.clear();
	inner2glo_.clear();

	for (size_t i = 0; i < n_v; i++)
	{
		if (seam_v_status[i] || mesh_.is_boundary(VH(i)))
		{
			seam2glo_.push_back(i);
		}
		else
		{
			inner2glo_.push_back(i);
			glo2inner_[i] = n_inner;
			++n_inner;
		}
	}
}

void GaussMin::SetInnerSmoothMatrix()
{
	int n_inner = inner2glo_.size();

	std::vector<T> triplet;
	triplet.reserve(n_inner * 10 * 2);
	for (size_t i = 0; i < n_inner; i++)
	{
		for (const VH& adj_v : mesh_.vv_range(VH(inner2glo_[i])))
		{
			triplet.push_back(T(i, adj_v.idx(), -1.0));
			triplet.push_back(T(i, inner2glo_[i], 1.0));
		}
	}

	m_sm_.resize(n_inner, mesh_.n_vertices());
	m_sm_.setFromTriplets(triplet.begin(), triplet.end());
}

void GaussMin::SetSeamSmoothMatrix()
{
	int n_seam = seam2glo_.size();
	const auto& seam_status = seg_mesh_.GetSeam();

	std::vector<T> triplet;
	triplet.reserve(n_seam * 10 * 2);
	for (int i = 0; i < n_seam; i++)
	{
		// boundary
		if (mesh_.is_boundary(VH(seam2glo_[i])))
		{
			for (const HEH& adj_he : mesh_.voh_range(VH(seam2glo_[i])))
			{
				if (!mesh_.is_boundary(mesh_.to_vertex_handle(adj_he))) continue;

				triplet.push_back(T(i, mesh_.to_vertex_handle(adj_he).idx(), -1.0));
				triplet.push_back(T(i, seam2glo_[i], 1.0));
			}
		}
		else
		{
			for (const HEH& adj_he : mesh_.voh_range(VH(seam2glo_[i])))
			{
				if (!seam_status[mesh_.edge_handle(adj_he).idx()]) continue;

				triplet.push_back(T(i, mesh_.to_vertex_handle(adj_he).idx(), -1.0));
				triplet.push_back(T(i, seam2glo_[i], 1.0));
			}
		}
	}

	m_se_.resize(n_seam, mesh_.n_vertices());
	m_se_.setFromTriplets(triplet.begin(), triplet.end());
}

void GaussMin::SetIdentityMatrix()
{
	m_id_.resize(mesh_.n_vertices(), mesh_.n_vertices());
	m_id_.setIdentity();
}

void GaussMin::SetInitTriplet()
{
	int n_v = mesh_.n_vertices();
	int n_inner = inner2glo_.size();
	int n_seam = seam2glo_.size();

	init_triplet_.clear();
	init_triplet_.reserve(3 * m_id_.nonZeros() + 3 * m_sm_.nonZeros() + 3 * m_se_.nonZeros());

	int start_idx = 0;

	for (int k = 0; k < m_id_.outerSize(); ++k)
	{
		for (SpMat::InnerIterator it(m_id_, k); it; ++it)
		{
			init_triplet_.push_back(T(start_idx + it.row(),
				it.col(), 
				lambda_close_ * it.value()));
			init_triplet_.push_back(T(start_idx + n_v + it.row(),
				n_v + it.col(),
				lambda_close_ * it.value()));
			init_triplet_.push_back(T(start_idx + 2 * n_v + it.row(),
				2 * n_v + it.col(), 
				lambda_close_ * it.value()));
		}
	}
	start_idx += 3 * n_v;

	for (int k = 0; k < m_sm_.outerSize(); ++k)
	{
		for (SpMat::InnerIterator it(m_sm_, k); it; ++it)
		{
			init_triplet_.push_back(T(start_idx + it.row(),
				it.col(),
				lambda_smooth_ * it.value()));
			init_triplet_.push_back(T(start_idx + n_inner + it.row(),
				n_v + it.col(),
				lambda_smooth_ * it.value()));
			init_triplet_.push_back(T(start_idx + 2 * n_inner + it.row(),
				2 * n_v + it.col(),
				lambda_smooth_ * it.value()));
		}
	}
	start_idx += 3 * n_inner;

	for (int k = 0; k < m_se_.outerSize(); ++k)
	{
		for (SpMat::InnerIterator it(m_se_, k); it; ++it)
		{
			init_triplet_.push_back(T(start_idx + it.row(),
				it.col(),
				lambda_seam_ * it.value()));
			init_triplet_.push_back(T(start_idx + n_seam + it.row(),
				n_v + it.col(),
				lambda_seam_ * it.value()));
			init_triplet_.push_back(T(start_idx + 2 * n_seam + it.row(),
				2 * n_v + it.col(),
				lambda_seam_ * it.value()));
		}
	}
}

void GaussMin::Opt()
{
	const int& n_v = mesh_.n_vertices();

	double mean_len = MeshTools::AverageEdgeLength(mesh_);

	norm_K_.resize(inner2glo_.size(), 1);
	m_direction_.resize(n_v, 3);
	r_.resize(6 * n_v + inner2glo_.size() + 2 * n_v);

	double prev_energy = DBL_MAX;
	double min_energy = DBL_MAX;
	Eigen::MatrixX3d min_pos = ori_pos_;

	double gamma;
	int update_cout = 0;
	bool search_status = false;
	for (size_t i = 0; i < iter_num_; i++)
	{
		UpdateEnergy();

		//std::cout << i << " th: "
		//	<< " total energy: " << total_energy_
		//	<< " close energy: " << close_energy_
		//	<< " seam energy: " << seam_energy_
		//	<< " smooth energy: " << smooth_energy_
		//	<< " dev energy: " << dev_energy_ << std::endl;

		//std::cout << "energy error: " << prev_energy - dev_energy_ << std::endl;
		if (abs(prev_energy - dev_energy_) < break_diff_)
		{
			break;
		}
		prev_energy = dev_energy_;

		//std::cout << cur_bound_ << std::endl;
		//std::cout << norm_K_ << std::endl;

		if (min_energy > dev_energy_)
		{
			min_energy = dev_energy_;
			min_pos = new_pos_;
		}

		ComputeJ();

		if (i % 100 == 0)
		{
			m_JTJ_ = J_.transpose() * J_;

			m_ldlt_.compute(m_JTJ_);
		}

		v_gradient_ = J_.transpose() * r_;

		v_direction_ = -m_ldlt_.solve(v_gradient_);
		//v_direction_ = -v_gradient_;

		m_direction_.col(0) = v_direction_.segment(0, n_v);
		m_direction_.col(1) = v_direction_.segment(n_v, n_v);
		m_direction_.col(2) = v_direction_.segment(2 * n_v, n_v);


		//std::ofstream file_g("gradient.txt");
		//file_g << m_direction_ << std::endl;
		//file_g.close();

		//std::cout << m_direction_ << std::endl;

		//double max_norm = 0;
		//for (size_t i = 0; i < m_direction_.rows(); i++)
		//{
		//	max_norm = std::max(m_direction_.row(i).norm(), max_norm);
		//}

		//if (max_norm < step_precent_ * m_dist_ * m_dist_)
		//{
		//	break;
		//}

		//std::cout << i << " th: "
		//	<< " total energy: " << total_energy_
		//	<< " close energy: " << close_energy_
		//	<< " seam energy: " << seam_energy_
		//	<< " smooth energy: " << smooth_energy_
		//	<< " dev energy: " << dev_energy_
		//	<< " max norm: " << max_norm << std::endl;

		bool search_status = LinearSearch(gamma, 0.9);


		//std::cout << "gamma: " << gamma
		//	<< "  search status: " << search_status << std::endl;

		//system("pause");

		if (search_status)
		{
			new_pos_ = new_pos_ + gamma * m_direction_;

			//Projection();
		}
		else
		{
			break;
		}

		//if (update_cout > update_lambda_num_ - 1)
		//{
		//	lambda_dev_ *= update_lambda_precent_;
		//}
	}

	new_pos_ = min_pos;
}

bool GaussMin::LinearSearch(double& alpha_i, double c2)
{
	// strong Wolfe condition
	int search_num = 200;
	double c1 = 0;
	double t = 2;
	//double alpha_max = 1e10;

	double alpha_i_1 = 0;
	alpha_i = 1.0;

	double phi_0 = total_energy_;
	double phi_0_d = v_direction_.transpose() * v_gradient_;

	//std::cout << phi_0_d << std::endl;

	Eigen::MatrixX3d init_pos = new_pos_;

	bool search_status = false;
	double alpha_lo = alpha_i_1, alpha_hi = alpha_i;
	{
		double phi_i_d;

		for (size_t i = 0; i < search_num; i++)
		{
			new_pos_ = init_pos + alpha_i * m_direction_;

			//Projection();
			UpdateEnergy();

			if (total_energy_ > phi_0 + c1 * alpha_i * phi_0_d)
			{
				alpha_lo = alpha_i_1;
				alpha_hi = alpha_i;
				break;
			}

			ComputeJ();
			v_gradient_ = J_.transpose() * r_;
			phi_i_d = v_direction_.transpose() * v_gradient_;

			if (abs(phi_i_d) < -c2 * phi_0_d)
			{
				search_status = true;
				break;
			}

			if (phi_i_d >= 0)
			{
				alpha_lo = alpha_i;
				alpha_hi = alpha_i_1;
				break;
			}

			alpha_i_1 = alpha_i;
			alpha_i = alpha_i_1 * t;
		}
	}

	if (!search_status)
	{
		double phi_i_d;
		
		new_pos_ = init_pos + alpha_lo * m_direction_;
		UpdateEnergy();
		double phi_lo = total_energy_;
		
		int cout = 0;
		do
		{
			alpha_i = (alpha_lo + alpha_hi) / 2;

			new_pos_ = init_pos + alpha_i * m_direction_;
			UpdateEnergy();

			//std::cout << "cout: " << cout
			//	<< " total energy: " << total_energy_
			//	<< " close energy: " << close_energy_
			//	<< " seam energy: " << seam_energy_
			//	<< " smooth energy: " << smooth_energy_
			//	<< " dev energy: " << dev_energy_ << std::endl;

			if (total_energy_ > phi_0 + c1 * alpha_i * phi_0_d 
				|| total_energy_ >= phi_lo)
			{
				alpha_hi = alpha_i;
			}
			else
			{
				ComputeJ();
				v_gradient_ = J_.transpose() * r_;
				phi_i_d = v_direction_.transpose() * v_gradient_;

				if (abs(phi_i_d) < -c2 * phi_0_d)
				{
					search_status = true;
					break;
				}

				if (phi_i_d * (alpha_hi - alpha_lo) >= 0)
				{
					alpha_hi = alpha_lo;
				}

				alpha_lo = alpha_i;
				
				new_pos_ = init_pos + alpha_lo * m_direction_;
				UpdateEnergy();
				phi_lo = total_energy_;
			}

			cout++;
		} while (cout < search_num);
	}

	new_pos_ = init_pos;
	return total_energy_ < phi_0;
}

std::vector<T> GaussMin::ComputeTempTriplet(const int& start_idx,
	const int& thread_id)
{
	int coor_idx = thread_id / 2;

	double sign;
	if (thread_id % 2 == 0)
	{
		sign = 1;
	}
	else
	{
		sign = -1;
	}

	int n_v = mesh_.n_vertices();
	std::vector<T> triplet;
	triplet.reserve(10 * n_v);

	Eigen::MatrixX3d change_pos = new_pos_;
	for (int v_id = 0; v_id < n_v; v_id++)
	{
		const double& temp_step = sign * step_h_[v_id];

		int val_id = v_id + coor_idx * n_v;

		change_pos(v_id, coor_idx) += temp_step;

		if (glo2inner_[v_id] != -1)
		{
			triplet.push_back(T(
				start_idx + glo2inner_[v_id],
				val_id,
				lambda_dev_ * 
				0.5 * AngleDefect(VH(v_id), change_pos) / temp_step	
				/ cur_bound_ * p_ 
				* pow(norm_K_(glo2inner_[v_id]), p_ - 1)));
		}

		for (const VH& adj_v : mesh_.vv_range(VH(v_id)))
		{
			if (glo2inner_[adj_v.idx()] != -1)
			{
				triplet.push_back(T(
					start_idx + glo2inner_[adj_v.idx()],
					val_id,
					lambda_dev_ *
					0.5 * AngleDefect(adj_v, change_pos) / temp_step
					/ cur_bound_ * p_ 
					* pow(norm_K_(glo2inner_[adj_v.idx()]), p_ - 1)));
			}
		}

		change_pos(v_id, coor_idx) -= temp_step;
	}

	return triplet;
}

template<typename Scalar>
inline Scalar GaussMin::AngleDefect(const VH& v_h, const Eigen::Matrix<Scalar, -1, 3>& temp_pos)
{
	Scalar scalar_pi = M_PI;
	Scalar angle_defect = 2 * scalar_pi;

	int cur_idx, prev_idx, next_idx;
	Eigen::Matrix<Scalar, -1, 1> cur_vec, prev_vec;
	for (const HEH& cur_he : mesh_.voh_range(v_h))
	{
		if (!mesh_.face_handle(cur_he).is_valid()) continue;

		cur_idx = v_h.idx();
		prev_idx = mesh_.opposite_vh(cur_he).idx();
		next_idx = mesh_.to_vertex_handle(cur_he).idx();

		cur_vec = (temp_pos.row(next_idx) - temp_pos.row(cur_idx)).normalized();
		prev_vec = (temp_pos.row(cur_idx) - temp_pos.row(prev_idx)).normalized();

		angle_defect -= scalar_pi - acos(cur_vec.dot(prev_vec));
	}

	return angle_defect;
}
