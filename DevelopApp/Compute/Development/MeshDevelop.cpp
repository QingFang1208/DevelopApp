#include "MeshDevelop.h"
#include <iostream>
#include <omp.h>
#include "svd/svd.hpp"
using OpenMesh::Vec3d;

MeshDevelop::MeshDevelop()
{
}

double meanEdgeLength(const Mesh& mesh)
{
	double meanEdgeLen = 0;
	for (auto eh : mesh.edges())
	{
		meanEdgeLen += mesh.calc_edge_length(eh);
	}
	meanEdgeLen /= mesh.n_edges();
	return meanEdgeLen;
}

double meanEdgeScale(const Mesh& mesh, Vec3d* pt, double* fA, double* eA, double* vA)
{
	double meanEdgeLen = 0;
	for (auto eh : mesh.edges())
	{
		meanEdgeLen += mesh.calc_edge_length(eh);
	}
	meanEdgeLen /= mesh.n_edges();

	for (auto vh : mesh.vertices())
	{
		pt[vh.idx()] = mesh.point(vh) / meanEdgeLen;
	}

	for (auto fh : mesh.faces())
	{
		auto he_h = mesh.halfedge_handle(fh);
		auto next_he = mesh.next_halfedge_handle(he_h);

		OpenMesh::Vec3d he_vector = mesh.calc_edge_vector(he_h);
		OpenMesh::Vec3d next_vector = mesh.calc_edge_vector(next_he);

		fA[fh.idx()] = 0.5 * cross(he_vector, next_vector).norm() / meanEdgeLen / meanEdgeLen;

		//fA[fh.idx()] = mesh.calc_face_area(fh) / meanEdgeLen / meanEdgeLen;
	}

	for (auto vh : mesh.vertices())
	{
		vA[vh.idx()] = 0;
	}

	for (auto eh : mesh.edges())
	{
		eA[eh.idx()] = 0;
	}

	for (auto fh : mesh.faces())
	{
		for (auto fhh : mesh.fh_range(fh))
		{
			eA[fhh.idx() >> 1] += fA[fh.idx()] / 3;
			vA[mesh.to_vertex_handle(fhh).idx()] += fA[fh.idx()] / 3;
		}
	}

	return meanEdgeLen;
}

double unifyAreaScale(const Mesh& mesh, Vec3d* pt, double* vA, double* fA)
{
	double Area = 0;
	for (auto fh : mesh.faces())
	{
		auto he_h = mesh.halfedge_handle(fh);
		auto next_he = mesh.next_halfedge_handle(he_h);

		OpenMesh::Vec3d he_vector = mesh.calc_edge_vector(he_h);
		OpenMesh::Vec3d next_vector = mesh.calc_edge_vector(next_he);

		fA[fh.idx()] = 0.5 * cross(he_vector, next_vector).norm();

		//fA[fh.idx()] = mesh.calc_face_area(fh);
		Area += fA[fh.idx()];
	}

	for (auto fh : mesh.faces())
	{
		fA[fh.idx()] /= Area;
	}

	for (auto vh : mesh.vertices())
	{
		vA[vh.idx()] = 0;
		for (auto vfh : mesh.vf_range(vh))
		{
			vA[vh.idx()] += fA[vfh.idx()];
		}
		vA[vh.idx()] /= 3;
	}

	double sqrtArea = sqrt(Area);

	for (auto vh : mesh.vertices())
	{
		pt[vh.idx()] = mesh.point(vh) / sqrtArea;
	}

	return sqrtArea;
}

void meshInfo(const Mesh& mesh, Vec3d* pt, double* fA, double* eA, double* vA)
{
	for (int i = 0; i < mesh.n_vertices(); ++i)
	{
		pt[i] = mesh.points()[i];
	}

	for (auto vh : mesh.vertices())
	{
		vA[vh.idx()] = 0;
	}

	for (auto eh : mesh.edges())
	{
		eA[eh.idx()] = 0;
	}

	for (auto fh : mesh.faces())
	{
		auto he_h = mesh.halfedge_handle(fh);
		auto next_he = mesh.next_halfedge_handle(he_h);

		OpenMesh::Vec3d he_vector = mesh.calc_edge_vector(he_h);
		OpenMesh::Vec3d next_vector = mesh.calc_edge_vector(next_he);

		fA[fh.idx()] = 0.5 * cross(he_vector, next_vector).norm();


		//fA[fh.idx()] = mesh.calc_face_area(fh);

		for (auto fhh : mesh.fh_range(fh))
		{
			eA[fhh.idx() >> 1] += fA[fh.idx()] / 3;
			vA[mesh.to_vertex_handle(fhh).idx()] += fA[fh.idx()] / 3;
		}
	}
}

double unifyScale(const Mesh& mesh, Vec3d* pt, double* fA, double* eA, double* vA)
{
	Vec3d boxMin(DBL_MAX), boxMax(-DBL_MAX);

	for (auto vh : mesh.vertices())
	{
		boxMin.minimize(mesh.point(vh));
		boxMax.maximize(mesh.point(vh));
	}

	double halfBoxLen = (boxMax - boxMin).max() / 2;

	for (auto vh : mesh.vertices())
	{
		pt[vh.idx()] = mesh.point(vh) / halfBoxLen;
	}

	for (auto fh : mesh.faces())
	{
		auto he_h = mesh.halfedge_handle(fh);
		auto next_he = mesh.next_halfedge_handle(he_h);

		OpenMesh::Vec3d he_vector = mesh.calc_edge_vector(he_h);
		OpenMesh::Vec3d next_vector = mesh.calc_edge_vector(next_he);

		fA[fh.idx()] = 0.5 * cross(he_vector, next_vector).norm() / halfBoxLen / halfBoxLen;


		//fA[fh.idx()] = mesh.calc_face_area(fh) / halfBoxLen / halfBoxLen;
	}

	for (auto vh : mesh.vertices())
	{
		vA[vh.idx()] = 0;
	}

	for (auto eh : mesh.edges())
	{
		eA[eh.idx()] = 0;
	}

	for (auto fh : mesh.faces())
	{
		for (auto fhh : mesh.fh_range(fh))
		{
			eA[fhh.idx() >> 1] += fA[fh.idx()] / 2;
			vA[mesh.to_vertex_handle(fhh).idx()] += fA[fh.idx()] / 3;
		}
	}

	return halfBoxLen;
}

double boxDiagonalLen(Vec3d* pt, int vnum)
{
	Vec3d boxMin(DBL_MAX), boxMax(-DBL_MAX);

	for (int i = 0; i < vnum; ++i)
	{
		boxMin.minimize(pt[i]);
		boxMax.maximize(pt[i]);
	}

	return (boxMax - boxMin).norm();
}

void computeNormals(const MeshTopo& topo, MatrixXd& vpfnen)
{
	auto vp = Eigen::Map<Eigen::MatrixXd>(&vpfnen(0, 0), 3, topo.vN);
	auto fn = Eigen::Map<Eigen::MatrixXd>(&vpfnen(0, topo.vN), 3, topo.fN);
	auto en = Eigen::Map<Eigen::MatrixXd>(&vpfnen(0, topo.vN + topo.fN), 3, topo.eN);

	for (int i = 0; i < topo.fN; ++i)
	{
		int v[3];
		v[0] = topo.f2v[i][0];
		Vector3d ea = vp.col(topo.f2v[i][1]) - vp.col(topo.f2v[i][0]);
		Vector3d eb = vp.col(topo.f2v[i][2]) - vp.col(topo.f2v[i][1]);
		fn.col(i) = ea.cross(eb);
		fn.col(i).normalize();
	}

	en.setZero();
	for (int i = 0; i < topo.hN; ++i)
	{
		if (topo.h2f[i] < 0) continue;
		en.col(i >> 1) += fn.col(topo.h2f[i]);
	}

	for (int i = 0; i < topo.eN; ++i)
	{
		en.col(i).normalize();
	}
}

void MeshDevelop::oriMesh(const Mesh& mesh)
{
	tmp = mesh;

	m_topo.init(mesh);
	m_vpfnen0 = MatrixXd::Zero(3, m_topo.vN + m_topo.fN + m_topo.eN + m_topo.vN);

	m_fixV.resize(m_topo.vN);
	for (int i = 0; i < m_topo.vN; ++i) m_fixV[i] = false;
	for (int i = 0; i < m_topo.vN; ++i)
	{
		if (m_topo.vb[i]) m_fixV[i] = true;
	}

	m_svdA.resize(m_topo.hN);
	m_svdM.resize(m_topo.hN);
	m_svdU = m_svdV = m_svdM;
	m_svd.resize(m_topo.hN);
	m_eig.resize(m_topo.hN);

	m_fA = VectorXd::Zero(m_topo.fN);
	m_eA = VectorXd::Zero(m_topo.eN);
	m_vA = VectorXd::Zero(m_topo.vN);

	m_scale = unifyScale(mesh, (Vec3d*)m_vpfnen0.data(), m_fA.data(), m_eA.data(), m_vA.data());
	computeNormals(m_topo, m_vpfnen0);
	m_dist = boxDiagonalLen((Vec3d*)m_vpfnen0.data(), m_topo.vN);

	m_fA = m_fA.mean() * VectorXd::Ones(m_fA.size());
	m_eA = m_fA.mean() * VectorXd::Ones(m_eA.size());

	m_el = m_ew = m_ed = VectorXd::Zero(m_topo.eN);
	m_eb0 = MatrixXd::Zero(3, m_topo.eN);
	for (int i = 0; i < m_topo.eN; ++i)
	{
		m_eb0.col(i) = m_vpfnen0.col(m_topo.h2v[i << 1]) - m_vpfnen0.col(m_topo.h2v[(i << 1) + 1]);
		m_el(i) = m_eb0.col(i).norm();
		assert(m_el(i) > 0);
		m_eb0.col(i) /= m_el(i);
	}

	std::vector<Vec3d> f2a, f2cot;
	faceAngle(mesh.points(), m_topo, f2a);

	f2cot.assign(f2a.size(), Vec3d(0));
	for (int i = 0; i < f2a.size(); ++i)
	{
		f2cot[i][0] = 0.5;//fabs(0.5 / tan(f2a[i][0]));
		f2cot[i][1] = 0.5;//fabs(0.5 / tan(f2a[i][1]));
		f2cot[i][2] = 0.5;//fabs(0.5 / tan(f2a[i][2]));
	}

	for (int i = 0; i < f2cot.size(); ++i)
	{
		m_ew(m_topo.f2h[i][0] >> 1) += f2cot[i][1];
		m_ew(m_topo.f2h[i][1] >> 1) += f2cot[i][2];
		m_ew(m_topo.f2h[i][2] >> 1) += f2cot[i][0];
	}

	int h0, h1;

	for (int i = 0; i < m_topo.eN; ++i)
	{
		h0 = i << 1;
		h1 = h0 + 1;

		if (m_topo.h2f[h0] < 0 || m_topo.h2f[h1] < 0) continue;

		m_ed(i) += 2 * m_eA(i);
		m_ed(m_topo.phf2e[h0][1]) += m_eA(i);
		m_ed(m_topo.phf2e[h0][2]) += m_eA(i);
		m_ed(m_topo.phf2e[h1][1]) += m_eA(i);
		m_ed(m_topo.phf2e[h1][2]) += m_eA(i);
	}

	for (int i = 0; i < m_topo.fN; ++i)
	{
		m_ed(m_topo.f2h[i][0] >> 1) += m_fA(i);
		m_ed(m_topo.f2h[i][1] >> 1) += m_fA(i);
		m_ed(m_topo.f2h[i][2] >> 1) += m_fA(i);
	}

	std::vector<Eigen::Triplet<double>> trips;
	trips.reserve(9 * f2a.size());
	for (int i = 0; i < f2a.size(); ++i)
	{
		trips.emplace_back(m_topo.f2v[i][0], m_topo.f2v[i][0], f2cot[i][1] + f2cot[i][2]);
		trips.emplace_back(m_topo.f2v[i][1], m_topo.f2v[i][1], f2cot[i][0] + f2cot[i][2]);
		trips.emplace_back(m_topo.f2v[i][2], m_topo.f2v[i][2], f2cot[i][0] + f2cot[i][1]);
		trips.emplace_back(m_topo.f2v[i][0], m_topo.f2v[i][1], -f2cot[i][2]);
		trips.emplace_back(m_topo.f2v[i][1], m_topo.f2v[i][0], -f2cot[i][2]);
		trips.emplace_back(m_topo.f2v[i][1], m_topo.f2v[i][2], -f2cot[i][0]);
		trips.emplace_back(m_topo.f2v[i][2], m_topo.f2v[i][1], -f2cot[i][0]);
		trips.emplace_back(m_topo.f2v[i][2], m_topo.f2v[i][0], -f2cot[i][1]);
		trips.emplace_back(m_topo.f2v[i][0], m_topo.f2v[i][2], -f2cot[i][1]);
	}

	m_lap.resize(m_topo.vN, m_topo.vN);
	m_lap.setFromTriplets(trips.begin(), trips.end());
	SPMatrixXd I(m_topo.vN, m_topo.vN);
	I.setIdentity();
	m_ldlt.compute(m_lap + m_wP * I * m_vA.asDiagonal());

	Eigen::MatrixXi faceE(3, m_topo.fN);
	for (int i = 0; i < m_topo.fN; i++)
	{
		faceE(0, i) = m_topo.f2h[i][0] >> 1;
		faceE(1, i) = m_topo.f2h[i][1] >> 1;
		faceE(2, i) = m_topo.f2h[i][2] >> 1;
	}

	Eigen::VectorXd faceSL = DBL_MAX * Eigen::VectorXd::Ones(m_topo.fN);
	for (int i = 0; i < m_topo.fN; i++)
	{
		faceSL(i) = faceSL(i) < m_el(faceE(0, i)) ? faceSL(i) : m_el(faceE(0, i));
		faceSL(i) = faceSL(i) < m_el(faceE(1, i)) ? faceSL(i) : m_el(faceE(1, i));
		faceSL(i) = faceSL(i) < m_el(faceE(2, i)) ? faceSL(i) : m_el(faceE(2, i));
	}

	triNeq.init(m_el, faceE, 1e-4, faceSL * 0.1);

	EN0_i_00.resize(m_topo.eN, 0);
	EN0_i_01.resize(m_topo.eN, 0);
	EN0_i_02.resize(m_topo.eN, 0);
	EN0_i_03.resize(m_topo.eN, 0);
	EN0_i_10.resize(m_topo.eN, 0);
	EN0_i_11.resize(m_topo.eN, 0);
	EN0_i_12.resize(m_topo.eN, 0);
	EN0_i_13.resize(m_topo.eN, 0);
	EN0_i_20.resize(m_topo.eN, 0);
	EN0_i_21.resize(m_topo.eN, 0);
	EN0_i_22.resize(m_topo.eN, 0);
	EN0_i_23.resize(m_topo.eN, 0);

	EN_i_00.resize(m_topo.eN, 0);
	EN_i_01.resize(m_topo.eN, 0);
	EN_i_02.resize(m_topo.eN, 0);
	EN_i_03.resize(m_topo.eN, 0);
	EN_i_10.resize(m_topo.eN, 0);
	EN_i_11.resize(m_topo.eN, 0);
	EN_i_12.resize(m_topo.eN, 0);
	EN_i_13.resize(m_topo.eN, 0);
	EN_i_20.resize(m_topo.eN, 0);
	EN_i_21.resize(m_topo.eN, 0);
	EN_i_22.resize(m_topo.eN, 0);
	EN_i_23.resize(m_topo.eN, 0);

	e0fn0_i_0.resize(m_topo.eN, 0);
	e0fn0_i_1.resize(m_topo.eN, 0);
	e0fn0_i_2.resize(m_topo.eN, 0);
	e0fn_i_0.resize(m_topo.eN, 0);
	e0fn_i_1.resize(m_topo.eN, 0);
	e0fn_i_2.resize(m_topo.eN, 0);
	e1fn0_i_0.resize(m_topo.eN, 0);
	e1fn0_i_1.resize(m_topo.eN, 0);
	e1fn0_i_2.resize(m_topo.eN, 0);
	e1fn_i_0.resize(m_topo.eN, 0);
	e1fn_i_1.resize(m_topo.eN, 0);
	e1fn_i_2.resize(m_topo.eN, 0);

	diag_i_0.resize(m_topo.eN, 0);
	diag_i_1.resize(m_topo.eN, 0);
	diag_i_2.resize(m_topo.eN, 0);
	diag_i_3.resize(m_topo.eN, 0);
}

void MeshDevelop::tarMesh(const Mesh& mesh, double dist)
{
	m_dist *= dist;

	VectorXd fA = VectorXd::Zero(m_topo.fN);
	VectorXd eA = VectorXd::Zero(m_topo.eN);
	VectorXd vA = VectorXd::Zero(m_topo.vN);

	m_vpfnen = MatrixXd::Zero(3, m_topo.vN + m_topo.fN + m_topo.eN + m_topo.vN);
	meshInfo(mesh, (Vec3d*)m_vpfnen.data(), fA.data(), eA.data(), vA.data());
	m_vpfnen.block(0, 0, 3, m_topo.vN) /= m_scale;
	computeNormals(m_topo, m_vpfnen);

	MatrixXd vp_init = m_vpfnen.block(0, 0, 3, m_topo.vN);
	AA = std::make_unique<AndersonAcceleration>(5, vp_init.size(), vp_init.size());
	AA->init(vp_init);

	m_vt = m_vpfnen.block(0, 0, 3, m_topo.vN);
	updateVT();

	m_es = VectorXd::Ones(m_topo.eN);

	std::vector<Matrix3d> fDhDeReU(m_topo.fN + 2 * m_topo.hN, Matrix3d::Identity(3, 3));
	m_fDhDeReU = Eigen::Map<MatrixXd>(&fDhDeReU[0](0, 0), 9, fDhDeReU.size());

	for (int i = 0; i < 10; ++i)
	{
		updateES();
		updateER();
		updateEU();
	}
}

double MeshDevelop::develop(int maxIter)
{
	//double energy_last = DBL_MAX;
	//double energy_new = DBL_MAX;
	//MatrixXd vp_last, vp_new;
	//double err_last = DBL_MAX, err_new = DBL_MAX;
	//auto vp = m_vpfnen.block(0, 0, 3, m_topo.vN);
	//vp_last = vp_new = vp;

	//localStep();
	//int rejectionNum = 0;

	//for (int i = 0; i < maxIter; ++i)
	//{
	//	globalStep();

	//	vp_last = vp;
	//	err_last = energy();
	//	AA->compute(vp_last, vp_new);
	//	vp = vp_new;
	//	localStep();
	//	err_new = energy();

	//	if (err_new > err_last)
	//	{
	//		++rejectionNum;
	//		vp = vp_last;
	//		AA->reset(vp_last);
	//		localStep();
	//	}

	//	// if (energyRP() >= 1e-16) break;

	//	if ((i + 1) % 100 == 0)
	//	{
	//		energy_last = energy_new;
	//		energy_new = energyD();

	//		std::cout << i << "th : " << energy() << '\t' << m_wD * energyD()
	//			<< '\t' << energyDP() + energyFN() << '\t' << m_wP * energyRP() << '\t' << m_wD
	//			<< "\tRejection: " << rejectionNum << "/100" << std::endl;

	//		rejectionNum = 0;

	//		if (energy_last < DBL_MAX && energy_last - energy_new < 5e-3 * energy_last)
	//		{
	//			m_wD += 1e-4;
	//			vp_new = vp;
	//			AA->reset(vp_new);
	//			localStep();
	//		}
	//	}


	//	/*if ((i + 1) % 500 == 0)
	//	{
	//		m_wD = fmin(5e-3, m_wD + 1e-4);
	//		vp_new = vp;
	//		AA->reset(vp_new);
	//	}*/
	//}

	for (int i = 0; i < maxIter; ++i)
	{
		updateD();
		energyD();

		calcEN();
		updateER();
		updateVP();
		updateVT();
		updateES();

		calcFN();
		updateEU();

		updateER();
		updateVP();
		updateVT();
		updateES();

		if ((i + 1) % 100 == 0)
		{
			std::cout << i << "th : " << energy() << '\t' << energyD()
				<< '\t' << energyDP() + energyFN() << '\t' << energyRP() << '\t' << m_wD << std::endl;
		}
		updatePara(i);

		//if ((i + 1) % 1000 == 0)
		//{
		//	outMesh(tmp);
		//	MeshTools::WriteMesh(tmp, "./tmp" + std::to_string(i) + ".obj");
		//	//MeshTools::WriteMesh(tmp, "./MED/tmp" + std::to_string((i + 1) / 2000) + ".obj");
		//}
	}

	updateD();
	return energyD();
}

void MeshDevelop::updatePara(int iter)
{
	if ((iter + 1) % 1000 == 0)
	{
		m_wD = fmin(5, m_wD + 0.1);
	}

	//if ((iter + 1) % 1000 == 0)
	//{
	//	m_wD = fmin(1, m_wD + 0.1);
	//}
}

void MeshDevelop::outMesh(Mesh& mesh)
{
	Vec3d pt;
	for (auto vh : mesh.vertices())
	{
		pt[0] = m_vpfnen(0, vh.idx());
		pt[1] = m_vpfnen(1, vh.idx());
		pt[2] = m_vpfnen(2, vh.idx());
		pt *= m_scale;

		mesh.set_point(vh, pt);
	}
}

void MeshDevelop::localStep()
{
	updateD();
	updateER();
	updateEU();
	updateVT();
	updateES();
}

void MeshDevelop::globalStep()
{
	updateVP();
}

void MeshDevelop::updateN()
{
	updateD();	//std::cout << energy() << '\t' << m_wD * energyD() << '\t' << energyDP() + energyFN() << '\t' << m_wP * energyRP() << std::endl;

	calcEN();
	updateER();	//std::cout << energy() << '\t' << m_wD * energyD() << '\t' << energyDP() + energyFN() << '\t' << m_wP * energyRP() << std::endl;

	calcFN();
	updateEU();	//std::cout << energy() << '\t' << m_wD * energyD() << '\t' << energyDP() + energyFN() << '\t' << m_wP * energyRP() << std::endl;

}

void MeshDevelop::updateV()
{
	updateVP();	//std::cout << energy() << '\t' << m_wD * energyD() << '\t' << energyDP() + energyFN() << '\t' << m_wP * energyRP() << std::endl;
	updateVT();	//std::cout << energy() << '\t' << m_wD * energyD() << '\t' << energyDP() + energyFN() << '\t' << m_wP * energyRP() << std::endl;
	updateES();	//std::cout << energy() << '\t' << m_wD * energyD() << '\t' << energyDP() + energyFN() << '\t' << m_wP * energyRP() << std::endl;
}

//void MeshDevelop::updateD()
//{
//	auto en0 = m_vpfnen0.block(0, m_topo.vN + m_topo.fN, 3, m_topo.eN);
//	auto fD = (Matrix3d*)m_fDhDeReU.data();
//	auto hD = (Matrix3d*)m_fDhDeReU.data() + m_topo.fN;
//	auto eR = (Matrix3d*)m_fDhDeReU.data() + m_topo.fN + m_topo.hN;
//	MatrixXd eleRen = en0;
//
//	for (int i = 0; i < m_topo.eN; ++i)
//	{
//		eleRen.col(i) = eR[i] * en0.col(i);
//	}
//
//#pragma omp parallel for
//	for (int i = 0; i < m_topo.fN; ++i)
//	{
//		fD[i].col(0) = eleRen.col(m_topo.f2h[i][0] >> 1);
//		fD[i].col(1) = eleRen.col(m_topo.f2h[i][1] >> 1);
//		fD[i].col(2) = eleRen.col(m_topo.f2h[i][2] >> 1);
//
//		if (m_topo.h2f[m_topo.f2h[i][0] ^ 1] < 0 || m_topo.h2f[m_topo.f2h[i][1] ^ 1] < 0
//			|| m_topo.h2f[m_topo.f2h[i][2] ^ 1] < 0) continue;
//
//		//m_eig[i].compute(fD[i] * fD[i].transpose());
//		//m_svdA[i] = m_eig[i].eigenvectors().col(0);
//		//fD[i].col(0) -= m_svdA[i].dot(fD[i].col(0)) * m_svdA[i];
//		//fD[i].col(1) -= m_svdA[i].dot(fD[i].col(1)) * m_svdA[i];
//		//fD[i].col(2) -= m_svdA[i].dot(fD[i].col(2)) * m_svdA[i];
//
//		//if (fabs(fD[i].determinant()) > 1e-4)
//		//{
//#ifdef USE_FAST_SVD
//			computeMat3SVD(fD[i].data(), m_svdU[i].data(), m_svdV[i].data(), m_svdA[i].data());
//#else
//			m_svd[i].compute(fD[i], Eigen::ComputeThinU | Eigen::ComputeThinV);
//			m_svdU[i] = m_svd[i].matrixU();	m_svdV[i] = m_svd[i].matrixV();	m_svdA[i] = m_svd[i].singularValues();
//#endif // USE_FAST_SVD
//
//			m_svdA[i](2) = 0;	
//			fD[i] = m_svdU[i] * m_svdA[i].asDiagonal() * m_svdV[i].transpose();
//
//		//}
//	}
//
//#pragma omp parallel for
//	for (int i = 0; i < m_topo.hN; ++i)
//	{
//		if (m_topo.h2f[i] < 0 || m_topo.h2f[i ^ 1] < 0) continue;
//
//		hD[i].col(0) = eleRen.col(m_topo.phf2e[i][0]);
//		hD[i].col(1) = eleRen.col(m_topo.phf2e[i][1]);
//		hD[i].col(2) = eleRen.col(m_topo.phf2e[i][2]);
//
//		if (m_topo.h2f[m_topo.phf2e[i][1] << 1] < 0 || m_topo.h2f[m_topo.phf2e[i][1] << 1 ^ 1] < 0 ||
//			m_topo.h2f[m_topo.phf2e[i][2] << 1] < 0 || m_topo.h2f[m_topo.phf2e[i][2] << 1 ^ 1] < 0) continue;
//
//		m_eig[i].compute(hD[i] * hD[i].transpose());
//		m_svdA[i] = m_eig[i].eigenvectors().col(0);
//		hD[i].col(0) -= m_svdA[i].dot(hD[i].col(0)) * m_svdA[i];
//		hD[i].col(1) -= m_svdA[i].dot(hD[i].col(1)) * m_svdA[i];
//		hD[i].col(2) -= m_svdA[i].dot(hD[i].col(2)) * m_svdA[i];
//
//		//if (fabs(hD[i].determinant()) > 1e-4)
//		//{
//#ifdef USE_FAST_SVD
//			computeMat3SVD(hD[i].data(), m_svdU[i].data(), m_svdV[i].data(), m_svdA[i].data());
//#else
//			m_svd[i].compute(hD[i], Eigen::ComputeThinU | Eigen::ComputeThinV);
//			m_svdU[i] = m_svd[i].matrixU();	m_svdV[i] = m_svd[i].matrixV();	m_svdA[i] = m_svd[i].singularValues();
//#endif
//			m_svdA[i](2) = 0;	
//			hD[i] = m_svdU[i] * m_svdA[i].asDiagonal() * m_svdV[i].transpose();
//		//}
//	}
//}

void MeshDevelop::updateD()
{
	auto en0 = m_vpfnen0.block(0, m_topo.vN + m_topo.fN, 3, m_topo.eN);
	auto fD = (Matrix3d*)m_fDhDeReU.data();
	auto hD = (Matrix3d*)m_fDhDeReU.data() + m_topo.fN;
	auto eR = (Matrix3d*)m_fDhDeReU.data() + m_topo.fN + m_topo.hN;
	MatrixXd eleRen = en0;

	for (int i = 0; i < m_topo.eN; ++i)
	{
		eleRen.col(i) = eR[i] * en0.col(i);
	}

#pragma omp parallel for
	for (int i = 0; i < m_topo.fN; ++i)
	{
		fD[i].col(0) = eleRen.col(m_topo.f2h[i][0] >> 1);
		fD[i].col(1) = eleRen.col(m_topo.f2h[i][1] >> 1);
		fD[i].col(2) = eleRen.col(m_topo.f2h[i][2] >> 1);

		if (m_topo.h2f[m_topo.f2h[i][0] ^ 1] < 0 || m_topo.h2f[m_topo.f2h[i][1] ^ 1] < 0
			|| m_topo.h2f[m_topo.f2h[i][2] ^ 1] < 0) continue;

		//m_eig[i].compute(fD[i] * fD[i].transpose());
		//m_svdA[i] = m_eig[i].eigenvectors().col(0);
		//fD[i].col(0) -= m_svdA[i].dot(fD[i].col(0)) * m_svdA[i];
		//fD[i].col(1) -= m_svdA[i].dot(fD[i].col(1)) * m_svdA[i];
		//fD[i].col(2) -= m_svdA[i].dot(fD[i].col(2)) * m_svdA[i];

		//if (fabs(fD[i].determinant()) > 1e-4)
		//{
#ifdef USE_FAST_SVD
			computeMat3SVD(fD[i].data(), m_svdU[i].data(), m_svdV[i].data(), m_svdA[i].data());
#else
			m_svd[i].compute(fD[i], Eigen::ComputeThinU | Eigen::ComputeThinV);
			m_svdU[i] = m_svd[i].matrixU();	m_svdV[i] = m_svd[i].matrixV();	m_svdA[i] = m_svd[i].singularValues();
#endif // USE_FAST_SVD

			m_svdA[i](2) = 0;
			//fD[i] = m_svdU[i] * m_svdA[i].asDiagonal() * m_svdV[i].transpose();

			const double& U_00 = m_svdU[i](0, 0);
			const double& U_01 = m_svdU[i](0, 1);
			const double& U_02 = m_svdU[i](0, 2);
			const double& U_10 = m_svdU[i](1, 0);
			const double& U_11 = m_svdU[i](1, 1);
			const double& U_12 = m_svdU[i](1, 2);
			const double& U_20 = m_svdU[i](2, 0);
			const double& U_21 = m_svdU[i](2, 1);
			const double& U_22 = m_svdU[i](2, 2);

			const double& V_00 = m_svdV[i](0, 0);
			const double& V_01 = m_svdV[i](0, 1);
			const double& V_02 = m_svdV[i](0, 2);
			const double& V_10 = m_svdV[i](1, 0);
			const double& V_11 = m_svdV[i](1, 1);
			const double& V_12 = m_svdV[i](1, 2);
			const double& V_20 = m_svdV[i](2, 0);
			const double& V_21 = m_svdV[i](2, 1);
			const double& V_22 = m_svdV[i](2, 2);

			const double& A_0 = m_svdA[i](0);
			const double& A_1 = m_svdA[i](1);
			const double& A_2 = m_svdA[i](2);

			fD[i](0, 0)
				= U_00 * A_0 * V_00
				+ U_01 * A_1 * V_01
				+ U_02 * A_2 * V_02;

			fD[i](0, 1)
				= U_00 * A_0 * V_10
				+ U_01 * A_1 * V_11
				+ U_02 * A_2 * V_12;

			fD[i](0, 2)
				= U_00 * A_0 * V_20
				+ U_01 * A_1 * V_21
				+ U_02 * A_2 * V_22;

			fD[i](1, 0)
				= U_10 * A_0 * V_00
				+ U_11 * A_1 * V_01
				+ U_12 * A_2 * V_02;

			fD[i](1, 1)
				= U_10 * A_0 * V_10
				+ U_11 * A_1 * V_11
				+ U_12 * A_2 * V_12;

			fD[i](1, 2)
				= U_10 * A_0 * V_20
				+ U_11 * A_1 * V_21
				+ U_12 * A_2 * V_22;

			fD[i](2, 0)
				= U_20 * A_0 * V_00
				+ U_21 * A_1 * V_01
				+ U_22 * A_2 * V_02;

			fD[i](2, 1)
				= U_20 * A_0 * V_10
				+ U_21 * A_1 * V_11
				+ U_22 * A_2 * V_12;

			fD[i](2, 2)
				= U_20 * A_0 * V_20
				+ U_21 * A_1 * V_21
				+ U_22 * A_2 * V_22;
		//}
	}

#pragma omp parallel for
	for (int i = 0; i < m_topo.hN; ++i)
	{
		if (m_topo.h2f[i] < 0 || m_topo.h2f[i ^ 1] < 0) continue;

		hD[i].col(0) = eleRen.col(m_topo.phf2e[i][0]);
		hD[i].col(1) = eleRen.col(m_topo.phf2e[i][1]);
		hD[i].col(2) = eleRen.col(m_topo.phf2e[i][2]);

		if (m_topo.h2f[m_topo.phf2e[i][1] << 1] < 0 || m_topo.h2f[m_topo.phf2e[i][1] << 1 ^ 1] < 0 ||
			m_topo.h2f[m_topo.phf2e[i][2] << 1] < 0 || m_topo.h2f[m_topo.phf2e[i][2] << 1 ^ 1] < 0) continue;

		//m_eig[i].compute(hD[i] * hD[i].transpose());
		//m_svdA[i] = m_eig[i].eigenvectors().col(0);
		//hD[i].col(0) -= m_svdA[i].dot(hD[i].col(0)) * m_svdA[i];
		//hD[i].col(1) -= m_svdA[i].dot(hD[i].col(1)) * m_svdA[i];
		//hD[i].col(2) -= m_svdA[i].dot(hD[i].col(2)) * m_svdA[i];

		//if (fabs(hD[i].determinant()) > 1e-4)
		//{
#ifdef USE_FAST_SVD
			computeMat3SVD(hD[i].data(), m_svdU[i].data(), m_svdV[i].data(), m_svdA[i].data());
#else
			m_svd[i].compute(hD[i], Eigen::ComputeThinU | Eigen::ComputeThinV);
			m_svdU[i] = m_svd[i].matrixU();	m_svdV[i] = m_svd[i].matrixV();	m_svdA[i] = m_svd[i].singularValues();
#endif // USE_FAST_SVD

			m_svdA[i](2) = 0;
			//hD[i] = m_svdU[i] * m_svdA[i].asDiagonal() * m_svdV[i].transpose();

			const double& U_00 = m_svdU[i](0, 0);
			const double& U_01 = m_svdU[i](0, 1);
			const double& U_02 = m_svdU[i](0, 2);
			const double& U_10 = m_svdU[i](1, 0);
			const double& U_11 = m_svdU[i](1, 1);
			const double& U_12 = m_svdU[i](1, 2);
			const double& U_20 = m_svdU[i](2, 0);
			const double& U_21 = m_svdU[i](2, 1);
			const double& U_22 = m_svdU[i](2, 2);

			const double& V_00 = m_svdV[i](0, 0);
			const double& V_01 = m_svdV[i](0, 1);
			const double& V_02 = m_svdV[i](0, 2);
			const double& V_10 = m_svdV[i](1, 0);
			const double& V_11 = m_svdV[i](1, 1);
			const double& V_12 = m_svdV[i](1, 2);
			const double& V_20 = m_svdV[i](2, 0);
			const double& V_21 = m_svdV[i](2, 1);
			const double& V_22 = m_svdV[i](2, 2);

			const double& A_0 = m_svdA[i](0);
			const double& A_1 = m_svdA[i](1);
			const double& A_2 = m_svdA[i](2);

			hD[i](0, 0)
				= U_00 * A_0 * V_00
				+ U_01 * A_1 * V_01
				+ U_02 * A_2 * V_02;

			hD[i](0, 1)
				= U_00 * A_0 * V_10
				+ U_01 * A_1 * V_11
				+ U_02 * A_2 * V_12;

			hD[i](0, 2)
				= U_00 * A_0 * V_20
				+ U_01 * A_1 * V_21
				+ U_02 * A_2 * V_22;

			hD[i](1, 0)
				= U_10 * A_0 * V_00
				+ U_11 * A_1 * V_01
				+ U_12 * A_2 * V_02;

			hD[i](1, 1)
				= U_10 * A_0 * V_10
				+ U_11 * A_1 * V_11
				+ U_12 * A_2 * V_12;

			hD[i](1, 2)
				= U_10 * A_0 * V_20
				+ U_11 * A_1 * V_21
				+ U_12 * A_2 * V_22;

			hD[i](2, 0)
				= U_20 * A_0 * V_00
				+ U_21 * A_1 * V_01
				+ U_22 * A_2 * V_02;

			hD[i](2, 1)
				= U_20 * A_0 * V_10
				+ U_21 * A_1 * V_11
				+ U_22 * A_2 * V_12;

			hD[i](2, 2)
				= U_20 * A_0 * V_20
				+ U_21 * A_1 * V_21
				+ U_22 * A_2 * V_22;

		//}
	}
}

//void MeshDevelop::updateEU()
//{
//	auto fn = m_vpfnen.block(0, m_topo.vN, 3, m_topo.fN);
//	auto fn0 = m_vpfnen0.block(0, m_topo.vN, 3, m_topo.fN);
//	auto eR = (Matrix3d*)m_fDhDeReU.data() + m_topo.fN + m_topo.hN;
//	auto eU = (Matrix3d*)m_fDhDeReU.data() + m_topo.fN + m_topo.hN + m_topo.eN;
//
//	//calcFN();
//
//	std::vector<Matrix3d> EN, EN0, ERT;
//	EN.resize(m_topo.eN);
//	EN0.resize(m_topo.eN);
//	ERT.resize(m_topo.eN);
//
//#pragma omp parallel for
//	for (int i = 0; i < m_topo.eN; ++i)
//	{
//		int h0, h1;
//
//		h0 = i << 1;
//		h1 = h0 + 1;
//
//		if (m_topo.h2f[h0] < 0 || m_topo.h2f[h1] < 0) continue;
//		ERT[i] = eR[i].transpose();
//
//		EN[i].col(0) = EN0[i].col(0) = m_eb0.col(i);
//		EN0[i].col(1) = fn0.col(m_topo.h2f[h0]);
//		EN0[i].col(2) = ERT[i] * fn.col(m_topo.h2f[h1]);
//		EN0[i].col(2) = EN0[i].col(2) - EN0[i].col(2).dot(EN0[i].col(0)) * EN0[i].col(0);
//
//		EN[i].col(1) = ERT[i] * fn.col(m_topo.h2f[h0]);
//		EN[i].col(2) = fn0.col(m_topo.h2f[h1]);
//		EN[i].col(1) = EN[i].col(1) - EN[i].col(1).dot(EN[i].col(0)) * EN[i].col(0);
//
//		//if ((eU[i] * EN0[i] - EN[i]).squaredNorm() > 1e-8)
//		{
//			m_svdM[i] = EN0[i] * EN[i].transpose();
//
//#ifdef USE_FAST_SVD
//			computeMat3SVD(m_svdM[i].data(), m_svdU[i].data(), m_svdV[i].data(), m_svdA[i].data());
//#else
//			m_svd[i].compute(m_svdM[i], Eigen::ComputeThinU | Eigen::ComputeThinV);
//			m_svdU[i] = m_svd[i].matrixU();	m_svdV[i] = m_svd[i].matrixV();
//#endif // USE_FAST_SVD
// 
//			eU[i] = m_svdV[i] * m_svdU[i].transpose();
//
//			/*if (eU[i].determinant() < 0)
//			{
//				U.col(2) = -U.col(2);
//				eU[i] = V * U.transpose();
//			}*/
//		}
//	}
//}

void MeshDevelop::updateEU()
{
	auto fn = m_vpfnen.block(0, m_topo.vN, 3, m_topo.fN);
	auto fn0 = m_vpfnen0.block(0, m_topo.vN, 3, m_topo.fN);
	auto eR = (Matrix3d*)m_fDhDeReU.data() + m_topo.fN + m_topo.hN;
	auto eU = (Matrix3d*)m_fDhDeReU.data() + m_topo.fN + m_topo.hN + m_topo.eN;

	//calcFN();

#pragma omp parallel for
	for (int i = 0; i < m_topo.eN; ++i)
	{
		int h0, h1;

		h0 = i << 1;
		h1 = h0 + 1;

		if (m_topo.h2f[h0] < 0 || m_topo.h2f[h1] < 0) continue;

		//ERT[i] = eR[i].transpose();
		const Eigen::Matrix3d& eR_i = eR[i];
		const double& ERT_00 = eR_i(0, 0);
		const double& ERT_01 = eR_i(1, 0);
		const double& ERT_02 = eR_i(2, 0);
		const double& ERT_10 = eR_i(0, 1);
		const double& ERT_11 = eR_i(1, 1);
		const double& ERT_12 = eR_i(2, 1);
		const double& ERT_20 = eR_i(0, 2);
		const double& ERT_21 = eR_i(1, 2);
		const double& ERT_22 = eR_i(2, 2);

		double& EN0_00 = EN0_i_00[i];
		double& EN0_01 = EN0_i_01[i];
		double& EN0_02 = EN0_i_02[i];
		double& EN0_10 = EN0_i_10[i];
		double& EN0_11 = EN0_i_11[i];
		double& EN0_12 = EN0_i_12[i];
		double& EN0_20 = EN0_i_20[i];
		double& EN0_21 = EN0_i_21[i];
		double& EN0_22 = EN0_i_22[i];

		double& EN_00 = EN_i_00[i];
		double& EN_01 = EN_i_01[i];
		double& EN_02 = EN_i_02[i];
		double& EN_10 = EN_i_10[i];
		double& EN_11 = EN_i_11[i];
		double& EN_12 = EN_i_12[i];
		double& EN_20 = EN_i_20[i];
		double& EN_21 = EN_i_21[i];
		double& EN_22 = EN_i_22[i];

		//EN[i].col(0) = EN0[i].col(0) = m_eb0.col(i);
		const Vector3d& m_eb0_i = m_eb0.col(i);
		EN_00 = EN0_00 = m_eb0_i(0);
		EN_10 = EN0_10 = m_eb0_i(1);
		EN_20 = EN0_20 = m_eb0_i(2);

		//EN0[i].col(1) = fn0.col(m_topo.h2f[h0]);
		const Vector3d& fn0_h0 = fn0.col(m_topo.h2f[h0]);
		EN0_01 = fn0_h0(0);
		EN0_11 = fn0_h0(1);
		EN0_21 = fn0_h0(2);

		//EN0[i].col(2) = ERT[i] * fn.col(m_topo.h2f[h1]);
		const Vector3d& fn_h1 = fn.col(m_topo.h2f[h1]);
		EN0_02
			= ERT_00 * fn_h1(0)
			+ ERT_01 * fn_h1(1)
			+ ERT_02 * fn_h1(2);
		EN0_12
			= ERT_10 * fn_h1(0)
			+ ERT_11 * fn_h1(1)
			+ ERT_12 * fn_h1(2);
		EN0_22
			= ERT_20 * fn_h1(0)
			+ ERT_21 * fn_h1(1)
			+ ERT_22 * fn_h1(2);

		//EN0[i].col(2) = EN0[i].col(2) - EN0[i].col(2).dot(EN0[i].col(0)) * EN0[i].col(0);
		double EN0_2_dot_EN_0_0
			= EN0_02 * EN0_00
			+ EN0_12 * EN0_10
			+ EN0_22 * EN0_20;
		EN0_02 = EN0_02 - EN0_2_dot_EN_0_0 * EN0_00;
		EN0_12 = EN0_12 - EN0_2_dot_EN_0_0 * EN0_10;
		EN0_22 = EN0_22 - EN0_2_dot_EN_0_0 * EN0_20;

		//EN[i].col(1) = ERT[i] * fn.col(m_topo.h2f[h0]);
		const Vector3d& fn_h0 = fn.col(m_topo.h2f[h0]);
		EN_01
			= ERT_00 * fn_h0(0)
			+ ERT_01 * fn_h0(1)
			+ ERT_02 * fn_h0(2);
		EN_11
			= ERT_10 * fn_h0(0)
			+ ERT_11 * fn_h0(1)
			+ ERT_12 * fn_h0(2);
		EN_21
			= ERT_20 * fn_h0(0)
			+ ERT_21 * fn_h0(1)
			+ ERT_22 * fn_h0(2);

		//EN[i].col(2) = fn0.col(m_topo.h2f[h1]);
		const Vector3d& fn0_h1 = fn0.col(m_topo.h2f[h1]);
		EN_02 = fn0_h1(0);
		EN_12 = fn0_h1(1);
		EN_22 = fn0_h1(2);

		//EN[i].col(1) = EN[i].col(1) - EN[i].col(1).dot(EN[i].col(0)) * EN[i].col(0);
		double EN_1_dot_EN_0
			= EN_01 * EN_00
			+ EN_11 * EN_10
			+ EN_21 * EN_20;
		EN_01 = EN_01 - EN_1_dot_EN_0 * EN_00;
		EN_11 = EN_11 - EN_1_dot_EN_0 * EN_10;
		EN_21 = EN_21 - EN_1_dot_EN_0 * EN_20;

		//if ((eU[i] * EN0[i] - EN[i]).squaredNorm() > 1e-8)
		{
			//m_svdM[i] = EN0[i] * EN[i].transpose();

			// (0, 0)
			m_svdM[i](0, 0)
				= EN0_00 * EN_00
				+ EN0_01 * EN_01
				+ EN0_02 * EN_02;

			// (0, 1)
			m_svdM[i](0, 1)
				= EN0_00 * EN_10
				+ EN0_01 * EN_11
				+ EN0_02 * EN_12;

			// (0, 2)
			m_svdM[i](0, 2)
				= EN0_00 * EN_20
				+ EN0_01 * EN_21
				+ EN0_02 * EN_22;

			// (1, 0)
			m_svdM[i](1, 0)
				= EN0_10 * EN_00
				+ EN0_11 * EN_01
				+ EN0_12 * EN_02;

			// (1, 1)
			m_svdM[i](1, 1)
				= EN0_10 * EN_10
				+ EN0_11 * EN_11
				+ EN0_12 * EN_12;

			// (1, 2)
			m_svdM[i](1, 2)
				= EN0_10 * EN_20
				+ EN0_11 * EN_21
				+ EN0_12 * EN_22;

			// (2, 0)
			m_svdM[i](2, 0)
				= EN0_20 * EN_00
				+ EN0_21 * EN_01
				+ EN0_22 * EN_02;

			// (2, 1)
			m_svdM[i](2, 1)
				= EN0_20 * EN_10
				+ EN0_21 * EN_11
				+ EN0_22 * EN_12;

			// (2, 2)
			m_svdM[i](2, 2)
				= EN0_20 * EN_20
				+ EN0_21 * EN_21
				+ EN0_22 * EN_22;

#ifdef USE_FAST_SVD
			computeMat3SVD(m_svdM[i].data(), m_svdU[i].data(), m_svdV[i].data(), m_svdA[i].data());
#else
			m_svd[i].compute(m_svdM[i], Eigen::ComputeThinU | Eigen::ComputeThinV);
			m_svdU[i] = m_svd[i].matrixU();	m_svdV[i] = m_svd[i].matrixV();
#endif // USE_FAST_SVD

			//eU[i] = m_svdV[i] * m_svdU[i].transpose();
			const double& V_00 = m_svdV[i](0, 0);
			const double& V_01 = m_svdV[i](0, 1);
			const double& V_02 = m_svdV[i](0, 2);
			const double& V_10 = m_svdV[i](1, 0);
			const double& V_11 = m_svdV[i](1, 1);
			const double& V_12 = m_svdV[i](1, 2);
			const double& V_20 = m_svdV[i](2, 0);
			const double& V_21 = m_svdV[i](2, 1);
			const double& V_22 = m_svdV[i](2, 2);

			const double& U_00 = m_svdU[i](0, 0);
			const double& U_01 = m_svdU[i](0, 1);
			const double& U_02 = m_svdU[i](0, 2);
			const double& U_10 = m_svdU[i](1, 0);
			const double& U_11 = m_svdU[i](1, 1);
			const double& U_12 = m_svdU[i](1, 2);
			const double& U_20 = m_svdU[i](2, 0);
			const double& U_21 = m_svdU[i](2, 1);
			const double& U_22 = m_svdU[i](2, 2);

			// (0, 0)
			eU[i](0, 0)
				= V_00 * U_00
				+ V_01 * U_01
				+ V_02 * U_02;

			// (0, 1)
			eU[i](0, 1)
				= V_00 * U_10
				+ V_01 * U_11
				+ V_02 * U_12;

			// (0, 2)
			eU[i](0, 2)
				= V_00 * U_20
				+ V_01 * U_21
				+ V_02 * U_22;

			// (1, 0)
			eU[i](1, 0)
				= V_10 * U_00
				+ V_11 * U_01
				+ V_12 * U_02;

			// (1, 1)
			eU[i](1, 1)
				= V_10 * U_10
				+ V_11 * U_11
				+ V_12 * U_12;

			// (1, 2)
			eU[i](1, 2)
				= V_10 * U_20
				+ V_11 * U_21
				+ V_12 * U_22;

			// (2, 0)
			eU[i](2, 0)
				= V_20 * U_00
				+ V_21 * U_01
				+ V_22 * U_02;

			// (2, 1)
			eU[i](2, 1)
				= V_20 * U_10
				+ V_21 * U_11
				+ V_22 * U_12;

			// (2, 2)
			eU[i](2, 2)
				= V_20 * U_20
				+ V_21 * U_21
				+ V_22 * U_22;

			/*if (eU[i].determinant() < 0)
			{
				U.col(2) = -U.col(2);
				eU[i] = V * U.transpose();
			}*/
		}
	}
}


//void MeshDevelop::updateER()
//{
//	auto vp = m_vpfnen.block(0, 0, 3, m_topo.vN);
//	auto fn = m_vpfnen.block(0, m_topo.vN, 3, m_topo.fN);
//	auto en = m_vpfnen.block(0, m_topo.vN + m_topo.fN, 3, m_topo.eN);
//	auto fn0 = m_vpfnen0.block(0, m_topo.vN, 3, m_topo.fN);
//	auto en0 = m_vpfnen0.block(0, m_topo.vN + m_topo.fN, 3, m_topo.eN);
//	auto eR = (Matrix3d*)m_fDhDeReU.data() + m_topo.fN + m_topo.hN;
//	auto eU = (Matrix3d*)m_fDhDeReU.data() + m_topo.fN + m_topo.hN + m_topo.eN;
//
//	//calcEN();
//
//	std::vector<Eigen::Matrix<double, 3, 4>> EN, EN0;
//	std::vector<Vector3d> e0fn0, e0fn, e1fn0, e1fn;
//	std::vector<Eigen::Vector4d> diag;
//	EN.resize(m_topo.eN);
//	EN0.resize(m_topo.eN);
//	e0fn0.resize(m_topo.eN);
//	e0fn.resize(m_topo.eN);
//	e1fn0.resize(m_topo.eN);
//	e1fn.resize(m_topo.eN);
//	diag.resize(m_topo.eN);
//
//#pragma omp parallel for
//	for (int i = 0; i < m_topo.eN; ++i)
//	{
//		int h0, h1;
//		h0 = i << 1;
//		h1 = h0 + 1;
//
//		if (m_topo.h2f[h0] < 0)	e0fn0[i] = e0fn[i] = Vector3d::Zero();
//		else
//		{
//			e0fn0[i] = fn0.col(m_topo.h2f[h0]);
//			e0fn[i] = fn.col(m_topo.h2f[h0]);
//		}
//
//		if (m_topo.h2f[h1] < 0)	e1fn0[i] = e1fn[i] = Vector3d::Zero();
//		else
//		{
//			e1fn0[i] = fn0.col(m_topo.h2f[h1]);
//			e1fn[i] = fn.col(m_topo.h2f[h1]);
//		}
//
//		EN0[i].col(0) = m_es(i) * m_eb0.col(i);
//		EN0[i].col(1) = en0.col(i);
//		EN0[i].col(2) = eU[i] * e0fn0[i];
//		EN0[i].col(3) = eU[i].transpose() * e1fn0[i];
//
//		EN[i].col(0) = (vp.col(m_topo.h2v[h0]) - vp.col(m_topo.h2v[h1])) / m_el(i);
//		EN[i].col(1) = en.col(i);
//		EN[i].col(2) = e0fn[i];
//		EN[i].col(3) = e1fn[i];
//
//		diag[i](0) = m_ew(i) * pow(m_el(i), 2.0);
//		diag[i](1) = m_wD * m_ed(i); // related to eA and fA
//		diag[i](2) = diag[i](3) = m_ew(i) * pow(m_el(i), 2.0);
//
//
//		//if (((eR[i] * EN0[i] - EN[i]).cwiseAbs2() * diag[i].asDiagonal()).sum() > 1e-8 * pow(m_el(i), 2.0))
//		{
//			m_svdM[i] = EN0[i] * diag[i].asDiagonal() * EN[i].transpose();
//
//#ifdef USE_FAST_SVD
//			computeMat3SVD(m_svdM[i].data(), m_svdU[i].data(), m_svdV[i].data(), m_svdA[i].data());
//#else		
//			m_svd[i].compute(m_svdM[i], Eigen::ComputeThinU | Eigen::ComputeThinV);
//			m_svdU[i] = m_svd[i].matrixU();	m_svdV[i] = m_svd[i].matrixV();
//#endif // USE_FAST_SVD
//			
//			eR[i] = m_svdV[i] * m_svdU[i].transpose();
//
//			/*if (eR[i].determinant() < 0)
//			{
//				U.col(2) = -U.col(2);
//				eR[i] = V * U.transpose();
//			}*/
//		}
//	}
//}

void MeshDevelop::updateER()
{
	auto vp = m_vpfnen.block(0, 0, 3, m_topo.vN);
	auto fn = m_vpfnen.block(0, m_topo.vN, 3, m_topo.fN);
	auto en = m_vpfnen.block(0, m_topo.vN + m_topo.fN, 3, m_topo.eN);
	auto fn0 = m_vpfnen0.block(0, m_topo.vN, 3, m_topo.fN);
	auto en0 = m_vpfnen0.block(0, m_topo.vN + m_topo.fN, 3, m_topo.eN);
	auto eR = (Matrix3d*)m_fDhDeReU.data() + m_topo.fN + m_topo.hN;
	auto eU = (Matrix3d*)m_fDhDeReU.data() + m_topo.fN + m_topo.hN + m_topo.eN;

	//calcEN();

#pragma omp parallel for
	for (int i = 0; i < m_topo.eN; ++i)
	{
		int h0, h1;
		h0 = i << 1;
		h1 = h0 + 1;

		double& e0fn0_0 = e0fn0_i_0[i];
		double& e0fn0_1 = e0fn0_i_1[i];
		double& e0fn0_2 = e0fn0_i_2[i];
		double& e0fn_0 =  e0fn_i_0[i];
		double& e0fn_1 =  e0fn_i_1[i];
		double& e0fn_2 =  e0fn_i_2[i];
		double& e1fn0_0 =  e1fn0_i_0[i];
		double& e1fn0_1 =  e1fn0_i_1[i];
		double& e1fn0_2 =  e1fn0_i_2[i];
		double& e1fn_0 = e1fn_i_0[i];
		double& e1fn_1 = e1fn_i_1[i];
		double& e1fn_2 = e1fn_i_2[i];

		const int& h2f_h0 = m_topo.h2f[h0];
		if (h2f_h0 > 0)
		{
			const Vector3d& fn0_h0 = fn0.col(h2f_h0);
			const Vector3d& fn_h0 = fn.col(h2f_h0);
			e0fn0_0 = fn0_h0(0);
			e0fn0_1 = fn0_h0(1);
			e0fn0_2 = fn0_h0(2);
			e0fn_0 = fn_h0(0);
			e0fn_1 = fn_h0(1);
			e0fn_2 = fn_h0(2);
		}

		const int& h2f_h1 = m_topo.h2f[h1];
		if (h2f_h1 > 0)
		{
			const Vector3d& fn0_h1 = fn0.col(h2f_h1);
			const Vector3d& fn_h1 = fn.col(h2f_h1);
			e1fn0_0 = fn0_h1(0);
			e1fn0_1 = fn0_h1(1);
			e1fn0_2 = fn0_h1(2);
			e1fn_0 = fn_h1(0);
			e1fn_1 = fn_h1(1);
			e1fn_2 = fn_h1(2);
		}

		//EN0[i].col(0) = m_es(i) * m_eb0.col(i);
		//EN0[i].col(1) = en0.col(i);
		//EN0[i].col(2) = eU[i] * e0fn0[i];
		//EN0[i].col(3) = eU[i].transpose() * e1fn0[i];

		const double& m_es_i = m_es(i);

		double& EN0_00 = EN0_i_00[i];
		double& EN0_01 = EN0_i_01[i];
		double& EN0_02 = EN0_i_02[i];
		double& EN0_03 = EN0_i_03[i];
		double& EN0_10 = EN0_i_10[i];
		double& EN0_11 = EN0_i_11[i];
		double& EN0_12 = EN0_i_12[i];
		double& EN0_13 = EN0_i_13[i];
		double& EN0_20 = EN0_i_20[i];
		double& EN0_21 = EN0_i_21[i];
		double& EN0_22 = EN0_i_22[i];
		double& EN0_23 = EN0_i_23[i];

		const Vector3d& m_eb0_i = m_eb0.col(i);
		EN0_00 = m_es_i * m_eb0_i(0);
		EN0_10 = m_es_i * m_eb0_i(1);
		EN0_20 = m_es_i * m_eb0_i(2);

		const Vector3d& en0_i = en0.col(i);
		EN0_01 = en0_i(0);
		EN0_11 = en0_i(1);
		EN0_21 = en0_i(2);

		const Eigen::Matrix3d& eU_i = eU[i];
		const double& eU_00 = eU_i(0, 0);
		const double& eU_01 = eU_i(0, 1);
		const double& eU_02 = eU_i(0, 2);
		const double& eU_10 = eU_i(1, 0);
		const double& eU_11 = eU_i(1, 1);
		const double& eU_12 = eU_i(1, 2);
		const double& eU_20 = eU_i(2, 0);
		const double& eU_21 = eU_i(2, 1);
		const double& eU_22 = eU_i(2, 2);

		EN0_02
			= eU_00 * e0fn0_0
			+ eU_01 * e0fn0_1
			+ eU_02 * e0fn0_2;
		EN0_12
			= eU_10 * e0fn0_0
			+ eU_11 * e0fn0_1
			+ eU_12 * e0fn0_2;
		EN0_22
			= eU_20 * e0fn0_0
			+ eU_21 * e0fn0_1
			+ eU_22 * e0fn0_2;

		EN0_03
			= eU_00 * e1fn0_0
			+ eU_10 * e1fn0_1
			+ eU_20 * e1fn0_2;
		EN0_13
			= eU_01 * e1fn0_0
			+ eU_11 * e1fn0_1
			+ eU_21 * e1fn0_2;
		EN0_23
			= eU_02 * e1fn0_0
			+ eU_12 * e1fn0_1
			+ eU_22 * e1fn0_2;

		double& EN_00 = EN_i_00[i];
		double& EN_01 = EN_i_01[i];
		double& EN_02 = EN_i_02[i];
		double& EN_03 = EN_i_03[i];
		double& EN_10 = EN_i_10[i];
		double& EN_11 = EN_i_11[i];
		double& EN_12 = EN_i_12[i];
		double& EN_13 = EN_i_13[i];
		double& EN_20 = EN_i_20[i];
		double& EN_21 = EN_i_21[i];
		double& EN_22 = EN_i_22[i];
		double& EN_23 = EN_i_23[i];

		const Vector3d& vp_h0 = vp.col(m_topo.h2v[h0]);
		const Vector3d& vp_h1 = vp.col(m_topo.h2v[h1]);
		const double& el_i = m_el(i);

		const Vector3d& en_i = en.col(i);
		EN_00 = (vp_h0(0) - vp_h1(0)) / el_i;
		EN_10 = (vp_h0(1) - vp_h1(1)) / el_i;
		EN_20 = (vp_h0(2) - vp_h1(2)) / el_i;
		EN_01 = en_i(0);
		EN_11 = en_i(1);
		EN_21 = en_i(2);
		EN_02 = e0fn_0;
		EN_12 = e0fn_1;
		EN_22 = e0fn_2;
		EN_03 = e1fn_0;
		EN_13 = e1fn_1;
		EN_23 = e1fn_2;

		double& diag_0 = diag_i_0[i];
		double& diag_1 = diag_i_1[i];
		double& diag_2 = diag_i_2[i];
		double& diag_3 = diag_i_3[i];

		diag_0 = m_ew(i) * el_i * el_i;
		diag_1 = m_wD * m_ed(i); // related to eA and fA
		diag_2 = diag_3 = diag_0;


		//if (((eR[i] * EN0[i] - EN[i]).cwiseAbs2() * diag[i].asDiagonal()).sum() > 1e-8 * pow(m_el(i), 2.0))
		{
			//m_svdM[i] = EN0[i] * diag[i].asDiagonal() * EN[i].transpose();

			// (0, 0)
			m_svdM[i](0, 0)
				= EN0_00 * diag_0 * EN_00
				+ EN0_01 * diag_1 * EN_01
				+ EN0_02 * diag_2 * EN_02
				+ EN0_03 * diag_3 * EN_03;

			// (0, 1)
			m_svdM[i](0, 1)
				= EN0_00 * diag_0 * EN_10
				+ EN0_01 * diag_1 * EN_11
				+ EN0_02 * diag_2 * EN_12
				+ EN0_03 * diag_3 * EN_13;

			// (0, 2)
			m_svdM[i](0, 2)
				= EN0_00 * diag_0 * EN_20
				+ EN0_01 * diag_1 * EN_21
				+ EN0_02 * diag_2 * EN_22
				+ EN0_03 * diag_3 * EN_23;

			// (1, 0)
			m_svdM[i](1, 0)
				= EN0_10 * diag_0 * EN_00
				+ EN0_11 * diag_1 * EN_01
				+ EN0_12 * diag_2 * EN_02
				+ EN0_13 * diag_3 * EN_03;

			// (1, 1)
			m_svdM[i](1, 1)
				= EN0_10 * diag_0 * EN_10
				+ EN0_11 * diag_1 * EN_11
				+ EN0_12 * diag_2 * EN_12
				+ EN0_13 * diag_3 * EN_13;

			// (1, 2)
			m_svdM[i](1, 2)
				= EN0_10 * diag_0 * EN_20
				+ EN0_11 * diag_1 * EN_21
				+ EN0_12 * diag_2 * EN_22
				+ EN0_13 * diag_3 * EN_23;

			// (2, 0)
			m_svdM[i](2, 0)
				= EN0_20 * diag_0 * EN_00
				+ EN0_21 * diag_1 * EN_01
				+ EN0_22 * diag_2 * EN_02
				+ EN0_23 * diag_3 * EN_03;

			// (2, 1)
			m_svdM[i](2, 1)
				= EN0_20 * diag_0 * EN_10
				+ EN0_21 * diag_1 * EN_11
				+ EN0_22 * diag_2 * EN_12
				+ EN0_23 * diag_3 * EN_13;

			// (2, 2)
			m_svdM[i](2, 2)
				= EN0_20 * diag_0 * EN_20
				+ EN0_21 * diag_1 * EN_21
				+ EN0_22 * diag_2 * EN_22
				+ EN0_23 * diag_3 * EN_23;


#ifdef USE_FAST_SVD
			computeMat3SVD(m_svdM[i].data(), m_svdU[i].data(), m_svdV[i].data(), m_svdA[i].data());
#else		
			m_svd[i].compute(m_svdM[i], Eigen::ComputeThinU | Eigen::ComputeThinV);
			m_svdU[i] = m_svd[i].matrixU();	m_svdV[i] = m_svd[i].matrixV();
#endif // USE_FAST_SVD

			//eR[i] = m_svdV[i] * m_svdU[i].transpose();
			const double& V_00 = m_svdV[i](0, 0);
			const double& V_01 = m_svdV[i](0, 1);
			const double& V_02 = m_svdV[i](0, 2);
			const double& V_10 = m_svdV[i](1, 0);
			const double& V_11 = m_svdV[i](1, 1);
			const double& V_12 = m_svdV[i](1, 2);
			const double& V_20 = m_svdV[i](2, 0);
			const double& V_21 = m_svdV[i](2, 1);
			const double& V_22 = m_svdV[i](2, 2);

			const double& U_00 = m_svdU[i](0, 0);
			const double& U_01 = m_svdU[i](0, 1);
			const double& U_02 = m_svdU[i](0, 2);
			const double& U_10 = m_svdU[i](1, 0);
			const double& U_11 = m_svdU[i](1, 1);
			const double& U_12 = m_svdU[i](1, 2);
			const double& U_20 = m_svdU[i](2, 0);
			const double& U_21 = m_svdU[i](2, 1);
			const double& U_22 = m_svdU[i](2, 2);

			// (0, 0)
			eR[i](0, 0)
				= V_00 * U_00
				+ V_01 * U_01
				+ V_02 * U_02;

			// (0, 1)
			eR[i](0, 1)
				= V_00 * U_10
				+ V_01 * U_11
				+ V_02 * U_12;

			// (0, 2)
			eR[i](0, 2)
				= V_00 * U_20
				+ V_01 * U_21
				+ V_02 * U_22;

			// (1, 0)
			eR[i](1, 0)
				= V_10 * U_00
				+ V_11 * U_01
				+ V_12 * U_02;

			// (1, 1)
			eR[i](1, 1)
				= V_10 * U_10
				+ V_11 * U_11
				+ V_12 * U_12;

			// (1, 2)
			eR[i](1, 2)
				= V_10 * U_20
				+ V_11 * U_21
				+ V_12 * U_22;

			// (2, 0)
			eR[i](2, 0)
				= V_20 * U_00
				+ V_21 * U_01
				+ V_22 * U_02;

			// (2, 1)
			eR[i](2, 1)
				= V_20 * U_10
				+ V_21 * U_11
				+ V_22 * U_12;

			// (2, 2)
			eR[i](2, 2)
				= V_20 * U_20
				+ V_21 * U_21
				+ V_22 * U_22;

			/*if (eR[i].determinant() < 0)
			{
				U.col(2) = -U.col(2);
				eR[i] = V * U.transpose();
			}*/
		}
	}
}

void MeshDevelop::calcEN()
{
	auto en = m_vpfnen.block(0, m_topo.vN + m_topo.fN, 3, m_topo.eN);
	auto en0 = m_vpfnen0.block(0, m_topo.vN + m_topo.fN, 3, m_topo.eN);
	auto fD = (Matrix3d*)m_fDhDeReU.data();
	auto hD = (Matrix3d*)m_fDhDeReU.data() + m_topo.fN;
	auto eR = (Matrix3d*)m_fDhDeReU.data() + m_topo.fN + m_topo.hN;

	en.setZero();

	int h0, h1;

	for (int i = 0; i < m_topo.eN; ++i)
	{
		h0 = i << 1;
		h1 = h0 + 1;

		if (m_topo.h2f[h0] < 0 || m_topo.h2f[h1] < 0) continue;

		en.col(i) += (hD[h0].col(0) + hD[h1].col(0)) * m_eA(i);
		en.col(m_topo.phf2e[h0][1]) += hD[h0].col(1) * m_eA(i);
		en.col(m_topo.phf2e[h0][2]) += hD[h0].col(2) * m_eA(i);
		en.col(m_topo.phf2e[h1][1]) += hD[h1].col(1) * m_eA(i);
		en.col(m_topo.phf2e[h1][2]) += hD[h1].col(2) * m_eA(i);
	}

	for (int i = 0; i < m_topo.fN; ++i)
	{
		en.col(m_topo.f2h[i][0] >> 1) += fD[i].col(0) * m_fA(i);
		en.col(m_topo.f2h[i][1] >> 1) += fD[i].col(1) * m_fA(i);
		en.col(m_topo.f2h[i][2] >> 1) += fD[i].col(2) * m_fA(i);
	}

	for (int i = 0; i < m_topo.eN; ++i)
	{
		en.col(i) /= m_ed(i);
	}
}

void MeshDevelop::calcFN()
{
	VectorXd fw = VectorXd::Zero(m_topo.fN);
	auto fn = m_vpfnen.block(0, m_topo.vN, 3, m_topo.fN);
	auto fn0 = m_vpfnen0.block(0, m_topo.vN, 3, m_topo.fN);
	auto eR = (Matrix3d*)m_fDhDeReU.data() + m_topo.fN + m_topo.hN;
	auto eU = (Matrix3d*)m_fDhDeReU.data() + m_topo.fN + m_topo.hN + m_topo.eN;

	fn.setZero();
	int h0, h1;

	for (int i = 0; i < m_topo.eN; ++i)
	{
		h0 = i << 1;
		h1 = h0 + 1;

		if (m_topo.h2f[h0] >= 0)
		{
			fn.col(m_topo.h2f[h0]) += eR[i] * (eU[i] * fn0.col(m_topo.h2f[h0])) * pow(m_el(i), 2) * m_ew(i);
			fw(m_topo.h2f[h0]) += pow(m_el(i), 2) * m_ew(i);
		}

		if (m_topo.h2f[h1] >= 0)
		{
			fn.col(m_topo.h2f[h1]) += eR[i] * (eU[i].transpose() * fn0.col(m_topo.h2f[h1])) * pow(m_el(i), 2) * m_ew(i);
			fw(m_topo.h2f[h1]) += pow(m_el(i), 2) * m_ew(i);
		}
	}

	for (int i = 0; i < m_topo.fN; ++i)
	{
		fn.col(i) /= fw(i);
	}
}

void MeshDevelop::updateVP()
{
	auto vp = m_vpfnen.block(0, 0, 3, m_topo.vN);
	auto eR = (Matrix3d*)m_fDhDeReU.data() + m_topo.fN + m_topo.hN;
	//vp = m_wP * m_vt;
	vp = m_wP * (m_vA.asDiagonal() * m_vt.transpose()).transpose();

	Vector3d dp;

	int h0, h1;

	for (int i = 0; i < m_topo.eN; ++i)
	{
		h0 = i << 1;
		h1 = h0 + 1;

		dp = eR[i] * m_eb0.col(i) * (m_el(i) * m_es(i) * m_ew(i));
		vp.col(m_topo.h2v[h0]) += dp;
		vp.col(m_topo.h2v[h1]) -= dp;
	}

	vp = m_ldlt.solve(vp.transpose()).transpose();
}

void MeshDevelop::updateVT()
{
	auto vp = m_vpfnen.block(0, 0, 3, m_topo.vN);
	auto vp0 = m_vpfnen0.block(0, 0, 3, m_topo.vN);

	m_vt = vp0;

	for (int i = 0; i < m_topo.vN; ++i)
	{
		if (m_fixV[i]) continue;
		m_vt.col(i) = vp0.col(i) + (vp.col(i) - vp0.col(i)) / fmax(1, (vp.col(i) - vp0.col(i)).norm() / m_dist);
	}
}

void MeshDevelop::updateES()
{
	int h0, h1;

	auto vp = m_vpfnen.block(0, 0, 3, m_topo.vN);
	auto vp0 = m_vpfnen0.block(0, 0, 3, m_topo.vN);
	auto eR = (Matrix3d*)m_fDhDeReU.data() + m_topo.fN + m_topo.hN;

	for (int i = 0; i < m_topo.eN; ++i)
	{
		h0 = i << 1;
		h1 = h0 + 1;

		m_es(i) = (vp.col(m_topo.h2v[h0]) - vp.col(m_topo.h2v[h1])).dot(eR[i] * m_eb0.col(i));
		//if (m_es(i) < triNeq.m_eps) m_es(i) = triNeq.m_eps;
		//m_es(i) = (vp.col(m_topo.h2v[h0]) - vp.col(m_topo.h2v[h1])).dot(eR[i] * m_eb0.col(i)) * m_el(i) * m_ew(i);

		//m_es(i) = (vp.col(m_topo.h2v[h0]) - vp.col(m_topo.h2v[h1])).dot(eR[i] * m_eb0.col(i)) / m_el(i);
		//if (m_es(i) < 0.75) m_es(i) = 0.75;
		//if (m_es(i) > 1.25) m_es(i) = 1.25;
	}

	m_es = triNeq.opt(m_es).cwiseQuotient(m_el);

	/*m_es = m_eldlt.solve(m_es);
	for (int i = 0; i < m_topo.eN; ++i)
	{
		if (m_es(i) < 0.8) m_es(i) = 0.8;
		if (m_es(i) > 1.2) m_es(i) = 1.2;
	}*/
}

double MeshDevelop::energy()
{
	return m_wD * energyD() + energyDP() + energyFN() + m_wP * energyRP();
}

double MeshDevelop::energyD()
{
	dev_energy.assign(m_topo.eN, 0);

	double rt = 0;
	Matrix3d ens0, ens1;
	auto en0 = m_vpfnen0.block(0, m_topo.vN + m_topo.fN, 3, m_topo.eN);
	auto fD = (Matrix3d*)m_fDhDeReU.data();
	auto hD = (Matrix3d*)m_fDhDeReU.data() + m_topo.fN;
	auto eR = (Matrix3d*)m_fDhDeReU.data() + m_topo.fN + m_topo.hN;
	MatrixXd eRen = en0;

	for (int i = 0; i < m_topo.eN; ++i)
	{
		eRen.col(i) = eR[i] * en0.col(i);
	}

	int h0, h1;
	for (int i = 0; i < m_topo.eN; ++i)
	{
		h0 = i << 1;
		h1 = h0 + 1;

		if (m_topo.h2f[h0] < 0 || m_topo.h2f[h1] < 0) continue;

		ens0.col(0) = ens1.col(0) = eRen.col(i);
		ens0.col(1) = eRen.col(m_topo.phf2e[h0][1]);
		ens0.col(2) = eRen.col(m_topo.phf2e[h0][2]);
		ens1.col(1) = eRen.col(m_topo.phf2e[h1][1]);
		ens1.col(2) = eRen.col(m_topo.phf2e[h1][2]);

		double temp_value = ((ens0 - hD[h0]).squaredNorm() + (ens1 - hD[h1]).squaredNorm());
		rt += temp_value;// *m_eA(i);

		dev_energy[i] += temp_value;
	}

	for (int i = 0; i < m_topo.fN; ++i)
	{
		ens0.col(0) = eRen.col(m_topo.f2h[i][0] >> 1);
		ens0.col(1) = eRen.col(m_topo.f2h[i][1] >> 1);
		ens0.col(2) = eRen.col(m_topo.f2h[i][2] >> 1);


		double temp_value = (ens0 - fD[i]).squaredNorm();
		rt += temp_value;// *m_fA[i];

		dev_energy[m_topo.f2h[i][0] >> 1] += temp_value;
		dev_energy[m_topo.f2h[i][1] >> 1] += temp_value;
		dev_energy[m_topo.f2h[i][2] >> 1] += temp_value;
	}

	return 0.5 * rt;
}

double MeshDevelop::energyFN()
{
	double rt = 0;
	auto fn = m_vpfnen.block(0, m_topo.vN, 3, m_topo.fN);
	auto fn0 = m_vpfnen0.block(0, m_topo.vN, 3, m_topo.fN);
	auto eR = (Matrix3d*)m_fDhDeReU.data() + m_topo.fN + m_topo.hN;
	auto eU = (Matrix3d*)m_fDhDeReU.data() + m_topo.fN + m_topo.hN + m_topo.eN;

	int h0, h1;

	for (int i = 0; i < m_topo.eN; ++i)
	{
		h0 = i << 1;
		h1 = h0 + 1;

		if (m_topo.h2f[h0] >= 0)
		{
			rt += pow(m_el(i), 2) * m_ew(i) *
				(eR[i] * (eU[i] * fn0.col(m_topo.h2f[h0])) - fn.col(m_topo.h2f[h0])).squaredNorm();
		}

		if (m_topo.h2f[h1] >= 0)
		{
			rt += pow(m_el(i), 2) * m_ew(i) *
				(eR[i] * (eU[i].transpose() * fn0.col(m_topo.h2f[h1])) - fn.col(m_topo.h2f[h1])).squaredNorm();
		}
	}

	return 0.5 * rt;
}

double MeshDevelop::energyDP()
{
	double rt = 0;
	auto vp = m_vpfnen.block(0, 0, 3, m_topo.vN);
	auto eR = (Matrix3d*)m_fDhDeReU.data() + m_topo.fN + m_topo.hN;

	for (int i = 0; i < m_topo.fN; ++i)
	{
		int e0 = m_topo.f2h[i][0] >> 1;
		int e1 = m_topo.f2h[i][1] >> 1;
		int e2 = m_topo.f2h[i][2] >> 1;

		rt += m_fA(i) * pow(m_es(e0) - m_es(e1), 2.0);
		rt += m_fA(i) * pow(m_es(e1) - m_es(e2), 2.0);
		rt += m_fA(i) * pow(m_es(e2) - m_es(e0), 2.0);
	}
	rt *= m_similar;

	int h0, h1;
	for (int i = 0; i < m_topo.eN; ++i)
	{
		h0 = i << 1;
		h1 = h0 + 1;

		rt += pow(m_el(i), 2) * m_ew(i) * (eR[i] * m_eb0.col(i) * m_es(i) -
			(vp.col(m_topo.h2v[h0]) - vp.col(m_topo.h2v[h1])) / m_el(i)).squaredNorm();
	}

	return 0.5 * rt;
}

double MeshDevelop::energyRP()
{
	//return 0.5 * (m_vpfnen.block(0, 0, 3, m_topo.vN) - m_vt).squaredNorm();
	return 0.5 * (m_vA.asDiagonal() * (m_vpfnen.block(0, 0, 3, m_topo.vN) -
		m_vt).cwiseAbs2().transpose()).sum();
}

bool MeshDevelop::inBound()
{
	auto vp = m_vpfnen.block(0, 0, 3, m_topo.vN);
	auto vp0 = m_vpfnen0.block(0, 0, 3, m_topo.vN);

	double maxDist = 0;
	for (int i = 0; i < m_topo.vN; ++i)
	{
		double colNorm = (vp.col(i) - vp0.col(i)).norm();
		if (colNorm > maxDist)
		{
			maxDist = colNorm;
		}
	}

	return maxDist <= m_dist;
}