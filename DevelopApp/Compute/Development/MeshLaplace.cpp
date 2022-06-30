#include "MeshLaplace.h"
using OpenMesh::Vec3d;

void faceArea(const OpenMesh::Vec3d* points, const MeshTopo& topo, std::vector<double>& fA)
{
	fA.clear();
	fA.assign(topo.f2v.size(), 0);

	Vec3d e0, e1, e2;
	for (int i = 0; i < fA.size(); ++i)
	{
		e0 = points[topo.f2v[i][1]] - points[topo.f2v[i][0]];
		e1 = points[topo.f2v[i][2]] - points[topo.f2v[i][1]];

		fA[i] = (e0 % e1).norm() / 2;
	}
}

void faceArea(const Mesh& mesh, std::vector<double>& fA)
{
	MeshTopo topo(mesh);
	faceArea(mesh.points(), topo, fA);
}

void vertArea(const OpenMesh::Vec3d* points, const MeshTopo& topo, std::vector<double>& vA)
{
	vA.clear();
	vA.assign(topo.vN, 0);

	double fA;
	Vec3d e0, e1;
	for (int i = 0; i < topo.f2v.size(); ++i)
	{
		e0 = points[topo.f2v[i][1]] - points[topo.f2v[i][0]];
		e1 = points[topo.f2v[i][2]] - points[topo.f2v[i][1]];

		fA = (e0 % e1).norm() / 6;

		vA[topo.f2v[i][0]] += fA;
		vA[topo.f2v[i][1]] += fA;
		vA[topo.f2v[i][2]] += fA;
	}
}

void vertArea(const Mesh& mesh, std::vector<double>& vA)
{
	MeshTopo topo(mesh);
	vertArea(mesh.points(), topo, vA);
}

void vertNormal(const OpenMesh::Vec3d* points, const MeshTopo& topo, std::vector<OpenMesh::Vec3d>& vN)
{
	vN.clear();
	vN.assign(topo.vN, Vec3d(0));

	OpenMesh::Vec3d fn;
	Vec3d e0, e1;
	for (int i = 0; i < topo.f2v.size(); ++i)
	{
		e0 = points[topo.f2v[i][1]] - points[topo.f2v[i][0]];
		e1 = points[topo.f2v[i][2]] - points[topo.f2v[i][1]];

		fn = (e0 % e1);
		vN[topo.f2v[i][0]] += fn;
		vN[topo.f2v[i][1]] += fn;
		vN[topo.f2v[i][2]] += fn;
	}

	for (int i = 0; i < topo.vN; ++i)
	{
		vN[i].normalize();
	}
}

void vertNormal(const Mesh& mesh, std::vector<OpenMesh::Vec3d>& vN)
{
	MeshTopo topo(mesh);
	vertNormal(mesh.points(), topo, vN);
}

void vertGauss(const OpenMesh::Vec3d* points, const MeshTopo& topo, std::vector<double>& vK)
{
	vK.clear();
	vK.assign(topo.vN, 2 * M_PI);
	
	for (int i = 0; i < vK.size(); ++i)
	{
		if (topo.vb[i]) vK[i] = M_PI;
	}

	std::vector<Vec3d> f2a;
	faceAngle(points, topo, f2a);

	for (int i = 0; i < f2a.size(); ++i)
	{
		vK[topo.f2v[i][0]] -= f2a[i][0];
		vK[topo.f2v[i][1]] -= f2a[i][1];
		vK[topo.f2v[i][2]] -= f2a[i][2];
	}
}

void vertGauss(const Mesh& mesh, std::vector<double>& vK)
{
	MeshTopo topo(mesh);
	vertGauss(mesh.points(), topo, vK);
}

void faceAngle(const Vec3d* points, const MeshTopo& topo, std::vector<Vec3d>& f2a)
{
	f2a.clear();
	f2a.assign(topo.f2v.size(), Vec3d(0));

	Vec3d e0, e1, e2;
	for (int i = 0; i < f2a.size(); ++i)
	{
		e0 = points[topo.f2v[i][1]] - points[topo.f2v[i][0]];		e0.normalize();
		e1 = points[topo.f2v[i][2]] - points[topo.f2v[i][1]];		e1.normalize();
		e2 = points[topo.f2v[i][0]] - points[topo.f2v[i][2]];		e2.normalize();

		f2a[i][0] = acos(-e2 | e0);
		f2a[i][1] = acos(-e0 | e1);
		f2a[i][2] = acos(-e1 | e2);
	}
}

void faceAngle(const Mesh& mesh, std::vector<OpenMesh::Vec3d>& f2a)
{
	MeshTopo topo(mesh);
	faceAngle(mesh.points(), topo, f2a);
}

void laplaceTriplet(const OpenMesh::Vec3d* points, const MeshTopo& topo, std::vector<Eigen::Triplet<double>>& trips)
{
	std::vector<Vec3d> f2a, f2cot;
	faceAngle(points, topo, f2a);

	f2cot.assign(f2a.size(), Vec3d(0));
	for (int i = 0; i < f2a.size(); ++i)
	{
		f2cot[i][0] = 0.5 / tan(f2a[i][0]);
		f2cot[i][1] = 0.5 / tan(f2a[i][1]);
		f2cot[i][2] = 0.5 / tan(f2a[i][2]);
	}

	trips.clear();
	for (int i = 0; i < f2a.size(); ++i)
	{
		trips.emplace_back(topo.f2v[i][0], topo.f2v[i][0], f2cot[i][1] + f2cot[i][2]);
		trips.emplace_back(topo.f2v[i][1], topo.f2v[i][1], f2cot[i][0] + f2cot[i][2]);
		trips.emplace_back(topo.f2v[i][2], topo.f2v[i][2], f2cot[i][0] + f2cot[i][1]);
		trips.emplace_back(topo.f2v[i][0], topo.f2v[i][1], -f2cot[i][2]);
		trips.emplace_back(topo.f2v[i][1], topo.f2v[i][0], -f2cot[i][2]);
		trips.emplace_back(topo.f2v[i][1], topo.f2v[i][2], -f2cot[i][0]);
		trips.emplace_back(topo.f2v[i][2], topo.f2v[i][1], -f2cot[i][0]);
		trips.emplace_back(topo.f2v[i][2], topo.f2v[i][0], -f2cot[i][1]);
		trips.emplace_back(topo.f2v[i][0], topo.f2v[i][2], -f2cot[i][1]);
	}
}

void laplaceTriplet(const Mesh& mesh, std::vector<Eigen::Triplet<double>>& trips)
{
	MeshTopo topo(mesh);
	laplaceTriplet(mesh.points(), topo, trips);
}