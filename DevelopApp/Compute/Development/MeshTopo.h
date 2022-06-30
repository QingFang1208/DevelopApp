#pragma once
#include "../MeshDefinition.h"

struct MeshTopo
{
	int vN = 0, fN = 0, eN=0, hN = 0;
	std::vector<bool> vb;
	std::vector<int> h2f, h2v;
	std::vector<OpenMesh::Vec3i> f2v, f2f, f2h, phf2e;

	MeshTopo();
	MeshTopo(const Mesh& mesh);
	void init(const Mesh& mesh);
	void clear();

	void phfToE();
};

void phfToE(const MeshTopo& topo, std::vector<OpenMesh::Vec3i>& phf2e);
