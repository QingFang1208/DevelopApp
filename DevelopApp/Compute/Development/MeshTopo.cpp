#include "MeshTopo.h"
using OpenMesh::Vec3i;

MeshTopo::MeshTopo()
{

}

MeshTopo::MeshTopo(const Mesh& mesh)
{
	init(mesh);
}

void MeshTopo::init(const Mesh& mesh)
{
	vN = mesh.n_vertices();		fN = mesh.n_faces();	eN = mesh.n_edges();	hN = mesh.n_halfedges();

	vb.assign(vN, false);			h2f.assign(hN, -1);			h2v.assign(hN, -1);
	f2v.assign(fN, Vec3i(-1));		f2f.assign(fN, Vec3i(-1));	f2h.assign(fN, Vec3i(-1));

	for (auto hh : mesh.halfedges())
	{
		h2f[hh.idx()] = mesh.face_handle(hh).idx();
		h2v[hh.idx()] = mesh.to_vertex_handle(hh).idx();

		vb[h2v[hh.idx()]] = vb[h2v[hh.idx()]] || (h2f[hh.idx()] < 0);
	}

	for (auto fh : mesh.faces())
	{
		auto fh_it(mesh.cfh_iter(fh));
		f2h[fh.idx()][0] = fh_it->idx();	f2v[fh.idx()][0] = mesh.to_vertex_handle(*fh_it).idx();		++fh_it;
		f2h[fh.idx()][1] = fh_it->idx();	f2v[fh.idx()][1] = mesh.to_vertex_handle(*fh_it).idx();		++fh_it;
		f2h[fh.idx()][2] = fh_it->idx();	f2v[fh.idx()][2] = mesh.to_vertex_handle(*fh_it).idx();

		f2f[fh.idx()][0] = h2f[f2h[fh.idx()][0] ^ 1];
		f2f[fh.idx()][1] = h2f[f2h[fh.idx()][1] ^ 1];
		f2f[fh.idx()][2] = h2f[f2h[fh.idx()][2] ^ 1];
	}

	phfToE();
}

void MeshTopo::clear()
{
	vN = 0;		fN = 0;		hN = 0;

	vb.clear();		h2f.clear();	h2v.clear();
	f2v.clear();	f2f.clear();	f2h.clear();
}

void MeshTopo::phfToE()
{
	phf2e.clear();
	phf2e.assign(h2f.size(), Vec3i(-1));

	for (auto fh : f2h)
	{
		if (h2f[fh[0] ^ 1] >= 0)
		{
			phf2e[fh[0]][0] = fh[0] >> 1;
			phf2e[fh[0]][1] = fh[2] >> 1;
			phf2e[fh[0] ^ 1][2] = fh[1] >> 1;
		}

		if (h2f[fh[1] ^ 1] >= 0)
		{
			phf2e[fh[1]][0] = fh[1] >> 1;
			phf2e[fh[1]][1] = fh[0] >> 1;
			phf2e[fh[1] ^ 1][2] = fh[2] >> 1;
		}

		if (h2f[fh[2] ^ 1] >= 0)
		{
			phf2e[fh[2]][0] = fh[2] >> 1;
			phf2e[fh[2]][1] = fh[1] >> 1;
			phf2e[fh[2] ^ 1][2] = fh[0] >> 1;
		}
	}
}