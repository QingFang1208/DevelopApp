#pragma once
#include <queue>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#ifdef _DEBUG
#pragma comment(lib, "OpenMeshCored.lib")
#pragma comment(lib, "OpenMeshToolsd.lib")
#else
#pragma comment(lib, "OpenMeshCore.lib")
#pragma comment(lib, "OpenMeshTools.lib")
#endif
struct MeshTraits : public OpenMesh::DefaultTraits
{
	typedef OpenMesh::Vec3d Point;
	typedef OpenMesh::Vec3d Normal;
	typedef OpenMesh::Vec3d TexCoord3D;
	VertexAttributes(OpenMesh::Attributes::Status | OpenMesh::Attributes::Normal);
	FaceAttributes(OpenMesh::Attributes::Status | OpenMesh::Attributes::Normal);
	EdgeAttributes(OpenMesh::Attributes::Status);
	HalfedgeAttributes(OpenMesh::Attributes::Status);

	typedef typename std::pair<Point, Point> Link;               // for local links
	typedef typename std::vector<Link> Link_vector;                        // for out links
	typedef typename std::deque<Link> Link_deque;            // for in links				

	FaceTraits
	{
	public:
		int m_tag = 0;                        // general-purpose tag
		// links
		Link_vector m_face_out_links;      // point links from this facet
		Link_deque m_face_in_links;  // point links to this facet
		Link_deque m_edge_in_links;
		Link_deque m_vertex_in_links;
	};

	EdgeTraits
	{
	public:
		// links
		Link_vector m_edge_out_links;         // point link from this edge
	};

	HalfedgeTraits
	{
		HalfedgeT() :face_he_var(-1)
		{
		};
	private:
		int face_he_var;
	public:
		int get_face_he_var()const { return face_he_var; };
		void set_face_he_var(int fhe) { face_he_var = fhe; };
	};

	VertexTraits
	{
		VertexT() : new_pos_fixed(false)
		{
		};
	private:
		OpenMesh::Vec3d new_pos;//can be used for deformation and parameterization
		bool new_pos_fixed;
	public:
		void set_New_Pos(const OpenMesh::Vec3d& n_p) { new_pos = n_p; }
		const OpenMesh::Vec3d& get_New_Pos()const { return new_pos; }
		void set_new_pos_fixed(bool f) { new_pos_fixed = f; }
		bool get_new_pos_fixed()const { return new_pos_fixed; }

		// links
		Link m_vertex_out_link;                 // Links from this vertex
	};
};
typedef OpenMesh::TriMesh_ArrayKernelT<MeshTraits> Mesh;

typedef OpenMesh::EdgeHandle EH;
typedef OpenMesh::FaceHandle FH;
typedef OpenMesh::HalfedgeHandle HEH;
typedef OpenMesh::VertexHandle VH;

class MeshTools
{
public:
	static bool ReadMesh(Mesh & mesh, const std::string & filename);
	static bool ReadOBJ(Mesh & mesh, const std::string & filename);
	//static bool ReadOFF(Mesh & mesh, const std::string & filename);
	static bool WriteMesh(const Mesh & mesh, const std::string & filename, const std::streamsize & precision = 6);
	static bool WriteOBJ(const Mesh & mesh, const std::string & filename, const std::streamsize & precision = 6);
	static double Area(const Mesh & mesh);
	static double AverageEdgeLength(const Mesh & mesh);
	static bool HasBoundary(const Mesh & mesh);
	static bool HasOneComponent(const Mesh & mesh);
	static int Genus(const Mesh & mesh);
	static void BoundingBox(const Mesh & mesh, Mesh::Point & bmax, Mesh::Point & bmin);
	static void Reassign(const Mesh & mesh1, Mesh & mesh2);
	static void ComputeGaussianCurvature(const Mesh& mesh, std::vector<double>& v_gauss, const bool& bound_remove = true);
};
