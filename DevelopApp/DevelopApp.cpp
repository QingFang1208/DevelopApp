#include "Compute\Optimization.h"

int main(int argc, char* argv[])
{
	// Input
	std::string input_name = "model/Bunny.obj";
	//std::string input_name = "model/Tortoise.obj";
	//std::cout << "Read mesh: " << filename << std::endl;
	Mesh input_mesh;
	MeshTools::ReadMesh(input_mesh, input_name);

	// ============ Deformation ==============
	SegMesh seg_mesh(input_mesh);
	std::cout << "opt " << std::endl;
	Optimization opt(seg_mesh);
	opt.Run();

	//// ============ Deformation ==============
	//std::string filename = "model/input.obj";
	//std::cout << "Read mesh: " << filename << std::endl;
	//Mesh mesh;
	//MeshTools::ReadMesh(mesh, filename);
	//std::cout << "Face num: " << mesh.n_faces() << std::endl;

	//SegMesh seg_mesh(mesh);
	//std::cout << "opt " << std::endl;
	//Optimization opt(seg_mesh);
	//opt.Deformation();

	//// ============ Optimization with deformation ==============
	//std::string filename = "mesh_data/deform_9.obj";
	//std::cout << "Read mesh: " << filename << std::endl;
	//Mesh def_mesh;
	//MeshTools::ReadMesh(def_mesh, filename);

	//std::cout << "Face num: " << def_mesh.n_faces() << std::endl;

	//SegMesh seg_mesh(def_mesh);
	//std::cout << "opt " << std::endl;
	//Optimization opt(seg_mesh);
	//opt.Segmentation();


	//// ============ Optimization with deformation and segmentation ==============
	//// Deformed
	//std::string deformed_name = "mesh_data/first_deform.obj";
	//Mesh mesh;
	//MeshTools::ReadMesh(mesh, deformed_name);

	//std::string segname = "mesh_data/first_deform_seg.txt";
	//SegMesh seg_mesh(mesh);
	//seg_mesh.ReadSeg(segname);

	//std::cout << "opt " << std::endl;
	//Optimization opt(seg_mesh);
	//opt.SetOriMesh(input_mesh);
	//opt.Refinement();

	return 0;
}


