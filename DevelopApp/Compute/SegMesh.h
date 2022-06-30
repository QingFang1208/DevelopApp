#pragma once
#include <fstream>
#include <iostream>
#include "../MeshDefinition.h"

class SegMesh
{
public:
	SegMesh(Mesh& mesh);

	void ReadSeg(const std::string& file_name);
	
	void Init();

	Mesh& GetMesh();

	std::vector<bool>& GetSeam();
	void BoundToIdx();
	
	std::vector<int>& GetSegId();
	void IdxToBound();

	int& GetSegNum();

	void ComputeVertexStatus(std::vector<bool>& seam_v_status);

	void ComputeSegFaceNum(std::vector<int>& f_num);

	std::vector<double> GetAbsGauss(bool remove_seam = false);
	
	void SetModifiedGauss(std::vector<double> v_gauss);
	std::vector<double>& GetModifiedGauss();

	void ComputeVertexCount(std::vector<int>& v_cout);

	void ComputeVertexSeg(std::vector<int>& v_seg);

	void WriteSeg(const std::string& file_name);

	void WriteMesh(const std::string& file_name);


private:
	Mesh& mesh_;

	int seg_num_;
	std::vector<bool> seam_status_;
	std::vector<int> seg_id_;

	std::vector<double> abs_v_gauss_;
	std::vector<double> modified_gauss_;
};

