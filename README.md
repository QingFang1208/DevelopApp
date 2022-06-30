DevelopApp
===

This is a C++ implementation of the piecewise developable approximation in the following paper:

[Developability-Driven Piecewise Approximations for Triangular Meshes](https://rec.ustc.edu.cn/share/d5e47230-f81b-11ec-9466-2ba4385ffe49).
Zheng-Yu Zhao\*, [Qing Fang](https://qingfang1208.github.io/)\*, Wenqing Ouyang, Zheng Zhang, [Ligang Liu](http://staff.ustc.edu.cn/~lgliu/), [Xiao-Ming Fu](https://ustc-gcl-f.github.io/).
*ACM Transactions on Graphics (SIGGRAPH)*, 41(4), 2022. (*Joint first authors)

The code is written by Zheng-Yu Zhao and Qing Fang using Microsoft Visual Studio 2019 based on the [Surface Mesh Framework](http://staff.ustc.edu.cn/~fuxm/code/index.html#sec_surface_framework), with OpenMP support for parallel programming.

## External Libraries

* [Matlab](https://www.mathworks.com/) (2020b)
* [Eigen](http://eigen.tuxfamily.org/) (3.4)
* [OpenMesh](https://www.openmesh.org/) (8.0)
* OpenMP support

## Output

* final.obj --- the final piecewise approximated result.
* final_seg.txt --- the parts id for triangles in the final result.
