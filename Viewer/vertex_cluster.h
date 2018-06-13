#ifndef VERTEX_CLUSTER
#define VERTEX_CLUSTER

#include <eigen3/Eigen/Geometry>
#include <vector>
#include "glwidget.h"


class vertexCluster{

public:
    void Cluster(const std::vector<float> &verts, const std::vector<int> &faces, Eigen::Vector3f &min, Eigen::Vector3f &max);
    void QuadCluster(const std::vector<float> &verts, const std::vector<int> &faces, const std::vector<float> &normals);
    //send to buffer
    std::vector<float> LODvertices;
    std::vector<int> LODfaces;
    std::vector<float> LODnormals;

    std::vector<float> LODvertices1;
    std::vector<int> LODfaces1;
    std::vector<float> LODnormals1;

    std::vector<float> LODvertices2;
    std::vector<int> LODfaces2;
    std::vector<float> LODnormals2;

    std::vector<float> LODvertices3;
    std::vector<int> LODface3;
    std::vector<float> LODnormals3;
private:
    //rebuild face
    std::vector<double> x_coords;
    std::vector<double> y_coords;
    std::vector<double> z_coords;
    std::vector<unsigned int> per_cell;
    std::vector<int> new_faces;
    std::map <unsigned int, unsigned int> check_cells; //Cell ids (includes empties) to vert coords (no empties)

    };


#endif // VERTEX_CLUSTER
