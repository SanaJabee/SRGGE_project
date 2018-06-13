#ifndef VERTEX_CLUSTER
#define VERTEX_CLUSTER

#include <eigen3/Eigen/Geometry>
#include <vector>
#include "glwidget.h"


class vertClust {
    //generate on a regular grid
    //using the mean average vertex for each vert
    //generate 4 levels of detail
public:
    void Cluster(const std::vector<float> &verts, const std::vector<int> &faces, Eigen::Vector3f &min, Eigen::Vector3f &max);
    void QuadCluster(const std::vector<float> &verts, const std::vector<int> &faces, const std::vector<float> &normals);
    //send to buffer
    std::vector<float> LODverts;
    std::vector<int> LODfaces;
    std::vector<float> LODnormals;

    std::vector<float> LODverts1;
    std::vector<int> LODfaces1;
    std::vector<float> LODnormals1;

    std::vector<float> LODverts2;
    std::vector<int> LODfaces2;
    std::vector<float> LODnormals2;

    std::vector<float> LODverts3;
    std::vector<int> LODface3;
    std::vector<float> LODnormals3;

private:
    //rebuild faces
    std::vector<unsigned int> verts_per_cell;
    std::vector<int> old_to_new; //new verts. index = old verts index
    std::vector<double> new_coords_x;
    std::vector<double> new_coords_y;
    std::vector<double> new_coords_z; //coords. index = cell id (0s included)
    std::map <unsigned int, unsigned int> cellToVertIndex; //Cell ids (includes empties) to vert coords (no empties)

};


#endif // VERTEX_CLUSTER
