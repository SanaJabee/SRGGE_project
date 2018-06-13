#include <fstream>
#include <iostream>
#include <map>
#include <memory>
using namespace std;
#include "vertex_cluster.h"
#include "mesh_io.h"

void Add3Items(float i1, float i2, float i3, std::vector<float> *vector) {
  (*vector).push_back(i1);
  (*vector).push_back(i2);
  (*vector).push_back(i3);
}

void Add3Ints(int i1, int i2, int i3, std::vector<int> *vector) {
  (*vector).push_back(i1);
  (*vector).push_back(i2);
  (*vector).push_back(i3);
}

void ComputeVertexNormals(const std::vector<float> &vertices,
                          const std::vector<int> &faces,
                          std::vector<float> *normals) {
  const int kFaces = faces.size();
  std::vector<float> face_normals(kFaces, 0);

  for (int i = 0; i < kFaces; i += 3) {
    Eigen::Vector3d v1(vertices[faces[i] * 3], vertices[faces[i] * 3 + 1],
                       vertices[faces[i] * 3 + 2]);
    Eigen::Vector3d v2(vertices[faces[i + 1] * 3],
                       vertices[faces[i + 1] * 3 + 1],
                       vertices[faces[i + 1] * 3 + 2]);
    Eigen::Vector3d v3(vertices[faces[i + 2] * 3],
                       vertices[faces[i + 2] * 3 + 1],
                       vertices[faces[i + 2] * 3 + 2]);
    Eigen::Vector3d v1v2 = v2 - v1;
    Eigen::Vector3d v1v3 = v3 - v1;
    Eigen::Vector3d normal = v1v2.cross(v1v3);

    if (normal.norm() < 0.00001) {
      normal = Eigen::Vector3d(0.0, 0.0, 0.0);
    } else {
      normal.normalize();
    }

    for (int j = 0; j < 3; ++j) face_normals[i + j] = normal[j];
  }

  const int kVertices = vertices.size();
  normals->resize(kVertices, 0);
  for (int i = 0; i < kFaces; i += 3) {
    for (int j = 0; j < 3; ++j) {
      int idx = faces[i + j];
      Eigen::Vector3d v1(vertices[faces[i + j] * 3],
                         vertices[faces[i + j] * 3 + 1],
                         vertices[faces[i + j] * 3 + 2]);
      Eigen::Vector3d v2(vertices[faces[i + (j + 1) % 3] * 3],
                         vertices[faces[i + (j + 1) % 3] * 3 + 1],
                         vertices[faces[i + (j + 1) % 3] * 3 + 2]);
      Eigen::Vector3d v3(vertices[faces[i + (j + 2) % 3] * 3],
                         vertices[faces[i + (j + 2) % 3] * 3 + 1],
                         vertices[faces[i + (j + 2) % 3] * 3 + 2]);

      Eigen::Vector3d v1v2 = v2 - v1;
      Eigen::Vector3d v1v3 = v3 - v1;
      double angle = acos(v1v2.dot(v1v3) / (v1v2.norm() * v1v3.norm()));

      if (angle == angle) {
        for (int k = 0; k < 3; ++k) {
          (*normals)[idx * 3 + k] += face_normals[i + k] * angle;
        }
      }
    }
  }

  const int kNormals = normals->size();
  for (int i = 0; i < kNormals; i += 3) {
    Eigen::Vector3d normal((*normals)[i], (*normals)[i + 1], (*normals)[i + 2]);
    if (normal.norm() > 0.00001) {
      normal.normalize();
    } else {
      normal = Eigen::Vector3d(0, 0, 0);
    }

    for (int j = 0; j < 3; ++j) (*normals)[i + j] = normal[j];
  }
}

//..............................................

// Vertex cluster - Average all verts in cell + octree

//··············································

void vertexCluster::Cluster(const std::vector<float> &verts, const std::vector<int> &faces, Eigen::Vector3f &min, Eigen::Vector3f &max){
   //General :find the bounding box of mesh
    float step = (max[0] - min[0])/8; //set step to be based on bounding box
    int dimX = (int) ((max[0] - min[0])/step + 0.5) + 1;
    int dimY = (int) ((max[1] - min[1])/step + 0.5) + 1;
    int dimZ = (int) ((max[2] - min[2])/step + 0.5) + 1;

    int totalCells = dimX*dimY*dimZ;
    const unsigned int numVertices = verts.size() / 3;
    verts_per_cell.resize(totalCells, 0);
    new_coords_x.resize(totalCells, 0); new_coords_y.resize(totalCells, 0); new_coords_z.resize(totalCells, 0);
    //std::vector<float> cellMin; std::vector<float> cellMax; //store cell min amd max position (octree) - well see about this.
    //cellMin.resize(totalCells); cellMax.resize(totalCells);
   std::cout << "Total number of cells: " << totalCells << std::endl;
    std::cout << "Started vertex cluster. Initial grid dimmensions: " << dimX << " x " << dimY << " x " << dimZ << std::endl;
    std::cout << "Initial number of vertices: " << numVertices << std::endl;
//min max bounding boxes for each axis
    std::cout << "Bounding box x: " << min[0] << ", " << max[0] << std::endl;
    std::cout << "Bounding box y: " << min[1] << ", " << max[1] << std::endl;
    std::cout << "Bounding box z: " << min[2] << ", " << max[2] << std::endl;

    unsigned int curCell = 0; //track current cell
    for (unsigned int v = 0; v < numVertices; v++){ //For each vert
        curCell = 0;
        bool stop = 0;
        for(float x = min[0]; x <= (max[0] + step) && stop == 0; x+=step){
            for(float y = min[1]; y <= (max[1] + step) && stop == 0; y+=step){
                for(float z = min[2]; z <= (max[2] + step) && stop == 0; z+=step){
                if (verts[v*3] >= x && verts[v*3] < (x + step) //vert.x
                && verts[v*3+1] >= y && verts[v*3+1] < (y + step) //vert.y
                && verts[v*3+2] >= z && verts[v*3+2] < (z + step) && stop == 0) //vert.z
                    { //vertex is inside cell. Add to count. Do nothing on empty cells

                    old_to_new.push_back(curCell);
                    new_coords_x[curCell] += verts[v*3]; new_coords_y[curCell] += verts[v*3+1]; new_coords_z[curCell] += verts[(v*3)+2];
                    verts_per_cell[curCell] ++;
                    stop = 1; //cell was found. Stop cell loop.
                    } // End of veretices condition
                curCell ++; //look into each grid cell. - current cell
                }
            }
        } //End of cell loop (and stop condition)

    }//End of vert loop. Verts mapped to cells.
    unsigned int currentLevel;
    unsigned int curCellwVert = 0; //track cur cell - ignores empties
    for(unsigned int i = 0; i < verts_per_cell.size(); i++){ //cell loop
        if (verts_per_cell[i] > 0){ //compute current octree level verts
            cellToVertIndex[i] = curCellwVert;
            float vecX = 0; float vecY = 0; float vecZ = 0;
            vecX = new_coords_x[i]; vecY = new_coords_y[i]; vecZ = new_coords_z[i];
            vecX = vecX/verts_per_cell[i]; vecY = vecY/verts_per_cell[i]; vecZ = vecZ/verts_per_cell[i]; //this is the new vert
            Add3Items(vecX, vecY, vecZ, &LODvertices); // add it to the final list
            curCellwVert ++;
            if (verts_per_cell[i] > 4){ //compute lower octree level list
                   //
            }
            /*if (vecX >= min[0] && vecX < max[0] //vert.x
            && vecY >= min[1] && vecY < max[1] //vert.y
            && vecZ >= min[2] && vecZ < max[2]){} //vert.z
            else {std::cout << vecX << ", " << vecY << ", " << vecZ << ", " << i << ". " << std::endl;}*/
             //Check if any new vert ended outside of bounding box

        }
    }//End of cell loop. New positions defined.

    const unsigned int numFaces = faces.size() / 3;
    int faceVerts[4] = {0, 0, 0, 0}; //face verts. [3] = weather or not this face is kept.
    for (unsigned int f = 0; f < numFaces; f++){ //loop through faces. Find where their verts are. Change the indices to the new verts.
        faceVerts[0] = cellToVertIndex[old_to_new[faces[f*3+0]]]; faceVerts[1] = cellToVertIndex[old_to_new[faces[f*3+1]]]; faceVerts[2] = cellToVertIndex[old_to_new[faces[f*3+2]]]; faceVerts[3] = 0;
        if (faceVerts[0] != faceVerts[1]) {faceVerts[3] += 1;} if (faceVerts[0] != faceVerts[2]) {faceVerts[3] += 1;}
        if (faceVerts[1] != faceVerts[2]) {faceVerts[3] += 1;} if (faceVerts[3] == 3){ //this face is marked to be kept.
            Add3Ints(faceVerts[0], faceVerts[1], faceVerts[2], &LODfaces);
        }
    }// End of face loop. New faces generated.

    ComputeVertexNormals(LODvertices, LODfaces, &LODnormals);// Generate new normals.
    std::cout << "End of vertex clustering." << std::endl;
}



//..............................................


/*
ADV: shape preserving
Define 8 normal directions
build the reg grid/ octree
compute mean with the QEM
*/

/*
octrees:
Bounding boxes -> 1st level of the grid -> redivide cells based on amount of vers verts cell.

loop through all faces -> find the verts that have been collapsed -> change them to the new verts.
*/



