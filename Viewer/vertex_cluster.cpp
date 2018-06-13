#include <fstream>
#include <iostream>
#include <map>
#include <memory>
using namespace std;
#include "vertex_cluster.h"
#include "mesh_io.h"

//General Functions
void get_norms(const std::vector<float> &vertices,const std::vector<int> &faces,std::vector<float> *normals) {
  const int face_size = faces.size();
  std::vector<float> face_normals(face_size, 0);
  for (int i = 0; i < face_size; i += 3) {
    Eigen::Vector3d vertex1(vertices[faces[i] * 3], vertices[faces[i] * 3 + 1],vertices[faces[i] * 3 + 2]);
    Eigen::Vector3d vertex2(vertices[faces[i + 1] * 3], vertices[faces[i + 1] * 3 + 1],vertices[faces[i + 1] * 3 + 2]);
    Eigen::Vector3d vertex3(vertices[faces[i + 2] * 3],vertices[faces[i + 2] * 3 + 1],vertices[faces[i + 2] * 3 + 2]);
    Eigen::Vector3d vertex1vertex2 = vertex2 - vertex1;
    Eigen::Vector3d vertex1vertex3 = vertex3 - vertex1;
    Eigen::Vector3d normal = vertex1vertex2.cross(vertex1vertex3);

    if (normal.norm() < 0.00001) {
      normal = Eigen::Vector3d(0.0, 0.0, 0.0);
    } else {
      normal.normalize();
    }

    for (int j = 0; j < 3; ++j) face_normals[i + j] = normal[j];
  }

  const int vertex_size = vertices.size();
  normals->resize(vertex_size, 0);
  for (int i = 0; i < face_size; i += 3) {
    for (int j = 0; j < 3; ++j) {
      int idx = faces[i + j];
      Eigen::Vector3d vertex1(vertices[faces[i + j] * 3],vertices[faces[i + j] * 3 + 1],vertices[faces[i + j] * 3 + 2]);
      Eigen::Vector3d vertex2(vertices[faces[i + (j + 1) % 3] * 3],vertices[faces[i + (j + 1) % 3] * 3 + 1],vertices[faces[i + (j + 1) % 3] * 3 + 2]);
      Eigen::Vector3d vertex3(vertices[faces[i + (j + 2) % 3] * 3],vertices[faces[i + (j + 2) % 3] * 3 + 1],vertices[faces[i + (j + 2) % 3] * 3 + 2]);
      Eigen::Vector3d vertex1vertex2 = vertex2 - vertex1;
      Eigen::Vector3d vertex1vertex3 = vertex3 - vertex1;
      double vert_data = acos(vertex1vertex2.dot(vertex1vertex3) / (vertex1vertex2.norm() * vertex1vertex3.norm()));

      if (vert_data == vert_data) {
        for (int k = 0; k < 3; ++k) {
          (*normals)[idx * 3 + k] += face_normals[i + k] * vert_data;
        }
      }
    }
  }

  const int normals_size = normals->size();
  for (int i = 0; i < normals_size; i += 3) {
    Eigen::Vector3d normal((*normals)[i], (*normals)[i + 1], (*normals)[i + 2]);
    if (normal.norm() > 0.00001) {
      normal.normalize();
    } else {
      normal = Eigen::Vector3d(0, 0, 0);
    }

    for (int j = 0; j < 3; ++j) (*normals)[i + j] = normal[j];
  }
}

void float_vector(float i1, float i2, float i3, std::vector<float> *vector) {
  (*vector).push_back(i1);
  (*vector).push_back(i2);
  (*vector).push_back(i3);
}

void int_vector(int i1, int i2, int i3, std::vector<int> *vector) {
  (*vector).push_back(i1);
  (*vector).push_back(i2);
  (*vector).push_back(i3);
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
    per_cell.resize(totalCells, 0);
    x_coords.resize(totalCells, 0); y_coords.resize(totalCells, 0); z_coords.resize(totalCells, 0);
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

                    new_faces.push_back(curCell);
                    x_coords[curCell] += verts[v*3]; y_coords[curCell] += verts[v*3+1]; z_coords[curCell] += verts[(v*3)+2];
                    per_cell[curCell] ++;
                    stop = 1; //cell was found. Stop cell loop.
                    } // End of veretices condition
                curCell ++; //look into each grid cell. - current cell
                }
            }
        } //End of cell loop (and stop condition)

    }//End of vert loop. Verts mapped to cells.
    unsigned int currentLevel;
    unsigned int curCellwVert = 0; //track cur cell - ignores empties
    for(unsigned int i = 0; i < per_cell.size(); i++){ //cell loop
        if (per_cell[i] > 0){ //compute current octree level verts
            check_cells[i] = curCellwVert;
            float vecX = 0; float vecY = 0; float vecZ = 0;
            vecX = x_coords[i]; vecY = y_coords[i]; vecZ = z_coords[i];
            vecX = vecX/per_cell[i]; vecY = vecY/per_cell[i]; vecZ = vecZ/per_cell[i]; //this is the new vert
            float_vector(vecX, vecY, vecZ, &LODvertices); // add it to the final list
            curCellwVert ++;
            if (per_cell[i] > 4){ //compute lower octree level list
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
        faceVerts[0] = check_cells[new_faces[faces[f*3+0]]]; faceVerts[1] = check_cells[new_faces[faces[f*3+1]]]; faceVerts[2] = check_cells[new_faces[faces[f*3+2]]]; faceVerts[3] = 0;
        if (faceVerts[0] != faceVerts[1]) {faceVerts[3] += 1;} if (faceVerts[0] != faceVerts[2]) {faceVerts[3] += 1;}
        if (faceVerts[1] != faceVerts[2]) {faceVerts[3] += 1;} if (faceVerts[3] == 3){ //this face is marked to be kept.
            int_vector(faceVerts[0], faceVerts[1], faceVerts[2], &LODfaces);
        }
    }// End of face loop. New faces generated.

    get_norms(LODvertices, LODfaces, &LODnormals);// Generate new normals.
    std::cout << "End of vertex clustering." << std::endl;
}



//..............................................

// Vertex cluster - quadratic error matrix calculations

//··············································

void vertexCluster::QuadCluster(const std::vector<float> &verts, const std::vector<int> &faces, const std::vector<float> &normals) {
//Find the vertices planes: V->F struct, vertex per face count.
    std::vector<Eigen::Matrix4f> KMat; KMat.resize(faces.size()); //fundamental quadratics matrix - one per face
    std::vector<Eigen::Matrix4f> QMat; KMat.resize(verts.size()); //Q matrix - one per valid vertex
    std::vector<std::vector<unsigned int>> VeF; VeF.resize(verts.size()); //faces per vert struct - one per vert
    std::vector<std::pair<unsigned int, unsigned int>> validPairs; //valid vertex pairs

    float pC[4]; // (temp) store coeficients of plane of face
    unsigned int fV[3]; // (temp) stores the vertices of this face

    for (unsigned int i=0; i < faces.size()/3; i++){ //FACES loop

        for (unsigned int v=0; v <= 3; v++) {//run through vert indices in current face
            //Building struct:
            unsigned int curFaceIndex = i*3 + v;
            VeF[curFaceIndex].push_back(faces[i]); //add the VERT INDEX to the curr FACE vector VeF.
            //Build data for matrix computation
            pC[v] = normals [i*3 + v]; //get a, b, c of plane from normal.
            fV[v] = faces [i*3 + v];
        }

        //Compute Matrices:
        //Get plane d from ax + by + cz + d = 0
        float x = verts[faces[i*3]+0]; //coords of first vertex of this face
        float y = verts[faces[i*3]+1];
        float z = verts[faces[i*3]+2];
        pC[3] = -(pC[0]*x + pC[1]*y + pC[2]*z);
        KMat[i] << pC[0]*pC[0], pC[0]*pC[1], pC[0]*pC[2], pC[0]*pC[3],
                   pC[1]*pC[0], pC[1]*pC[1], pC[1]*pC[2], pC[1]*pC[3],
                   pC[2]*pC[0], pC[2]*pC[1], pC[2]*pC[2], pC[2]*pC[3],
                   pC[3]*pC[0], pC[3]*pC[1], pC[3]*pC[2], pC[3]*pC[3];

    }//end of face loop

    //
    for (unsigned int i=0; i < verts.size()/3; i++){//VERTS loop
        //1)Compute Q matrices as the sum of K matrices
        std::vector <unsigned int> myFaces = VeF[i];
        QMat[i] << 0, 0, 0, 0,
                   0, 0, 0, 0,
                   0, 0, 0, 0,
                   0, 0, 0, 0;
        for (unsigned int j = 0; j < myFaces.size(); j++){ //sum the matrices
        QMat[i] += KMat[myFaces [j]]; //add the face's K mat
        }
    }
}

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



