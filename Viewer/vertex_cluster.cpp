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
int vertsInLastLOD = 0;
void vertexCluster::Cluster(const std::vector<float> &verts, const std::vector<int> &faces, Eigen::Vector3f &min, Eigen::Vector3f &max){
   //General :find the bounding box of mesh
    float pos = (max[0] - min[0])/4;
    int x_axis = (int) ((max[0] - min[0])/pos + 0.5) + 1;
    int y_axis = (int) ((max[1] - min[1])/pos + 0.5) + 1;
    int z_axis = (int) ((max[2] - min[2])/pos + 0.5) + 1;

    int cells = x_axis*y_axis*z_axis;
    const unsigned int vert_count = verts.size() / 3;
    per_cell.resize(cells, 0);
    x_coords.resize(cells, 0); y_coords.resize(cells, 0); z_coords.resize(cells, 0);


    std::cout << "::::::General Initial Data::::: " << std::endl;
    std::cout << "Total number of vertices: " << vert_count << std::endl;
    std::cout << "Total number of cells: " << cells << std::endl;
    std::cout << "::::::::Started vertex cluster::::::: " << x_axis << " x " << y_axis << " x " << z_axis << std::endl;

    unsigned int arbitrary_cell = 0; //track current cell
    for (unsigned int v = 0; v < vert_count; v++){
        arbitrary_cell = 0;
        bool pos_end = 0;
        for(float x = min[0]; x <= (max[0] + pos) && pos_end == 0; x+=pos){
            for(float y = min[1]; y <= (max[1] + pos) && pos_end == 0; y+=pos){
                for(float z = min[2]; z <= (max[2] + pos) && pos_end == 0; z+=pos){
                if (verts[v*3] >= x && verts[v*3] < (x + pos) //vert.x
                && verts[v*3+1] >= y && verts[v*3+1] < (y + pos) //vert.y
                && verts[v*3+2] >= z && verts[v*3+2] < (z + pos) && pos_end == 0)
                    {
                    new_faces.push_back(arbitrary_cell);
                    cellMin[arbitrary_cell] = {x, y, z};
                    x_coords[arbitrary_cell] += verts[v*3]; y_coords[arbitrary_cell] += verts[v*3+1]; z_coords[arbitrary_cell] += verts[(v*3)+2];
                    verts_in_cell[arbitrary_cell].push_back(v);
                    per_cell[arbitrary_cell] ++;
                    pos_end = 1;
                    }
                arbitrary_cell ++;
                }
            }
        }
    }
    unsigned int Level;
    unsigned int empty_arbitrary_cell = 0; //track cur cell - ignores empties
    int new_level = 0;
    for(unsigned int i = 0; i < cells; i++){
        if (per_cell[i] > 0){ //compute current octree level verts
            check_cells[i] = empty_arbitrary_cell;
            float vecX = 0; float vecY = 0; float vecZ = 0;
            vecX = x_coords[i]; vecY = y_coords[i]; vecZ = z_coords[i];
            vecX = vecX/per_cell[i]; vecY = vecY/per_cell[i]; vecZ = vecZ/per_cell[i]; //this is the new vert

            if (i >= vertsInLastLOD){
                           vecX = vecX/per_cell[i]; vecY = vecY/per_cell[i]; vecZ = vecZ/per_cell[i]; //this is the new vert
                       }


            if (Level == 0) {float_vector(vecX, vecY, vecZ, &LODvertices);} // add it to the final list
            else if (Level == 1) {float_vector(vecX, vecY, vecZ, &LODvertices1);} // add it to the final list
            else if (Level == 2) {float_vector(vecX, vecY, vecZ, &LODvertices2);} // add it to the final list
            else if (Level == 3) {float_vector(vecX, vecY, vecZ, &LODvertices3);} // add it to the final list

            if (per_cell[i] > 36 && Level < 4){ //compute lower octree level if there are enough verts in this cell
                float lower_pos = pos/((Level +1)*2); //pos for next LOD
                for (int c = 0;c<8;c++){per_cell.push_back(0);
                    x_coords.push_back(0);
                    y_coords.push_back(0);
                    z_coords.push_back(0);
                    cellMin.push_back({0, 0, 0});} //add 8 new verts


                for (int v=0; v<per_cell[i]; v++){//loop through verts inside this cell
                    int pos_end = 0; //pos_end checking when cell is found
                    int subCellndex = 0; //index of subcell
                    int vertIndex = verts_in_cell[i][v]; //GETTING THE WRONG INFO!!!!!!!!
                                //std::cout << verts_in_cell[i][v] << std::endl;
                    float xV = verts[3*vertIndex+0]; float yV = verts[3*vertIndex+1]; float zV = verts[3*vertIndex+2];
                    for(float x = cellMin[i][0]; x <= (cellMin[i][0] + pos - 0.001); x+=lower_pos){
                    for(float y = cellMin[i][1]; y <= (cellMin[i][1] + pos - 0.001); y+=lower_pos){
                    for(float z = cellMin[i][2]; z <= (cellMin[i][2] + pos - 0.001); z+=lower_pos){
                        if (pos_end == 0 && xV >= x && xV < (x + lower_pos + 0.001) //vert.x
                        && yV >= y && yV < (y + lower_pos + 0.001) //vert.y
                        && zV >= z && zV < (z + lower_pos + 0.001)){
                                        int newCellInd = cells + new_level + subCellndex;
                                        x_coords[newCellInd] += xV;
                                        y_coords[newCellInd] += yV;
                                        z_coords[newCellInd] += zV;
                                        per_cell[newCellInd] ++;
                                        verts_in_cell[newCellInd].push_back(vertIndex);
                                        new_faces_again[vertIndex] = newCellInd;
                                        if(cellMin[newCellInd][0] == 0){cellMin[newCellInd] = {x, y, z};}

                                        //old_to_new[vertIndex] = cells+new_level+subCellndex;
                                        pos_end = 1;
                                    }

                                    subCellndex ++;
                                }}}
                            }
                            new_level ++;
                            per_cell[i] = 0;
                        }
                        else {
                           for (int v=0; v<per_cell[i]; v++){new_faces_again[verts_in_cell[i][v]] = i;}}
                        empty_arbitrary_cell ++;


                        }
                if(i == cells -1 && Level < 4)
                {
                    i=0;
                    vertsInLastLOD = cells;
                    cells = per_cell.size();
                    computeFacesNormals(Level, faces, new_faces);
                    Level++;
                    empty_arbitrary_cell = 0; //track cur cell - ignores empties
                    new_faces = new_faces_again; //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! probably a better way to do this.
                    } //change the LOD level
                }//End of cell loop. New positions defined.
            // Generate new normals.
                std::cout << "End of vertex clustering." << std::endl;
            }

void vertexCluster::computeFacesNormals(int LOD, const std::vector<int> &faces, const std::vector<int> &oldNew){
    const unsigned int numFaces = faces.size() / 3;
    int faceVerts[4] = {0, 0, 0, 0}; //face verts. [3] = weather or not this face is kept.
    for (unsigned int f = 0; f < numFaces; f++){ //loop through faces. Find where their verts are. Change the indices to the new verts.
        faceVerts[0] = check_cells[oldNew[faces[f*3+0]]]; faceVerts[1] = check_cells[oldNew[faces[f*3+1]]]; faceVerts[2] = check_cells[oldNew[faces[f*3+2]]]; faceVerts[3] = 0;
        if (faceVerts[0] != faceVerts[1]) {faceVerts[3] += 1;} if (faceVerts[0] != faceVerts[2]) {faceVerts[3] += 1;}
        if (faceVerts[1] != faceVerts[2]) {faceVerts[3] += 1;} if (faceVerts[3] == 3){ //this face is marked to be kept.
            if (LOD == 0){int_vector(faceVerts[0], faceVerts[1], faceVerts[2], &LODfaces);}
            else if (LOD == 1){int_vector(faceVerts[0], faceVerts[1], faceVerts[2], &LODfaces1);}
            else if (LOD == 2){int_vector(faceVerts[0], faceVerts[1], faceVerts[2], &LODfaces2);}
            else if (LOD == 3){int_vector(faceVerts[0], faceVerts[1], faceVerts[2], &LODfaces3);}
                  }
                }// End of face loop. New faces generated.
    if (LOD == 0){get_norms(LODvertices, LODfaces, &LODnormals);}
    else if (LOD == 1){get_norms(LODvertices1, LODfaces1, &LODnormals1);}
    else if (LOD == 2){get_norms(LODvertices2, LODfaces2, &LODnormals2);}
    else if (LOD == 3){get_norms(LODvertices3, LODfaces3, &LODnormals3);}
}



