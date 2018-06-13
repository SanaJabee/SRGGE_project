// Author: Marc Comino 2018

#include "./glwidget.h"

#include <fstream>
#include <iostream>
#include <memory>

#include "./triangle_mesh.h"
#include "./mesh_io.h"
#include "./vertex_cluster.h"

namespace {

const float kFieldOfView = 60;
const float kZNear = 0.0001;
const float kZFar = 10;

const char kVertexShaderFile[] = "../shaders/phong.vert";
const char kFragmentShaderFile[] = "../shaders/phong.frag";

const int kVertexAttributeIdx = 0;
const int kNormalAttributeIdx = 1;

bool ReadFile(const std::string filename, std::string *shader_source) {
  std::ifstream infile(filename.c_str());

  if (!infile.is_open() || !infile.good()) {
    std::cerr << "Error " + filename + " not found." << std::endl;
    return false;
  }

  std::stringstream stream;
  stream << infile.rdbuf();
  infile.close();
  *shader_source = stream.str();
  return true;
}

}  // namespace

GLWidget::GLWidget(QWidget *parent)
    : QGLWidget(parent), initialized_(false), width_(0.0), height_(0.0) {
  setFocusPolicy(Qt::StrongFocus);
}

GLWidget::~GLWidget() {}

bool GLWidget::LoadModel(QString filename) {
  std::string file = filename.toUtf8().constData();
  uint pos = file.find_last_of(".");
  std::string type = file.substr(pos + 1);

  std::unique_ptr<data_representation::TriangleMesh> mesh =
      std::unique_ptr<data_representation::TriangleMesh>(new data_representation::TriangleMesh);

  bool res = false;
  if (type.compare("ply") == 0) {
    res = data_representation::ReadFromPly(file, mesh.get());
  }

  if (res) {
    mesh_.reset(mesh.release());
    camera_.UpdateModel(mesh_->min_, mesh_->max_);

    GLuint coordBufferID; //VBO verts
    GLuint indexBuffersID; //VBO faces
    GLuint normBuffersID; //VBO norms

    //get new vertices from VERTEX CLUSTERING
    static vertClust vC;
    vC.Cluster(mesh_->vertices_, mesh_->faces_, mesh_->min_, mesh_->max_);
    //VERTEX CLUSTERING END
    std::vector<float> verts = vC.LODverts;
    std::vector<int> faces = vC.LODfaces;
    std::vector<float> normals = vC.LODnormals;

    emit SetFaces(QString(std::to_string(faces.size() / 3).c_str()));
    emit SetVertices(QString(std::to_string(verts.size() / 3).c_str()));

    // TODO: Create / Initialize buffers.
    //Vertex buffer objects
    //Create empty VAO/ Buffers
    glGenVertexArrays(1, &VAO); //generate VAO
    glBindVertexArray(VAO); //bind VAO
    //-----//
    glGenBuffers(1, &coordBufferID); //VBO verts
    glGenBuffers(1, &indexBuffersID); //VBO faces
    glGenBuffers(1, &normBuffersID); //VBO faces
    //Bind/ pass data verts
    glBindBuffer(GL_ARRAY_BUFFER, coordBufferID); //bind coord buffer
    glBufferData(GL_ARRAY_BUFFER, verts.size()*sizeof(float), // sizeof(boxVerts)/sizeof(boxVerts[0]),
                   &verts[0], GL_STATIC_DRAW); // ONCE - pass data to buffer
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0); // VAO - set id and size
    glEnableVertexAttribArray(0); // Enable 0
    //Bind / pass data normals
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, normBuffersID); // VBO
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, normals.size()*sizeof(float),
                   &normals[0], GL_STATIC_DRAW);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, 0);
    glEnableVertexAttribArray(1);
    //Bind / pass data faces
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indexBuffersID); // VBO
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, faces.size()*sizeof(int),
                   &faces[0], GL_STATIC_DRAW);
    glBindVertexArray(0);//unbind VAO (!!)
    //BUFFERS!
    // END.

    return true;
  }

  return false;
}

void GLWidget::initializeGL() {
  glewInit();

  glEnable(GL_NORMALIZE);
  glEnable(GL_CULL_FACE);
  glCullFace(GL_BACK);
  glEnable(GL_DEPTH_TEST);

  std::string vertex_shader, fragment_shader;
  bool res = ReadFile(kVertexShaderFile, &vertex_shader) &&
             ReadFile(kFragmentShaderFile, &fragment_shader);

  if (!res) exit(0);

  program_ = std::unique_ptr<QOpenGLShaderProgram>(new QOpenGLShaderProgram);
  program_->addShaderFromSourceCode(QOpenGLShader::Vertex,
                                    vertex_shader.c_str());
  program_->addShaderFromSourceCode(QOpenGLShader::Fragment,
                                    fragment_shader.c_str());
  program_->bindAttributeLocation("vertex", kVertexAttributeIdx);
  program_->bindAttributeLocation("normal", kNormalAttributeIdx);
  program_->link();

  initialized_ = true;
}

void GLWidget::resizeGL(int w, int h) {
  if (h == 0) h = 1;
  width_ = w;
  height_ = h;

  camera_.SetViewport(0, 0, w, h);
  camera_.SetProjection(kFieldOfView, kZNear, kZFar);
}

void GLWidget::mousePressEvent(QMouseEvent *event) {
  if (event->button() == Qt::LeftButton) {
    camera_.StartRotating(event->x(), event->y());
  }
  if (event->button() == Qt::RightButton) {
    camera_.StartZooming(event->x(), event->y());
  }
  updateGL();
}

void GLWidget::mouseMoveEvent(QMouseEvent *event) {
  camera_.SetRotationX(event->y());
  camera_.SetRotationY(event->x());
  camera_.SafeZoom(event->y());
  updateGL();
}

void GLWidget::mouseReleaseEvent(QMouseEvent *event) {
  if (event->button() == Qt::LeftButton) {
    camera_.StopRotating(event->x(), event->y());
  }
  if (event->button() == Qt::RightButton) {
    camera_.StopZooming(event->x(), event->y());
  }
  updateGL();
}

void GLWidget::keyPressEvent(QKeyEvent *event) {
  if (event->key() == Qt::Key_Up) camera_.Zoom(-1);
  if (event->key() == Qt::Key_Down) camera_.Zoom(1);

  if (event->key() == Qt::Key_Left) camera_.Rotate(-1);
  if (event->key() == Qt::Key_Right) camera_.Rotate(1);

  if (event->key() == Qt::Key_W) camera_.Zoom(-1);
  if (event->key() == Qt::Key_S) camera_.Zoom(1);

  if (event->key() == Qt::Key_A) camera_.Rotate(-1);
  if (event->key() == Qt::Key_D) camera_.Rotate(1);

  updateGL();
}

void GLWidget::paintGL() {
  glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  if (initialized_) {
    camera_.SetViewport();

    Eigen::Matrix4f projection = camera_.SetProjection();
    Eigen::Matrix4f view = camera_.SetView();

    if (mesh_ != nullptr) {


      //Rendering model grid
      const int n = 2;
      const float scaler = 1.5;
      for  (int i = -n+1; i < n; ++i) {
        for  (int j = -n+1; j < n; ++j) { // Draw multiple copies for of an object

          Eigen::Matrix4f model = camera_.SetModel();
          Eigen::Matrix4f inverseView = view.inverse();
          Eigen::Vector3f camera (inverseView(0, 3), inverseView(1, 3), inverseView(2, 3));
          Eigen::Vector3f translation (i*scaler, 0, j*scaler);
          Eigen::Vector3f dist = camera - translation;
          const Eigen::Affine3f kTranslation(Eigen::Translation3f(i*scaler, 0, j*scaler));
          model = kTranslation * model;
          distance = dist.norm()/5;

          Eigen::Matrix4f t = view * model;
          Eigen::Matrix3f normal;
          for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) normal(i, j) = t(i, j);

          normal = normal.inverse().transpose();

          program_->bind();

          GLuint projection_location = program_->uniformLocation("projection");
          glUniformMatrix4fv(projection_location, 1, GL_FALSE, projection.data());

          GLuint view_location = program_->uniformLocation("view");
          glUniformMatrix4fv(view_location, 1, GL_FALSE, view.data());

          GLuint model_location = program_->uniformLocation("model");
          glUniformMatrix4fv(model_location, 1, GL_FALSE, model.data());

          GLuint LOD_location = program_->uniformLocation("LOD");
          glUniform1f(LOD_location, distance);

          GLuint normal_matrix_location =
              program_->uniformLocation("normal_matrix");
          glUniformMatrix3fv(normal_matrix_location, 1, GL_FALSE, normal.data());

          glBindVertexArray(VAO); //Draw
          glDrawElements(GL_TRIANGLES, mesh_->faces_.size(), GL_UNSIGNED_INT, 0);
          glBindVertexArray(0); //unbind the VAO (!!)

            }

        //Rendering mode 2 - Museus rooms.
        /*tile based representation.
        -color coded floor plan
        */
          }
          glEnd();

          // END.
        }
      }



  // TODO: Implement framerate displaying.

  // emit SetFramerate(...);

  // END.
}
