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
       std::unique_ptr<data_representation::TriangleMesh>(new data_representation::TriangleMesh());
     //std::make_unique<data_representation::TriangleMesh>();

  bool res = false;
  if (type.compare("ply") == 0) {
    res = data_representation::ReadFromPly(file, mesh.get());
  }

  if (res) {
    mesh_.reset(mesh.release());
    camera_.UpdateModel(mesh_->min_, mesh_->max_);

    GLuint BufferCoords;
    GLuint BufferIndxs;
    GLuint Buffernorms;

    std::vector<float> verts = mesh_->vertices_;
    std::vector<int> faces = mesh_->faces_;
    std::vector<float> normals = mesh_->normals_;
    glGenVertexArrays(1, &VAO);
    glBindVertexArray(VAO);
    glGenBuffers(1, &BufferCoords);
    glGenBuffers(1, &BufferIndxs);
    glGenBuffers(1, &Buffernorms);
    //Bind/ pass data verts
    glBindBuffer(GL_ARRAY_BUFFER, BufferCoords);
    glBufferData(GL_ARRAY_BUFFER, verts.size()*sizeof(float),
    &verts[0], GL_STATIC_DRAW);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
    glEnableVertexAttribArray(0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, Buffernorms);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, normals.size()*sizeof(float),
    &normals[0], GL_STATIC_DRAW);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, 0);
    glEnableVertexAttribArray(1);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, BufferIndxs);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, faces.size()*sizeof(int),
    &faces[0], GL_STATIC_DRAW);
    glBindVertexArray(0);
    emit SetFaces(QString(std::to_string(mesh_->faces_.size() / 3).c_str()));
    emit SetVertices(
        QString(std::to_string(mesh_->vertices_.size() / 3).c_str()));
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

  program_ = std::unique_ptr<QOpenGLShaderProgram>(new(QOpenGLShaderProgram));
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
    Eigen::Matrix4f model = camera_.SetModel();

    Eigen::Matrix4f t = view * model;
    Eigen::Matrix3f normal;
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j) normal(i, j) = t(i, j);

    normal = normal.inverse().transpose();

    if (mesh_ != nullptr) {
      program_->bind();

      GLuint MyShader = program_->programId();

      glBindAttribLocation(MyShader, 0, "vert");
      glBindAttribLocation(MyShader, 1, "normal");

      GLuint projection_location = program_->uniformLocation("projection");
      glUniformMatrix4fv(projection_location, 1, GL_FALSE, projection.data());

      GLuint view_location = program_->uniformLocation("view");
      glUniformMatrix4fv(view_location, 1, GL_FALSE, view.data());

      GLuint model_location = program_->uniformLocation("model");
      glUniformMatrix4fv(model_location, 1, GL_FALSE, model.data());

      GLuint normal_matrix_location =
          program_->uniformLocation("normal_matrix");
      glUniformMatrix3fv(normal_matrix_location, 1, GL_FALSE, normal.data());

      // TODO: Implement model rendering.

      glBindVertexArray(VAO); //Draw
                glDrawElements(GL_TRIANGLES, mesh_->faces_.size(), GL_UNSIGNED_INT, 0);
                glBindVertexArray(0); //unbind the VAO (!!)

        }
      }
      glEnd();

      // END.
    }

  // TODO: Implement framerate displaying.

  // emit SetFramerate(...);

  // END.
