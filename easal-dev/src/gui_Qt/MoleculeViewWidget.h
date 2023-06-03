#ifndef MOLECULEVIEWWIDGET_H
#define MOLECULEVIEWWIDGET_H

#include <Source/gui/AtlasNode.h>
#include <Source/gui/SaveLoader.h>
#include <Source/gui/Settings.h>
#include <Source/gui/Utils.h>
#include <Source/gui_Qt/BondViewWidget.h>
#include <Source/gui_Qt/Mesh3D.h>
#include <Source/gui_Qt/SharedDataGUI.h>

#include <QKeyEvent>
#include <QMouseEvent>
#include <QOpenGLBuffer>
#include <QOpenGLContext>
#include <QOpenGLExtraFunctions>
#include <QOpenGLShaderProgram>
#include <QOpenGLVertexArrayObject>
#include <QOpenGLWidget>
#include <QTimer>
#include <cmath>
#include <iostream>

#define PI 3.14159265

using namespace std;

class MoleculeViewWidget : public QOpenGLWidget {
  Q_OBJECT
 public:
  MoleculeViewWidget(SharedDataGUI *sharedData, BondViewWidget *bondViewWidget);
  /** @brief to render some mlecules we we to load nodes using save loader
   * but save loader is not included in constructor o we need to set it in
   * setSaveLoader */
  void setSaveLoader(SaveLoader *saveLoader);
  void setAtlasNode(AtlasNode *atlasNode);
  void setMolecularUnits(MolecularUnit *molecularUnitA,
                         MolecularUnit *molecularUnitB);
  void initializeGL();
  void resizeGL(int w, int h);
  void paintGL();

 private:
  BondViewWidget *bondViewWidget;
  //[ CAMERA]
  /** @brief vars above take care of camera
   *  cameraPosition represents current coordinates of camera
   *  cameraDirection represents unit vector located at origin, camera looks in
   * direction of this vector */
  QVector3D m_cameraPosition;
  QVector3D m_cameraDirection;
  QVector3D m_cameraUp;
  QVector3D m_cameraLeft;
  double m_cameraRotationalSpeed;
  double m_cameraMoveSpeed;
  double angleYaw;
  double anglePitch;
  /** @brief vars above take care of mouse */
  //[~CAMERA]
  //[ MOUSE]
  bool m_mouseIsPressed;
  QPoint m_mouseLastPosition;
  //[~MOUSE]
  //[ TIMER]
  /** @brief Timer that triggers screen updates once each m_updateTimeInMS*/
  QTimer m_timer;
  double m_updateTimeInMS;
  //[ TIMER]
  //[ UNIFORMS]
  /** @brief Matrices below is part of model view projection concept used in 3D
   * graphics (I will crify if they are row or column major (it determines
   * multiplication order) normal matrix is used for flat shading */
  QMatrix4x4 m_modelMatrix;
  QMatrix4x4 m_viewMatrix;
  QMatrix4x4 m_projectionMatrix;
  QMatrix4x4 m_normalMatrix;
  QMatrix4x4 m_rotationAroundOriginMatrix;
  QVector3D m_lightPositionVector;
  //[ UNIFORMS]
  //[ UNIFORM LOCATIONS]
  /** @brief Uniforms are OpenGL "pointers" to location of matrices above in
   * shader program */
  int m_modelMatrixLocation;
  int m_viewMatrixLocation;
  int m_projectionMatrixLocation;
  int m_normalMatrixLocation;
  int m_lightPositionVectorLocation;
  int m_colorVectorLocation;
  int m_rotationAroundOriginMatrixLocation;
  //[ UNIFORM LOCATIONS]
  //[ MESHES3D]
  Mesh3D m_sphere;
  Mesh3D m_cylinder;
  //[~MESHES3D]
  //[ ATLAS NODE]
  bool atlasNodeIsSet;
  AtlasNode *atlasNode;
  QMatrix4x4 moleculeARotation;
  QMatrix4x4 moleculeBTransform;
  MolecularUnit *molecularUnitA;
  MolecularUnit *molecularUnitB;
  vector<Atom *> atomsA;
  vector<Atom *> atomsB;
  vector<std::pair<int, int>> bondSolid;
  vector<std::pair<int, int>> bondVariable;
  ActiveConstraintRegion *currentACR;
  int currentNodeID;
  int currentSpaceID;
  int currentFlipID;
  //[~ATLAS NODE]
  //[ VBO AND VAO FOR MESHES3D]
  QOpenGLBuffer *m_vbo_sphere;
  QOpenGLBuffer *m_vbo_cylinder;
  QOpenGLVertexArrayObject *m_vao_sphere;
  QOpenGLVertexArrayObject *m_vao_cylinder;
  //[~VBO AND VAO FOR MESHES3D]
  //[ MISC]
  QOpenGLShaderProgram *m_program;
  SaveLoader *saveLoader;
  SharedDataGUI *sharedData;
  //[~MISC]
  /** @brief initMesh3D function fills given vbo and vao with mseh vertex data
   */
  void initMesh3D(QOpenGLBuffer *VBO, QOpenGLVertexArrayObject *VAO,
                  Mesh3D *mesh3D);
  /** @brief renderSphere function renders shpere with given size at position */
  void renderSpere(QVector3D position, float radious, QVector4D color);
  /** @brief draw cylinder between two points (z buffer is set to GL_ALWAYS to
   * render cylender on top of eveything) */
  void renderCylinder(QVector3D positionStart, QVector3D positionEnd,
                      float radious, QVector4D color);
  void renderAtoms();
  void renderBonds();
  void updateSpaceAndFlip();
  void keyPressEvent(QKeyEvent *keyEvent) Q_DECL_OVERRIDE;
  void mousePressEvent(QMouseEvent *mouseEvent) Q_DECL_OVERRIDE;
  void mouseReleaseEvent(QMouseEvent *mouseEvent) Q_DECL_OVERRIDE;
  void mouseMoveEvent(QMouseEvent *mouseEvent) Q_DECL_OVERRIDE;
  //    QOpenGLShaderProgram *m_program;
  QByteArray versionedShaderCode(const char *src);
};

#endif  // MOLECULEVIEWWIDGET_H
