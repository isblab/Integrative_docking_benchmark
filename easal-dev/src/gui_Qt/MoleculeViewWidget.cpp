#include "MoleculeViewWidget.h"

MoleculeViewWidget::MoleculeViewWidget(SharedDataGUI *sharedData,
                                       BondViewWidget *bondViewWidget) {
  this->sharedData = sharedData;
  this->bondViewWidget = bondViewWidget;
  //[ CAMERA SETTINGS]
  m_cameraPosition = QVector3D(0, 0, 15);
  m_cameraDirection = QVector3D(0, 0, -1);
  m_cameraRotationalSpeed = 2;
  m_cameraMoveSpeed = 0.1;
  m_rotationAroundOriginMatrix.setToIdentity();
  anglePitch = 0;
  angleYaw = 0;
  //[~CAMERA SETTINGS]
  //[ MOUSE SETTINGS]
  m_mouseIsPressed = false;
  m_mouseLastPosition = QPoint(0, 0);
  //[~MOUSE SETTINGS]
  //[ TIMER SETTINGS]
  // start argument determins update rate
  m_updateTimeInMS = 33;
  m_timer.start(m_updateTimeInMS);
  // repaint is function of QWidget that triggers paintEvent
  connect(&m_timer, SIGNAL(timeout()), this, SLOT(update()));
  //[ TIMER SETTINGS]
  //[ WIDGET SETTINGS]
  setFocusPolicy(Qt::StrongFocus);
  setUpdateBehavior(UpdateBehavior::NoPartialUpdate);
  //[~WIDGET SETTINGS]
  //[ MISC]
  currentNodeID = -1;
  currentSpaceID = -1;
  currentFlipID = -1;
  //[~MISC]
}

void MoleculeViewWidget::setSaveLoader(SaveLoader *saveLoader) {
  this->saveLoader = saveLoader;
}

QByteArray MoleculeViewWidget::versionedShaderCode(const char *src) {
  QByteArray versionedSrc;

  if (QOpenGLContext::currentContext()->isOpenGLES())
    versionedSrc.append(QByteArrayLiteral("#version 300 es\n"));
  else
    versionedSrc.append(QByteArrayLiteral("#version 330\n"));

  versionedSrc.append(src);
  return versionedSrc;
}
static const char *vertexShaderSource =
    "layout(location = 0) in vec4 vertex;\n"
    "layout(location = 1) in vec3 normal;\n"
    "uniform mat4 modelMatrix;\n"
    "uniform mat4 viewMatrix;\n"
    "uniform mat4 projectionMatrix;\n"
    "uniform mat4 normalMatrix;\n"
    "uniform mat4 rotationAroundOriginMatrix;\n"
    "out vec3 vertexNormal;\n"
    "out vec3 vertexCoordinate;\n"
    "void main() {\n"
    "   vertexNormal = mat3(transpose(inverse( rotationAroundOriginMatrix * "
    "modelMatrix))) * normal;\n"
    "   vertexCoordinate = vec3(modelMatrix * vertex);\n"
    "   gl_Position =   projectionMatrix * viewMatrix * "
    "rotationAroundOriginMatrix * modelMatrix * vertex;\n"
    "}\n";

static const char *fragmentShaderSource =
    "in highp vec3 vertexNormal;\n"
    "in highp vec3 vertexCoordinate;\n"
    "out highp vec4 fragColor;\n"
    "uniform highp vec3 lightPositionVector;\n"
    "uniform highp vec4 color;\n"
    "void main() {\n"
    "   highp vec3 lightVector = normalize(lightPositionVector - "
    "vertexCoordinate);\n"
    "   highp float lightFactor = max(dot(normalize(vertexNormal), "
    "lightVector), 0.0);\n"
    "   fragColor = clamp(color * 0.5 + color * 0.5 * lightFactor, 0.0, 1.0);\n"
    "}\n";

void MoleculeViewWidget::initializeGL() {
  QOpenGLFunctions *f = QOpenGLContext::currentContext()->functions();

  //[ SET SHADING PROGRAM]
  m_program = new QOpenGLShaderProgram;
  m_program->addShaderFromSourceCode(QOpenGLShader::Vertex,
                                     versionedShaderCode(vertexShaderSource));
  m_program->addShaderFromSourceCode(QOpenGLShader::Fragment,
                                     versionedShaderCode(fragmentShaderSource));
  m_program->link();
  //[ SET SHADING PROGRAM]

  //[ GET UNIFORM LOCATIONS]
  m_modelMatrixLocation = m_program->uniformLocation("modelMatrix");
  m_viewMatrixLocation = m_program->uniformLocation("viewMatrix");
  m_projectionMatrixLocation = m_program->uniformLocation("projectionMatrix");
  m_normalMatrixLocation = m_program->uniformLocation("normalMatrix");
  m_lightPositionVectorLocation =
      m_program->uniformLocation("lightPositionVector");
  m_rotationAroundOriginMatrixLocation =
      m_program->uniformLocation("rotationAroundOriginMatrix");
  m_colorVectorLocation = m_program->uniformLocation("color");
  //[~GET UNIFORM LOCATIONS]

  //[ MESH3D SETTINGS]
  // Spere
  m_sphere.getVertexDataFromFile("12by12UVsphere.obj");
  m_vao_sphere = new QOpenGLVertexArrayObject;
  ;
  m_vao_sphere->create();
  m_vbo_sphere = new QOpenGLBuffer;
  m_vbo_sphere->create();
  initMesh3D(m_vbo_sphere, m_vao_sphere, &m_sphere);
  // Cylinder
  m_cylinder.getVertexDataFromFile("12Cylinder.obj");
  m_vao_cylinder = new QOpenGLVertexArrayObject;
  ;
  m_vao_cylinder->create();
  m_vbo_cylinder = new QOpenGLBuffer;
  m_vbo_cylinder->create();
  initMesh3D(m_vbo_cylinder, m_vao_cylinder, &m_cylinder);
  //[~MESH3D SETTINGS]

  //[ OPENGL SETTINGS]
  f->glEnable(GL_DEPTH_TEST);
  f->glEnable(GL_CULL_FACE);
  f->glClearColor(1, 1, 1, 1);
  //[~OPENGL SETTINGS]
}

void MoleculeViewWidget::resizeGL(int w, int h) {}

void MoleculeViewWidget::paintGL() {
  // Now use QOpenGLExtraFunctions instead of QOpenGLFunctions as we want to
  // do more than what GL(ES) 2.0 offers.
  QOpenGLExtraFunctions *f = QOpenGLContext::currentContext()->extraFunctions();

  f->glClearColor(1, 1, 1, 1);
  f->glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  m_program->bind();

  //[ SET MVP] set values for model view transform matrices
  // Model matrix - responsible for translation, rotation, scale of object that
  // render View matrix - responsible for camera position and direction (set by
  // lookAt function) Projection matrix - responsible for field of view and
  // clipping planes (determine how far and close we can see)(set by perspective
  // function)
  m_modelMatrix.setToIdentity();
  m_viewMatrix.setToIdentity();
  m_viewMatrix.lookAt(m_cameraPosition, m_cameraPosition + m_cameraDirection,
                      QVector3D(0, 1, 0));
  m_projectionMatrix.setToIdentity();
  m_projectionMatrix.perspective(
      45.0f, GLfloat(this->size().width()) / this->size().height(), 0.01f,
      100.0f);
  m_normalMatrix.setToIdentity();
  static double rotationAngle = 0;
  rotationAngle += 2;
  // m_rotationAroundOriginMatrix.setToIdentity();
  // m_rotationAroundOriginMatrix.rotate(rotationAngle,QVector3D(0,1,0));
  m_lightPositionVector = QVector3D(100, 100, 300);

  m_program->setUniformValue(m_modelMatrixLocation, m_modelMatrix);
  m_program->setUniformValue(m_viewMatrixLocation, m_viewMatrix);
  m_program->setUniformValue(m_projectionMatrixLocation, m_projectionMatrix);
  m_program->setUniformValue(m_normalMatrixLocation, m_normalMatrix);
  m_program->setUniformValue(m_lightPositionVectorLocation,
                             m_lightPositionVector);
  m_program->setUniformValue(m_rotationAroundOriginMatrixLocation,
                             m_rotationAroundOriginMatrix);
  //[~SET MVP]

  if (sharedData->atlasNodeWasSet) {
    if (currentNodeID != sharedData->currentAtlasNode->getID()) {
      setAtlasNode(sharedData->currentAtlasNode);
    }
    renderAtoms();
    renderBonds();
  }
}

void MoleculeViewWidget::initMesh3D(QOpenGLBuffer *VBO,
                                    QOpenGLVertexArrayObject *VAO,
                                    Mesh3D *mesh3D) {
  QOpenGLFunctions *f = QOpenGLContext::currentContext()->functions();
  VAO->bind();
  VBO->bind();
  VBO->allocate(mesh3D->constData(), mesh3D->count() * sizeof(GLfloat));
  f->glEnableVertexAttribArray(0);
  f->glEnableVertexAttribArray(1);
  f->glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(GLfloat), 0);
  f->glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(GLfloat),
                           reinterpret_cast<void *>(3 * sizeof(GLfloat)));
  VBO->release();
  VAO->release();
}

void MoleculeViewWidget::renderSpere(QVector3D position, float radious,
                                     QVector4D color) {
  QOpenGLExtraFunctions *f = QOpenGLContext::currentContext()->extraFunctions();
  m_vao_sphere->bind();
  m_vbo_sphere->bind();

  m_modelMatrix.setToIdentity();
  m_modelMatrix.translate(position);
  m_modelMatrix.scale(QVector3D(radious, radious, radious));
  m_program->setUniformValue(m_modelMatrixLocation, m_modelMatrix);
  m_program->setUniformValue(m_colorVectorLocation, color);

  f->glCullFace(GL_FRONT);
  f->glDrawArrays(GL_TRIANGLES, 0, m_sphere.vertexCount());
  f->glCullFace(GL_BACK);
  f->glDrawArrays(GL_TRIANGLES, 0, m_sphere.vertexCount());

  m_vbo_sphere->release();
  m_vao_sphere->release();
}

void MoleculeViewWidget::renderCylinder(QVector3D positionStart,
                                        QVector3D positionEnd, float radious,
                                        QVector4D color) {
  //[ FIND ROTATION]
  QVector3D directionNotNormalized = QVector3D(positionEnd - positionStart);
  QVector3D direction = QVector3D(positionEnd - positionStart);
  QVector3D directionProjection = QVector3D(direction.x(), 0, direction.z());
  direction.normalize();
  QVector3D axisOfRotation =
      QVector3D::crossProduct(directionProjection, direction);
  float angle = acos(QVector3D::dotProduct(direction, directionProjection) /
                     (direction.length() * directionProjection.length()));
  angle = angle * 180 / PI + 90;
  //[~FIND ROTATION]
  QOpenGLExtraFunctions *f = QOpenGLContext::currentContext()->extraFunctions();
  m_vao_cylinder->bind();
  m_vbo_cylinder->bind();

  m_modelMatrix.setToIdentity();
  m_modelMatrix.translate(positionStart + directionNotNormalized / 2);
  //    static double addAngle=0;
  m_modelMatrix.rotate(angle, axisOfRotation);
  //    cout<<"Add Angle = "<<addAngle<<endl;
  m_modelMatrix.scale(
      QVector3D(radious, directionNotNormalized.length() / 2, radious));
  m_program->setUniformValue(m_modelMatrixLocation, m_modelMatrix);
  m_program->setUniformValue(m_colorVectorLocation, color);

  // Cylinder is alway rendered on tope of everything
  glDepthFunc(GL_ALWAYS);
  f->glCullFace(GL_FRONT);
  f->glDrawArrays(GL_TRIANGLES, 0, m_cylinder.vertexCount());
  f->glCullFace(GL_BACK);
  f->glDrawArrays(GL_TRIANGLES, 0, m_cylinder.vertexCount());
  glDepthFunc(GL_LESS);

  m_vbo_cylinder->release();
  m_vao_cylinder->release();
}

void MoleculeViewWidget::renderAtoms() {
  for (int i = 0; i < atomsA.size(); i++) {
    double *atomPosition = atomsA[i]->getLocation();
    double atomRadious = atomsA[i]->getRadius();
    renderSpere(QVector3D(atomPosition[0], atomPosition[1], atomPosition[2]),
                atomRadious, QVector4D(1, 0, 0, 1));
  }
  for (int i = 0; i < atomsB.size(); i++) {
    double *atomPosition = atomsB[i]->getLocation();
    double atomRadious = atomsB[i]->getRadius();
    renderSpere(QVector3D(atomPosition[0], atomPosition[1], atomPosition[2]),
                atomRadious, QVector4D(1, 1, 0, 1));
  }
}

void MoleculeViewWidget::renderBonds() {
  for (int i = 0; i < bondSolid.size(); i++) {
    double *bondStart, *bondEnd;
    bondStart = atomsA[bondSolid[i].second]->getLocation();
    bondEnd = atomsB[bondSolid[i].first]->getLocation();
    renderCylinder(QVector3D(bondStart[0], bondStart[1], bondStart[2]),
                   QVector3D(bondEnd[0], bondEnd[1], bondEnd[2]), 0.1,
                   QVector4D(0.1, 0.1, 0.1, 1));
  }
  for (int i = 0; i < bondVariable.size(); i++) {
    double *bondStart, *bondEnd;
    bondStart = atomsA[bondVariable[i].second]->getLocation();
    bondEnd = atomsB[bondVariable[i].first]->getLocation();
    renderCylinder(QVector3D(bondStart[0], bondStart[1], bondStart[2]),
                   QVector3D(bondEnd[0], bondEnd[1], bondEnd[2]), 0.1,
                   QVector4D(0.7, 0.7, 0.7, 1));
  }
}

void MoleculeViewWidget::updateSpaceAndFlip() {
  double fb[3][3], tb[3][3];
  currentACR->space[sharedData->currentSpaceID]
      ->getOrientations()[sharedData->currentFlipID]
      ->getFromTo(fb, tb);
  vector<double> trans_b = Utils::getTransMatrix(fb, tb);
  moleculeBTransform = QMatrix4x4(
      trans_b[0], trans_b[1], trans_b[2], trans_b[3], trans_b[4], trans_b[5],
      trans_b[6], trans_b[7], trans_b[8], trans_b[9], trans_b[10], trans_b[11],
      trans_b[12], trans_b[13], trans_b[14], trans_b[15]);

  atomsB = molecularUnitB->getXFAtoms(
      trans_b);  // be sure to properly delete XFAtoms2
  atomsA = molecularUnitA->getAtoms();
  bondSolid = atlasNode->getCG()->getParticipants();
  bondVariable = atlasNode->getCG()->getParamLines();
  atlasNodeIsSet = true;
  // Set bondData for bondViewWidget
  //    AtomBondData bondData;
  //    bondData.bondSolid = this->bondSolid;
  //    bondData.bondVariable = this->bondVariable;
  //    bondViewWidget->setBondData(&bondData);
  // sharedData->currentSpaceID = i;
  // sharedData->currentFlipID = orientations[0]->getFlipNum();
}

void MoleculeViewWidget::setAtlasNode(AtlasNode *atlasNode) {
  currentNodeID = atlasNode->getID();
  delete currentACR;
  currentACR = new ActiveConstraintRegion();
  saveLoader->loadNode(atlasNode->getID(), currentACR);
  sharedData->currentACR = currentACR;
  sharedData->ACRWasSet = true;
  if (currentACR->space.size() > 0) {
    // Iterate through this spaces to find the one that has orientation
    for (int i = 0; i < currentACR->space.size(); i++) {
      // Get first available orientation and return
      if (currentACR->space[i]->getOrientations().size() > 0) {
        //                this->atlasNode = currentAtlasNode;

        std::vector<Orientation *> orientations =
            currentACR->space[i]->getOrientations();
        double fb[3][3], tb[3][3];
        orientations[0]->getFromTo(fb, tb);
        vector<double> trans_b = Utils::getTransMatrix(fb, tb);
        moleculeBTransform =
            QMatrix4x4(trans_b[0], trans_b[1], trans_b[2], trans_b[3],
                       trans_b[4], trans_b[5], trans_b[6], trans_b[7],
                       trans_b[8], trans_b[9], trans_b[10], trans_b[11],
                       trans_b[12], trans_b[13], trans_b[14], trans_b[15]);

        atomsB = molecularUnitB->getXFAtoms(
            trans_b);  // be sure to properly delete XFAtoms2
        atomsA = molecularUnitA->getAtoms();
        bondSolid = atlasNode->getCG()->getParticipants();
        bondVariable = atlasNode->getCG()->getParamLines();
        atlasNodeIsSet = true;
        // Set bondData for bondViewWidget
        AtomBondData bondData;
        bondData.bondSolid = this->bondSolid;
        bondData.bondVariable = this->bondVariable;
        bondViewWidget->setBondData(&bondData);
        sharedData->currentSpaceID = i;
        sharedData->currentFlipID = orientations[0]->getFlipNum();
        return;
      }
    }
  }
}

void MoleculeViewWidget::setMolecularUnits(MolecularUnit *molecularUnitA,
                                           MolecularUnit *molecularUnitB) {
  //[ ATLAS NODE]
  atlasNodeIsSet = false;
  this->molecularUnitA = molecularUnitA;
  this->molecularUnitB = molecularUnitB;
  //[~ATLAS NODE]
}

void MoleculeViewWidget::keyPressEvent(QKeyEvent *keyEvent) {
  bool D = true;
  if (keyEvent->key() == Qt::Key_Q) {
    m_cameraPosition.setZ(m_cameraPosition.z() + m_cameraMoveSpeed);
  }
  if (keyEvent->key() == Qt::Key_E) {
    m_cameraPosition.setZ(m_cameraPosition.z() - m_cameraMoveSpeed);
  }
  if (keyEvent->key() == Qt::Key_W) {
    m_cameraPosition.setY(m_cameraPosition.y() + m_cameraMoveSpeed);
  }
  if (keyEvent->key() == Qt::Key_S) {
    m_cameraPosition.setY(m_cameraPosition.y() - m_cameraMoveSpeed);
  }
  if (keyEvent->key() == Qt::Key_A) {
    m_cameraPosition.setX(m_cameraPosition.x() - m_cameraMoveSpeed);
  }
  if (keyEvent->key() == Qt::Key_D) {
    m_cameraPosition.setX(m_cameraPosition.x() + m_cameraMoveSpeed);
  }
  if (D)
    cout << "Camera Position: " << m_cameraPosition.x() << " "
         << m_cameraPosition.y() << " " << m_cameraPosition.z() << endl;
}

void MoleculeViewWidget::mousePressEvent(QMouseEvent *mouseEvent) {
  bool D = true;
  if (D) cout << "Mouse Was Pressed" << endl;
  if (mouseEvent->button() == Qt::LeftButton) {
    m_mouseIsPressed = true;
    m_mouseLastPosition = mouseEvent->pos();
  }
}

void MoleculeViewWidget::mouseReleaseEvent(QMouseEvent *mouseEvent) {
  bool D = true;
  if (D) cout << "Mouse Was Released" << endl;
  if (mouseEvent->button() == Qt::LeftButton) {
    m_mouseIsPressed = false;
  }
}

void MoleculeViewWidget::mouseMoveEvent(QMouseEvent *mouseEvent) {
  bool D = false;
  if (m_mouseIsPressed) {
    //        //[ HANDLE MOUSE]Get change in x and y coords of mouse since last
    //        call of this function QPoint mouseDelta = m_mouseLastPosition -
    //        mouseEvent->pos(); m_mouseLastPosition = mouseEvent->pos();
    //        //[~HANDLE MOUSE]
    //        //[ MAKE QUATERNION]
    //        // Get rotation according to X (roll) and Y axis (yaw) using mouse
    //        data double angleYaw = ((double) mouseDelta.x() /
    //        (double)this->size().width())*360*m_cameraRotationalSpeed; double
    //        anglePitch = ((double) mouseDelta.y() /
    //        (double)this->size().height())*360*m_cameraRotationalSpeed; double
    //        angleRoll = 0;
    //        // Create quaternion correnspondent to the rotation
    //        QQuaternion quaternionCameraRotation;
    //        quaternionCameraRotation =
    //        QQuaternion::fromEulerAngles(anglePitch,angleYaw,angleRoll);
    //        //[~MAKE QUATERNION]
    //        //[ ROTATE USING MATRIX]
    //        // Get rotational matrix4x4
    //        QMatrix4x4 cameraDirectionRotationMatrix;
    //        cameraDirectionRotationMatrix.setToIdentity();
    //        cameraDirectionRotationMatrix.rotate(quaternionCameraRotation);
    //        // Get new camera direction by multiplying rotational matrxi by
    //        old camera direction vector m_cameraDirection =
    //        cameraDirectionRotationMatrix * m_cameraDirection;
    //        m_cameraDirection.normalize();
    //        //[~ROTATE USING MATRIX]
    //        if(D)cout<<"Camera direction "<<m_cameraDirection.x()<<"
    //        "<<m_cameraDirection.y()<<" "<<m_cameraDirection.z()<<endl;

    //[ HANDLE MOUSE]Get change in x and y coords of mouse since last call of
    //this function
    QPoint mouseDelta = m_mouseLastPosition - mouseEvent->pos();
    m_mouseLastPosition = mouseEvent->pos();
    //[~HANDLE MOUSE]
    //[ MAKE QUATERNION]
    // Get rotation according to X (roll) and Y axis (yaw) using mouse data
    angleYaw += ((double)mouseDelta.x() / (double)this->size().width()) * 360 *
                m_cameraRotationalSpeed;
    anglePitch += ((double)mouseDelta.y() / (double)this->size().height()) *
                  360 * m_cameraRotationalSpeed;
    double angleRoll = 0;
    // Create quaternion correnspondent to the rotation
    QQuaternion quaternionCameraRotation;
    quaternionCameraRotation =
        QQuaternion::fromEulerAngles(anglePitch, angleYaw, angleRoll);
    //[~MAKE QUATERNION]
    //[ ROTATE USING MATRIX]
    m_rotationAroundOriginMatrix.setToIdentity();
    m_rotationAroundOriginMatrix.rotate(quaternionCameraRotation);
  }
}
