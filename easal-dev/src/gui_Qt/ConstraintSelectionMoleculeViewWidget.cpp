/*
 This file is part of EASAL.

 EASAL is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 EASAL is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/*
 *  Created on: 2016-2017
 *      Author: Chkhaidze Giorgio
 */
#include "ConstraintSelectionMoleculeViewWidget.h"

ConstraintSelectionMoleculeViewWidget::ConstraintSelectionMoleculeViewWidget() {
  //[ CAMERA SETTINGS]
  camera.init(this->size(), QVector3D(0, 0, 0), QVector3D(0, 0, -1));
  //[~CAMERA SETTINGS]
  //[ TIMER SETTINGS]
  m_updateTimeInMS = 33;
  m_timer.start(m_updateTimeInMS);
  connect(&m_timer, SIGNAL(timeout()), this, SLOT(update()));
  //[ TIMER SETTINGS]
  setMinimumSize(600, 400);
  setFocusPolicy(Qt::StrongFocus);
}

void ConstraintSelectionMoleculeViewWidget::setData(
    MolecularUnit *molecularUnitA, MolecularUnit *molecularUnitB) {
  for (int i = 0; i < MOLECULAR_UNITS_NUMBER; i++) {
    translationMatrices[i].clear();
    rotationMatrices[i].clear();
    scaleMatrices[i].clear();
    colors[i].clear();
  }
  atomsA = molecularUnitA->getAtoms();
  atomsB = molecularUnitB->getAtoms();
  for (int i = 0; i < atomsA.size(); i++) {
    double *atomPosition = atomsA[i]->getLocation();
    double atomRadious = atomsA[i]->getRadius();

    QMatrix4x4 translationTmp;
    QMatrix4x4 rotationTmp;
    QMatrix4x4 scaleTmp;
    QVector4D colorTmp;

    translationTmp.translate(
        QVector3D(atomPosition[0], atomPosition[1], atomPosition[2]) +
        molecularUnitOffset3D[MOLECULAR_UNIT_A]);
    rotationTmp.setToIdentity();
    scaleTmp.scale(atomRadious);
    colorTmp = QVector4D(1, 0, 0, 1);

    translationMatrices[MOLECULAR_UNIT_A].push_back(translationTmp);
    rotationMatrices[MOLECULAR_UNIT_A].push_back(rotationTmp);
    scaleMatrices[MOLECULAR_UNIT_A].push_back(scaleTmp);
    colors[MOLECULAR_UNIT_A].push_back(colorTmp);
  }

  for (int i = 0; i < atomsB.size(); i++) {
    double *atomPosition = atomsB[i]->getLocation();
    double atomRadious = atomsB[i]->getRadius();

    QMatrix4x4 translationTmp;
    QMatrix4x4 rotationTmp;
    QMatrix4x4 scaleTmp;
    QVector4D colorTmp;

    translationTmp.translate(
        QVector3D(atomPosition[0], atomPosition[1], atomPosition[2]) +
        molecularUnitOffset3D[MOLECULAR_UNIT_B]);
    rotationTmp.setToIdentity();
    scaleTmp.scale(atomRadious);
    colorTmp = QVector4D(1, 1, 0, 1);

    translationMatrices[MOLECULAR_UNIT_B].push_back(translationTmp);
    rotationMatrices[MOLECULAR_UNIT_B].push_back(rotationTmp);
    scaleMatrices[MOLECULAR_UNIT_B].push_back(scaleTmp);
    colors[MOLECULAR_UNIT_B].push_back(colorTmp);
  }
}

void ConstraintSelectionMoleculeViewWidget::initializeGL() {
  renderer.initOpenGL();
  sphere.createAndBind("12by12UVsphere.obj");
}

void ConstraintSelectionMoleculeViewWidget::paintGL() {
  if (LMBisPressed || wheelWasRotated) {
    updateTransformationMatrices();
  }

  renderer.renderingProgram->setUniformValue(
      "lightPositionVector",
      camera.getEyeLocation() + camera.getRight() * 10 + camera.getUp() * 30);
  renderer.renderingProgram->setUniformValue(
      "cameraMatrix", camera.getCameraMatrix(mapFromGlobal(cursor().pos())));

  renderer.bindRenderProgramAndClear();
  // Render spheres that represent atoms
  renderer.render(&sphere, translationMatrices[MOLECULAR_UNIT_A],
                  rotationMatrices[MOLECULAR_UNIT_A],
                  scaleMatrices[MOLECULAR_UNIT_A], colors[MOLECULAR_UNIT_A]);
  renderer.render(&sphere, translationMatrices[MOLECULAR_UNIT_B],
                  rotationMatrices[MOLECULAR_UNIT_B],
                  scaleMatrices[MOLECULAR_UNIT_B], colors[MOLECULAR_UNIT_B]);
}

void ConstraintSelectionMoleculeViewWidget::resizeGL(int w, int h) {
  camera.resizeEvent(w, h);
}

void ConstraintSelectionMoleculeViewWidget::updateTransformationMatrices() {
  QPoint mouseCoordinates = mapFromGlobal(cursor().pos());
  QPoint mouseDelta = mouseCoordinates - lastMouseCoordinates;
  QVector2D rotationAngle =
      QVector2D((double)mouseDelta.x() / width() * rotationSpeed,
                (double)mouseDelta.y() / height() * rotationSpeed);

  double anglePitch = rotationAngle.y();  // around X
  double angleYaw = rotationAngle.x();    // around Y
  double angleRoll = 0;                   // around Z
  if (mouseRelatesToMolecularUnit != NONE) {
    for (int i = 0; i < translationMatrices[mouseRelatesToMolecularUnit].size();
         i++) {
      QMatrix4x4 rotationAroundOriginMatrix;
      QMatrix4x4 translationTmp;
      QVector4D tmpPointCoords(
          translationMatrices[mouseRelatesToMolecularUnit][i].data()[12],
          translationMatrices[mouseRelatesToMolecularUnit][i].data()[13],
          translationMatrices[mouseRelatesToMolecularUnit][i].data()[14],
          translationMatrices[mouseRelatesToMolecularUnit][i].data()[15]);
      // Translate to origine
      tmpPointCoords =
          tmpPointCoords -
          QVector4D(molecularUnitOffset3D[mouseRelatesToMolecularUnit]);
      rotationAroundOriginMatrix.rotate(anglePitch, QVector3D(1, 0, 0));
      rotationAroundOriginMatrix.rotate(angleYaw, QVector3D(0, 1, 0));

      if (wheelDelta > 0) {
        rotationAroundOriginMatrix.scale((wheelDelta / 120.0) * scaleSpeed);
        scaleMatrices[mouseRelatesToMolecularUnit][i].scale(
            (wheelDelta / 120.0) * scaleSpeed);
      } else if (wheelDelta < 0) {
        rotationAroundOriginMatrix.scale(1.0 /
                                         ((-wheelDelta / 120.0) * scaleSpeed));
        scaleMatrices[mouseRelatesToMolecularUnit][i].scale(
            1.0 / ((-wheelDelta / 120.0) * scaleSpeed));
      }

      tmpPointCoords = tmpPointCoords * rotationAroundOriginMatrix;
      // Translate back
      tmpPointCoords =
          tmpPointCoords +
          QVector4D(newMolecularUnitOffset3D[mouseRelatesToMolecularUnit]);

      translationMatrices[mouseRelatesToMolecularUnit][i].setToIdentity();
      translationMatrices[mouseRelatesToMolecularUnit][i].translate(
          tmpPointCoords.x(), tmpPointCoords.y(), tmpPointCoords.z());
    }
    molecularUnitOffset3D[mouseRelatesToMolecularUnit] =
        newMolecularUnitOffset3D[mouseRelatesToMolecularUnit];
  }

  lastMouseCoordinates = mapFromGlobal(cursor().pos());
  wheelWasRotated = false;
  wheelDelta = 0;
}

void ConstraintSelectionMoleculeViewWidget::keyPressEvent(QKeyEvent *event) {
  if (mouseRelatesToMolecularUnit != NONE) {
    if (event->key() == Qt::Key_A) {
      newMolecularUnitOffset3D[mouseRelatesToMolecularUnit].setX(
          newMolecularUnitOffset3D[mouseRelatesToMolecularUnit].x() -
          translationSpeed);
    }
    if (event->key() == Qt::Key_D) {
      newMolecularUnitOffset3D[mouseRelatesToMolecularUnit].setX(
          newMolecularUnitOffset3D[mouseRelatesToMolecularUnit].x() +
          translationSpeed);
    }
    if (event->key() == Qt::Key_S) {
      newMolecularUnitOffset3D[mouseRelatesToMolecularUnit].setY(
          newMolecularUnitOffset3D[mouseRelatesToMolecularUnit].y() -
          translationSpeed);
    }
    if (event->key() == Qt::Key_W) {
      newMolecularUnitOffset3D[mouseRelatesToMolecularUnit].setY(
          newMolecularUnitOffset3D[mouseRelatesToMolecularUnit].y() +
          translationSpeed);
    }
    if (event->key() == Qt::Key_C) {
      if (mouseRelatesToMolecularUnit == MOLECULAR_UNIT_A) {
        newMolecularUnitOffset3D[mouseRelatesToMolecularUnit] =
            QVector3D(-5, 0, 0);
      } else if (mouseRelatesToMolecularUnit == MOLECULAR_UNIT_B) {
        newMolecularUnitOffset3D[mouseRelatesToMolecularUnit] =
            QVector3D(5, 0, 0);
      }
    }
    updateTransformationMatrices();
  }
}

void ConstraintSelectionMoleculeViewWidget::mousePressEvent(
    QMouseEvent *mouseEvent) {
  if (mouseEvent->button() == Qt::LeftButton) {
    if (mapFromGlobal(cursor().pos()).x() < this->width() / 2) {
      mouseRelatesToMolecularUnit = MOLECULAR_UNIT_A;
    } else {
      mouseRelatesToMolecularUnit = MOLECULAR_UNIT_B;
    }
    LMBisPressed = true;
    lastMouseCoordinates = mapFromGlobal(cursor().pos());
  }
}

void ConstraintSelectionMoleculeViewWidget::mouseReleaseEvent(
    QMouseEvent *mouseEvent) {
  if (mouseEvent->button() == Qt::LeftButton) {
    LMBisPressed = false;
  }
}

void ConstraintSelectionMoleculeViewWidget::wheelEvent(
    QWheelEvent *wheelEvent) {
  if (mapFromGlobal(cursor().pos()).x() < this->width() / 2) {
    mouseRelatesToMolecularUnit = MOLECULAR_UNIT_A;
  } else {
    mouseRelatesToMolecularUnit = MOLECULAR_UNIT_B;
  }
  wheelWasRotated = true;
  wheelDelta = wheelEvent->delta();
}
