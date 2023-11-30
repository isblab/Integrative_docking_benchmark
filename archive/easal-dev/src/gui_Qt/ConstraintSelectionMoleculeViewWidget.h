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
#ifndef CONSTRAINTSELECTIONMOLECULEVIEWWIDGET_H
#define CONSTRAINTSELECTIONMOLECULEVIEWWIDGET_H
#include <easalcore/MolecularUnit.h>

#include <QtWidgets>
#include <vector>

#include "Camera.h"
#include "Mesh3D.h"
#include "Renderer.h"

#define MOLECULAR_UNITS_NUMBER 2

enum MolecularUnitType {
  NONE = -1,
  MOLECULAR_UNIT_A = 0,
  MOLECULAR_UNIT_B = 1
};

using namespace std;

class ConstraintSelectionMoleculeViewWidget : public QOpenGLWidget {
  Q_OBJECT
 public:
  ConstraintSelectionMoleculeViewWidget();
  void setData(MolecularUnit *molecularUnitA, MolecularUnit *molecularUnitB);
  void initializeGL() override;
  void paintGL() override;
  void resizeGL(int w, int h) override;

  vector<QMatrix4x4> translationMatrices[MOLECULAR_UNITS_NUMBER];
  vector<QMatrix4x4> rotationMatrices[MOLECULAR_UNITS_NUMBER];
  vector<QMatrix4x4> scaleMatrices[MOLECULAR_UNITS_NUMBER];
  vector<QVector4D> colors[MOLECULAR_UNITS_NUMBER];

  Camera camera;

 private:
  void updateTransformationMatrices();

  void keyPressEvent(QKeyEvent *event) Q_DECL_OVERRIDE;
  void mousePressEvent(QMouseEvent *mouseEvent) Q_DECL_OVERRIDE;
  void mouseReleaseEvent(QMouseEvent *mouseEvent) Q_DECL_OVERRIDE;
  void wheelEvent(QWheelEvent *wheelEvent) Q_DECL_OVERRIDE;

  MolecularUnit *molecularUnitA;
  MolecularUnit *molecularUnitB;

  vector<Atom *> atomsA;
  vector<Atom *> atomsB;

  Renderer renderer;
  Mesh3D sphere;
  //[ TIMER]
  /** @brief Timer that triggers screen updates once each m_updateTimeInMS*/
  QTimer m_timer;
  double m_updateTimeInMS;
  //[ TIMER]
  QVector3D molecularUnitOffset3D[2] = {QVector3D(-5, 0, 0),
                                        QVector3D(5, 0, 0)};
  QVector3D newMolecularUnitOffset3D[2] = {
      QVector3D(-5, 0, 0),
      QVector3D(
          5, 0,
          0)};  ///< we need it to translate molec,ules using keyborad (WASD)

  MolecularUnitType mouseRelatesToMolecularUnit = NONE;
  bool LMBisPressed = false;
  bool wheelWasRotated = false;
  int wheelDelta = 0;
  QPoint lastMouseCoordinates;
  double rotationSpeed = 360;
  double scaleSpeed = 1.1;
  double translationSpeed = 0.5;
};

#endif  // CONSTRAINTSELECTIONMOLECULEVIEWWIDGET_H
