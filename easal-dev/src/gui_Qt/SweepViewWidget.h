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
#ifndef SWEEPVIEWWIDGET_H
#define SWEEPVIEWWIDGET_H

#include <easalcore/AtlasNode.h>
#include <easalcore/SaveLoader.h>
#include <easalcore/Settings.h>
#include <easalcore/Utils.h>
#include <gui_Qt/BondViewWidget.h>
#include <gui_Qt/Camera.h>
#include <gui_Qt/Mesh3D.h>
#include <gui_Qt/Renderer.h>
#include <gui_Qt/SharedDataGUI.h>

#include <QCheckBox>
#include <QGridLayout>
#include <QKeyEvent>
#include <QMouseEvent>
#include <QOpenGLBuffer>
#include <QOpenGLContext>
#include <QOpenGLExtraFunctions>
#include <QOpenGLShaderProgram>
#include <QOpenGLVertexArrayObject>
#include <QOpenGLWidget>
#include <QPushButton>
#include <QTimer>
#include <cmath>
#include <iostream>

#define PI 3.14159265

// debug includes
#include "ctime"

using namespace std;

class SweepViewWidget : public QOpenGLWidget {
  Q_OBJECT
 public:
  SweepViewWidget(SharedDataGUI *sharedData);
  /** @brief to render some mlecules we we to load nodes using save loader
   * but save loader is not included in constructor o we need to set it in
   * setSaveLoader */
  void setSaveLoader(SaveLoader *saveLoader);
  void setAtlasNode(AtlasNode *atlasNode);
  void setMolecularUnits(MolecularUnit *molecularUnitA,
                         MolecularUnit *molecularUnitB);
  void initializeGL() override;
  void paintGL() override;
  void resizeGL(int w, int h) override;
 private Q_SLOTS:
  void playButtonSlot();
  void stopButtonSlot();
  void nextButtonSlot();
  void previousButtonSlot();
  void updateNextAtomB();
  void updateTransformationMatrices();

 private:
  //[ COMPONENTS]
  Camera camera;
  Renderer renderer;
  //[~COMPONENTS]
  //[ DEPENDENCIES]
  SaveLoader *saveLoader;
  SharedDataGUI *sharedData;
  MolecularUnit *molecularUnitA;
  MolecularUnit *molecularUnitB;
  //[~DEPENDENCIES]
  //[ TIMER]
  /** @brief Timer that triggers screen updates once each m_updateTimeInMS*/
  QTimer m_timer;
  double m_updateTimeInMS;
  //[ TIMER]
  //[ UNIFORMS]
  vector<QMatrix4x4> translationMatricesSphere;
  vector<QMatrix4x4> rotationMatricesSphere;
  vector<QMatrix4x4> scaleMatricesSphere;
  vector<QVector4D> colorVectorsSphere;
  //[ UNIFORMS]
  //[ MESHES3D]
  Mesh3D sphere;
  //[~MESHES3D]
  //[ ATLAS NODE]
  bool atlasNodeIsSet;
  AtlasNode *atlasNode;
  QMatrix4x4 moleculeARotation;
  QMatrix4x4 moleculeBTransform;
  vector<Atom *> atomsA;
  vector<vector<Atom *>> atomsB[8];
  ActiveConstraintRegion *currentACR = 0;
  bool noGoodOrientation;
  int currentNodeID = -1;
  //[~ATLAS NODE]
  //[ GUI RELATED ELEMENTS]
  QCheckBox checkBoxes[8];
  QPushButton playButton;
  QPushButton stopButton;
  QPushButton previousButton;
  QPushButton nextButton;
  bool playOrientationsOneByOne;
  bool stopPlayingOrientations;
  bool nextOrientation;
  bool previousOrientation;
  QTimer orientationTimer;
  double msPerOrientation;
  QVector2D currentAtomBIndicies;
  //[~GUI RELATED ELEMENTS]
  //[ MISC]
  QVector4D colors[8];
  QGridLayout gridLayout;
  //[~MISC]

  void clearTarnsformationMatrices();
  void keyPressEvent(QKeyEvent *keyEvent) Q_DECL_OVERRIDE;
  void mousePressEvent(QMouseEvent *mouseEvent) Q_DECL_OVERRIDE;
  void mouseReleaseEvent(QMouseEvent *mouseEvent) Q_DECL_OVERRIDE;
  void wheelEvent(QWheelEvent *wheelEvent) Q_DECL_OVERRIDE;
};

#endif  // SWEEPVIEWWIDGET_H
