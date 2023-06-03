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
#include "SweepViewWidget.h"

SweepViewWidget::SweepViewWidget(SharedDataGUI *sharedData) {
  this->sharedData = sharedData;
  //[ CAMERA SETTINGS]
  camera.init(this->size(), QVector3D(0, 0, 0), QVector3D(0, 0, -1));
  //[~CAMERA SETTINGS]
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
  //[ COLORS]
  colors[0] = QVector4D(0.75, 0, 0.75, 1);
  colors[1] = QVector4D(0.75, 0, 0.5, 1);
  colors[2] = QVector4D(0.75, 1, 0, 1);
  colors[3] = QVector4D(0.25, 1, 0.25, 1);

  colors[4] = QVector4D(0, 1, 0.75, 1);
  colors[5] = QVector4D(0, 1, 1, 1);
  colors[6] = QVector4D(0, 0.5, 0.25, 1);
  colors[7] = QVector4D(0, 0.5, 0.75, 1);
  //[~COLORS]
  //[ PLAY CONTROL VARS]
  playOrientationsOneByOne = false;
  stopPlayingOrientations = false;
  nextOrientation = false;
  previousOrientation = false;
  msPerOrientation = 1000;
  // Connect timer that will allow to play orientations one by one
  connect(&orientationTimer, SIGNAL(timeout()), this, SLOT(updateNextAtomB()));
  //[~PLAY CONTROL VARS]
  //[ MISC]
  noGoodOrientation = true;
  gridLayout.setAlignment(Qt::AlignBottom);
  for (int i = 0; i < 8; i++) {
    checkBoxes[i].setText("Flip " + QString::number(i));
    checkBoxes[i].setTristate(true);
    connect(&checkBoxes[i], SIGNAL(stateChanged(int)), this,
            SLOT(updateTransformationMatrices()));
    gridLayout.addWidget(&checkBoxes[i], 1, i, 1, 2);
  }
  playButton.setText("Play");
  connect(&playButton, SIGNAL(pressed()), this, SLOT(playButtonSlot()));
  gridLayout.addWidget(&playButton, 0, 1, 1, 1);
  stopButton.setText("Stop");
  connect(&stopButton, SIGNAL(pressed()), this, SLOT(stopButtonSlot()));
  gridLayout.addWidget(&stopButton, 0, 2, 1, 1);
  nextButton.setText("Next");
  connect(&nextButton, SIGNAL(pressed()), this, SLOT(nextButtonSlot()));
  gridLayout.addWidget(&nextButton, 0, 3, 1, 1);
  previousButton.setText("Prev");
  connect(&previousButton, SIGNAL(pressed()), this, SLOT(previousButtonSlot()));
  gridLayout.addWidget(&previousButton, 0, 0, 1, 1);
  setLayout(&gridLayout);
  //[~MISC]
}

void SweepViewWidget::setSaveLoader(SaveLoader *saveLoader) {
  this->saveLoader = saveLoader;
}

void SweepViewWidget::initializeGL() {
  renderer.initOpenGL();
  sphere.createAndBind("12by12UVsphere.obj");
}

void SweepViewWidget::paintGL() {
  if (sharedData->atlasNodeWasSet) {
    if (currentNodeID != sharedData->currentAtlasNode->getID()) {
      setAtlasNode(sharedData->currentAtlasNode);
    }
    renderer.renderingProgram->setUniformValue(
        "lightPositionVector",
        camera.getEyeLocation() + camera.getRight() * 10 + camera.getUp() * 30);
    renderer.renderingProgram->setUniformValue(
        "cameraMatrix", camera.getCameraMatrix(mapFromGlobal(cursor().pos())));

    renderer.bindRenderProgramAndClear();
    renderer.render(&sphere, translationMatricesSphere, rotationMatricesSphere,
                    scaleMatricesSphere, colorVectorsSphere);
  }
}

void SweepViewWidget::resizeGL(int w, int h) { camera.resizeEvent(w, h); }

void SweepViewWidget::playButtonSlot() {
  playOrientationsOneByOne = !playOrientationsOneByOne;
  if (playOrientationsOneByOne) {
    orientationTimer.start(msPerOrientation);
    // Find orientation to start (if non is available then assgne -1,-1 value
    currentAtomBIndicies = QVector2D(-1, -1);
    for (int i = 0; i < 8; i++) {
      if (checkBoxes[i].isEnabled() && checkBoxes[i].isChecked()) {
        currentAtomBIndicies = QVector2D(i, 0);
        break;
      }
    }
  } else {
    orientationTimer.stop();
  }
  updateTransformationMatrices();
}

void SweepViewWidget::stopButtonSlot() {
  if (playOrientationsOneByOne && !stopPlayingOrientations) {
    orientationTimer.stop();
    stopPlayingOrientations = true;
  } else if (playOrientationsOneByOne && stopPlayingOrientations) {
    orientationTimer.start(msPerOrientation);
    stopPlayingOrientations = false;
  }
  updateTransformationMatrices();
}

void SweepViewWidget::nextButtonSlot() {
  if (currentAtomBIndicies == QVector2D(-1, -1)) {
    return;
  }
  for (int i = 0; i < 9; i++) {
    int flipNumber = ((int)currentAtomBIndicies.x() + i) % 8;
    if (!(checkBoxes[flipNumber].isEnabled() &&
          checkBoxes[flipNumber].isChecked())) {
      continue;
    }
    if (i == 0) {
      for (int j = currentAtomBIndicies.y() + 1; j < atomsB[flipNumber].size();
           j++) {
        currentAtomBIndicies = QVector2D(flipNumber, j);
        updateTransformationMatrices();
        return;
      }
    } else {
      for (int j = 0; j < atomsB[flipNumber].size(); j++) {
        currentAtomBIndicies = QVector2D(flipNumber, j);
        updateTransformationMatrices();
        return;
      }
    }
  }
  updateTransformationMatrices();
}

void SweepViewWidget::previousButtonSlot() {
  if (currentAtomBIndicies == QVector2D(-1, -1)) {
    return;
  }
  for (int i = 0; i < 9; i++) {
    int flipNumber = ((int)currentAtomBIndicies.x() - i + 16) % 8;
    if (!(checkBoxes[flipNumber].isEnabled() &&
          checkBoxes[flipNumber].isChecked())) {
      continue;
    }
    if (i == 0) {
      for (int j = currentAtomBIndicies.y() - 1; j >= 0; j--) {
        currentAtomBIndicies = QVector2D(flipNumber, j);
        updateTransformationMatrices();
        return;
      }
    } else {
      for (int j = atomsB[flipNumber].size() - 1; j >= 0; j--) {
        currentAtomBIndicies = QVector2D(flipNumber, j);
        updateTransformationMatrices();
        return;
      }
    }
  }
  updateTransformationMatrices();
}

void SweepViewWidget::updateNextAtomB() {
  if (currentAtomBIndicies == QVector2D(-1, -1)) {
    return;
  }
  for (int i = 0; i < 9; i++) {
    int flipNumber = ((int)currentAtomBIndicies.x() + i) % 8;
    if (!(checkBoxes[flipNumber].isEnabled() &&
          checkBoxes[flipNumber].isChecked())) {
      continue;
    }
    if (i == 0) {
      for (int j = currentAtomBIndicies.y() + 1; j < atomsB[flipNumber].size();
           j++) {
        currentAtomBIndicies = QVector2D(flipNumber, j);
        updateTransformationMatrices();
        return;
      }
    } else {
      for (int j = 0; j < atomsB[flipNumber].size(); j++) {
        currentAtomBIndicies = QVector2D(flipNumber, j);
        updateTransformationMatrices();
        return;
      }
    }
  }
  updateTransformationMatrices();
}

void SweepViewWidget::updateTransformationMatrices() {
  clearTarnsformationMatrices();

  for (int i = 0; i < atomsA.size(); i++) {
    double *atomPosition = atomsA[i]->getLocation();
    double atomRadious = atomsA[i]->getRadius();

    QMatrix4x4 translationTmp;
    QMatrix4x4 rotationTmp;
    QMatrix4x4 scaleTmp;
    QVector4D colorTmp;

    translationTmp.translate(
        QVector3D(atomPosition[0], atomPosition[1], atomPosition[2]));
    rotationTmp.setToIdentity();
    scaleTmp.scale(atomRadious);
    colorTmp = QVector4D(1, 0, 0, 1);

    translationMatricesSphere.push_back(translationTmp);
    rotationMatricesSphere.push_back(rotationTmp);
    scaleMatricesSphere.push_back(scaleTmp);
    colorVectorsSphere.push_back(colorTmp);
  }

  // Render all orientations with flips that user selected using checkboxes
  if (!playOrientationsOneByOne) {
    for (int h = 0; h < 8; h++) {
      if (checkBoxes[h].isEnabled() && checkBoxes[h].isChecked()) {
        for (int i = 0; i < atomsB[h].size(); i++) {
          for (int j = 0; j < atomsB[h][i].size(); j++) {
            double *atomPosition = atomsB[h][i][j]->getLocation();
            double atomRadious = atomsB[h][i][j]->getRadius();

            QMatrix4x4 translationTmp;
            QMatrix4x4 rotationTmp;
            QMatrix4x4 scaleTmp;
            QVector4D colorTmp;

            translationTmp.translate(
                QVector3D(atomPosition[0], atomPosition[1], atomPosition[2]));
            rotationTmp.setToIdentity();
            scaleTmp.scale(atomRadious);
            colorTmp = colors[h];

            translationMatricesSphere.push_back(translationTmp);
            rotationMatricesSphere.push_back(rotationTmp);
            scaleMatricesSphere.push_back(scaleTmp);
            colorVectorsSphere.push_back(colorTmp);
          }
        }
      }
    }
  }
  // Render only one particular orientation because now one by one mode is
  // selected
  else {
    if (currentAtomBIndicies != QVector2D(-1, -1)) {
      int h = currentAtomBIndicies.x();
      int i = currentAtomBIndicies.y();
      for (int j = 0; j < atomsB[h][i].size(); j++) {
        double *atomPosition = atomsB[h][i][j]->getLocation();
        double atomRadious = atomsB[h][i][j]->getRadius();

        QMatrix4x4 translationTmp;
        QMatrix4x4 rotationTmp;
        QMatrix4x4 scaleTmp;
        QVector4D colorTmp;

        translationTmp.translate(
            QVector3D(atomPosition[0], atomPosition[1], atomPosition[2]));
        rotationTmp.setToIdentity();
        scaleTmp.scale(atomRadious);
        colorTmp = colors[h];

        translationMatricesSphere.push_back(translationTmp);
        rotationMatricesSphere.push_back(rotationTmp);
        scaleMatricesSphere.push_back(scaleTmp);
        colorVectorsSphere.push_back(colorTmp);
      }
    }
  }
}

void SweepViewWidget::clearTarnsformationMatrices() {
  translationMatricesSphere.clear();
  rotationMatricesSphere.clear();
  scaleMatricesSphere.clear();
  colorVectorsSphere.clear();
}

void SweepViewWidget::setAtlasNode(AtlasNode *atlasNode) {
  currentNodeID = atlasNode->getID();
  // Delete old ACR and create new to feed it to saveloader that will allow to
  // get orientation data about atlasNode
  if (currentACR != 0) {
    delete currentACR;
  }
  currentACR = new ActiveConstraintRegion();
  saveLoader->loadNode(atlasNode->getID(), currentACR);

  atomsA = molecularUnitA->getAtoms();
  atlasNodeIsSet = true;
  noGoodOrientation = !atlasNode->hasAnyGoodOrientation();
  currentAtomBIndicies = QVector2D(-1, -1);

  // Clear previous atomB data that is used for rendering
  for (int h = 0; h < 8; h++) {
    for (int i = 0; i < atomsB[h].size(); i++) {
      for (int j = 0; j > atomsB[h][i].size(); j++) {
        atomsB[h][i].clear();
      }
    }
    atomsB[h].clear();
  }

  vector<CayleyPoint*> sampledPoints = currentACR->GetSamplePoints();

  if (sampledPoints.size() > 0) {
    // Iterate through this spaces to find the one that has orientation
    for (int i = 0; i < sampledPoints.size(); i++) {
      if (sampledPoints[i]->getOrientations().size() > 0) {
        //                this->atlasNode = currentAtlasNode;
        std::vector<Orientation *> orientations =
            sampledPoints[i]->getOrientations();
        for (int j = 0; j < orientations.size(); j++) {
          double fb[3][3], tb[3][3];
          orientations[j]->getFromTo(fb, tb);
          vector<double> trans_b = Utils::getTransMatrix(fb, tb);
          moleculeBTransform =
              QMatrix4x4(trans_b[0], trans_b[1], trans_b[2], trans_b[3],
                         trans_b[4], trans_b[5], trans_b[6], trans_b[7],
                         trans_b[8], trans_b[9], trans_b[10], trans_b[11],
                         trans_b[12], trans_b[13], trans_b[14], trans_b[15]);
          // Apparently orientations can repeat many times that leads to
          // rendering dozens of same molecules that introduces overhead To avoid
          // it we check if atomsB[flip] already has molecule with current
          // orientation
          std::vector<Atom *> tmpAtoms = molecularUnitB->getXFAtoms(trans_b);
          bool orientationsAreSame = false;
          for (int k = 0; k < atomsB[orientations[j]->getFlipNum()].size();
               k++) {
            for (int l = 0; l < tmpAtoms.size(); l++) {
              if (tmpAtoms[l]->getLocation()[0] ==
                      atomsB[orientations[j]->getFlipNum()][k][l]
                          ->getLocation()[0] &&
                  tmpAtoms[l]->getLocation()[1] ==
                      atomsB[orientations[j]->getFlipNum()][k][l]
                          ->getLocation()[1] &&
                  tmpAtoms[l]->getLocation()[2] ==
                      atomsB[orientations[j]->getFlipNum()][k][l]
                          ->getLocation()[2]) {
                orientationsAreSame = true;
                goto loopEnd1;
              }
            }
          }
        loopEnd1:
          if (!orientationsAreSame) {
            atomsB[orientations[j]->getFlipNum()].push_back(
                tmpAtoms);  // be sure to properly delete XFAtoms2
          }
        }
      }
    }
  }
  // Take care of witness points in the same way as it was done with normal
  // points above
  vector<CayleyPoint*> witnessPoints = currentACR->GetWitnessPoints();
  if (witnessPoints.size() > 0) {
    // Iterate through this spaces to find the one that has orientation
    for (int i = 0; i < witnessPoints.size(); i++) {
      if (witnessPoints[i]->getOrientations().size() > 0) {
        //                this->atlasNode = currentAtlasNode;
        std::vector<Orientation *> orientations =
            witnessPoints[i]->getOrientations();
        for (int j = 0; j < orientations.size(); j++) {
          double fb[3][3], tb[3][3];
          orientations[j]->getFromTo(fb, tb);
          vector<double> trans_b = Utils::getTransMatrix(fb, tb);
          moleculeBTransform =
              QMatrix4x4(trans_b[0], trans_b[1], trans_b[2], trans_b[3],
                         trans_b[4], trans_b[5], trans_b[6], trans_b[7],
                         trans_b[8], trans_b[9], trans_b[10], trans_b[11],
                         trans_b[12], trans_b[13], trans_b[14], trans_b[15]);
          // Apparently orientations can repeat many times that leads to
          // rendering dozens of same molecules that introduces overhead To avoid
          // it we check if atomsB[flip] already has molecule with current
          // orientation
          std::vector<Atom *> tmpAtoms = molecularUnitB->getXFAtoms(trans_b);
          bool orientationsAreSame = false;
          for (int k = 0; k < atomsB[orientations[j]->getFlipNum()].size();
               k++) {
            for (int l = 0; l < tmpAtoms.size(); l++) {
              if (tmpAtoms[l]->getLocation()[0] ==
                      atomsB[orientations[j]->getFlipNum()][k][l]
                          ->getLocation()[0] &&
                  tmpAtoms[l]->getLocation()[1] ==
                      atomsB[orientations[j]->getFlipNum()][k][l]
                          ->getLocation()[1] &&
                  tmpAtoms[l]->getLocation()[2] ==
                      atomsB[orientations[j]->getFlipNum()][k][l]
                          ->getLocation()[2]) {
                orientationsAreSame = true;
                goto loopEnd2;
              }
            }
          }
        loopEnd2:
          if (!orientationsAreSame) {
            atomsB[orientations[j]->getFlipNum()].push_back(
                tmpAtoms);  // be sure to properly delete XFAtoms2
          }
        }
      }
    }
  }
  for (int i = 0; i < 8; i++) {
    if (atomsB[i].size() > 0) {
      sharedData->currentFlips[i] = true;
      sharedData->currentFlipID = i;
      checkBoxes[i].setTristate(false);
      checkBoxes[i].setEnabled(true);
    } else {
      sharedData->currentFlips[i] = false;
      checkBoxes[i].setTristate(true);
      checkBoxes[i].setCheckState(Qt::PartiallyChecked);
      checkBoxes[i].setEnabled(false);
    }
  }

  updateTransformationMatrices();

  return;
}

void SweepViewWidget::setMolecularUnits(MolecularUnit *molecularUnitA,
                                        MolecularUnit *molecularUnitB) {
  //[ ATLAS NODE]
  atlasNodeIsSet = false;
  this->molecularUnitA = molecularUnitA;
  this->molecularUnitB = molecularUnitB;
  //[~ATLAS NODE]
}

void SweepViewWidget::keyPressEvent(QKeyEvent *keyEvent) {
  camera.keyPressEvent(keyEvent);
}

void SweepViewWidget::mousePressEvent(QMouseEvent *mouseEvent) {
  camera.mousePressEvent(mouseEvent);
}

void SweepViewWidget::mouseReleaseEvent(QMouseEvent *mouseEvent) {
  camera.mouseReleaseEvent(mouseEvent);
}

void SweepViewWidget::wheelEvent(QWheelEvent *wheelEvent) {
  camera.wheelEvent(wheelEvent);
}
