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
#include "RealisationViewWidget.h"

RealisationViewWidget::RealisationViewWidget(SharedDataGUI *sharedData) {
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
  gridLayout.setAlignment(Qt::AlignBottom);
  gridLayout.addWidget(&showAllOrientationsCheckBox, 0, 0, 1, 1);
  gridLayout.addWidget(&orientationSelectionLabel, 1, 0, 1, 1);
  gridLayout.addWidget(&orientationSelectionComboBox, 1, 1, 1, 1);
  connect(&showAllOrientationsCheckBox, SIGNAL(stateChanged(int)), this,
          SLOT(showAllOrintationsCheckBoxSlot()));
  connect(&orientationSelectionComboBox, SIGNAL(currentIndexChanged(int)), this,
          SLOT(orientationSelectionComboBoxSlot()));
  setLayout(&gridLayout);
  setFocusPolicy(Qt::StrongFocus);
  setUpdateBehavior(UpdateBehavior::NoPartialUpdate);
  //[~WIDGET SETTINGS]
  //[ COLORS]
  colors[0] = QVector4D(0.75, 0, 0.75, 1);
  colors[1] = QVector4D(0.75, 0, 0.25, 1);
  colors[2] = QVector4D(0.75, 1, 0, 1);
  colors[3] = QVector4D(0.25, 1, 0.25, 1);

  colors[4] = QVector4D(0, 1, 0.75, 1);
  colors[5] = QVector4D(0, 1, 1, 1);
  colors[6] = QVector4D(0, 0.5, 0.25, 1);
  colors[7] = QVector4D(0, 0.5, 0.75, 1);
  //[~COLORS]
  //[ MISC]
  showAllOrientationsCheckBox.setText("Show all orientations ");
  orientationSelectionLabel.setText("Select flip ");
  currentNodeID = -1;
  currentCayleyPointID = -1;
  currentFlipID = -1;
  currentACR = 0;
  //[~MISC]
}

void RealisationViewWidget::setSaveLoader(SaveLoader *saveLoader) {
  this->saveLoader = saveLoader;
}

void RealisationViewWidget::initializeGL() {
  renderer.initOpenGL();
  sphere.createAndBind("12by12UVsphere.obj");
  cylinder.createAndBind("12Cylinder.obj");
}

void RealisationViewWidget::paintGL() {
  // We need to block it to prevent crashes
  orientationSelectionComboBox.blockSignals(true);
  if (sharedData->atlasNodeWasSet) {
    updateAtlasNode();
    updateCayleyPoint();

    renderer.renderingProgram->setUniformValue(
        "lightPositionVector",
        camera.getEyeLocation() + camera.getRight() * 10 + camera.getUp() * 30);
    renderer.renderingProgram->setUniformValue(
        "cameraMatrix", camera.getCameraMatrix(mapFromGlobal(cursor().pos())));

    renderer.bindRenderProgramAndClear();
    // Render spheres that represent atoms
    renderer.render(&sphere, translationMatricesSphere, rotationMatricesSphere,
                    scaleMatricesSphere, colorVectorsSphere);
    // Render cylinders that represent bonds
    renderer.disableDepthTest();
    renderer.render(&cylinder, translationMatricesCylinder,
                    rotationMatricesCylinder, scaleMatricesCylinder,
                    colorVectorsCylinder);
    renderer.enableDepthTest();
  }
  orientationSelectionComboBox.blockSignals(false);
}

void RealisationViewWidget::resizeGL(int w, int h) { camera.resizeEvent(w, h); }

void RealisationViewWidget::showAllOrintationsCheckBoxSlot() {
  if (sharedData->atlasNodeWasSet) {
    if (currentNodeID == sharedData->currentAtlasNode->getID()) {
      updateTransformationMatrices();
    }
  }
}

void RealisationViewWidget::orientationSelectionComboBoxSlot() {
  if (sharedData->atlasNodeWasSet) {
    if (currentNodeID == sharedData->currentAtlasNode->getID()) {
      updateTransformationMatrices();
    }
  }
}

void RealisationViewWidget::updateAtlasNode() {
  if (currentNodeID != sharedData->currentAtlasNode->getID()) {
    setAtlasNode(sharedData->currentAtlasNode);
  }
}

void RealisationViewWidget::updateCayleyPoint() {
  if (currentCayleyPointID != sharedData->currentCayleyPointID &&
      sharedData->currentCayleyPointID >= 0 &&
      sharedData->currentCayleyPointID < currentACR->GetSamplePoints().size()) {
    currentCayleyPointID = sharedData->currentCayleyPointID;
    currentCayleyPointType = sharedData->currentCayleyPointType;
    setCayleyPoint(currentCayleyPointID, currentCayleyPointType);
  }
}

void RealisationViewWidget::updateTransformationMatrices() {
  clearTarnsformationMatrices();

  //[ UPDATE ATOMS]
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
  for (int i = 0; i < atomsB.size(); i++) {
    for (int j = 0; j < atomsB[i].size(); j++) {
      if (orientationSelectionComboBox.currentIndex() == i ||
          showAllOrientationsCheckBox.isChecked()) {
        double *atomPosition = atomsB[i][j]->getLocation();
        double atomRadious = atomsB[i][j]->getRadius();

        QMatrix4x4 translationTmp;
        QMatrix4x4 rotationTmp;
        QMatrix4x4 scaleTmp;
        QVector4D colorTmp;

        translationTmp.translate(
            QVector3D(atomPosition[0], atomPosition[1], atomPosition[2]));
        rotationTmp.setToIdentity();
        scaleTmp.scale(atomRadious);
        colorTmp = colors[flipNumberCorrespondentToID[i]];

        translationMatricesSphere.push_back(translationTmp);
        rotationMatricesSphere.push_back(rotationTmp);
        scaleMatricesSphere.push_back(scaleTmp);
        colorVectorsSphere.push_back(colorTmp);
      }
    }
    if (!showAllOrientationsCheckBox.isChecked() &&
        orientationSelectionComboBox.currentIndex() == i) {
      break;
    }
  }
  //[ UPDATE ATOMS]
  //[ UPDATE BONDS]
  for (int i = 0; i < atomsB.size(); i++) {
    for (int j = 0; j < bondSolid.size(); j++) {
      if (orientationSelectionComboBox.currentIndex() == i ||
          showAllOrientationsCheckBox.isChecked()) {
        double *bondStart, *bondEnd;
        bondStart = atomsA[bondSolid[j].first]->getLocation();
        bondEnd = atomsB[i][bondSolid[j].second]->getLocation();

        QMatrix4x4 translationTmp;
        QMatrix4x4 rotationTmp;
        QMatrix4x4 scaleTmp;
        QVector4D colorTmp = QVector4D(0.1, 0.1, 0.1, 1);

        transformCylinderAccordingToData(
            QVector3D(bondStart[0], bondStart[1], bondStart[2]),
            QVector3D(bondEnd[0], bondEnd[1], bondEnd[2]), 0.1, &translationTmp,
            &rotationTmp, &scaleTmp);

        translationMatricesCylinder.push_back(translationTmp);
        rotationMatricesCylinder.push_back(rotationTmp);
        scaleMatricesCylinder.push_back(scaleTmp);
        colorVectorsCylinder.push_back(colorTmp);
      }
    }
    if (!showAllOrientationsCheckBox.isChecked() &&
        orientationSelectionComboBox.currentIndex() == i) {
      break;
    }
  }
  for (int i = 0; i < atomsB.size(); i++) {
    for (int j = 0; j < bondVariable.size(); j++) {
      if (orientationSelectionComboBox.currentIndex() == i ||
          showAllOrientationsCheckBox.isChecked()) {
        double *bondStart, *bondEnd;
        bondStart = atomsA[bondVariable[j].first]->getLocation();
        bondEnd = atomsB[i][bondVariable[j].second]->getLocation();

        QMatrix4x4 translationTmp;
        QMatrix4x4 rotationTmp;
        QMatrix4x4 scaleTmp;
        QVector4D colorTmp = QVector4D(0.7, 0.7, 0.7, 1);

        transformCylinderAccordingToData(
            QVector3D(bondStart[0], bondStart[1], bondStart[2]),
            QVector3D(bondEnd[0], bondEnd[1], bondEnd[2]), 0.1, &translationTmp,
            &rotationTmp, &scaleTmp);

        translationMatricesCylinder.push_back(translationTmp);
        rotationMatricesCylinder.push_back(rotationTmp);
        scaleMatricesCylinder.push_back(scaleTmp);
        colorVectorsCylinder.push_back(colorTmp);
      }
    }

    if (!showAllOrientationsCheckBox.isChecked() &&
        orientationSelectionComboBox.currentIndex() == i) {
      break;
    }
  }
  //[~UPDATE BONDS]
}

void RealisationViewWidget::updateBondData() {
  bondSolid = sharedData->currentAtlasNode->getCG()->getParticipants();
  bondVariable = sharedData->currentAtlasNode->getCG()->getParamLines();

  // Set bondData for bondViewWidget
  AtomBondData bondData;
  bondData.bondSolid = this->bondSolid;
  bondData.bondVariable = this->bondVariable;
}

void RealisationViewWidget::updateOrientationSelectionComboBoxValues(
    std::vector<Orientation *> orientations) {
  orientationSelectionComboBox.clear();
  flipNumberCorrespondentToID.clear();
  orientationSelectionComboBox.blockSignals(true);
  for (int i = 0; i < orientations.size(); i++) {
    orientationSelectionComboBox.addItem(
        QString::number(orientations[i]->getFlipNum()));
    flipNumberCorrespondentToID.push_back(orientations[i]->getFlipNum());
  }
  orientationSelectionComboBox.blockSignals(false);
}

void RealisationViewWidget::clearTarnsformationMatrices() {
  translationMatricesSphere.clear();
  rotationMatricesSphere.clear();
  scaleMatricesSphere.clear();
  colorVectorsSphere.clear();

  translationMatricesCylinder.clear();
  rotationMatricesCylinder.clear();
  scaleMatricesCylinder.clear();
  colorVectorsCylinder.clear();
}

void RealisationViewWidget::setAtlasNode(AtlasNode *atlasNode) {
  currentNodeID = atlasNode->getID();
  if (currentACR != 0) {
    delete currentACR;
  }
  currentACR = new ActiveConstraintRegion();
  saveLoader->loadNode(atlasNode->getID(), currentACR);
  //    sharedData->currentACR = currentACR;
  //    sharedData->ACRWasSet = true;
  vector<CayleyPoint*> space = currentACR->GetSamplePoints();
  vector<CayleyPoint*> witspace = currentACR->GetWitnessPoints();
  if (space.size() > 0) {
    // Iterate through this spaces to find the one that has orientation
    for (int i = 0; i < space.size(); i++) {
      // Get first available orientation and return
      if (space[i]->getOrientations().size() > 0) {
        setCayleyPoint(i, PointType::SampledPoint);
        break;
      }
    }
  } else if (witspace.size() > 0) {
    // Iterate through this spaces to find the one that has orientation
    for (int i = 0; i < witspace.size(); i++) {
      if (witspace[i]->getOrientations().size() > 0) {
        setCayleyPoint(i, PointType::WitnessPoint);
        break;
      }
    }
  }

  updateBondData();

  atlasNodeIsSet = true;
}

void RealisationViewWidget::setCayleyPoint(int cayleyPointID,
                                           PointType pointType) {
  atomsB.clear();

  //    sharedData->currentCayleyPointID = cayleyPointID;
  //    sharedData->currentFlipID = orientations[0]->getFlipNum();

  atomsA = molecularUnitA->getAtoms();

  if (pointType == PointType::SampledPoint) { 
    vector<CayleyPoint*> space = currentACR->GetSamplePoints();
    if (space[cayleyPointID]->getOrientations().size() > 0) {
      std::vector<Orientation *> orientations =
          space[cayleyPointID]->getOrientations();
      for (int i = 0; i < orientations.size(); i++) {
        double fb[3][3], tb[3][3];
        orientations[i]->getFromTo(fb, tb);
        vector<double> trans_b = Utils::getTransMatrix(fb, tb);
        moleculeBTransform =
            QMatrix4x4(trans_b[0], trans_b[1], trans_b[2], trans_b[3],
                       trans_b[4], trans_b[5], trans_b[6], trans_b[7],
                       trans_b[8], trans_b[9], trans_b[10], trans_b[11],
                       trans_b[12], trans_b[13], trans_b[14], trans_b[15]);

        atomsB.push_back(molecularUnitB->getXFAtoms(
            trans_b));  // be sure to properly delete XFAtoms2
      }
      updateOrientationSelectionComboBoxValues(orientations);
    }
  } else if (pointType == PointType::WitnessPoint) {
	vector<CayleyPoint*> witspace = currentACR->GetWitnessPoints();
    if (witspace[cayleyPointID]->getOrientations().size() > 0) {
      std::vector<Orientation *> orientations =
          witspace[cayleyPointID]->getOrientations();
      for (int i = 0; i < orientations.size(); i++) {
        double fb[3][3], tb[3][3];
        orientations[0]->getFromTo(fb, tb);
        vector<double> trans_b = Utils::getTransMatrix(fb, tb);
        moleculeBTransform =
            QMatrix4x4(trans_b[0], trans_b[1], trans_b[2], trans_b[3],
                       trans_b[4], trans_b[5], trans_b[6], trans_b[7],
                       trans_b[8], trans_b[9], trans_b[10], trans_b[11],
                       trans_b[12], trans_b[13], trans_b[14], trans_b[15]);

        atomsB.push_back(molecularUnitB->getXFAtoms(
            trans_b));  // be sure to properly delete XFAtoms2
      }
      updateOrientationSelectionComboBoxValues(orientations);
    }
  }
  updateTransformationMatrices();
}

void RealisationViewWidget::setMolecularUnits(MolecularUnit *molecularUnitA,
                                              MolecularUnit *molecularUnitB) {
  //[ ATLAS NODE]
  atlasNodeIsSet = false;
  this->molecularUnitA = molecularUnitA;
  this->molecularUnitB = molecularUnitB;
  //[~ATLAS NODE]
}

void RealisationViewWidget::keyPressEvent(QKeyEvent *keyEvent) {
  camera.keyPressEvent(keyEvent);
}

void RealisationViewWidget::mousePressEvent(QMouseEvent *mouseEvent) {
  camera.mousePressEvent(mouseEvent);
}

void RealisationViewWidget::mouseReleaseEvent(QMouseEvent *mouseEvent) {
  camera.mouseReleaseEvent(mouseEvent);
}

void RealisationViewWidget::wheelEvent(QWheelEvent *wheelEvent) {
  camera.wheelEvent(wheelEvent);
}

void RealisationViewWidget::transformCylinderAccordingToData(
    QVector3D positionStart, QVector3D positionEnd, float radious,
    QMatrix4x4 *translation, QMatrix4x4 *rotation, QMatrix4x4 *scale) {
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

  translation->setToIdentity();
  rotation->setToIdentity();
  scale->setToIdentity();

  translation->translate(positionStart + directionNotNormalized / 2);
  rotation->rotate(angle, axisOfRotation);
  scale->scale(
      QVector3D(radious, directionNotNormalized.length() / 2, radious));
}
