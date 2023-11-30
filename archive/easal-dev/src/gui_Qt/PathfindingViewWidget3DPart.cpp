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
#include "PathfindingViewWidget3DPart.h"

PathfindingViewWidget3DPart::PathfindingViewWidget3DPart(
    SharedDataGUI *sharedData) {
  this->sharedData = sharedData;
  //[ TIMER SETTINGS]
  updateTimeInMS = 30;
  timer.start(updateTimeInMS);
  connect(&timer, SIGNAL(timeout()), this, SLOT(update()));
  //[ TIMER SETTINGS]
  //[ WIDGET SETTINGS]
  setFocusPolicy(Qt::StrongFocus);
  setUpdateBehavior(UpdateBehavior::PartialUpdate);
  //[~WIDGET SETTINGS]
  //[ MANAGE LAYOUTS]
  int row = 0;
  mainGridLayout.addWidget(&bondViewWidget, row++, 0, 1, 1);
  mainGridLayout.addWidget(&selectedMoleculeOnTopCheckBox, row++, 0, 1, 1);
  mainGridLayout.addWidget(&slider, row++, 0, 1, 1);
  mainGridLayout.setAlignment(Qt::AlignBottom);
  connect(&selectedMoleculeOnTopCheckBox, SIGNAL(toggled(bool)), this,
          SLOT(selectedMoleculeOnTopCheckBoxSlot()));
  selectedMoleculeOnTopCheckBox.setText("Selected node realisations on top");
  setLayout(&mainGridLayout);
  //[~MANAGE LAYOUTS]
}

void PathfindingViewWidget3DPart::setData(Atlas *atlas, SaveLoader *saveLoader,
                                          MolecularUnit *molecularUnitA,
                                          MolecularUnit *molecularUnitB) {
  this->atlas = atlas;
  this->saveLoader = saveLoader;
  this->molecularUnitA = molecularUnitA;
  this->molecularUnitB = molecularUnitB;
  atomsA = molecularUnitA->getAtoms();
}

void PathfindingViewWidget3DPart::setPath(vector<PathfindingReturnTuple> path) {
  for (int i = 0; i < path.size(); i++)  {
    // read values in path vector and store them in
	// groupedOrientationMarices grouped by atlasNodeID
    auto atlasNodeIDToGroupedAtomsBFirstIndexIterator =
        atlasNodeIDToGroupedAtomsBFirstIndex.find(path[i].atlasNodeID);
    if (atlasNodeIDToGroupedAtomsBFirstIndexIterator ==
        atlasNodeIDToGroupedAtomsBFirstIndex.end())  {// Such index wasn't found so we add it 
      atlasNodeIDToGroupedAtomsBFirstIndex[path[i].atlasNodeID] =
          groupedAtomsB.size();  // Added index to map
      vector<vector<Atom *>> moleculeOfSameAtlasNodeTmp;  // create new vector tha will be added to
                                       // groupedOrientationMarices
      double fb[3][3], tb[3][3];
      path[i].orientation->getFromTo(fb, tb);
      vector<double> trans_b = Utils::getTransMatrix(fb, tb);
      vector<Atom *> atomsBTmp = molecularUnitB->getXFAtoms(
          trans_b);  // get all atoms of current molecule B
      moleculeOfSameAtlasNodeTmp.push_back(atomsBTmp);
      groupedAtomsB.push_back(moleculeOfSameAtlasNodeTmp);
      // add node ID to nodeIDs
      nodeIDs.push_back(path[i].atlasNodeID);
    } else  // index exists so we push current molecules B atoms there
    {
      int groupedAtomsBFirstIndex =
          atlasNodeIDToGroupedAtomsBFirstIndexIterator->second;
      double fb[3][3], tb[3][3];
      path[i].orientation->getFromTo(fb, tb);
      vector<double> trans_b = Utils::getTransMatrix(fb, tb);
      vector<Atom *> atomsBTmp = molecularUnitB->getXFAtoms(
          trans_b);  // get all atoms of current molecule B
      groupedAtomsB[groupedAtomsBFirstIndex].push_back(atomsBTmp);
    }
  }
  // Free space by deleting orientations manually
  /*for(int i=0;i<path.size();i++)
  {
      delete path[i].orientation;
  }*/
  // Set total number of molecules B
  totalNumberOfMoleculesB = path.size();
  isPathSet = true;
  currentAtlasNodeID = nodeIDs[0];
  updateTransformationMatrices();
  updateBondData();
  slider.setPathfindingData(nodeIDs);
}

void PathfindingViewWidget3DPart::selectedMoleculeOnTopCheckBoxSlot() {
  if (sharedData->atlasNodeWasSet && isPathSet) {
    updateTransformationMatrices();
  }
}

void PathfindingViewWidget3DPart::initializeGL() {
  renderer.initOpenGL();
  sphere.createAndBind("12by12UVsphere.obj");
  cylinder.createAndBind("12Cylinder.obj");
}

void PathfindingViewWidget3DPart::paintGL() {
  if (sharedData->atlasNodeWasSet && isPathSet) {
    if (slider.getCurrentValue() != currentAtlasNodeID) {
      currentAtlasNodeID = slider.getCurrentValue();
      updateTransformationMatrices();
      updateBondData();
    }
    renderer.renderingProgram->setUniformValue(
        "lightPositionVector",
        camera.getEyeLocation() + camera.getRight() * 10 + camera.getUp() * 30);
    renderer.renderingProgram->setUniformValue(
        "cameraMatrix", camera.getCameraMatrix(mapFromGlobal(cursor().pos())));

    renderer.bindRenderProgramAndClear();
    renderer.renderWithTransperansy(&sphere, translationMatrices,
                                    rotationMatrices, scaleMatrices,
                                    colorVectors);
  }
}

void PathfindingViewWidget3DPart::resizeGL(int w, int h) {
  camera.resizeEvent(w, h);
}

void PathfindingViewWidget3DPart::updateTransformationMatrices() {
  clearTarnsformationMatrices();
  transperancyStrength = 0.3;
  // Prepare matrices for molecule A
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

    translationMatrices.push_back(translationTmp);
    rotationMatrices.push_back(rotationTmp);
    scaleMatrices.push_back(scaleTmp);
    colorVectors.push_back(colorTmp);
  }
  for (int i = 0; i < groupedAtomsB.size(); i++) {
    for (int j = 0; j < groupedAtomsB[i].size(); j++) {
      for (int k = 0; k < groupedAtomsB[i][j].size(); k++) {
        double *atomPosition = groupedAtomsB[i][j][k]->getLocation();
        double atomRadious = groupedAtomsB[i][j][k]->getRadius();

        QMatrix4x4 translationTmp;
        QMatrix4x4 rotationTmp;
        QMatrix4x4 scaleTmp;
        QVector4D colorTmp;

        translationTmp.translate(
            QVector3D(atomPosition[0], atomPosition[1], atomPosition[2]));
        rotationTmp.setToIdentity();
        scaleTmp.scale(atomRadious);
        // Check if current index corresponds to current atlas node ID
        if (atlasNodeIDToGroupedAtomsBFirstIndex[currentAtlasNodeID] == i) {
          colorTmp = QVector4D(0, 0, 1, 1);
        } else {
          if (selectedMoleculeOnTopCheckBox.isChecked()) {
            colorTmp = QVector4D(
                0.8, 0.8, 0.8, 0.99);  // if alpha value is less than 1 renderer
                                       // will render without depth test
          } else {
            colorTmp = QVector4D(0.8, 0.8, 0.8, 1);
          }
        }

        translationMatrices.push_back(translationTmp);
        rotationMatrices.push_back(rotationTmp);
        scaleMatrices.push_back(scaleTmp);
        colorVectors.push_back(colorTmp);
      }
    }
  }
}

void PathfindingViewWidget3DPart::updateBondData() {
  AtomBondData bondData;
  bondData.bondSolid = (*atlas)[currentAtlasNodeID]->getCG()->getParticipants();
  bondData.bondVariable =
      (*atlas)[currentAtlasNodeID]->getCG()->getParamLines();
  bondViewWidget.setBondData(&bondData);
}

void PathfindingViewWidget3DPart::clearTarnsformationMatrices() {
  translationMatrices.clear();
  rotationMatrices.clear();
  scaleMatrices.clear();
  colorVectors.clear();
}

void PathfindingViewWidget3DPart::keyPressEvent(QKeyEvent *keyEvent) {
  camera.keyPressEvent(keyEvent);
}

void PathfindingViewWidget3DPart::mousePressEvent(QMouseEvent *mouseEvent) {
  camera.mousePressEvent(mouseEvent);
}

void PathfindingViewWidget3DPart::mouseReleaseEvent(QMouseEvent *mouseEvent) {
  camera.mouseReleaseEvent(mouseEvent);
}

void PathfindingViewWidget3DPart::wheelEvent(QWheelEvent *wheelEvent) {
  camera.wheelEvent(wheelEvent);
}
