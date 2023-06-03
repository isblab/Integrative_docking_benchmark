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
#include "ConstraintSelectionWindow.h"

ConstraintSelectionWindow::ConstraintSelectionWindow() {
  preloadNodeNumberLable.setText("Preload node number:");
  loadButton.setText("Load");
  acceptButton.setText("Accept");
  cancelButton.setText("Cancel");
  nodeValuesRangeLable.setText("Range");
  for (int i = 0; i < 6; i++) {
    connectCheckBox[i].setText("Connect");
    connect(&connectCheckBox[i], SIGNAL(stateChanged(int)), this,
            SLOT(updateConnectionData()));
    connect(&comboBoxesA[i], SIGNAL(currentIndexChanged(int)), this,
            SLOT(updateConnectionData()));
    connect(&comboBoxesB[i], SIGNAL(currentIndexChanged(int)), this,
            SLOT(updateConnectionData()));
  }
  connectCheckBox[0].setChecked(
      true);  ///< First connection is always set to true to avoid situations
              ///< when no connection is seleceted by user
  constraintSelectionPanelWidget.setMinimumWidth(240);
  constraintSelectionPanelWidget.setMaximumWidth(240);
  preloadNodeNumberLable.setMaximumHeight(30);
  int columns0and1Width = 100;
  nodeValuesRangeLable.setMinimumWidth(columns0and1Width);
  nodeValuesRangeLable.setMaximumWidth(columns0and1Width);
  preloadNodeNumberLineEdit.setMinimumWidth(columns0and1Width);
  preloadNodeNumberLineEdit.setMaximumWidth(columns0and1Width);

  constraintSelectionBondViewWidget.setConstraintSelectionMoleculeViewWidget(
      &constraintSelectionMoleculeViewWidget);  ///< we need it to get
                                                ///< transfomation matrices to
                                                ///< reder atoms IDs correctly
  constraintSelectionBondViewWidget.setConnectionData(
      &connectionData);  ///< we need it to know where to render bond
  setLayouts();
}

void ConstraintSelectionWindow::setData(Atlas *atlas,
                                        MolecularUnit *molecularUnitA,
                                        MolecularUnit *molecularUnitB) {
  this->atlas = atlas;
  this->molecularUnitA = molecularUnitA;
  this->molecularUnitB = molecularUnitB;
  constraintSelectionMoleculeViewWidget.setData(molecularUnitA, molecularUnitB);

  // Set all possible values for dropdowns
  QStringList dropDownValues[2];
  for (int i = 0; i < molecularUnitA->getAtoms().size(); i++) {
    dropDownValues[MOLECULAR_UNIT_A] << QString::number(i);
  }
  for (int i = 0; i < molecularUnitB->getAtoms().size(); i++) {
    dropDownValues[MOLECULAR_UNIT_B] << QString::number(i);
  }
  for (int i = 0; i < 6; i++) {
    comboBoxesA[i].addItems(dropDownValues[MOLECULAR_UNIT_A]);
    comboBoxesB[i].addItems(dropDownValues[MOLECULAR_UNIT_B]);
  }
}

void ConstraintSelectionWindow::updateConnectionData() {
  for (int i = 0; i < 6; i++) {
    connectionData.atomIDFromMUA[i] = comboBoxesA[i].currentIndex();
    connectionData.atomIDFromMUB[i] = comboBoxesB[i].currentIndex();
    connectionData.isConnected[i] = connectCheckBox[i].isChecked();
  }
}

void ConstraintSelectionWindow::setLayouts() {
  int row = 0;
  mainGridLayout.addWidget(&constraintSelectionMoleculeViewWidget, 0, 0, 1, 1);
  mainGridLayout.addWidget(&constraintSelectionBondViewWidget, 0, 0, 1, 1);
  mainGridLayout.addWidget(&constraintSelectionPanelWidget, 0, 1, 1, 1);

  // Rows 0
  constraintSelectionPanelGridLayout.addWidget(&preloadNodeNumberLable, row++,
                                               0, 1, 3);
  constraintSelectionPanelGridLayout.addWidget(&nodeValuesRangeLable, row, 0, 1,
                                               1);
  constraintSelectionPanelGridLayout.addWidget(&preloadNodeNumberLineEdit, row,
                                               1, 1, 1);
  constraintSelectionPanelGridLayout.addWidget(&loadButton, row++, 2, 1, 1);
  // Rows 1-6
  for (int i = 0; i < 6; i++) {
    if (i == 0) {
      constraintSelectionPanelGridLayout.addWidget(&comboBoxesA[i], i + row, 0,
                                                   1, 1);
      constraintSelectionPanelGridLayout.addWidget(&comboBoxesB[i], i + row, 1,
                                                   1, 1);
    } else {
      constraintSelectionPanelGridLayout.addWidget(&comboBoxesA[i], i + row, 0,
                                                   1, 1);
      constraintSelectionPanelGridLayout.addWidget(&comboBoxesB[i], i + row, 1,
                                                   1, 1);
      constraintSelectionPanelGridLayout.addWidget(&connectCheckBox[i], i + row,
                                                   2, 1, 2);
    }
  }
  row += 6;
  // Last row
  constraintSelectionPanelGridLayout.addWidget(&acceptButton, row, 0, 1, 1);
  constraintSelectionPanelGridLayout.addWidget(&cancelButton, row, 1, 1, 1);

  constraintSelectionPanelWidget.setLayout(&constraintSelectionPanelGridLayout);

  setLayout(&mainGridLayout);

  connect(&acceptButton, SIGNAL(clicked(bool)), this,
          SLOT(acceptButtonClicked()));
}

ActiveConstraintGraph *ConstraintSelectionWindow::getConstraintGraph() {
  return graph;
}

void ConstraintSelectionWindow::acceptButtonClicked() {
  vector<pair<int, int> > parts;
  for (int i = 0; i < 6; i++) {
    if (connectCheckBox[i].isChecked()) {
      parts.push_back(make_pair(comboBoxesA[i].currentIndex(),
                                comboBoxesB[i].currentIndex()));
    }
  }
  graph = new ActiveConstraintGraph(parts);
  this->close();
}
