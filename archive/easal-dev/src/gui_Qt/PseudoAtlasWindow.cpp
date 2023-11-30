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
 *  Created on: 2019
 *  Author: Rahul Prabhu
 */
#include "PseudoAtlasWindow.h"

PseudoAtlasWindow::PseudoAtlasWindow() {
  this->setWindowTitle("Pseudo Atlas Statistics Window");

  this->resize(500, 500);

  // Create the layout Widgets
  centralWidget = new QWidget(this);
  verticalLayoutWidget_1 = new QWidget(centralWidget);
  verticalLayoutWidget_2 = new QWidget(centralWidget);
  verticalLayoutWidget_3 = new QWidget(centralWidget);
  horizontalLayoutWidget_1 = new QWidget(centralWidget);
  horizontalLayoutWidget_2 = new QWidget(centralWidget);

  // Create the layouts
  verticalLayout_1 = new QVBoxLayout(verticalLayoutWidget_1);
  verticalLayout_2 = new QVBoxLayout(verticalLayoutWidget_2);
  verticalLayout_3 = new QVBoxLayout(verticalLayoutWidget_3);
  horizontalLayout_1 = new QHBoxLayout(horizontalLayoutWidget_1);
  horizontalLayout_2 = new QHBoxLayout(horizontalLayoutWidget_2);

  // Set the Geometry of the Layout Widgets
  verticalLayoutWidget_1->setGeometry(QRect(20, 30, 210, 321));
  verticalLayoutWidget_2->setGeometry(QRect(180, 30, 50, 321));
  verticalLayoutWidget_3->setGeometry(QRect(380, 30, 210, 321));
  horizontalLayoutWidget_1->setGeometry(QRect(0, 350, 470, 80));
  horizontalLayoutWidget_2->setGeometry(QRect(30, 430, 470, 80));

  // Crate the labels
  dimensionHeader = new QLabel(verticalLayoutWidget_1);
  atlasHeader = new QLabel(verticalLayoutWidget_2);
  pseudoAtlasHeader = new QLabel(verticalLayoutWidget_3);

  for (int i = 0; i < 6; i++) {
    dimension[i] = new QLabel(verticalLayoutWidget_1);
    atlasNodes[i] = new QLabel(verticalLayoutWidget_2);
    pseudoAtlasNodes[i] = new QLabel(verticalLayoutWidget_3);
  }

  // Create the buttons
  generatePseudoAtlasButton = new QPushButton(horizontalLayoutWidget_1);
  exitButton = new QPushButton(horizontalLayoutWidget_1);

  // Create a status bar
  statusBar = new QLabel(horizontalLayoutWidget_2);

  // Add the labels to the widgets

  for (int i = 0; i < 6; i++) {
    verticalLayout_1->addWidget(dimension[i]);
    verticalLayout_2->addWidget(atlasNodes[i]);
    verticalLayout_3->addWidget(pseudoAtlasNodes[i]);
  }

  // Add the buttons to the widgets
  horizontalLayout_1->addWidget(generatePseudoAtlasButton);
  horizontalLayout_1->addWidget(exitButton);

  // Set text in labels
  dimensionHeader->setText("Dimension");
  atlasHeader->setText("Atlas");
  pseudoAtlasHeader->setText("PseudoAtlas");
  for (int i = 0; i < 6; i++) {
    string dim = to_string(i) + "D Nodes";
    dimension[i]->setText(dim.c_str());
    atlasNodes[i]->setText("0");
    pseudoAtlasNodes[i]->setText("0");
  }

  // Set text on Buttons
  generatePseudoAtlasButton->setText("Generate Pseudo Atlas");
  exitButton->setText("Exit");

  generatePseudoAtlasButton->resize(QSize(100, 20));
  exitButton->resize(QSize(100, 20));

  statusBar->resize(QSize(450, 60));

  // Connect the Slots
  connect(exitButton, SIGNAL(pressed()), this, SLOT(exit()));
  connect(generatePseudoAtlasButton, SIGNAL(pressed()), this,
          SLOT(generatePseudoAtlas()));
}

void PseudoAtlasWindow::generatePseudoAtlas() {
  BasicStatistics *bs = new BasicStatistics();

  statusBar->setStyleSheet("QLabel {color : black;}");
  statusBar->setText("Counting Number of Regions.");

  int status = bs->countRegions();

  if (status == -1) {
    statusBar->setStyleSheet("QLabel {color : red;}");
    statusBar->setText(
        "Atlas not found. Please load a previously generated atlas \nor create "
        "a new one from the File->Create menu.");
    return;
  }

  int *atlasRegions = bs->getNumRegions();

  for (int i = 0; i < 6; i++) {
    atlasNodes[i]->setText(to_string(atlasRegions[i]).c_str());
  }

  statusBar->setText("Finished counting Atlas.");
  statusBar->setText("Generating Pseudo Atlas...");

  // pseudoAtlas *pA = new pseudoAtlas(&df, psA, psB, snl);
  status = bs->generatePseudoAtlas();

  if (status == -1) {
    statusBar->setStyleSheet("QLabel {color : red;}");
    statusBar->setText("Atlas Not present!");
    return;
  }

  statusBar->setText("Finished Generating Pseudo Atlas.");
  int *PARegions = bs->getNumPARegions();

  for (int i = 0; i < 6; i++) {
    pseudoAtlasNodes[i]->setText(to_string(PARegions[i]).c_str());
  }

  delete bs;
}

void PseudoAtlasWindow::exit() { this->close(); }
