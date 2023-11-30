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
#include "CayleyStatisticsWindow.h"

CayleyStatisticsWindow::CayleyStatisticsWindow(SharedDataGUI *sharedData) {
  this->sharedData = sharedData;
  this->setWindowTitle("Cayley Statistics Window");

  this->resize(750, 500);

  // Create the layout Widgets
  centralWidget = new QWidget(this);
  verticalLayoutWidget_1 = new QWidget(centralWidget);
  verticalLayoutWidget_2 = new QWidget(centralWidget);
  verticalLayoutWidget_3 = new QWidget(centralWidget);
  verticalLayoutWidget_4 = new QWidget(centralWidget);
  horizontalLayoutWidget_1 = new QWidget(centralWidget);
  horizontalLayoutWidget_2 = new QWidget(centralWidget);

  // Create the layouts
  verticalLayout_1 = new QVBoxLayout(verticalLayoutWidget_1);
  verticalLayout_2 = new QVBoxLayout(verticalLayoutWidget_2);
  verticalLayout_3 = new QVBoxLayout(verticalLayoutWidget_3);
  verticalLayout_4 = new QVBoxLayout(verticalLayoutWidget_4);
  horizontalLayout_1 = new QHBoxLayout(horizontalLayoutWidget_1);
  horizontalLayout_2 = new QHBoxLayout(horizontalLayoutWidget_2);

  // Set the Geometry of the Layout Widgets
  verticalLayoutWidget_1->setGeometry(QRect(20, 30, 210, 321));
  verticalLayoutWidget_2->setGeometry(QRect(180, 30, 50, 321));
  verticalLayoutWidget_3->setGeometry(QRect(380, 30, 210, 321));
  verticalLayoutWidget_4->setGeometry(QRect(580, 30, 210, 321));
  horizontalLayoutWidget_1->setGeometry(QRect(0, 350, 470, 80));
  horizontalLayoutWidget_2->setGeometry(QRect(30, 430, 470, 80));

  // Crate the labels
  dimensionHeader = new QLabel(verticalLayoutWidget_1);
  minValuesHeader = new QLabel(verticalLayoutWidget_2);
  maxValuesHeader = new QLabel(verticalLayoutWidget_3);
  averageValuesHeader = new QLabel(verticalLayoutWidget_4);

  for (int i = 0; i < 6; i++) {
    dimension[i] = new QLabel(verticalLayoutWidget_1);
    minValues[i] = new QLabel(verticalLayoutWidget_2);
    maxValues[i] = new QLabel(verticalLayoutWidget_3);
    averageValues[i] = new QLabel(verticalLayoutWidget_4);
  }

  // Create the buttons
  generateCayleyStatsButton = new QPushButton(horizontalLayoutWidget_1);
  exitButton = new QPushButton(horizontalLayoutWidget_1);

  // Create a status bar
  statusBar = new QLabel(horizontalLayoutWidget_2);

  // Add the labels to the widgets

  for (int i = 0; i < 6; i++) {
    verticalLayout_1->addWidget(dimension[i]);
    verticalLayout_2->addWidget(minValues[i]);
    verticalLayout_3->addWidget(maxValues[i]);
    verticalLayout_4->addWidget(averageValues[i]);
  }

  // Add the buttons to the widgets
  horizontalLayout_1->addWidget(generateCayleyStatsButton);
  horizontalLayout_1->addWidget(exitButton);

  // Set text in labels
  dimensionHeader->setText("Dimension");
  minValuesHeader->setText("Min Cayley Param");
  maxValuesHeader->setText("Max Cayley Param");
  averageValuesHeader->setText("Average Cayley Param");
  for (int i = 0; i < 6; i++) {
    string dim = to_string(i) + "D Nodes";
    dimension[i]->setText(dim.c_str());
    minValues[i]->setText("0");
    maxValues[i]->setText("0");
    averageValues[i]->setText("0");
  }

  // Set text on Buttons
  generateCayleyStatsButton->setText("Generate Cayley Statistics");
  exitButton->setText("Exit");

  generateCayleyStatsButton->resize(QSize(100, 20));
  exitButton->resize(QSize(100, 20));

  statusBar->resize(QSize(450, 60));

  // Connect the Slots
  connect(exitButton, SIGNAL(pressed()), this, SLOT(exit()));
  connect(generateCayleyStatsButton, SIGNAL(pressed()), this,
          SLOT(generateCayleyStatistics()));
}

void CayleyStatisticsWindow::generateCayleyStatistics() {
  BasicStatistics *bs = new BasicStatistics();

  statusBar->setStyleSheet("QLabel {color : black;}");
  statusBar->setText("Computing Cayley statistics.");

  // int status =
  // bs->generateCayleyStats(sharedData->currentAtlasNode->getID());
  int status = bs->generateAtlasStats();

  if (status == -1) {
    statusBar->setStyleSheet("QLabel {color : red;}");
    statusBar->setText(
        "Atlas not found. Please load a previously generated atlas \nor create "
        "a new one from the File->Create menu.");
    return;
  }

  double *cayleyMin = bs->getCayleyMin();
  double *cayleyMax = bs->getCayleyMax();
  double *cayleyAve = bs->getOverallGeometricAverage();

  for (int i = 0; i < 6; i++) {
    minValues[i]->setText(to_string(cayleyMin[i]).c_str());
    maxValues[i]->setText(to_string(cayleyMax[i]).c_str());
    averageValues[i]->setText(to_string(cayleyAve[i]).c_str());
  }

  statusBar->setText("Finished computing Cayley statistics.");

  delete bs;
}

void CayleyStatisticsWindow::exit() { this->close(); }
