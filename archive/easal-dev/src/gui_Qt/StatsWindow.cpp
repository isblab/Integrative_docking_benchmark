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
#include "StatsWindow.h"

StatsWindow::StatsWindow() {
  this->setWindowTitle("Atlas Statistics Window");

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
  verticalLayoutWidget_1->setGeometry(QRect(10, 30, 230, 320));
  verticalLayoutWidget_2->setGeometry(QRect(240, 30, 50, 320));
  verticalLayoutWidget_3->setGeometry(QRect(290, 30, 210, 320));
  horizontalLayoutWidget_1->setGeometry(QRect(0, 350, 470, 80));
  horizontalLayoutWidget_2->setGeometry(QRect(30, 430, 470, 80));

  // Crate the labels
  numNodes = new QLabel(verticalLayoutWidget_1);
  numGoodSamples = new QLabel(verticalLayoutWidget_1);
  numCollisionSamples = new QLabel(verticalLayoutWidget_1);
  numTotalSamples = new QLabel(verticalLayoutWidget_1);

  eq1 = new QLabel(verticalLayoutWidget_2);
  eq2 = new QLabel(verticalLayoutWidget_2);
  eq3 = new QLabel(verticalLayoutWidget_2);
  eq4 = new QLabel(verticalLayoutWidget_2);

  numNodesResults = new QLabel(verticalLayoutWidget_3);
  numGoodSamplesResults = new QLabel(verticalLayoutWidget_3);
  numCollisionSamplesResults = new QLabel(verticalLayoutWidget_3);
  numTotalSamplesResults = new QLabel(verticalLayoutWidget_3);

  // Create the buttons
  generateStatsButton = new QPushButton(horizontalLayoutWidget_1);
  exitButton = new QPushButton(horizontalLayoutWidget_1);

  // Create a status bar
  statusBar = new QLabel(horizontalLayoutWidget_2);

  // Add the labels to the widgets
  verticalLayout_1->addWidget(numNodes);
  verticalLayout_1->addWidget(numGoodSamples);
  verticalLayout_1->addWidget(numCollisionSamples);
  verticalLayout_1->addWidget(numTotalSamples);

  verticalLayout_2->addWidget(eq1);
  verticalLayout_2->addWidget(eq2);
  verticalLayout_2->addWidget(eq3);
  verticalLayout_2->addWidget(eq4);

  verticalLayout_3->addWidget(numNodesResults);
  verticalLayout_3->addWidget(numGoodSamplesResults);
  verticalLayout_3->addWidget(numCollisionSamplesResults);
  verticalLayout_3->addWidget(numTotalSamplesResults);

  // Add the buttons to the widgets
  horizontalLayout_1->addWidget(generateStatsButton);
  horizontalLayout_1->addWidget(exitButton);

  // Set text in labels
  numNodes->setText("Number of Nodes");
  numGoodSamples->setText("Number of Good Samples");
  numCollisionSamples->setText("Number of Collision Samples");
  numTotalSamples->setText("Number of Total Samples");

  eq1->setText("=");
  eq2->setText("=");
  eq3->setText("=");
  eq4->setText("=");

  numNodesResults->setText("0");
  numGoodSamplesResults->setText("0");
  numCollisionSamplesResults->setText("0");
  numTotalSamplesResults->setText("0");

  // Set text on Buttons
  generateStatsButton->setText("Generate Statistics");
  exitButton->setText("Exit");

  generateStatsButton->resize(QSize(100, 20));
  exitButton->resize(QSize(100, 20));

  statusBar->resize(QSize(450, 60));

  // Connect the Slots
  connect(exitButton, SIGNAL(pressed()), this, SLOT(exit()));
  connect(generateStatsButton, SIGNAL(pressed()), this, SLOT(generateStats()));
}

void StatsWindow::generateStats() {
  BasicStatistics *bs = new BasicStatistics();
  int status = bs->generateAtlasStats();
  if (status == -1) {
    statusBar->setStyleSheet("QLabel {color : red;}");
    statusBar->setText(
        "Atlas not found. Please load a previously generated atlas \nor create "
        "a new one from the File->Create menu.");
    return;
  }
  statusBar->setStyleSheet("QLabel {color : black;}");
  statusBar->setText("Generating Statistics. Please wait.");
  int *nra = bs->getNumRegions();
  int *nsa = bs->getNumSamples();
  int *ngsa = bs->getNumGoodSamples();
  int *ncsa = bs->getNumCollisionSamples();

  int nr = 0;
  int ns = 0;
  int ngs = 0;
  int ncs = 0;

  for (int i = 0; i < 6; i++) {
    nr += nra[i];
    ns += nsa[i];
    ngs += ngsa[i];
    ncs += ncsa[i];
  }

  numNodesResults->setText(std::to_string(nr).c_str());
  numGoodSamplesResults->setText(std::to_string(ngs).c_str());
  numCollisionSamplesResults->setText(std::to_string(ncs).c_str());
  numTotalSamplesResults->setText(std::to_string(ns).c_str());

  statusBar->setText("Finished Generating Statistics.");

  delete bs;
}

void StatsWindow::exit() { this->close(); }
