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
#ifndef STATSWINDOW_H
#define STATSWINDOW_H

#include <QtWidgets>
#include <cstdio>
#include <iostream>

#include "easalcore/Atlas.h"
#include "easalcore/BasicStatistics.h"
#include "easalcore/SaveLoader.h"
#include "easalcore/Settings.h"

using namespace std;

class MainWindow;

class StatsWindow : public QDialog {
  Q_OBJECT

 public:
  StatsWindow();

 private Q_SLOTS:
  void exit();
  void generateStats();

 private:
  Atlas *atlas;
  SaveLoader *snl;

  QWidget *centralWidget;
  QWidget *verticalLayoutWidget_1;
  QWidget *verticalLayoutWidget_2;
  QWidget *verticalLayoutWidget_3;
  QWidget *horizontalLayoutWidget_1;
  QWidget *horizontalLayoutWidget_2;

  QVBoxLayout *verticalLayout_1;
  QVBoxLayout *verticalLayout_2;
  QVBoxLayout *verticalLayout_3;
  QHBoxLayout *horizontalLayout_1;
  QHBoxLayout *horizontalLayout_2;

  QLabel *numNodes;
  QLabel *numGoodSamples;
  QLabel *numCollisionSamples;
  QLabel *numTotalSamples;

  QLabel *eq1;
  QLabel *eq2;
  QLabel *eq3;
  QLabel *eq4;

  QLabel *numNodesResults;
  QLabel *numGoodSamplesResults;
  QLabel *numCollisionSamplesResults;
  QLabel *numTotalSamplesResults;

  QPushButton *generateStatsButton;
  QPushButton *exitButton;

  QLabel *statusBar;
};

#endif  // STATSWINDOW_H
