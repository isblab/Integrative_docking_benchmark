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
#ifndef ABOUT_WINDOW_H
#define ABOUT_WINDOW_H

#include <easalcore/Atlas.h>
#include <easalcore/BasicStatistics.h>
#include <easalcore/SaveLoader.h>
#include <easalcore/Settings.h>
#include <gui_Qt/SharedDataGUI.h>

#include <QtWidgets/QTextBrowser>
#include <QtWidgets>
#include <cstdio>
#include <iostream>

using namespace std;

class AboutWindow : public QDialog {
  Q_OBJECT

 public:
  AboutWindow();

 private Q_SLOTS:
  void exit();

 private:
  // QTextEdit *License;
  QTextBrowser *License;

  QWidget *centralWidget;
  QWidget *verticalLayoutWidget_1;
  QWidget *verticalLayoutWidget_2;
  QWidget *verticalLayoutWidget_3;
  QWidget *verticalLayoutWidget_4;
  QWidget *horizontalLayoutWidget_1;
  QWidget *horizontalLayoutWidget_2;
  QWidget *horizontalLayoutWidget_3;

  QVBoxLayout *verticalLayout_1;
  QVBoxLayout *verticalLayout_2;
  QVBoxLayout *verticalLayout_3;
  QVBoxLayout *verticalLayout_4;
  QHBoxLayout *horizontalLayout_1;
  QHBoxLayout *horizontalLayout_2;
  QHBoxLayout *horizontalLayout_3;

  QLabel *introduction;
  QLabel *dimensionHeader;

  QLabel *minValuesHeader;

  QLabel *maxValuesHeader;

  QLabel *averageValuesHeader;

  QLabel *dimension[6];

  QLabel *minValues[6];

  QLabel *maxValues[6];

  QLabel *averageValues[6];

  QLabel *statusBar;

  QPushButton *generateCayleyStatsButton;

  SharedDataGUI *sharedData;

  QPushButton *exitButton;
};

#endif
