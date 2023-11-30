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
#ifndef INPUTWINDOW_H
#define INPUTWINDOW_H

#include <easalcore/Settings.h>

#include <QDir>
#include <QtWidgets>
#include <cstdio>
#include <iostream>

#include "AdvancedOptionWindow.h"
#include "MessageWindow.h"

using namespace std;

class MainWindow;

class InputWindow : public QDialog {
  Q_OBJECT

 public:
  InputWindow();

 private Q_SLOTS:

  void exit();
  void accept();
  void browse1();
  void browse2();
  void browse3();
  void browse4();
  void advancedOptions();

 private:
  MessageWindow warningWindow{"Errors:"};

  QVBoxLayout *mainLayout;

  QGroupBox *data;
  QGroupBox *bounds;
  QGroupBox *stepSize;
  QGroupBox *acceptExit;

  QLineEdit *lineEditDataForMoleculeA;
  QLineEdit *lineEditDataForMoleculeB;
  QLineEdit *lineEditPredefinedInteractions;
  QLineEdit *lineEditDataDirectory;

  QPushButton *buttonBrowse1;
  QPushButton *buttonSetData1;
  QPushButton *buttonBrowse2;
  QPushButton *buttonSetData2;
  QPushButton *buttonBrowse3;
  QPushButton *buttonSetData3;
  QPushButton *buttonBrowse4;

  QPushButton *advancedOptionsButton;

  QPushButton *buttonAccept;
  QPushButton *buttonExit;

  QLineEdit *boundingTreshholdLowerBoundLambdaLineEdit2_1;
  QLineEdit *boundingTreshholdLowerBoundDeltaLineEdit2_2;
  QLineEdit *boundingTreshholdUpperBoundLambdaLineEdit3_1;
  QLineEdit *boundingTreshholdUpperBoundDeltaLineEdit3_2;
  QLineEdit *collisionTreshholdLowerBoundLambdaLineEdit5_1;
  QLineEdit *collisionTreshholdLowerBoundDeltaLineEdit5_2;
  QLineEdit *thetaLowLineEdit7_1;
  QLineEdit *thetaHighLineEdit7_2;

  QLineEdit *lineEditStepSize;

  QString fileFormatString;

  AdvancedOptionWindow *advancedOptionsWindowInstance;

  int numberOfIncorrectInputs = 0;

  void setDataGroupBox();
  void setBoundsGroupBox();
  void setStepSizeGroupBox();
  void setAcceptExitGroupBox();
  bool isInputCorrect();
  bool doesFileExist(string fileName);
  bool doesDirectoryExist(string directoryName);
};

#endif  // INPUTWINDOW_H
