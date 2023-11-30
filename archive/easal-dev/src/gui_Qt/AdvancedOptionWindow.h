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
#ifndef ADVANCEDOPTIONSWINDOW_H
#define ADVANCEDOPTIONSWINDOW_H
#include <QDialog>
#include <QInputDialog>
#include <QtWidgets>

class AdvancedOptionWindow : public QDialog {
  Q_OBJECT

 public:
  AdvancedOptionWindow();

 private Q_SLOTS:

  void accept();
  void cancel();

 private:
  QVBoxLayout *mainLayout;
  QGridLayout *checkBoxGroupLayout;
  QGridLayout *lineEditGroupLayout;
  QGridLayout *acceptExitGroupLayout;

  QGroupBox *checkBoxGroup;
  QGroupBox *lineEditGroup;
  QGroupBox *acceptExitGroup;

  QCheckBox *reverseWitnessCheckBox;
  QCheckBox *dynamicStepSizeAmongCheckBox;
  QCheckBox *useParticipatingAtomsZDistanceCheckBox;
  QCheckBox *ShortRangeSamplingCheckBox;
  QCheckBox *wholeCollisionsCheckBox;
  QCheckBox *dynamicStepSizeWithinCheckBox;
  QCheckBox *fourDRootNodeCheckBox;
  QCheckBox *reversePairDumbbellsCheckBox;

  QLabel *reverseWitnessLable;
  QLabel *dynamicStepSizeAmongLable;
  QLabel *useParticipatingAtomsZDistanceLable;
  QLabel *shortRangeSamplingLable;
  QLabel *wholeCollisionsLable;
  QLabel *dynamicStepSizeWithinLable;
  QLabel *fourDRootNodeLable;
  QLabel *reversePairDumbbellsLable;
  QLabel *participatingAtomsZDistanceLable;
  QLabel *participatingAtomIndexLowLable;
  QLabel *initial4DContactSeperationLowLable;
  QLabel *savePointsFrequencyLable;
  QLabel *participatingAtomIndexHighLable;
  QLabel *initial4DContactSeperationHighLable;

  QLineEdit *participatingAtomsZDistanceLineEdit;
  QLineEdit *participatingAtomIndexLowLineEdit;
  QLineEdit *initial4DContactSeperationLowLineEdit;
  QLineEdit *savePointsFrequencyLineEdit;
  QLineEdit *participatingAtomIndexHighLineEdit;
  QLineEdit *initial4DContactSeperationHighLineEdit;

  QPushButton *buttonAccept;
  QPushButton *buttonCancel;

  void setCheckBoxGroup();
  void setlineEditGroup();
  void setAcceptExitGroup();
};

#endif  // ADVANCEDOPTIONSWINDOW_H
