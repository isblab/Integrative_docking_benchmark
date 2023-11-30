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
#ifndef CUSTOMSLIDER_H
#define CUSTOMSLIDER_H

#include <QGridLayout>
#include <QKeyEvent>
#include <QLabel>
#include <QPainter>
#include <QTimer>
#include <QWidget>
#include <iostream>
#include <vector>

using namespace std;

// Slider for CayleySpaceViewWidget to select values of 4th and 5th constrain
class CustomSlider : public QWidget {
  Q_OBJECT
 public:
  CustomSlider();
  void setCayleySpaceViewData(double max, double min, int nodeID,
                              std::vector<double> validPoints);
  void setPathfindingData(vector<int> atlasNodesIDs);
  double getCurrentValue();

 private:
  enum SliderType { CAYLEY_SAPCE_VIEW, PATHFINDING };
  SliderType sliderType;
  //[ TIMER]
  /** @brief Timer that triggers screen updates once each m_updateTimeInMS*/
  QTimer m_timer;
  double m_updateTimeInMS;
  //[ TIMER]
  int xOffset = 5;
  int sliderMarkHight = 3;
  int sliderHight = 5;
  int currentValueIndex;
  int currentNodeID = -1;
  double max, min;
  bool mouseWasJustClicked = false;
  std::vector<double> values;
  std::vector<int> displayValues;

  double currentValue;
  bool dataWasSet = false;
  QLabel minValueLineEdit;
  QLabel maxValueLineEdit;
  QLabel currentValueLineEdit;
  QGridLayout gridLayout;
  void updateCurrentValue();
  void paintEvent(QPaintEvent *event) Q_DECL_OVERRIDE;
  void mousePressEvent(QMouseEvent *event) Q_DECL_OVERRIDE;
  void mouseReleaseEvent(QMouseEvent *event) Q_DECL_OVERRIDE;
  void keyPressEvent(QKeyEvent *event) Q_DECL_OVERRIDE;
};

#endif  // CUSTOMSLIDER_H
