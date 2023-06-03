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
#include "CustomSlider.h"

CustomSlider::CustomSlider() {
  //[ TIMER SETTINGS]
  // start argument determins update rate
  m_updateTimeInMS = 33;
  m_timer.start(m_updateTimeInMS);
  // repaint or update is function of QWidget that triggers paintEvent
  connect(&m_timer, SIGNAL(timeout()), this, SLOT(update()));
  //[ TIMER SETTINGS]
  //[ HANDEL FOCUS]
  // without focus policy mouse and keyboard events cannot work
  this->setFocusPolicy(Qt::StrongFocus);
  //[~HANDLE FOCUS]

  minValueLineEdit.setText("NaN");
  currentValueLineEdit.setText("NaN");
  maxValueLineEdit.setText("NaN");
  setFixedHeight(55);
  gridLayout.addWidget(&minValueLineEdit, 0, 0, 1, 1);
  gridLayout.addWidget(&currentValueLineEdit, 0, 1, 1, 1);
  gridLayout.addWidget(&maxValueLineEdit, 0, 2, 1, 1);
  setLayout(&gridLayout);
}

void CustomSlider::setCayleySpaceViewData(double max, double min, int nodeID,
                                          std::vector<double> validPoints) {
  if (validPoints.size() > 1) {
    this->max = max;
    this->min = min;
    dataWasSet = true;
    currentValue = validPoints[0];
    currentValueIndex = 0;
    currentNodeID = nodeID;
    values = validPoints;
    minValueLineEdit.setText("Min = " + QString::number(min));
    maxValueLineEdit.setText("Max = " + QString::number(max));
    currentValueLineEdit.setText("Current = " + QString::number(currentValue));
  } else if (validPoints.size() == 1) {
    this->max = validPoints[0] + 1;
    this->min = validPoints[0] - 1;
    dataWasSet = true;
    currentValue = validPoints[0];
    currentValueIndex = 0;
    currentNodeID = nodeID;
    values = validPoints;
    minValueLineEdit.setText("Min = " + QString::number(this->min));
    maxValueLineEdit.setText("Max = " + QString::number(this->max));
    currentValueLineEdit.setText("Current = " + QString::number(currentValue));
  }
  sliderType = CAYLEY_SAPCE_VIEW;
}

void CustomSlider::setPathfindingData(vector<int> atlasNodesIDs) {
  this->max = atlasNodesIDs.size() - 1;
  this->min = 0;
  dataWasSet = true;
  currentValue = 0;
  vector<double> tmpValidPoints;
  for (int i = 0; i < atlasNodesIDs.size(); i++) {
    tmpValidPoints.push_back(i);
  }
  values = tmpValidPoints;
  displayValues = atlasNodesIDs;

  minValueLineEdit.setText("Min = " + QString::number(atlasNodesIDs[0]));
  maxValueLineEdit.setText(
      "Max = " + QString::number(atlasNodesIDs[atlasNodesIDs.size() - 1]));
  currentValueLineEdit.setText("Current = " +
                               QString::number(atlasNodesIDs[0]));

  sliderType = PATHFINDING;
}

double CustomSlider::getCurrentValue() {
  if (sliderType == CAYLEY_SAPCE_VIEW) {
    return currentValue;
  } else if (sliderType == PATHFINDING) {
    // Returns atlas node ID with index of current value
    return displayValues[currentValue];
  }
}

void CustomSlider::updateCurrentValue() {
  QPoint mousePosition = mapFromGlobal(QCursor::pos());
  if (mousePosition.y() > minValueLineEdit.height() && values.size() > 0) {
    int xVertLineCoord =
        (this->width() - 2 * xOffset) / (max - min) * (values[0] - min) +
        xOffset;
    double delta = abs(xVertLineCoord - mousePosition.x());
    int minDeltaIndex = 0;
    for (int i = 1; i < values.size(); i++) {
      xVertLineCoord =
          (this->width() - 2 * xOffset) / (max - min) * (values[i] - min) +
          xOffset;
      if (abs(xVertLineCoord - mousePosition.x()) < delta) {
        delta = abs(xVertLineCoord - mousePosition.x());
        minDeltaIndex = i;
      }
    }
    currentValue = values[minDeltaIndex];
    if (sliderType == CAYLEY_SAPCE_VIEW) {
      currentValueLineEdit.setText("Current = " +
                                   QString::number(currentValue));
    } else if (sliderType == PATHFINDING) {
      currentValueLineEdit.setText(
          "Current = " + QString::number(displayValues[currentValue]));
    }
  }
}

void CustomSlider::paintEvent(QPaintEvent *event) {
  if (dataWasSet) {
    QPainter painter(this);
    double yLineCoord = 45;
    //        double yLineCoord = this->height()/2+minValueLineEdit.height();
    painter.drawLine(xOffset, yLineCoord, this->width() - xOffset, yLineCoord);
    painter.drawLine(xOffset, yLineCoord, this->width() - xOffset, yLineCoord);
    for (int i = 0; i < values.size(); i++) {
      double xVertLineCoord =
          (this->width() - 2 * xOffset) / (max - min) * (values[i] - min) +
          xOffset;
      // Draw vertical marks
      painter.drawLine(xVertLineCoord, yLineCoord + sliderMarkHight,
                       xVertLineCoord, yLineCoord - sliderMarkHight);
    }
    // Handle mouse input
    if (mouseWasJustClicked) {
      updateCurrentValue();
      mouseWasJustClicked = false;
    }
    // Draw slider
    double xSliderCoord =
        (this->width() - 2 * xOffset) / (max - min) * (currentValue - min) +
        xOffset;
    painter.drawRect(xSliderCoord - 2, yLineCoord - sliderHight, 4,
                     sliderHight * 2);
  }
}

void CustomSlider::mousePressEvent(QMouseEvent *event) {
  mouseWasJustClicked = true;
}

void CustomSlider::mouseReleaseEvent(QMouseEvent *event) {}

void CustomSlider::keyPressEvent(QKeyEvent *event) {
  if (dataWasSet) {
    if (event->key() == Qt::Key_Right) {
      int newIndex = ((currentValueIndex - 1) + values.size()) % values.size();
      currentValue = values[newIndex];
      currentValueIndex = newIndex;
      if (sliderType == CAYLEY_SAPCE_VIEW) {
        currentValueLineEdit.setText("Current = " +
                                     QString::number(currentValue));
      } else if (sliderType == PATHFINDING) {
        currentValueLineEdit.setText(
            "Current = " + QString::number(displayValues[currentValue]));
      }
    }
    if (event->key() == Qt::Key_Left) {
      int newIndex = ((currentValueIndex + 1) + values.size()) % values.size();
      currentValue = values[newIndex];
      currentValueIndex = newIndex;
      if (sliderType == CAYLEY_SAPCE_VIEW) {
        currentValueLineEdit.setText("Current = " +
                                     QString::number(currentValue));
      } else if (sliderType == PATHFINDING) {
        currentValueLineEdit.setText(
            "Current = " + QString::number(displayValues[currentValue]));
      }
    }
  }
}
