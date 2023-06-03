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
#include "BondViewWidget.h"

BondViewWidget::BondViewWidget() {
  //[ TIMER SETTINGS]
  // start argument determins update rate
  m_updateTimeInMS = 33;
  m_timer.start(m_updateTimeInMS);
  // repaint is function of QWidget that triggers paintEvent
  connect(&m_timer, SIGNAL(timeout()), this, SLOT(update()));
  //[~TIMER SETTINGS]
  //[ DRAWING SIZES]
  xAtomDistance = 100;
  yAtomDistance = 37;
  atomCircleRadius = 8;
  this->setMinimumHeight(225);
  this->setMinimumWidth(130);
  //[~DRAWING SIZES]
  dataWasSet = false;
}

void BondViewWidget::setBondData(AtomBondData *atlasBondData) {
  dataWasSet = true;

  bondSolid.clear();
  bondSolid.resize(atlasBondData->bondSolid.size());
  bondVariable.clear();
  bondVariable.resize(atlasBondData->bondVariable.size());
  // Get all solid connections
  for (int i = 0; i < atlasBondData->bondSolid.size(); i++) {
    bondSolid[i] = atlasBondData->bondSolid[i];
  }
  // Get all variable size connections
  for (int i = 0; i < atlasBondData->bondVariable.size(); i++) {
    bondVariable[i] = atlasBondData->bondVariable[i];
  }
  IDA.clear();
  IDB.clear();
  // Find which atoms participate in these connections
  for (int i = 0; i < bondSolid.size(); i++) {
    if (std::find(IDA.begin(), IDA.end(), bondSolid[i].first) == IDA.end()) {
      IDA.push_back(bondSolid[i].first);
    }
    if (std::find(IDB.begin(), IDB.end(), bondSolid[i].second) == IDB.end()) {
      IDB.push_back(bondSolid[i].second);
    }
  }
  for (int i = 0; i < bondVariable.size(); i++) {
    if (std::find(IDA.begin(), IDA.end(), bondVariable[i].first) == IDA.end()) {
      IDA.push_back(bondVariable[i].first);
    }
    if (std::find(IDB.begin(), IDB.end(), bondVariable[i].second) ==
        IDB.end()) {
      IDB.push_back(bondVariable[i].second);
    }
  }
  // And sort them to render
  // in order from smallest index to highest
  std::sort(IDA.begin(), IDA.end());
  std::sort(IDB.begin(), IDB.end());
}

void BondViewWidget::paintEvent(QPaintEvent *event) {
  if (dataWasSet) {
    int numNodesA = 0;
    int numNodesB = 0;
    QPainter painter(this);
    // First we render bonds because they are one layer below cirles that
    // represent atoms
    for (int i = 0; i < bondSolid.size(); i++) {
      // We already know x coordinates of atoms that will render because they
      // positioned in two rows So below we find y coordinates (y coordinate
      // depends on position of atoms in IDA or IDB) And these atom indecies we
      // get from bonds that we iterate trough
      int y1 =
          yAtomDistance *
          (0.5 + std::distance(IDA.begin(), std::find(IDA.begin(), IDA.end(),
                                                      bondSolid[i].first)));
      int y2 =
          yAtomDistance *
          (0.5 + std::distance(IDB.begin(), std::find(IDB.begin(), IDB.end(),
                                                      bondSolid[i].second)));
      painter.setPen(
          QPen(Qt::darkGray, 3, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
      painter.drawLine(atomCircleRadius * 3, y1,
                       xAtomDistance + atomCircleRadius, y2);
    }
    for (int i = 0; i < bondVariable.size(); i++) {
      int y1 =
          yAtomDistance *
          (0.5 + std::distance(IDA.begin(), std::find(IDA.begin(), IDA.end(),
                                                      bondVariable[i].first)));
      int y2 =
          yAtomDistance *
          (0.5 + std::distance(IDB.begin(), std::find(IDB.begin(), IDB.end(),
                                                      bondVariable[i].second)));
      painter.setPen(
          QPen(Qt::darkGray, 3, Qt::DashLine, Qt::RoundCap, Qt::RoundJoin));
      painter.drawLine(atomCircleRadius * 3, y1,
                       xAtomDistance + atomCircleRadius, y2);
    }
    // Now we render cirles that represent atoms
    for (int i = 0; i < 6; i++) {
      if (i < IDA.size()) {
        painter.setPen(
            QPen(Qt::black, 3, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
        painter.setBrush(Qt::yellow);
        painter.drawEllipse(
            QPoint(atomCircleRadius * 3, yAtomDistance * (0.5 + i)),
            atomCircleRadius * 2, atomCircleRadius * 2);
        painter.drawText(atomCircleRadius * 3, yAtomDistance * (0.5 + i),
                         QString::number(IDA[i]));
      } else {
        painter.setPen(
            QPen(Qt::black, 1, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
        painter.setBrush(Qt::lightGray);
        painter.drawEllipse(
            QPoint(atomCircleRadius * 3, yAtomDistance * (0.5 + i)),
            atomCircleRadius * 2, atomCircleRadius * 2);
      }
    }
    for (int i = 0; i < 6; i++) {
      if (i < IDB.size()) {
        painter.setPen(
            QPen(Qt::black, 3, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
        painter.setBrush(Qt::red);
        painter.drawEllipse(
            QPoint(xAtomDistance + atomCircleRadius, yAtomDistance * (0.5 + i)),
            atomCircleRadius * 2, atomCircleRadius * 2);
        painter.drawText(xAtomDistance + atomCircleRadius,
                         yAtomDistance * (0.5 + i), QString::number(IDB[i]));
      } else {
        painter.setPen(
            QPen(Qt::black, 1, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
        painter.setBrush(Qt::lightGray);
        painter.drawEllipse(
            QPoint(xAtomDistance + atomCircleRadius, yAtomDistance * (0.5 + i)),
            atomCircleRadius * 2, atomCircleRadius * 2);
      }
    }
  }
}
