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
#include "ConstraintSelectionBondViewWidget.h"

ConstraintSelectionBondViewWidget::ConstraintSelectionBondViewWidget() {
  setAttribute(Qt::WA_TransparentForMouseEvents);
  //[ TIMER SETTINGS]
  m_updateTimeInMS = 33;
  m_timer.start(m_updateTimeInMS);
  connect(&m_timer, SIGNAL(timeout()), this, SLOT(update()));
  //[ TIMER SETTINGS]
}

void ConstraintSelectionBondViewWidget::
    setConstraintSelectionMoleculeViewWidget(
        ConstraintSelectionMoleculeViewWidget *value) {
  constraintSelectionMoleculeViewWidget = value;
}

void ConstraintSelectionBondViewWidget::paintEvent(QPaintEvent *e) {
  QPainter painter(this);
  QFont font;
  font.setPixelSize(fontSize);
  painter.setFont(font);
  QPen pen;
  pen.setWidth(4);
  painter.setPen(pen);

  vector<QVector4D> pointsProjectionCoords;

  for (int i = 0; i < MOLECULAR_UNITS_NUMBER; i++) {
    for (int j = 0;
         j <
         constraintSelectionMoleculeViewWidget->translationMatrices[i].size();
         j++) {
      QVector4D pointCoords(
          constraintSelectionMoleculeViewWidget->translationMatrices[i][j]
              .data()[12],
          constraintSelectionMoleculeViewWidget->translationMatrices[i][j]
              .data()[13],
          constraintSelectionMoleculeViewWidget->translationMatrices[i][j]
              .data()[14],
          constraintSelectionMoleculeViewWidget->translationMatrices[i][j]
              .data()[15]);
      QVector4D pointProjectionCoords =
          constraintSelectionMoleculeViewWidget->camera.getLastCameraMatrix() *
          pointCoords;
      pointProjectionCoords.setW(j);
      pointsProjectionCoords.push_back(pointProjectionCoords);
    }
  }

  //    sort(pointsProjectionCoords.begin(),pointsProjectionCoords.end(),[](const
  //    QVector4D &a,QVector4D &b) -> bool
  //    {
  //        return a.z()>b.z();
  //    });

  for (int i = 0; i < pointsProjectionCoords.size(); i++) {
    int xCoord = (1 + pointsProjectionCoords[i].x()) / 2.0 * width();
    int yCoord =
        height() - (1 + pointsProjectionCoords[i].y()) / 2.0 * height();
    painter.drawText(xCoord, yCoord,
                     QString::number(pointsProjectionCoords[i].w()));
  }

  for (int i = 0; i < 6; i++) {
    if (connectionData->isConnected[i]) {
      int numberOfAtomsInMUA = constraintSelectionMoleculeViewWidget
                                   ->translationMatrices[MOLECULAR_UNIT_A]
                                   .size();

      int xCoord1 =
          (1 + pointsProjectionCoords[connectionData->atomIDFromMUA[i]].x()) /
          2.0 * width();
      int yCoord1 =
          height() -
          (1 + pointsProjectionCoords[connectionData->atomIDFromMUA[i]].y()) /
              2.0 * height();
      int xCoord2 =
          (1 + pointsProjectionCoords[numberOfAtomsInMUA +
                                      connectionData->atomIDFromMUB[i]]
                   .x()) /
          2.0 * width();
      int yCoord2 =
          height() -
          (1 + pointsProjectionCoords[numberOfAtomsInMUA +
                                      connectionData->atomIDFromMUB[i]]
                   .y()) /
              2.0 * height();
      painter.drawLine(xCoord1, yCoord1, xCoord2, yCoord2);
    }
  }
}

void ConstraintSelectionBondViewWidget::setConnectionData(
    ConnectionData *value) {
  connectionData = value;
}
