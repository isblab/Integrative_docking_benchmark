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
#ifndef CONSTRAINTSELECTIONBONDVIEWWIDGET_H
#define CONSTRAINTSELECTIONBONDVIEWWIDGET_H

#include <QtWidgets>
#include <algorithm>
#include <vector>

#include "ConstraintSelectionMoleculeViewWidget.h"

using namespace std;

struct ConnectionData {
  int atomIDFromMUA[6];
  int atomIDFromMUB[6];
  bool isConnected[6];
};

class ConstraintSelectionBondViewWidget : public QWidget {
 public:
  ConstraintSelectionBondViewWidget();
  void setConstraintSelectionMoleculeViewWidget(
      ConstraintSelectionMoleculeViewWidget *value);
  void setConnectionData(ConnectionData *value);

 private:
  void paintEvent(QPaintEvent *e) Q_DECL_OVERRIDE;
  //[ TIMER]
  /** @brief Timer that triggers screen updates once each m_updateTimeInMS*/
  QTimer m_timer;
  double m_updateTimeInMS;
  //[ TIMER]
  ConstraintSelectionMoleculeViewWidget *constraintSelectionMoleculeViewWidget;
  int fontSize = 18;
  ConnectionData *connectionData;
};

#endif  // CONSTRAINTSELECTIONBONDVIEWWIDGET_H
