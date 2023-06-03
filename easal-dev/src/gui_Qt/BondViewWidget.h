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
#ifndef BONDVIEWWIDGET_H
#define BONDVIEWWIDGET_H

#include <QPainter>
#include <QTimer>
#include <QVector4D>
#include <QWidget>
#include <algorithm>
#include <vector>

using namespace std;

struct AtomBondViewWidgetData {
  /** @brief position is coordinate of circle that repsents atom relative to
   * widget */
  QPoint position;
  QVector4D color;
  /** @brief ID si number of atom in molecule */
  int ID;
};

struct AtomBondData {
  /**
   * @brief bondSolid: vector of pair of int that holds ID of atom where fixed
   * size bond starts and ends
   */
  vector<std::pair<int, int>> bondSolid;
  /**
   * @brief bondSolid: vector of pair of int that holds ID of atom where variale
   * size bond starts and ends
   */
  vector<std::pair<int, int>> bondVariable;
};

class BondViewWidget : public QWidget {
  Q_OBJECT
 public:
  BondViewWidget();
  void setBondData(AtomBondData *atlasBondData);

 public Q_SLOTS:

 private:
  /** @brief xAtomDistance determines horizontal distance between atoms of
   * different molecules */
  int xAtomDistance;
  /** @brief yAtomDistance determines vertical distance between atoms of same
   * molecule */
  int yAtomDistance;
  int atomCircleRadius;
  //[ TIMER]
  /** @brief Timer that triggers screen updates once each m_updateTimeInMS*/
  QTimer m_timer;
  double m_updateTimeInMS;
  //[~TIMER]
  //[ BONDS]
  /** bonds hold pairs of indicies of atoms from molecule A and B */
  vector<std::pair<int, int>> bondSolid;
  vector<std::pair<int, int>> bondVariable;
  //[~BONDS]
  //[ IDs OF ATOMS]
  /** @brief IDA and IDB hold sorted IDs of atoms that participae in connection
   * in molecule A and B */
  vector<int> IDA;
  vector<int> IDB;
  //[~IDs OF ATOMS]
  bool dataWasSet;
  vector<AtomBondViewWidgetData> atomsData;
  void paintEvent(QPaintEvent *event) Q_DECL_OVERRIDE;
};

#endif  // BONDVIEWWIDGET_H
