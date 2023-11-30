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
#ifndef SHAREDDATAGUI_H
#define SHAREDDATAGUI_H

#include "easalcore/ActiveConstraintRegion.h"
#include "easalcore/Atlas.h"
#include "easalcore/AtlasBuilder.h"
#include "easalcore/AtlasNode.h"
#include "easalcore/SaveLoader.h"
#include "easalcore/Settings.h"
#include "gui_Qt/CayleyPointRenderingWrapper.h"
#include "iostream"
#include "vector"

enum PathfindingPointType { START = 0, END = 1 };

using namespace std;

class SharedDataGUI {
 public:
  SharedDataGUI();
  AtlasNode *currentAtlasNode;
  int currentAtlasNodeID;
  bool atlasNodeWasSet;
  bool isSamplingStarted =
      false;  ///< bool that indicates if user opened setting window and
              ///< succesfuly started sampling
  ActiveConstraintRegion *currentACR;
  bool ACRWasSet;
  int currentCayleyPointID;
  PointType currentCayleyPointType;
  int currentFlipID;
  bool currentFlips[8];
  bool pathfindingSelectCayleyPointNow[2];
  bool pathfindingCayleyPointValueWasChanged[2];
  int pathfindingSelectedAtlasNodeID[2];
  PointType pathfindingSelectedCayleyPointType[2];
  int pathfindingSelectedCayleyPointID[2];
  vector<int> pathfindingFlipVector[2];
};

#endif  // SHAREDDATAGUI_H
