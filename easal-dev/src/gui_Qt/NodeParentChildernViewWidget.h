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
#ifndef NODEPARENTCHILDERNVIEWWIDGET_H
#define NODEPARENTCHILDERNVIEWWIDGET_H

#include <easalcore/Atlas.h>
#include <easalcore/SaveLoader.h>

#include <QColor>
#include <QFont>
#include <QKeyEvent>
#include <QMouseEvent>
#include <QPainter>
#include <QTimer>
#include <QVector2D>
#include <QtWidgets>
#include <algorithm>
#include <iostream>
#include <unordered_map>
#include <vector>

#include "RealizationViewWidgetCornerVersion.h"
#include "SharedDataGUI.h"

using namespace std;

enum RelationToRoot { ROOT, CHILD, PARENT };

/**
 * @brief The DimensionData struct: holds data about nodes of
 * NodeParentChildernViewWidget with same dimension (0 to 5)
 */
struct DimensionData {
  RelationToRoot relationToRoot;
  vector<int> nodesIDs;
  int distanceToParent;
  int distanceToChild;
  int verticalOffset;
  int verticalDistanceBetweenBoxes;
  vector<vector<int>> nodesChildren;
  vector<vector<int>> nodesParents;
  unordered_map<int, int>
      mapIDToIndex;  ///< maps IDs of nodes into correspondent indecies of
                     ///< vector<int> nodesIDs, used to render lines between
                     ///< nodes fast
  QSize boxSize;
};

class NodeParentChildernViewWidget : public QWidget {
  Q_OBJECT
 public:
  NodeParentChildernViewWidget(SharedDataGUI *sharedData);
  void setAtlasNode(AtlasNode *currentAtlasNode);
  void setComponents(Atlas *atlas, SaveLoader *saveLoader);
  void updatePositions();

 private:
  //[ ATLAS NODE]
  AtlasNode *currentAtlasNode;
  int currentAtlasNodeID;
  //[~ATLAS NODE]
  //[ MOUSE]
  bool LMBisPressed;
  bool RMBisPressed;
  bool mouseWasJustPressed;
  bool mouseWasReleased;
  /** @brief We use this vector to move viewport with mouse*/
  QVector2D lastMouseCoords;
  /** @brief This point folds coordinates of last mouse click */
  QPoint mouseLastClickCoordinates;
  //[~MOUSE]
  //[ TIMER]
  QTimer *timer;
  int updateTimeInMS;
  //[~TIMER]
  //[ VIEWPORT]
  QVector2D viewTranslation;
  QVector2D viewScale;
  //[~VIEWPORT]
  //[ WIDGET]
  bool widgetIsActive;
  int activeWidgetID;
  bool nodeWasSet;
  //[~WIDGET]
  //[ DATA]
  Atlas *atlas;
  SaveLoader *saveLoader;
  //[~DATA]
  //[ SIZES]
  int xSize;
  int ySize;
  //[~SIZES]
  //[ MISC]
  Qt::GlobalColor colors[6];
  SharedDataGUI *sharedData;
  //[~MISC]
  DimensionData dimensionData[6];
  int rootDimension;

  void paintEvent(QPaintEvent *event);
  /** @brief We override mousePressEvent to handle mouse input */
  void mousePressEvent(QMouseEvent *event);
  void mouseReleaseEvent(QMouseEvent *event);
  /** @brief We override keyPressEvent to handle keyboard input */
  void keyPressEvent(QKeyEvent *event);
  void renderRootDimension(QPainter *painter);
  void renderChildrenDimensions(QPainter *painter);
  void renderParentsDimensions(QPainter *painter);
  void renderBox(int x, int y, int width, int height, int ID,
                 Qt::GlobalColor color, QPainter *painter);
  void renderIDString(int x, int y, int ID, QPainter *painter);
  void renderDimensions(QPainter *painter);
  void mouseNodePickingCheck(int x, int y, int ID);
  void clearDimensionData();
  vector<int> getChildernIDs(AtlasNode *node, int currentDim);
  vector<int> getParentsIDs(AtlasNode *node, int currentDim);
  void updateViewTranslationAndScale(QPainter *painter);
  void updateDimensionData(vector<int> nodesIDs, RelationToRoot relationToRoot);
};

#endif  // NODEPARENTCHILDERNVIEWWIDGET_H
