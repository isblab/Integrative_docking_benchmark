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
#ifndef ATLASVIEW_H
#define ATLASVIEW_H

#include <easalcore/Atlas.h>
#include <easalcore/SaveLoader.h>

#include <QKeyEvent>
#include <QMouseEvent>
#include <QPainter>
#include <QTimer>
#include <QVector2D>
#include <QWidget>
#include <algorithm>
#include <iostream>

#include "RealizationViewWidgetCornerVersion.h"
#include "SharedDataGUI.h"

using namespace std;

/**
 * @brief The FoldData struct: holds data about atlas view widget node
 * (basically rectangle that represents atlas node)
 */
struct FoldData {
  bool isFolded;
  bool wasSet;
  int parents[6];
};

/**
 * @brief The SelectedNode struct: holds data about currently selected atlas
 * view widget node (basically rectangle that represents atlas node)
 */
struct SelectedNode {
  int ID;
  int parents[6];
};
/**
 * @brief The LastAction enum: holds names of actions that we need to fold node
 * in atlasViewWidget
 */
enum LastAction { NODE_PICKED, FOLD_DUPLICATES, FOLD_DIMENSION };

class AtlasViewWidget : public QWidget {
  Q_OBJECT
 public:
  AtlasViewWidget(SharedDataGUI *sharedData);
  Atlas *atlas;
  SaveLoader *saveLoader;
  /** @brief Showall allow to render all nodes (including those with no good
     orientation
      TODO add function togglShowAll */
  bool showall;
  void setComponents(Atlas *atlas, SaveLoader *saveLoader);
  /**
   * @brief getCurrentNode: allows to get pointer to node that was selected by
   * mouse
   * @return
   */
  AtlasNode *getCurrentNode();

 private Q_SLOTS:
  void updateAtlas();

 private:
  //[ SIZE]
  // Stays for size of rendered node
  int rectXSize;
  int rectYSize;
  // Stays for horisontal distance between rect
  int xDistance;
  // Stays for vertical distance between rect
  int yDistance;
  // Stays for horizontal distance between block of 5dim nodes
  int xBlockDistance;
  // Stays for vertical distance between blocks of 5dim nodes
  int yBlockDistance;
  // Stays for horisontal offset due to presence of root node
  int xRootNodeOffset;
  // Stays for vertical offset due to presence of root node
  int yRootNodeOffset;
  // Tracks X coord of last root node
  int xRootNode;
  // Tracks Y coord of last root node
  int yRootNode;
  // Increases yDistance to render text for node id propely when nodes are
  // folded
  int yFoldBonus;
  //[~SIZE]
  //[ OFFSETS]
  // We check if setData funtion was called and atlas and saveloader were
  // initialized
  int maxNodeID;
  // Keeps track of max Y value
  int maxY;
  int maxYCurrentBlock;
  int y4LastRow;
  //[~OFFSETS]
  //[ VIEWPORT]
  QVector2D viewTranslation;
  QVector2D viewScale;
  //[~VIEWPORT]
  //[ FOLD]
  bool foldArray[6];
  // folds nodes that have index smaller then current max nodeID (that means
  // that node was rendered somewhere else befor
  bool foldDuplicatedNodes;
  bool unfold;
  int currentParents[6];
  vector<bool> foldBoolData;
  vector<FoldData> foldData;
  //[~FOLD]
  //[ MOUSE]
  bool LMBisPressed;
  bool RMBisPressed;
  bool mouseWasJustPressed;
  bool mouseWasReleased;
  /** @brief We use this vector to move viewport with mouse*/
  QVector2D lastMouseCoords;
  /** @brief This point folds coordinates of last mouse click */
  QPoint mouseLastClickCoordinates;
  LastAction lastAction;
  //[~MOUSE]
  //[ TIMER]
  QTimer *timer;
  QTimer *timer2;
  int updateTimeInMS;
  //[~TIMER]
  //[ MISC]
  SelectedNode selectedNode;
  /** @brief Each time we renrder atlas we increase currentFrameID by 1 */
  int currentFrameID;
  /** @brief variable above holds ID of node that we cliked on last recently */
  int lastSelectedNodeID;
  /** @brief We use updateTimeInMS together with timer to update screen */
  bool dataWasSet;
  bool isWidgetActive;
  int activeWidgetID;
  Qt::GlobalColor colors[6];
  /** @brief This bool indicates if we clicked on node recently */
  bool nodeWasSelected;
  /** @brief These vectors used to "move" and "zoom" our viewport */

  //[~MISC]
  SharedDataGUI *sharedData;
  /** @brief For all nodes we calculate if we fold them or not */
  void updateFoldData();
  void paintEvent(QPaintEvent *event);
  /** @brief We override mousePressEvent to handle mouse input */
  void mousePressEvent(QMouseEvent *event);
  void mouseReleaseEvent(QMouseEvent *event);
  /** @brief We override keyPressEvent to handle keyboard input */
  void keyPressEvent(QKeyEvent *event);
  /** @brief We override wheelEvent to zoom viewport */
  void wheelEvent(QWheelEvent *event);
  /** @brief Render rect or 3 rect (if isFolded) at coords x y, color and string
     id printe on top
      @arg painter and function itself can only be called with painter declared
     in paintEvent */
  void renderNodeAndCheckIfNodeWasSelected(int x, int y, Qt::GlobalColor color,
                                           int id, int dimension,
                                           QPainter *painter, bool isFolded);

  void renderLine(int x1, int y1, int x2, int y2, QPainter *painter);

  /**
   * @brief getChildren: allows to select childrens from neighbors vector
   * "nodes"
   * @param nodes: vector that hold indecies of neighbors
   * @param currentDimension: int that represents dimension of node that have
   * all "nodes" as neighbors
   * @return vector of ID that contains only children
   */
  std::vector<int> getChildren(std::vector<int> nodes, int currentDimension);
  /**
   * @brief getRootIndices: allows to get first nodes from tlas
   * @return
   */
  std::vector<int> getRootIndices();
  /**
   * @brief renderNodesRecursively: allow to go trough all root indecies and
   * render them and all tree of their children, grand children etc.
   * @param painter
   */
  void renderNodesRecursively(QPainter *painter);
  /**
   * @brief updateViewTranslationAndScale: allows to move viewport around atlas
   * that we render
   * @param painter
   */
  void updateViewTranslationAndScale(QPainter *painter);
  /**
   * @brief isIndexBiggerThanMaxID: checks if "index" is bigger than index of
   * all nodes tha we rendered before (for some folding reasons)
   * @param index: int that represnt index of node in atlas
   * @return
   */
  bool isIndexBiggerThanMaxID(int index);
  /**
   * @brief traverseNodes: renders all nodes that correspondent to indicies and
   * recursivly calls itself on their childern
   * @param indicies: int vector that represent indecies of nodes in atlas
   * @param dimension: dimension of these nodes
   * @param xPreviousDimension: x coordinate of rectangle that represented node
   * of previous demension
   * @param yPreviousDimension: y coordinate of rectangle that represented node
   * of previous demension
   * @param painter
   */
  void traverseNodes(vector<int> indicies, int dimension,
                     int xPreviousDimension, int yPreviousDimension,
                     QPainter *painter);
  /**
   * @brief compareParents: checks if currently selected node have same parents
   * as in data
   * @param data: struct that hold indecies of parent nodes
   * @param dimension: int that represent dimension of node that data refers to
   * @return
   */
  bool compareParents(SelectedNode data, int dimension);
  /**
   * @brief resetFoldData: cleans all fold related data
   */
  void resetFoldData();
  bool compareParents(FoldData data, int dimension);
};

#endif  // ATLASVIEW_H
