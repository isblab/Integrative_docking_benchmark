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
#include "AtlasViewWidget.h"

#include <vector>
#define FRAME_WASNT_RENDERED -100

AtlasViewWidget::AtlasViewWidget(SharedDataGUI *sharedData) {
  this->sharedData = sharedData;
  //[ HANDEL TIMER]
  updateTimeInMS = 1000.0 / 30.0;
  timer = new QTimer(this);
  timer2 = new QTimer(this);
  // repaint is function of QWidget that triggers paintEvent
  connect(timer, SIGNAL(timeout()), this, SLOT(repaint()));
  connect(timer2, SIGNAL(timeout()), this, SLOT(updateAtlas()));

  // start argument determins update rate
  timer->start(updateTimeInMS);
  timer2->start(3000);
  //[~HANDEL TIMER]
  //[ HANDEL FOCUS]
  // without focus policy mouse and keyboard events cannot work
  this->setFocusPolicy(Qt::StrongFocus);
  //[~HANDLE FOCUS]
  //[ HANDLE COLORS]
  colors[5] = Qt::yellow;
  colors[4] = Qt::red;
  colors[3] = Qt::green;
  colors[2] = Qt::cyan;
  colors[1] = Qt::darkBlue;
  colors[0] = Qt::darkGray;
  //[~HANDLE COLORS]
  //[ HANDLE MOUSE]
  LMBisPressed = false;
  RMBisPressed = false;
  mouseWasJustPressed = false;
  mouseWasReleased = true;
  //[~HANDLE MOUSE]
  //[ HANDLE WIDGET]
  isWidgetActive = false;
  //[ HANDLE WIDGET]
  //[ HANDLE MISC VARS]
  rectXSize = 10;
  rectYSize = 5;
  maxNodeID = 0;
  maxY = -yBlockDistance;
  xDistance = 20;
  yDistance = 15;
  xBlockDistance = 6 * xDistance + 5 * rectXSize;
  yBlockDistance = 30;
  xRootNodeOffset = rectXSize + 3 * xDistance;
  yRootNodeOffset = 40;
  xRootNode = 10;
  yRootNode = 0;
  yFoldBonus = 8;

  currentFrameID = 0;

  showall = true;
  isWidgetActive = true;
  activeWidgetID = 1;

  dataWasSet = false;

  unfold = false;
  nodeWasSelected = false;

  viewTranslation = QVector2D(0, 0);
  viewScale = QVector2D(1, 1);

  for (int i = 0; i < 6; i++) {
    foldArray[i] = false;
  }

  foldDuplicatedNodes = false;
  //[~HANDLE MISC VARS]
}

void AtlasViewWidget::setComponents(Atlas *atlas, SaveLoader *saveLoader) {
  // this->atlas=atlas;
  this->atlas = new Atlas();
  this->saveLoader = saveLoader;
  dataWasSet = true;
}

AtlasNode *AtlasViewWidget::getCurrentNode() {
  if (nodeWasSelected) {
    return atlas->getNodes()[lastSelectedNodeID];
  } else {
    return NULL;
  }
}

void AtlasViewWidget::renderNodesRecursively(QPainter *painter) {
  Settings *sett = Settings::getInstance();
  bool D = false;
  maxNodeID = -1;
  maxY = 0;
  updateViewTranslationAndScale(painter);
  if (sett->RootNodeCreation.dimension_of_rootNodes == 5) {
    traverseNodes(getRootIndices(), 5, 0, 0, painter);
  } else {
    traverseNodes(getRootIndices(), 4, xRootNode, yBlockDistance, painter);
  }
}

void AtlasViewWidget::updateViewTranslationAndScale(QPainter *painter) {
  bool D = false;
  int mouseX = mapFromGlobal(QCursor::pos()).x();
  int mouseY = mapFromGlobal(QCursor::pos()).y();
  if (mouseWasReleased && LMBisPressed) {
    lastMouseCoords.setX(mouseX);
    lastMouseCoords.setY(mouseY);
    mouseWasReleased = false;
  } else if (LMBisPressed) {
    //        viewTranslation.setX(viewTranslation.x()+(mouseX-lastMouseCoords.x()));
    viewTranslation.setY(
        min(viewTranslation.y() + (mouseY - lastMouseCoords.y()), 0.0f));
    lastMouseCoords.setX(mouseX);
    lastMouseCoords.setY(mouseY);
	mouseWasJustPressed = false;
  }
  painter->translate(viewTranslation.x(), viewTranslation.y());
  painter->scale(viewScale.x(), viewScale.y());
  if (D)
    cout << " Scale x " << viewScale.x() << " Scale y " << viewScale.y()
         << endl;
  if (D)
    cout << " view translation x " << viewTranslation.x() << " y "
         << viewTranslation.y() << endl;
}

void AtlasViewWidget::traverseNodes(vector<int> indicies, int dimension,
                                    int xPreviousDimension,
                                    int yPreviousDimension, QPainter *painter) {
  Settings *sett = Settings::getInstance();
  // it is hard to comment this function because it contains a lot of different
  // offsets that hard do describe using language so good luck editig it
  bool D = false;
  int xNodeCoordinate;
  int yNodeCoordinate;
  bool fold = false;
  bool unfoldSet = false;

  updateFoldData();

  for (int i = 0; i < indicies.size(); i++) {
    int nodeID = (*atlas)[indicies[i]]->getID();
    currentParents[dimension] = nodeID;
    vector<int> children =
        getChildren((*atlas)[indicies[i]]->getConnection(), dimension);
    bool sameParents = compareParents(foldData[nodeID], dimension);
    bool foldDataBool = foldData[nodeID].isFolded && sameParents;
    if (unfold) {
      fold = false;
    } else {
      if (foldData[nodeID].wasSet && sameParents) {
        fold = foldDataBool;
        unfold = true;
        unfoldSet = true;
      } else
        fold = (foldDuplicatedNodes && !isIndexBiggerThanMaxID(nodeID)) ||
               foldArray[dimension];
    }

    if (dimension == 0 || children.size() == 0) fold = false;

    switch (dimension) {
      case 5:
        xNodeCoordinate = xRootNode;
        yNodeCoordinate = maxY + yBlockDistance;
        y4LastRow = 0;
        if (yNodeCoordinate > maxY) maxY = yNodeCoordinate;
        break;
      case 4:
        xNodeCoordinate =
            xRootNodeOffset + xDistance +
            xBlockDistance *
                (i %
                 max((int)(this->width() / xBlockDistance / viewScale.x()), 1));
        // yNodeCoordinate;
        // If this block is first on current row we change its yCoord
        if (i % max((int)(this->width() / xBlockDistance / viewScale.x()), 1) ==
            0) {
          yNodeCoordinate = maxY + yBlockDistance + rectYSize / 2;
          if (i != 0) yNodeCoordinate += yRootNodeOffset;
          maxY = yNodeCoordinate;
          y4LastRow = yNodeCoordinate;
          if (sett->RootNodeCreation.dimension_of_rootNodes == 5) {
            renderLine(xPreviousDimension + rectXSize, yPreviousDimension,
                       xPreviousDimension + rectXSize,
                       y4LastRow - yBlockDistance, painter);
          }
        } else {
          yNodeCoordinate = y4LastRow;
          if (maxY < yNodeCoordinate) maxY = yNodeCoordinate;
        }
        maxYCurrentBlock = yNodeCoordinate;
        if (sett->RootNodeCreation.dimension_of_rootNodes == 5) {
          renderLine(xPreviousDimension + rectXSize, y4LastRow - yBlockDistance,
                     xNodeCoordinate, y4LastRow - yBlockDistance, painter);
          renderLine(xNodeCoordinate, y4LastRow - yBlockDistance,
                     xNodeCoordinate, yNodeCoordinate, painter);
        }
        break;
      default:
        xNodeCoordinate = xPreviousDimension + xDistance + rectXSize;
        // yNodeCoordinate;
        if (i != 0)
          yNodeCoordinate = maxYCurrentBlock + yDistance + rectYSize;
        else
          yNodeCoordinate = maxYCurrentBlock;

        maxYCurrentBlock = yNodeCoordinate;
        if (maxYCurrentBlock > maxY) maxY = maxYCurrentBlock;
        renderLine(xPreviousDimension + rectXSize,
                   yPreviousDimension + rectYSize / 2,
                   xPreviousDimension + rectXSize,
                   yNodeCoordinate + rectYSize / 2, painter);
        renderLine(xPreviousDimension + rectXSize,
                   yNodeCoordinate + rectYSize / 2, xNodeCoordinate,
                   yNodeCoordinate + rectYSize / 2, painter);
        break;
    }
    if (D) cout << "X " << xNodeCoordinate << " Y " << yNodeCoordinate << endl;
    renderNodeAndCheckIfNodeWasSelected(xNodeCoordinate, yNodeCoordinate,
                                        colors[dimension], nodeID, dimension,
                                        painter, fold);
    if (dimension != 0) {
      if (!fold) {
        traverseNodes(children, dimension - 1, xNodeCoordinate, yNodeCoordinate,
                      painter);
      }
    }
    if (unfoldSet) unfold = false;
  }
}

bool AtlasViewWidget::compareParents(FoldData data, int dimension) {
  for (int i = dimension; i < 6; i++) {
    if (data.parents[i] != currentParents[i]) return false;
  }
  return true;
}
bool AtlasViewWidget::compareParents(SelectedNode data, int dimension) {
  for (int i = dimension; i < 6; i++) {
    if (data.parents[i] != currentParents[i]) return false;
  }
  return true;
}

void AtlasViewWidget::resetFoldData() {
  for (int i = 0; i < foldData.size(); i++) {
    foldData[i].wasSet = false;
    foldData[i].isFolded = false;
  }
}

bool AtlasViewWidget::isIndexBiggerThanMaxID(int index) {
  if (maxNodeID < index) {
    maxNodeID = index;
    return true;
  } else {
    return false;
  }
}

void AtlasViewWidget::updateFoldData() {
  while (foldData.size() < atlas->number_of_nodes()) {
    FoldData tmp;
    tmp.isFolded = false;
    tmp.wasSet = false;
    for (int i = 0; i < 6; i++) {
      tmp.parents[i] = -1;
    }
    foldData.push_back(tmp);
  }
}

void AtlasViewWidget::paintEvent(QPaintEvent *event) {
  // D is debug var that enables console output
  bool D = false;
  // We check if setData funtion was called and atlas and saveloader were
  // initialized
  if (dataWasSet && isWidgetActive && atlas != NULL) {
    QPainter painter(this);
    renderNodesRecursively(&painter);
  }
}

void AtlasViewWidget::renderNodeAndCheckIfNodeWasSelected(int x, int y,
                                                          Qt::GlobalColor color,
                                                          int id, int dimension,
                                                          QPainter *painter,
                                                          bool isFolded) {
  bool D = true;
  bool isHighLighted = false;
  int highlightWidth = 5;

  // Check if node is in viewport
  if (!(y > -viewTranslation.y() / viewScale.y() &&
        y < -viewTranslation.y() / viewScale.y() +
                this->height() / viewScale.y())) {
    return;
  }

  if (selectedNode.ID == id) {
    if (compareParents(selectedNode, dimension)) {
      isHighLighted = true;
      int frameWidth = 1;
      // Set pen to determine rect bordr width
      QPen pen = painter->pen();
      pen.setWidth(frameWidth);
      pen.setColor(Qt::black);
      painter->setPen(pen);
      // Set brush thsat determines rect fill color
      QBrush brush = painter->brush();
      brush.setColor(Qt::lightGray);
      brush.setStyle(Qt::SolidPattern);
      painter->setBrush(brush);
      // Render rect with current pen and brush (highlight rect)
      painter->drawRect(QRect(x - highlightWidth, y - highlightWidth,
                              rectXSize + highlightWidth * 2,
                              rectYSize + highlightWidth * 2));
      brush.setColor(color);
      brush.setStyle(Qt::SolidPattern);
      painter->setBrush(brush);
      // Render rect with current pen and brush
      painter->drawRect(QRect(x, y, rectXSize, rectYSize));
    } else {
      isHighLighted = true;
      int frameWidth = 1;
      // Set pen to determine rect bordr width
      QPen pen = painter->pen();
      pen.setWidth(frameWidth);
      pen.setColor(Qt::black);
      painter->setPen(pen);
      // Set brush thsat determines rect fill color
      QBrush brush = painter->brush();
      brush.setColor(Qt::yellow);
      brush.setStyle(Qt::SolidPattern);
      painter->setBrush(brush);
      // Render rect with current pen and brush (highlight rect)
      painter->drawRect(QRect(x - highlightWidth, y - highlightWidth,
                              rectXSize + highlightWidth * 2,
                              rectYSize + highlightWidth * 2));
      brush.setColor(color);
      brush.setStyle(Qt::SolidPattern);
      painter->setBrush(brush);
      // Render rect with current pen and brush
      painter->drawRect(QRect(x, y, rectXSize, rectYSize));
    }
  } else {
    int frameWidth = 1;
    // Set pen to determine rect bordr width
    QPen pen = painter->pen();
    pen.setWidth(frameWidth);
    pen.setColor(Qt::black);
    painter->setPen(pen);
    // Set brush thsat determines rect fill color
    QBrush brush = painter->brush();
    brush.setColor(color);
    brush.setStyle(Qt::SolidPattern);
    painter->setBrush(brush);
    // Render rect with current pen and brush
    painter->drawRect(QRect(x, y, rectXSize, rectYSize));
  }
  if (isFolded) {
    painter->drawRect(QRect(x + 3, y + 3, rectXSize, rectYSize));
    painter->drawRect(QRect(x + 6, y + 6, rectXSize, rectYSize));
  }

  // Render node ID
  painter->drawText(x, y, QString::number(id));

  if (mouseWasJustPressed) {
    mouseLastClickCoordinates = mapFromGlobal(QCursor::pos());
  }
  // Check if this node was previously clicked
  if (mouseWasJustPressed &&
      mouseLastClickCoordinates.x() - viewTranslation.x() >=
          x * viewScale.x() &&
      mouseLastClickCoordinates.x() - viewTranslation.x() <=
          (x + rectXSize) * viewScale.x() &&
      mouseLastClickCoordinates.y() - viewTranslation.y() >=
          (y - 6) * viewScale.y() &&
      mouseLastClickCoordinates.y() - viewTranslation.y() <=
          (y + rectYSize) * viewScale.y()) {
    nodeWasSelected = true;
    mouseWasJustPressed = false;
    lastSelectedNodeID = id;
    foldData[id].isFolded = !foldData[id].isFolded;
    foldData[id].wasSet = true;
    selectedNode.ID = id;
    for (int i = dimension; i < 6; i++) {
      foldData[id].parents[i] = currentParents[i];
      selectedNode.parents[i] = currentParents[i];
    }
    lastAction = LastAction::NODE_PICKED;

    painter->fillRect(QRect(x - 4, y - 4, rectXSize + 8, rectYSize + 8),
                      Qt::white);
    // if(RMBisPressed) {
    sharedData->currentAtlasNode = atlas->getNodes()[id];
    sharedData->currentAtlasNodeID = id;
    sharedData->atlasNodeWasSet = true;
    //}
    if (D) cout << "You have clicked on node " << lastSelectedNodeID << endl;
  }
}

void AtlasViewWidget::renderLine(int x1, int y1, int x2, int y2,
                                 QPainter *painter) {
  // Check if node is in viewport
  if (!(y1 > -viewTranslation.y() / viewScale.y() &&
        y1 < -viewTranslation.y() / viewScale.y() +
                 this->height() / viewScale.y()) &&
      !(y2 > -viewTranslation.y() / viewScale.y() &&
        y2 < -viewTranslation.y() / viewScale.y() +
                 this->height() / viewScale.y())) {
    return;
  }
  painter->drawLine(x1, y1, x2, y2);
}

void AtlasViewWidget::mousePressEvent(QMouseEvent *event) {
  if (event->button() == Qt::LeftButton) {
    LMBisPressed = true;
    mouseWasJustPressed = true;
  } else if (event->button() == Qt::RightButton) {
    RMBisPressed = true;
    mouseWasJustPressed = true;
  }
}

void AtlasViewWidget::mouseReleaseEvent(QMouseEvent *event) {
  LMBisPressed = false;
  RMBisPressed = false;
  mouseWasReleased = true;
}

void AtlasViewWidget::keyPressEvent(QKeyEvent *event) {
  //[ MOVE VIEWPORT KEYS]
  if (event->key() == Qt::Key_A) {
    viewTranslation.setX(viewTranslation.x() - 10);
  }
  if (event->key() == Qt::Key_D) {
    viewTranslation.setX(viewTranslation.x() + 10);
  }
  if (event->key() == Qt::Key_S) {
    viewTranslation.setY(viewTranslation.y() - 10);
  }
  if (event->key() == Qt::Key_W) {
    viewTranslation.setY(viewTranslation.y() + 10);
  }
  if (event->key() == Qt::Key_0) {
    lastAction = LastAction::FOLD_DIMENSION;
  }
  if (event->key() == Qt::Key_1) {
    foldArray[1] = !foldArray[1];
    resetFoldData();
  }
  if (event->key() == Qt::Key_2) {
    foldArray[2] = !foldArray[2];
    resetFoldData();
  }
  if (event->key() == Qt::Key_3) {
    foldArray[3] = !foldArray[3];
    resetFoldData();
  }
  if (event->key() == Qt::Key_4) {
    foldArray[4] = !foldArray[4];
    resetFoldData();
  }
  if (event->key() == Qt::Key_5) {
    foldArray[5] = !foldArray[5];
    resetFoldData();
  }
  if (event->key() == Qt::Key_F) {
    foldDuplicatedNodes = !foldDuplicatedNodes;
    lastAction = LastAction::FOLD_DUPLICATES;
    resetFoldData();
  }
  if (event->key() == Qt::Key_C) {
    viewTranslation = QVector2D(0, 0);
    viewScale = QVector2D(1, 1);
  }
  if (event->key() == Qt::Key_Plus || event->key() == Qt::Key_Equal) {
    viewScale *= 1.2;
  }
  if (event->key() == Qt::Key_Minus) {
    viewScale /= 1.2;
  }

  //[~MOVE VIEWPORT KEYS]
}

void AtlasViewWidget::wheelEvent(QWheelEvent *event) {
  // delta represtents wheel rotation in some random units (120 per 1 wheel
  // "sector"
  double delta = event->delta();
  if (delta > 0) {
    viewScale *= 1.2;
  } else {
    viewScale /= 1.2;
  }
}

std::vector<int> AtlasViewWidget::getChildren(std::vector<int> nodes,
                                              int currentDimension) {
  // Atlas node has only getConnections() method that contains both connection
  // to parent and child We filter all children looking at parent
  // dimension(currentDimension and child dimension
  std::vector<int> result;
  for (int i = 0; i < nodes.size(); i++) {
    if ((*atlas)[nodes[i]]->getDim() < currentDimension &&
        (*atlas)[nodes[i]]->hasAnyGoodOrientation()) {
      result.push_back(nodes[i]);
    }
  }
  return result;
}

std::vector<int> AtlasViewWidget::getRootIndices() {
  std::vector<long unsigned int> tmp = atlas->getRootIndicies();
  std::vector<int> rootNodeIndices;
  rootNodeIndices.reserve(tmp.size());
  for (int i = 0; i < tmp.size(); i++) {
    rootNodeIndices.push_back(tmp[i]);
  }
  return rootNodeIndices;
}

void AtlasViewWidget::updateAtlas() {
  if (dataWasSet && isWidgetActive) {
    saveLoader->loadMapView(this->atlas);
    if (atlas->getNodes().size() == 0) {
      return;
    }
    if (sharedData->atlasNodeWasSet) {
      sharedData->currentAtlasNode =
          atlas->getNodes()[sharedData->currentAtlasNodeID];
    }
  }
}
