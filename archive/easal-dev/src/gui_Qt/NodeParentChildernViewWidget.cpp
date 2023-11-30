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
#include "NodeParentChildernViewWidget.h"

NodeParentChildernViewWidget::NodeParentChildernViewWidget(
    SharedDataGUI *sharedData) {
  this->sharedData = sharedData;
  //[ HANDEL TIMER]
  updateTimeInMS = 1000.0 / 30.0;
  timer = new QTimer(this);
  // start argument determins update rate
  timer->start(updateTimeInMS);
  // repaint is function of QWidget that triggers paintEvent
  connect(timer, SIGNAL(timeout()), this, SLOT(repaint()));
  //[~HANDEL TIMER]
  //[ HANDEL FOCUS]
  // without focus policy mouse and keyboard events cannot work
  this->setFocusPolicy(Qt::StrongFocus);
  //[~HANDLE FOCUS]
  //[ HANDLE VIEWPORT]
  viewTranslation = QVector2D(0, 0);
  viewScale = QVector2D(1, 1);
  //[~HANDLE VIEWPORT]
  //[ HANDLE MOUSE]
  LMBisPressed = false;
  RMBisPressed = false;
  mouseWasJustPressed = false;
  mouseWasReleased = true;
  //[~HANDLE MOUSE]
  //[ HANDLE WIDGET]
  widgetIsActive = false;
  nodeWasSet = false;
  //[ HANDLE WIDGET]
  //[ HANDLE SIZES]
  xSize = 10;
  ySize = 5;
  //    fontSize = ySize-fontOffset*2;
  //[~HANDLE SIZES]
  //[ HANDLE COLORS]
  colors[5] = Qt::yellow;
  colors[4] = Qt::red;
  colors[3] = Qt::green;
  colors[2] = Qt::cyan;
  colors[1] = Qt::darkCyan;
  colors[0] = Qt::darkGray;
  //[~HANDLE COLORS]

  //[ HANDLE MISC]
  currentAtlasNodeID = -1;
  //[~HANDLE MISC]
}

void NodeParentChildernViewWidget::paintEvent(QPaintEvent *event) {
  if (mouseWasJustPressed) {
    mouseLastClickCoordinates = mapFromGlobal(QCursor::pos());
  }
  if (sharedData->atlasNodeWasSet) {
    if (currentAtlasNodeID != sharedData->currentAtlasNode->getID()) {
      setAtlasNode(sharedData->currentAtlasNode);
    }
    QPainter painter(this);
    // Set pen that determines line width and rect boundary
    QPen pen = painter.pen();
    pen.setWidth(1);
    painter.setPen(pen);

    vector<int> tmp;
    tmp.push_back(sharedData->currentAtlasNode->getID());
    clearDimensionData();
    updateDimensionData(tmp, ROOT);
    renderDimensions(&painter);
    updateViewTranslationAndScale(
        &painter);  ///< must be called after renderDimensions(&painter) because
                    ///< of LMBisPressed
  }
  mouseWasJustPressed = false;
}

void NodeParentChildernViewWidget::setAtlasNode(AtlasNode *atlasNode) {
  if (atlasNode != 0) {
    this->currentAtlasNode = atlasNode;
    currentAtlasNodeID = atlasNode->getID();
    nodeWasSet = true;
  }
}

void NodeParentChildernViewWidget::setComponents(Atlas *atlas,
                                                 SaveLoader *saveLoader) {
  this->atlas = atlas;
  this->saveLoader = saveLoader;
}

/*
void NodeParentChildernViewWidget::updatePositions()
{
    if(sharedData->atlasNodeWasSet)
    {
        sharedData->currentAtlasNode;
    }
}
*/

void NodeParentChildernViewWidget::mousePressEvent(QMouseEvent *event) {
  if (event->button() == Qt::LeftButton) {
    LMBisPressed = true;
    mouseWasJustPressed = true;
  } else if (event->button() == Qt::RightButton) {
    RMBisPressed = true;
    mouseWasJustPressed = true;
  }
}

void NodeParentChildernViewWidget::mouseReleaseEvent(QMouseEvent *event) {
  LMBisPressed = false;
  RMBisPressed = false;
  mouseWasReleased = true;
}

void NodeParentChildernViewWidget::keyPressEvent(QKeyEvent *event) {
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

void NodeParentChildernViewWidget::renderRootDimension(QPainter *painter) {
  // Render root node
  int x = width() / 12 + (5 - rootDimension) * width() / 6;
  int y = height() / 2;
  if (rootDimension != 0) {
    for (int k = 0; k < dimensionData[rootDimension].nodesChildren[0].size();
         k++) {
      int childID = dimensionData[rootDimension].nodesChildren[0][k];
      int childIndex = dimensionData[rootDimension - 1].mapIDToIndex[childID];
      // Find child coords
      int childX = x + width() / 6;
      int childY =
          y -
          dimensionData[rootDimension - 1].verticalDistanceBetweenBoxes *
              dimensionData[rootDimension - 1].nodesIDs.size() / 2 +
          childIndex *
              dimensionData[rootDimension - 1].verticalDistanceBetweenBoxes;
      painter->drawLine(x + dimensionData[rootDimension].boxSize.width(),
                        y + dimensionData[rootDimension].boxSize.height() / 2,
                        childX, childY);
    }
  }

  if (rootDimension != 5) {
    for (int k = 0; k < dimensionData[rootDimension].nodesParents[0].size();
         k++) {
      int parentID = dimensionData[rootDimension].nodesParents[0][k];
      int parentIndex = dimensionData[rootDimension + 1].mapIDToIndex[parentID];
      // Find child coords
      int parentX =
          x - width() / 6 + dimensionData[rootDimension + 1].boxSize.width();
      int parentY =
          y -
          dimensionData[rootDimension + 1].verticalDistanceBetweenBoxes *
              dimensionData[rootDimension + 1].nodesIDs.size() / 2 +
          parentIndex *
              dimensionData[rootDimension + 1].verticalDistanceBetweenBoxes;
      painter->drawLine(x,
                        y + dimensionData[rootDimension].boxSize.height() / 2,
                        parentX, parentY);
    }
  }
  int highlightSize = 5;
  painter->fillRect(
      QRect(x - highlightSize, y - highlightSize,
            dimensionData[rootDimension].boxSize.width() + highlightSize * 2,
            dimensionData[rootDimension].boxSize.height() + highlightSize * 2),
      Qt::lightGray);  ///< highlight rect
  renderBox(x, y, dimensionData[rootDimension].boxSize.width(),
            dimensionData[rootDimension].boxSize.height(),
            dimensionData[rootDimension].nodesIDs[0], colors[rootDimension],
            painter);
}

void NodeParentChildernViewWidget::renderChildrenDimensions(QPainter *painter) {
  int x = width() / 12 + (5 - rootDimension) * width() / 6;
  int y = height() / 2;
  // Render children
  int childrenHorizontalOffset = x;
  int childrenVerticalOffset = 0;
  for (int i = rootDimension - 1; i >= 0; i--) {
    childrenHorizontalOffset += width() / 6;
    //        childrenHorizontalOffset += dimensionData[i].distanceToChild;
    childrenVerticalOffset = y - dimensionData[i].verticalDistanceBetweenBoxes *
                                     dimensionData[i].nodesIDs.size() / 2;

    for (int j = 0; j < dimensionData[i].nodesIDs.size(); j++) {
      renderBox(childrenHorizontalOffset,
                childrenVerticalOffset +
                    j * dimensionData[i].verticalDistanceBetweenBoxes,
                dimensionData[i].boxSize.width(),
                dimensionData[i].boxSize.height(), dimensionData[i].nodesIDs[j],
                colors[i], painter);

      if (i != 0) {
        for (int k = 0; k < dimensionData[i].nodesChildren[j].size(); k++) {
          int childID = dimensionData[i].nodesChildren[j][k];
          int childIndex = dimensionData[i - 1].mapIDToIndex[childID];
          // Find child coords
          int childVerticalOffset =
              y - dimensionData[i - 1].verticalDistanceBetweenBoxes *
                      dimensionData[i - 1].nodesIDs.size() / 2;
          int childX = childrenHorizontalOffset + width() / 6;
          int childY =
              childVerticalOffset +
              childIndex * dimensionData[i - 1].verticalDistanceBetweenBoxes;
          painter->drawLine(
              childrenHorizontalOffset + dimensionData[i].boxSize.width(),
              childrenVerticalOffset +
                  j * dimensionData[i].verticalDistanceBetweenBoxes,
              childX, childY);
        }
      }
    }
  }
}

void NodeParentChildernViewWidget::renderParentsDimensions(QPainter *painter) {
  int x = width() / 12 + (5 - rootDimension) * width() / 6;
  int y = height() / 2;
  // Render parents
  int parentsHorizontalOffset = x;
  int parentsVerticalOffset = 0;
  for (int i = rootDimension + 1; i < 6; i++) {
    parentsHorizontalOffset -= width() / 6;
    //        parentHorizontalOffset -= dimensionData[i].distanceToParent;
    parentsVerticalOffset = y - dimensionData[i].verticalDistanceBetweenBoxes *
                                    dimensionData[i].nodesIDs.size() / 2;
    for (int j = 0; j < dimensionData[i].nodesIDs.size(); j++) {
      renderBox(parentsHorizontalOffset,
                parentsVerticalOffset +
                    j * dimensionData[i].verticalDistanceBetweenBoxes,
                dimensionData[i].boxSize.width(),
                dimensionData[i].boxSize.height(), dimensionData[i].nodesIDs[j],
                colors[i], painter);
      if (i != 5) {
        for (int k = 0; k < dimensionData[i].nodesParents[j].size(); k++) {
          int parentID = dimensionData[i].nodesParents[j][k];
          int parentIndex = dimensionData[i + 1].mapIDToIndex[parentID];
          // Find child coords
          int parentVerticalOffset =
              y - dimensionData[i + 1].verticalDistanceBetweenBoxes *
                      dimensionData[i + 1].nodesIDs.size() / 2;
          int parentX = parentsHorizontalOffset - width() / 6 +
                        dimensionData[i + 1].boxSize.width();
          int parentY =
              parentVerticalOffset +
              parentIndex * dimensionData[i + 1].verticalDistanceBetweenBoxes;
          painter->drawLine(
              parentsHorizontalOffset,
              parentsVerticalOffset +
                  j * dimensionData[i].verticalDistanceBetweenBoxes,
              parentX, parentY);
        }
      }
    }
  }
}

void NodeParentChildernViewWidget::renderBox(int x, int y, int width,
                                             int height, int ID,
                                             Qt::GlobalColor color,
                                             QPainter *painter) {
  QBrush brush = painter->brush();
  brush.setColor(color);
  brush.setStyle(Qt::SolidPattern);
  painter->setBrush(brush);
  QRect tmp = QRect(x, y, width, height);
  painter->drawRect(tmp);
  renderIDString(x, y, ID, painter);
  mouseNodePickingCheck(x, y, ID);
}

void NodeParentChildernViewWidget::renderIDString(int x, int y, int ID,
                                                  QPainter *painter) {
  painter->drawText(x, y, QString::number(ID));
}

void NodeParentChildernViewWidget::renderDimensions(QPainter *painter) {
  renderRootDimension(painter);
  renderChildrenDimensions(painter);
  renderParentsDimensions(painter);
}

void NodeParentChildernViewWidget::mouseNodePickingCheck(int x, int y, int ID) {
  if (mouseLastClickCoordinates.x() - viewTranslation.x() >=
          x * viewScale.x() &&
      mouseLastClickCoordinates.x() - viewTranslation.x() <=
          (x + xSize) * viewScale.x() &&
      mouseLastClickCoordinates.y() - viewTranslation.y() >=
          (y - 6) * viewScale.y() &&
      mouseLastClickCoordinates.y() - viewTranslation.y() <=
          (y + ySize) * viewScale.y() &&
      LMBisPressed) {
    sharedData->currentAtlasNode = (*atlas)[ID];
    viewTranslation = QVector2D(0, 0);
  }
}

void NodeParentChildernViewWidget::clearDimensionData() {
  for (int i = 0; i < 6; i++) {
    for (int j = 0; j < dimensionData[i].nodesChildren.size(); j++) {
      dimensionData[i].nodesChildren[j].clear();
    }
    for (int j = 0; j < dimensionData[i].nodesParents.size(); j++) {
      dimensionData[i].nodesParents[j].clear();
    }
    dimensionData[i].nodesChildren.clear();
    dimensionData[i].nodesParents.clear();
    dimensionData[i].nodesIDs.clear();
    dimensionData[i].mapIDToIndex.clear();
  }
}

vector<int> NodeParentChildernViewWidget::getChildernIDs(AtlasNode *node,
                                                         int currentDim) {
  vector<int> IDs = node->getConnection();
  vector<int> children;
  for (int i = 0; i < IDs.size(); i++) {
    if ((*atlas)[IDs[i]]->getDim() < currentDim) {
      children.push_back(IDs[i]);
    }
  }
  return children;
}

vector<int> NodeParentChildernViewWidget::getParentsIDs(AtlasNode *node,
                                                        int currentDim) {
  vector<int> IDs = node->getConnection();
  vector<int> parents;
  for (int i = 0; i < IDs.size(); i++) {
    if ((*atlas)[IDs[i]]->getDim() > currentDim) {
      parents.push_back(IDs[i]);
    }
  }
  return parents;
}

void NodeParentChildernViewWidget::updateViewTranslationAndScale(
    QPainter *painter) {
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

void NodeParentChildernViewWidget::updateDimensionData(
    vector<int> nodesIDs, RelationToRoot relationToRoot) {
  if (nodesIDs.size() == 0) {
    return;
  }
  int dimension = (*atlas)[nodesIDs[0]]->getDim();
  dimensionData[dimension].nodesIDs = nodesIDs;
  dimensionData[dimension].relationToRoot = ROOT;
  dimensionData[dimension].boxSize = QSize(xSize, ySize);
  dimensionData[dimension].distanceToChild = 100;
  dimensionData[dimension].distanceToParent = 100;
  dimensionData[dimension].verticalDistanceBetweenBoxes = 30;
  dimensionData[dimension].verticalOffset = 0;

  if (relationToRoot == ROOT) {
    rootDimension = dimension;
  }

  if (relationToRoot == ROOT || relationToRoot == CHILD) {
    // Get all children
    vector<int> children;
    for (int i = 0; i < nodesIDs.size(); i++) {
      vector<int> tmpChildren =
          getChildernIDs((*atlas)[nodesIDs[i]], dimension);
      // Save children
      dimensionData[dimension].nodesChildren.push_back(tmpChildren);
      children.insert(children.end(), tmpChildren.begin(), tmpChildren.end());
    }
    // Sort and get rid of duplicates
    sort(children.begin(), children.end());
    children.erase(unique(children.begin(), children.end()), children.end());

    updateDimensionData(children, CHILD);
  }

  if (relationToRoot == ROOT || relationToRoot == PARENT) {
    // Get all children
    vector<int> parents;
    for (int i = 0; i < nodesIDs.size(); i++) {
      vector<int> tmpParents = getParentsIDs((*atlas)[nodesIDs[i]], dimension);
      // Save parents
      dimensionData[dimension].nodesParents.push_back(tmpParents);
      parents.insert(parents.end(), tmpParents.begin(), tmpParents.end());
    }
    // Sort and get rid of duplicates
    sort(parents.begin(), parents.end());
    parents.erase(unique(parents.begin(), parents.end()), parents.end());

    updateDimensionData(parents, PARENT);
  }

  dimensionData[dimension].mapIDToIndex.clear();

  for (int i = 0; i < nodesIDs.size(); i++) {
    dimensionData[dimension].mapIDToIndex.insert(
        {dimensionData[dimension].nodesIDs[i], i});
  }
}
