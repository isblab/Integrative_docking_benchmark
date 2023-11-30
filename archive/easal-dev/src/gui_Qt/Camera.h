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
#ifndef CAMERA_H
#define CAMERA_H

#include <QKeyEvent>
#include <QMatrix4x4>
#include <QMouseEvent>
#include <QPoint>
#include <QResizeEvent>
#include <QSize>
#include <QVector3D>
#include <QWheelEvent>
#include <iostream>

#define PI 3.14159265

using namespace std;

class Camera {
 public:
  Camera();
  //[ INPUT]
  void mousePressEvent(QMouseEvent *mouseEvent);
  void mouseReleaseEvent(QMouseEvent *mouseEvent);
  void wheelEvent(QWheelEvent *wheelEvent);
  void keyPressEvent(QKeyEvent *keyEvent);
  void resizeEvent(QResizeEvent *resizeEvent);
  void resizeEvent(int width, int height);
  //[~INPUT]
  //[ GETTERS SETTERS]
  /**
   * @brief getLastCameraMatrix: return camera matrix with last used mouse
   * coordinates
   * @param mouseCoordinates
   * @return
   */
  QMatrix4x4 getLastCameraMatrix();
  QMatrix4x4 getCameraMatrix(QPoint mouseCoordinates);
  /**
   * @brief getEyeLocation: allows to get current position from which camera
   * "looks" from
   * @return
   */
  QVector3D getEyeLocation();
  /**
   * @brief getRight: allows to get vector that points right from eye
   * @return
   */
  QVector3D getRight();
  /**
   * @brief getUp: allows to get vector that points up from eye
   * @return
   */
  QVector3D getUp();
  void setPosition(QVector3D position);
  /**
   * @brief init: sets initail data, works as constructor
   * @param screenSize
   * @param position: posioton of camera target
   * @param direction: direction where camera looks
   */
  void init(QSize screenSize, QVector3D position, QVector3D direction);
  //[~GETTERS SETTERS]

 private:
  //[ CAMERA SETTINGS]
  double rotationSpeed = 360;
  double moveSpeed = 0.1;
  double fov = 60;
  double nearPlaneDistance = -300;
  double farPlaneDistance = 300;
  double distanceToCamera = -10;
  QSize screenSize;
  QMatrix4x4 cameraMatrix;
  QMatrix4x4 viewMatrix;
  QMatrix4x4 projectionMatrix;
  QVector3D position = QVector3D(0, 0, 0);
  QVector3D direction = QVector3D(0, 0, -1);
  //[~CAMERA SETTINGS]

  double wheelSensitivity = 1.0 / 120.0;

  bool mouseIsPressed = false;
  bool mouseWasJustPressed = false;

  //[ HELPER VECTORS]
  QVector3D up = QVector3D(0, 1, 0);
  QVector3D right = QVector3D(1, 0, 0);
  //[~HELPER VECTORS]
  QPoint lastMouseCoordinates;

  void updateViewMatrix(QPoint mouseCoordinates);
};

#endif  // CAMERA_H
