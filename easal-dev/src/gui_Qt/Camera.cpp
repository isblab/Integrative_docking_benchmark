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
#include "Camera.h"

Camera::Camera() {}

void Camera::mousePressEvent(QMouseEvent *mouseEvent) {
  if (mouseEvent->button() == Qt::LeftButton) {
    mouseIsPressed = true;
    mouseWasJustPressed = true;
  }
}

void Camera::mouseReleaseEvent(QMouseEvent *mouseEvent) {
  if (mouseEvent->button() == Qt::LeftButton) {
    mouseIsPressed = false;
  }
}

void Camera::wheelEvent(QWheelEvent *wheelEvent) {
  distanceToCamera += double(wheelEvent->delta()) * wheelSensitivity;
  updateViewMatrix(lastMouseCoordinates);
  //    cout<<"Delta "<<wheelEvent->delta()<<" distanceToCAmera
  //    "<<distanceToCamera<<endl;
}

void Camera::keyPressEvent(QKeyEvent *keyEvent) {
  if (keyEvent->key() == Qt::Key_Q) {
    position += direction * moveSpeed;
  }
  if (keyEvent->key() == Qt::Key_E) {
    position -= direction * moveSpeed;
  }
  if (keyEvent->key() == Qt::Key_W) {
    position += up * moveSpeed;
  }
  if (keyEvent->key() == Qt::Key_S) {
    position -= up * moveSpeed;
  }
  if (keyEvent->key() == Qt::Key_A) {
    position -= right * moveSpeed;
  }
  if (keyEvent->key() == Qt::Key_D) {
    position += right * moveSpeed;
  }
  if (keyEvent->key() == Qt::Key_J) {
    direction = QVector3D(0, 1, 0);
    right = QVector3D(1, 0, 0);
    up = QVector3D(0, 0, -1);
    updateViewMatrix(lastMouseCoordinates);
  }
  if (keyEvent->key() == Qt::Key_K) {
    direction = QVector3D(0, 0, -1);
    right = QVector3D(1, 0, 0);
    up = QVector3D(0, -1, 0);
    updateViewMatrix(lastMouseCoordinates);
  }
  if (keyEvent->key() == Qt::Key_I) {
    direction = QVector3D(1, 0, 0);
    right = QVector3D(0, 0, 1);
    up = QVector3D(0, -1, 0);
    updateViewMatrix(lastMouseCoordinates);
  }
  updateViewMatrix(lastMouseCoordinates);
}

void Camera::resizeEvent(QResizeEvent *resizeEvent) {
  screenSize = resizeEvent->size();
  updateViewMatrix(lastMouseCoordinates);
}

void Camera::resizeEvent(int width, int height) {
  screenSize = QSize(width, height);
  updateViewMatrix(lastMouseCoordinates);
}

QMatrix4x4 Camera::getLastCameraMatrix() { return cameraMatrix; }

QMatrix4x4 Camera::getCameraMatrix(QPoint mouseCoordinates) {
  if (mouseWasJustPressed) {
    lastMouseCoordinates = mouseCoordinates;
    mouseWasJustPressed = false;
  }
  if (lastMouseCoordinates != mouseCoordinates && mouseIsPressed) {
    updateViewMatrix(mouseCoordinates);
  }
  return cameraMatrix;
}

QVector3D Camera::getEyeLocation() {
  return position + direction * distanceToCamera;
}

QVector3D Camera::getRight() { return right; }

QVector3D Camera::getUp() { return up; }

void Camera::setPosition(QVector3D position) {
  this->position = position;
  updateViewMatrix(lastMouseCoordinates);
}

void Camera::init(QSize screenSize, QVector3D position, QVector3D direction) {
  this->screenSize = screenSize;
  this->position = position;
  this->direction = direction;
}

void Camera::updateViewMatrix(QPoint mouseCoordinates) {
  //[ GET ANGLE]
  QPoint mouseDelta = mouseCoordinates - lastMouseCoordinates;
  QVector2D rotationAngle =
      QVector2D((double)mouseDelta.x() / screenSize.width() * rotationSpeed,
                (double)mouseDelta.y() / screenSize.height() * rotationSpeed);
  //    cout<<"Rotational angle "<<rotationAngle.x()<<"
  //    "<<rotationAngle.y()<<endl;
  //[~GET ANGLE]
  //[ ROTATE DIRECTION]
  double anglePitch = rotationAngle.y();  // around X
  double angleYaw = rotationAngle.x();    // around Y
  double angleRoll = 0;                   // around Z

  QMatrix4x4 directionRotationMatrix;
  directionRotationMatrix.rotate(angleYaw, up);
  right = (right * directionRotationMatrix).normalized();
  directionRotationMatrix.rotate(anglePitch, right);
  direction = (direction * directionRotationMatrix).normalized();
  up = (QVector3D::crossProduct(right, direction)).normalized();
  //    cout<<"Direction "<<direction.x()<<" "<<direction.y()<<"
  //    "<<direction.z()<<endl; cout<<"Postition "<<position.x()<<"
  //    "<<position.y()<<" "<<position.z()<<endl;
  //[~ROTATE DIRECTION]
  //[ SET MATRICES]
  viewMatrix.setToIdentity();
  viewMatrix.lookAt(position + direction * distanceToCamera, position, up);
  projectionMatrix.setToIdentity();
  double aspectRatio = (float)screenSize.width() / (float)screenSize.height();
  if (aspectRatio > 1) {
    projectionMatrix.ortho(-1.0f, 1.0f, -1.0f / aspectRatio, 1.0f / aspectRatio,
                           -10, 10000.0f);
  } else {
    projectionMatrix.ortho(-1.0f * aspectRatio, 1.0f * aspectRatio, -1.0f, 1.0f,
                           -10, 10000.0f);
  }

  if (distanceToCamera > 0) {
    projectionMatrix.scale(distanceToCamera);
  } else if (distanceToCamera < 0) {
    projectionMatrix.scale(1 / (-distanceToCamera));
  }

  //    projectionMatrix.perspective(fov,
  //    screenSize.width()/screenSize.height(), nearPlaneDistance,
  //    farPlaneDistance);
  cameraMatrix = projectionMatrix * viewMatrix;
  //[~SET MATRICES]

  lastMouseCoordinates = mouseCoordinates;
}
