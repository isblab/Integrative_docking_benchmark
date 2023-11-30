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
#ifndef RENDERER_H
#define RENDERER_H

#include <QMatrix4x4>
#include <QOpenGLBuffer>
#include <QOpenGLContext>
#include <QOpenGLExtraFunctions>
#include <QOpenGLFramebufferObject>
#include <QOpenGLFunctions>
#include <QOpenGLShaderProgram>
#include <QOpenGLVertexArrayObject>
#include <QPoint>
#include <QVector3D>
#include <cmath>
#include <iostream>
#include <vector>

#include "Mesh3D.h"

using namespace std;

class Renderer {
 public:
  Renderer();
  /**
   * @brief initOpenGL: prepares class for work, must be called before first
   * render call
   */
  void initOpenGL();
  /**
   * @brief bindRenderProgramAndClear: binds Renderer program and clears color
   * and depth buffer
   */
  void bindRenderProgramAndClear();
  /**
   * @brief bindPickingProgramAndClear: binds Picking program and clears color
   * and depth buffer
   */
  void bindPickingProgramAndClear();
  /**
   * @brief render: renders object made out of vertices in mesh using 3
   * transformation matrices and color(optional)
   * @param mesh: class that contains vertex data (coordinates of vertex and
   * normals) of object that we render
   * @param translationMatrix: vector of QMatrix4x4 that determines object
   * position
   * @param rotationMatrix: vector of QMatrix4x4 that determines object rotation
   * @param scaleMatrix: vector of QMatrix4x4 that determines object size
   * @param colorVector: vector of QVector4D that determines RGBA color of
   * object
   * @param colorIDOffset: int that is used for mousepicking reasons (it allows
   * to destinguish objects rendered at different render calls
   */
  void render(Mesh3D *mesh, vector<QMatrix4x4> translationMatrix,
              vector<QMatrix4x4> rotationMatrix, vector<QMatrix4x4> scaleMatrix,
              vector<QVector4D> colorVector = vector<QVector4D>(),
              int colorIDOffset = 0);
  /**
   * @brief enableDepthTest: enables depth testing to render everything properly
   */
  void enableDepthTest();
  /**
   * @brief disableDepthTest: disables depth testing to render stuff on top of
   * everything
   */
  void disableDepthTest();
  /**
   * @brief enableBlend: allow to take alpha color in consideration to deal with
   * transparency
   */
  void enableBlend();
  /**
   * @brief disableBlend: turns off transparency ignoring alpha channel
   */
  void disableBlend();
  /**
   * @brief getPixelValueAt: allows to access pixel at particular position and
   * get ID of object there or -1
   * @param mousePosition: QPoint that holds mouse coords
   * @param screenSize: QSize that hold width and height
   * @return int that represent object ID that mouse was clicked on(-1 if no
   * object)
   */
  void disableAll();
  void enableAll();
  int getPixelValueAt(QPoint mousePosition, QSize screenSize);
  /**
   * @brief renderingProgram: simple program with flat shading
   */
  QOpenGLShaderProgram *renderingProgram;
  /**
   * @brief pickingProgram: program that is used only for mouse picking
   */
  QOpenGLShaderProgram *pickingProgram;
  //[ UNIFORMS]
  /**
   * @brief cameraMatrix: QMatrix4x4 that determines where camera is located and
   * how it works
   */
  QMatrix4x4 cameraMatrix;
  /**
   * @brief lightPositionVector: QVector3D that hold coordinate of single light
   * source that shading program supports
   */
  QVector3D lightPositionVector;
  //[~UNIFORMS]

  void renderWithTransperansy(Mesh3D *mesh,
                              vector<QMatrix4x4> translationMatrix,
                              vector<QMatrix4x4> rotationMatrix,
                              vector<QMatrix4x4> scaleMatrix,
                              vector<QVector4D> colorVector);

 private:
  /**
   * @brief instancesPerDrawCall: int that deternies max number of same objects
   * that can be rendered per single render call. This value is hardcoded in
   * shader
   */
  double instancesPerDrawCall = 128;
  /**
   * @brief versionedShaderCode: chooses which version of OpenGL to use
   * @param src: pointer to shader program text
   * @return appended text of src
   */
  QByteArray versionedShaderCode(const char *src);
  bool renderingProgramIsBind = false;
  bool pickingProgramIsBind = false;

  //[ UNIFORM LOCATIONS]
  int translationMatrixLocation;
  int rotationMatrixLocation;
  int scaleMatrixLocation;
  int cameraMatrixLocation;
  int lightPositionVectorLocation;
  int colorVectorLocation;

  int translationMatrixLocationPicking;
  int rotationMatrixLocationPicking;
  int scaleMatrixLocationPicking;
  int colorIDLocationPicking;
  //[ UNIFORM LOCATIONS]
};

#endif  // RENDERER_H
