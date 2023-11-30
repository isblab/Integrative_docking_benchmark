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
#include "Renderer.h"

Renderer::Renderer() {}

static const char *vertexShaderSourceRendering =
    "layout(location = 0) in highp vec4 vertex;\n"
    "layout(location = 1) in highp vec3 normal;\n"
    "uniform highp mat4 translationMatrix[32];\n"
    "uniform highp mat4 rotationMatrix[32];\n"
    "uniform highp mat4 scaleMatrix[32];\n"
    "uniform highp mat4 cameraMatrix;\n"
    "out highp vec3 vertexNormal;\n"
    "out highp vec3 vertexCoordinate;\n"
    "flat out int instanceID;\n"
    "void main() {\n"
    "   instanceID = gl_InstanceID;\n"
    "   mat4 modelMatrix = translationMatrix[gl_InstanceID] * "
    "rotationMatrix[gl_InstanceID] * scaleMatrix[gl_InstanceID];\n"
    "   vertexNormal = mat3(transpose(inverse(modelMatrix))) * normal;\n"
    "   vertexCoordinate = vec3(modelMatrix * vertex);\n"
    "   gl_Position = cameraMatrix * modelMatrix * vertex;\n"
    "}\n";

static const char *fragmentShaderSourceRendering =
    "in highp vec3 vertexNormal;\n"
    "in highp vec3 vertexCoordinate;\n"
    "flat in int instanceID;\n"
    "out highp vec4 fragColor;\n"
    "uniform highp vec3 lightPositionVector;\n"
    "uniform highp vec4 colorVector[32];\n"
    "void main() {\n"
    "   highp vec3 lightVector = normalize(lightPositionVector - "
    "vertexCoordinate);\n"
    "   highp float lightFactor = dot(normalize(vertexNormal), lightVector);\n"
    "   fragColor = clamp(colorVector[instanceID] * 0.6 + "
    "colorVector[instanceID] * 0.25 * lightFactor, 0.0, 1.0);\n"
    "}\n";

static const char *vertexShaderSourcePicking =
    "layout(location = 0) in highp vec4 vertex;\n"
    "layout(location = 1) in highp vec3 normal;\n"
    "uniform highp mat4 translationMatrix[32];\n"
    "uniform highp mat4 rotationMatrix[32];\n"
    "uniform highp mat4 scaleMatrix[32];\n"
    "uniform highp mat4 cameraMatrix;\n"
    "flat out int instanceID;\n"
    "void main() {\n"
    "   instanceID = gl_InstanceID;\n"
    "   mat4 modelMatrix = translationMatrix[gl_InstanceID] * "
    "rotationMatrix[gl_InstanceID] * scaleMatrix[gl_InstanceID];\n"
    "   gl_Position = cameraMatrix * modelMatrix * vertex;\n"
    "}\n";

static const char *fragmentShaderSourcePicking =
    "in highp vec3 vertexNormal;\n"
    "in highp vec3 vertexCoordinate;\n"
    "flat in int instanceID;\n"
    "uniform int colorID;\n"
    "out highp vec4 fragColor;\n"
    "void main() {\n"
    "   fragColor = vec4(instanceID/32.0f,colorID/32.0f,0.345,1);\n"
    "}\n";

void Renderer::initOpenGL() {
  QOpenGLFunctions *openGLFunctions =
      QOpenGLContext::currentContext()->functions();
  //[ SET SHADING PROGRAM]
  renderingProgram = new QOpenGLShaderProgram;
  renderingProgram->addShaderFromSourceCode(
      QOpenGLShader::Vertex, versionedShaderCode(vertexShaderSourceRendering));
  renderingProgram->addShaderFromSourceCode(
      QOpenGLShader::Fragment,
      versionedShaderCode(fragmentShaderSourceRendering));
  renderingProgram->link();
  //[ SET SHADING PROGRAM]

  //[ SET MOUSEPICKING PROGRAM]
  pickingProgram = new QOpenGLShaderProgram;
  pickingProgram->addShaderFromSourceCode(
      QOpenGLShader::Vertex, versionedShaderCode(vertexShaderSourcePicking));
  pickingProgram->addShaderFromSourceCode(
      QOpenGLShader::Fragment,
      versionedShaderCode(fragmentShaderSourcePicking));
  pickingProgram->link();
  //[ SET MOUSEPICKING PROGRAM]

  //[ GET UNIFORM LOCATIONS]
  translationMatrixLocation =
      renderingProgram->uniformLocation("translationMatrix");
  rotationMatrixLocation = renderingProgram->uniformLocation("rotationMatrix");
  scaleMatrixLocation = renderingProgram->uniformLocation("scaleMatrix");
  cameraMatrixLocation = renderingProgram->uniformLocation("cameraMatrix");
  lightPositionVectorLocation =
      renderingProgram->uniformLocation("lightPositionVector");
  colorVectorLocation = renderingProgram->uniformLocation("colorVector");

  translationMatrixLocationPicking =
      pickingProgram->uniformLocation("translationMatrix");
  rotationMatrixLocationPicking =
      pickingProgram->uniformLocation("rotationMatrix");
  scaleMatrixLocationPicking = pickingProgram->uniformLocation("scaleMatrix");
  colorIDLocationPicking = pickingProgram->uniformLocation("colorID");
  //[~GET UNIFORM LOCATIONS]

  //[ OPENGL SETTINGS]
  openGLFunctions->glEnable(GL_DEPTH_TEST);
  openGLFunctions->glEnable(GL_CULL_FACE);
  openGLFunctions->glClearColor(1, 1, 1, 1);

  //[~OPENGL SETTINGS]
}

void Renderer::bindRenderProgramAndClear() {
  QOpenGLExtraFunctions *openGLFunctions =
      QOpenGLContext::currentContext()->extraFunctions();
  openGLFunctions->glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  renderingProgram->bind();
  renderingProgramIsBind = true;
  pickingProgramIsBind = false;
}

void Renderer::bindPickingProgramAndClear() {
  QOpenGLExtraFunctions *openGLFunctions =
      QOpenGLContext::currentContext()->extraFunctions();
  openGLFunctions->glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  pickingProgram->bind();
  pickingProgramIsBind = true;
  renderingProgramIsBind = false;
}

void Renderer::render(Mesh3D *mesh, vector<QMatrix4x4> translationMatrix,
                      vector<QMatrix4x4> rotationMatrix,
                      vector<QMatrix4x4> scaleMatrix,
                      vector<QVector4D> colorVector, int colorIDOffset) {
  QOpenGLExtraFunctions *openGLFunctions =
      QOpenGLContext::currentContext()->extraFunctions();
  mesh->vao.bind();
  mesh->vbo.bind();

  int numberOfObjectsToRender = translationMatrix.size();
  for (int i = 0; i < numberOfObjectsToRender; i += 32) {
    int numberOfObjectsToRenderNow = numberOfObjectsToRender - i;
    if (numberOfObjectsToRenderNow > 32) {
      numberOfObjectsToRenderNow = 32;
    }
    if (renderingProgramIsBind) {
      renderingProgram->setUniformValueArray(translationMatrixLocation,
                                             &translationMatrix[i],
                                             numberOfObjectsToRenderNow);
      renderingProgram->setUniformValueArray(rotationMatrixLocation,
                                             &rotationMatrix[i],
                                             numberOfObjectsToRenderNow);
      renderingProgram->setUniformValueArray(
          scaleMatrixLocation, &scaleMatrix[i], numberOfObjectsToRenderNow);
      renderingProgram->setUniformValueArray(
          colorVectorLocation, &colorVector[i], numberOfObjectsToRenderNow);
    } else if (pickingProgramIsBind) {
      float colorID = i;

      pickingProgram->setUniformValueArray(translationMatrixLocationPicking,
                                           &translationMatrix[i],
                                           numberOfObjectsToRenderNow);
      pickingProgram->setUniformValueArray(rotationMatrixLocationPicking,
                                           &rotationMatrix[i],
                                           numberOfObjectsToRenderNow);
      pickingProgram->setUniformValueArray(scaleMatrixLocationPicking,
                                           &scaleMatrix[i],
                                           numberOfObjectsToRenderNow);
      pickingProgram->setUniformValue(colorIDLocationPicking,
                                      i / 32 + colorIDOffset);
    }

    openGLFunctions->glDrawArraysInstanced(
        GL_TRIANGLES, 0, mesh->polygonCount(), numberOfObjectsToRenderNow);
  }

  mesh->vbo.release();
  mesh->vao.release();
}

void Renderer::renderWithTransperansy(Mesh3D *mesh,
                                      vector<QMatrix4x4> translationMatrix,
                                      vector<QMatrix4x4> rotationMatrix,
                                      vector<QMatrix4x4> scaleMatrix,
                                      vector<QVector4D> colorVector) {
  vector<QMatrix4x4> translationMatrixTransparent;
  vector<QMatrix4x4> rotationMatrixTransparent;
  vector<QMatrix4x4> scaleMatrixTransparent;
  vector<QVector4D> colorVectorTransparent;

  vector<QMatrix4x4> translationMatrixOpaque;
  vector<QMatrix4x4> rotationMatrixOpaque;
  vector<QMatrix4x4> scaleMatrixOpaque;
  vector<QVector4D> colorVectorOpaque;

  QOpenGLExtraFunctions *openGLFunctions =
      QOpenGLContext::currentContext()->extraFunctions();
  mesh->vao.bind();
  mesh->vbo.bind();

  int numberOfObjectsToRender = translationMatrix.size();
  for (int i = 0; i < numberOfObjectsToRender; i++) {
    if (colorVector[i].w() != 1)  // Transparent
    {
      translationMatrixTransparent.push_back(translationMatrix[i]);
      rotationMatrixTransparent.push_back(rotationMatrix[i]);
      scaleMatrixTransparent.push_back(scaleMatrix[i]);
      colorVectorTransparent.push_back(colorVector[i]);
    } else {
      translationMatrixOpaque.push_back(translationMatrix[i]);
      rotationMatrixOpaque.push_back(rotationMatrix[i]);
      scaleMatrixOpaque.push_back(scaleMatrix[i]);
      colorVectorOpaque.push_back(colorVector[i]);
    }
  }
  // Render transparent stuff
  disableDepthTest();
  enableBlend();
  int numberOfTransparentObjectsToRender = translationMatrixTransparent.size();
  for (int i = 0; i < numberOfTransparentObjectsToRender; i += 32) {
    int numberOfTransparentObjectsToRenderNow =
        numberOfTransparentObjectsToRender - i;
    if (numberOfTransparentObjectsToRenderNow > 32) {
      numberOfTransparentObjectsToRenderNow = 32;
    }
    if (renderingProgramIsBind) {
      renderingProgram->setUniformValueArray(
          translationMatrixLocation, &translationMatrixTransparent[i],
          numberOfTransparentObjectsToRenderNow);
      renderingProgram->setUniformValueArray(
          rotationMatrixLocation, &rotationMatrixTransparent[i],
          numberOfTransparentObjectsToRenderNow);
      renderingProgram->setUniformValueArray(
          scaleMatrixLocation, &scaleMatrixTransparent[i],
          numberOfTransparentObjectsToRenderNow);
      renderingProgram->setUniformValueArray(
          colorVectorLocation, &colorVectorTransparent[i],
          numberOfTransparentObjectsToRenderNow);
    }

    openGLFunctions->glDrawArraysInstanced(
        GL_TRIANGLES, 0, mesh->polygonCount(),
        numberOfTransparentObjectsToRenderNow);
  }
  disableBlend();
  enableDepthTest();
  // Render opaque stuff
  int numberOfOpaqueObjectsToRender = translationMatrixOpaque.size();
  for (int i = 0; i < numberOfOpaqueObjectsToRender; i += 32) {
    int numberOfOpaqueObjectsToRenderNow = numberOfOpaqueObjectsToRender - i;
    if (numberOfOpaqueObjectsToRenderNow > 32) {
      numberOfOpaqueObjectsToRenderNow = 32;
    }
    if (renderingProgramIsBind) {
      renderingProgram->setUniformValueArray(translationMatrixLocation,
                                             &translationMatrixOpaque[i],
                                             numberOfOpaqueObjectsToRenderNow);
      renderingProgram->setUniformValueArray(rotationMatrixLocation,
                                             &rotationMatrixOpaque[i],
                                             numberOfOpaqueObjectsToRenderNow);
      renderingProgram->setUniformValueArray(scaleMatrixLocation,
                                             &scaleMatrixOpaque[i],
                                             numberOfOpaqueObjectsToRenderNow);
      renderingProgram->setUniformValueArray(colorVectorLocation,
                                             &colorVectorOpaque[i],
                                             numberOfOpaqueObjectsToRenderNow);
    }

    openGLFunctions->glDrawArraysInstanced(GL_TRIANGLES, 0,
                                           mesh->polygonCount(),
                                           numberOfOpaqueObjectsToRenderNow);
  }

  mesh->vbo.release();
  mesh->vao.release();
}

void Renderer::enableDepthTest() {
  QOpenGLExtraFunctions *openGLFunctions =
      QOpenGLContext::currentContext()->extraFunctions();
  openGLFunctions->glEnable(GL_DEPTH_TEST);
}

void Renderer::disableDepthTest() {
  QOpenGLExtraFunctions *openGLFunctions =
      QOpenGLContext::currentContext()->extraFunctions();
  openGLFunctions->glDisable(GL_DEPTH_TEST);
}

void Renderer::enableBlend() {
  QOpenGLExtraFunctions *openGLFunctions =
      QOpenGLContext::currentContext()->extraFunctions();
  openGLFunctions->glEnable(GL_BLEND);
  openGLFunctions->glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
}

void Renderer::disableBlend() {
  QOpenGLExtraFunctions *openGLFunctions =
      QOpenGLContext::currentContext()->extraFunctions();
  openGLFunctions->glDisable(GL_BLEND);
}

void Renderer::disableAll() {
  QOpenGLExtraFunctions *openGLFunctions =
      QOpenGLContext::currentContext()->extraFunctions();
  openGLFunctions->glDisable(GL_DEPTH_TEST);
  openGLFunctions->glDisable(GL_CULL_FACE);
}

void Renderer::enableAll() {
  QOpenGLExtraFunctions *openGLFunctions =
      QOpenGLContext::currentContext()->extraFunctions();
  openGLFunctions->glEnable(GL_DEPTH_TEST);
  openGLFunctions->glEnable(GL_CULL_FACE);
}

int Renderer::getPixelValueAt(QPoint mousePosition, QSize screenSize) {
  QOpenGLExtraFunctions *openGLFunctions =
      QOpenGLContext::currentContext()->extraFunctions();
  openGLFunctions->glFinish();

  float *pixelValues = new float[3];
  openGLFunctions->glReadPixels(mousePosition.x(),
                                screenSize.height() - mousePosition.y(), 1, 1,
                                GL_RGB, GL_FLOAT, pixelValues);
  // cout<<"pixel values = "<<pixelValues[0]<<" "<<pixelValues[1]<<"
  // "<<pixelValues[2]<<endl;
  int ID;
  if ((int)(pixelValues[2] * 1000) != 345) {
    ID = -1;
  } else {
    ID = round(pixelValues[0] * 32) + round(pixelValues[1] * 32) * 32;
  }
  // cout<<ID<<endl;
  return ID;
}

QByteArray Renderer::versionedShaderCode(const char *src) {
  QByteArray versionedSrc;

  if (QOpenGLContext::currentContext()->isOpenGLES())
    versionedSrc.append(QByteArrayLiteral("#version 300 es\n"));
  else
    versionedSrc.append(QByteArrayLiteral("#version 330\n"));

  versionedSrc.append(src);
  return versionedSrc;
}
