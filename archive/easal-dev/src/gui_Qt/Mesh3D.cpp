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
/****************************************************************************
**
** Copyright (C) 2016 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of the examples of the Qt Toolkit.
**
** $QT_BEGIN_LICENSE:BSD$
** Commercial License Usage
** Licensees holding valid commercial Qt licenses may use this file in
** accordance with the commercial license agreement provided with the
** Software or, alternatively, in accordance with the terms contained in
** a written agreement between you and The Qt Company. For licensing terms
** and conditions see https://www.qt.io/terms-conditions. For further
** information use the contact form at https://www.qt.io/contact-us.
**
** BSD License Usage
** Alternatively, you may use this file under the terms of the BSD license
** as follows:
**
** "Redistribution and use in source and binary forms, with or without
** modification, are permitted provided that the following conditions are
** met:
**   * Redistributions of source code must retain the above copyright
**     notice, this list of conditions and the following disclaimer.
**   * Redistributions in binary form must reproduce the above copyright
**     notice, this list of conditions and the following disclaimer in
**     the documentation and/or other materials provided with the
**     distribution.
**   * Neither the name of The Qt Company Ltd nor the names of its
**     contributors may be used to endorse or promote products derived
**     from this software without specific prior written permission.
**
**
** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
** "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
** LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
** A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
** OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
** SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
** LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
** DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
** THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
** (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
** OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE."
**
** $QT_END_LICENSE$
**
****************************************************************************/

#include "Mesh3D.h"

#include <qmath.h>

Mesh3D::Mesh3D() {}

void Mesh3D::createAndBind(string filename) {
  getVertexDataFromFile(filename);
  bindVboAndVao();
}
void Mesh3D::getVertexDataFromFile(string filename) {
  ifstream file;
  string line;
  vector<float> vertexData;
  vector<float> normalData;
  vector<int> vertexIndicies;
  vector<float> normalIndicies;

  data.clear();
  file.open(filename);
  if (!file.is_open()) {
    cout << "Wasn't able to open file " << filename
         << " that hold vertex data for rendering" << endl;
  }
  //[ OBJ PARSER]
  while (getline(file, line)) {
    istringstream stringStream(line);
    string word;
    vector<string> words;
    while (getline(stringStream, word, ' ')) {
      words.push_back(word);
    }
    if (words.size() > 0) {
      if (words[0] == "v") {
        vertexData.push_back(atof(words[1].c_str()));
        vertexData.push_back(atof(words[2].c_str()));
        vertexData.push_back(atof(words[3].c_str()));
      } else if (words[0] == "vn") {
        normalData.push_back(atof(words[1].c_str()));
        normalData.push_back(atof(words[2].c_str()));
        normalData.push_back(atof(words[3].c_str()));
      } else if (words[0] == "f") {
        for (int i = 1; i < 4; i++) {
          char wordTmp[32];
          strcpy(wordTmp, words[i].c_str());
          char *token = strtok(wordTmp, "//");
          vertexIndicies.push_back(atoi(token) - 1);
          token = strtok(NULL, "//");
          normalIndicies.push_back(atoi(token) - 1);
        }
      }
    }
  }
  //[~OBJ PARSER]
  for (int i = 0; i < vertexIndicies.size(); i++) {
    // vetex
    // x
    data.push_back(vertexData[vertexIndicies[i] * 3]);
    // y
    data.push_back(vertexData[vertexIndicies[i] * 3 + 1]);
    // z
    data.push_back(vertexData[vertexIndicies[i] * 3 + 2]);

    // normal for vertex
    // x
    data.push_back(normalData[normalIndicies[i] * 3]);
    // y
    data.push_back(normalData[normalIndicies[i] * 3 + 1]);
    // z
    data.push_back(normalData[normalIndicies[i] * 3 + 2]);
  }
  numberOfFloats = vertexIndicies.size() * 6;

  dataWasSet = true;
}

void Mesh3D::bindVboAndVao() {
  if (dataWasSet) {
    QOpenGLFunctions *openGLFunctions =
        QOpenGLContext::currentContext()->functions();

    // create
    vao.create();
    vbo.create();
    // bind
    vao.bind();
    vbo.bind();

    // fill with data
    vbo.allocate(data.constData(), numberOfFloats * sizeof(GLfloat));
    openGLFunctions->glEnableVertexAttribArray(0);
    openGLFunctions->glEnableVertexAttribArray(1);
    openGLFunctions->glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE,
                                           6 * sizeof(GLfloat), 0);
    openGLFunctions->glVertexAttribPointer(
        1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(GLfloat),
        reinterpret_cast<void *>(3 * sizeof(GLfloat)));

    //        openGLFunctions->glDisableVertexAttribArray(0);
    //        openGLFunctions->glDisableVertexAttribArray(1);
    // release
    vbo.release();
    vao.release();
  }
}
