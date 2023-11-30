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
 *  Created on: 2019
 *  Author: Rahul Prabhu
 */

#include "AboutWindow.h"

AboutWindow::AboutWindow() {
  this->resize(750, 500);
  this->setWindowTitle("About EASAL");

  // Create layout Widgets
  centralWidget = new QWidget(this);
  // verticalLayoutWidget_1 = new QWidget(centralWidget);
  horizontalLayoutWidget_1 = new QWidget(centralWidget);
  horizontalLayoutWidget_2 = new QWidget(centralWidget);
  horizontalLayoutWidget_3 = new QWidget(centralWidget);

  // Create the layout
  // verticalLayout_1 = new QVBoxLayout(verticalLayoutWidget_1);
  horizontalLayout_1 = new QHBoxLayout(horizontalLayoutWidget_1);
  horizontalLayout_2 = new QHBoxLayout(horizontalLayoutWidget_2);
  horizontalLayout_3 = new QHBoxLayout(horizontalLayoutWidget_3);

  // Create Geometry
  // verticalLayoutWidget_1->setGeometry(QRect(0, 0, 700, 450));
  horizontalLayoutWidget_1->setGeometry(QRect(50, 0, 470, 80));
  horizontalLayoutWidget_2->setGeometry(QRect(50, 80, 625, 500));
  horizontalLayoutWidget_3->setGeometry(QRect(50, 430, 470, 80));

  // Create Labels
  introduction = new QLabel(horizontalLayoutWidget_1);

  // Create Text box
  License = new QTextBrowser(horizontalLayoutWidget_2);

  // Create Button
  exitButton = new QPushButton(horizontalLayoutWidget_3);

  // Set text
  introduction->setText(
      "EASAL - Efficient Atlasing and Search of Assembly Landscapes");
  introduction->resize(470, 80);

  // textBrowser->setGeometry(QRect(10, 480, 521, 201));
  License->resize(650, 300);
  License->setHtml(QApplication::translate(
      "MainWindow",
      "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" "
      "\"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
      "<html><head><meta name=\"qrichtext\" content=\"1\" /><style "
      "type=\"text/css\">\n"
      "p, li { white-space: pre-wrap; }\n"
      "</style></head><body style=\" font-family:'Ubuntu'; font-size:11pt; "
      "font-weight:400; font-style:normal;\">\n"
      "<p style=\" margin-top:0px; margin-bottom:0px; margin-left:2px; "
      "margin-right:2px; -qt-block-indent:0; text-indent:0px;\"><span style=\" "
      "font-size:12pt;\">EASAL is free software: you can redistribute it "
      "and/or modify it under the terms of the GNU General Public License as "
      "published by the Free Software Foundation, either version 3 of the "
      "License, or (at your option) any later version.</span></p>\n"
      "<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; "
      "margin-left:2px; margin-right:2px; -qt-block-indent:0; text-indent:0px; "
      "font-size:12pt;\"><br /></p>\n"
      "<p style=\" margin-top:0px; margin-bottom:0px; margin-left:2px; "
      "margin-right:2px; -qt-block-indent:0; text-inden"
      "t:0px;\"><span style=\" font-size:12pt;\">EASAL is distributed in the "
      "hope that it will be useful, but WITHOUT ANY WARRANTY; without even the "
      "implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR "
      "PURPOSE.  See the GNU General Public License for more "
      "details.</span></p>\n"
      "<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; "
      "margin-left:2px; margin-right:2px; -qt-block-indent:0; text-indent:0px; "
      "font-size:12pt;\"><br /></p>\n"
      "<p style=\" margin-top:0px; margin-bottom:0px; margin-left:2px; "
      "margin-right:2px; -qt-block-indent:0; text-indent:0px;\"><span style=\" "
      "font-size:12pt;\">You should have received a copy of the GNU General "
      "Public License along with this program.  If not, see "
      "&lt;http://www.gnu.org/licenses/&gt;.</span></p></body></html>",
      Q_NULLPTR));

  exitButton->setText("Exit");
  exitButton->resize(QSize(100, 20));
  connect(exitButton, SIGNAL(pressed()), this, SLOT(exit()));
}

void AboutWindow::exit() { this->close(); }
