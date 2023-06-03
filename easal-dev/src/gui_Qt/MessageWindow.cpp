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
#include "MessageWindow.h"

MessageWindow::MessageWindow(QString title) {
  this->setWindowTitle(title);
  this->setMaximumWidth(700);
  this->setMinimumWidth(200);

  //    messageLable.setWordWrap(true);
  messageLable.setText("Warning: ");

  acceptButton.setText("Accept");
  connect(&acceptButton, SIGNAL(pressed()), this, SLOT(acceptButtonSlot()));

  verticalBoxLayout.addWidget(&messageLable);
  verticalBoxLayout.addWidget(&acceptButton);
  this->setLayout(&verticalBoxLayout);
}

void MessageWindow::clearText() {
  messageLable.setText("");
  numberOfLines = 0;
}

void MessageWindow::addLineToMessage(QString text) {
  numberOfLines++;
  messageLable.setText(messageLable.text() + "\n" +
                       QString::number(numberOfLines) + ") " + text);
}

void MessageWindow::addLineToMessage(string text) {
  numberOfLines++;
  messageLable.setText(messageLable.text() + "\n" +
                       QString::number(numberOfLines) + ") " +
                       QString::fromStdString(text));
}

bool MessageWindow::isEmpty() { return numberOfLines == 0; }

void MessageWindow::acceptButtonSlot() { this->close(); }
