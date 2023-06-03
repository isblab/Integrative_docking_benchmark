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
#ifndef MESSAGEWINDOW_H
#define MESSAGEWINDOW_H

#include <QtWidgets>
#include <string>

using namespace std;

class MessageWindow : public QDialog {
  Q_OBJECT
 public:
  MessageWindow(QString title);
  void clearText();
  void addLineToMessage(QString text);
  void addLineToMessage(string text);
  bool isEmpty();
 private Q_SLOTS:
  void acceptButtonSlot();

 private:
  int numberOfLines = 0;
  QPushButton acceptButton;
  QLabel messageLable;
  QVBoxLayout verticalBoxLayout;
};

#endif  // MESSAGEWINDOW_H
