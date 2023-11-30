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
#include "Atom.h"

#include <stdlib.h>
//#include <iostream.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

using namespace std;

Atom::Atom() {
  for (int i = 0; i < 3; i++) this->loc[i] = 0.0f;
  this->radius = 0.0f;
}

Atom::Atom(Atom* copy) {
  for (int i = 0; i < 3; i++) {
    this->loc[i] = copy->loc[i];
  }
  this->radius = copy->radius;

  this->name = copy->name;
  this->atomID = copy->atomID;
  this->row = copy->row;
  this->xCol = copy->xCol;
  this->yCol = copy->yCol;
  this->zCol = copy->zCol;
}

Atom::Atom(double location[], double radius) {
  for (int i = 0; i < 3; i++) this->loc[i] = location[i];
  this->radius = radius;
  this->name = "";
  this->atomID = "";
}

Atom::Atom(double location[], double radius, string name, string atomID) {
  for (int i = 0; i < 3; i++) this->loc[i] = location[i];
  this->radius = radius;
  this->name = name;
  this->atomID = atomID;
}

Atom::~Atom() {}

void Atom::setLocation(double x, double y, double z) {
  this->loc[0] = x;
  this->loc[1] = y;
  this->loc[2] = z;
}

void Atom::setRadius(double radius) { this->radius = radius; }

void Atom::setName(string name) { this->name = name; }

void Atom::setAtomID(string atomID) { this->atomID = atomID; }

double Atom::getRadius() { return this->radius; }

string Atom::getName() { return this->name; }

vector<string> Atom::getRow() { return this->row; }

string Atom::getAtomID() { return this->atomID; }

void Atom::setLine(vector<string> linedata, size_t xcol, size_t ycol,
                   size_t zcol) {
  this->row.assign(linedata.begin(), linedata.end());
  this->xCol = xcol;
  this->yCol = ycol;
  this->zCol = zcol;
}

string Atom::getLine() {
  stringstream ss;

  bool pdb_format = false;

  if (!pdb_format) {
    for (size_t i = 0; i < this->row.size(); i++) {
      if (i == this->xCol)
        ss << this->loc[0] << "\t";
      else if (i == this->yCol)
        ss << this->loc[1] << "\t";
      else if (i == this->zCol)
        ss << this->loc[2] << "\t";
      else
        ss << this->row[i] << "\t";
    }

  } else {
    ss << setw(6) << left << this->row[0];
    ss << setw(5) << right << this->row[1];
    ss << "  ";
    ss << setw(3) << left << this->row[2];
    ss << setw(4) << right << this->row[3];
    ss << setw(2) << right << this->row[4];
    ss << setw(4) << right << this->row[5];
    ss << setw(4) << "";
    ss << setw(8) << setiosflags(ios::fixed) << setprecision(3) << right
       << this->loc[0];
    ss << setw(8) << setiosflags(ios::fixed) << setprecision(3) << right
       << this->loc[1];
    ss << setw(8) << setiosflags(ios::fixed) << setprecision(3) << right
       << this->loc[2];
    ss << setw(6) << right << this->row[9];
    ss << setw(6) << right << this->row[10];
    ss << setw(6) << right << this->row[11];
  }

  ss.flush();
  return ss.str();
}

double* Atom::getLocation() { return this->loc; }

ostream& operator<<(ostream& os, Atom& a) {
  os << "x= " << a.loc[0];
  os << "\ty= " << a.loc[1];
  os << "\tz= " << a.loc[2];
  os << "\name= " << a.name;
  os << "\trad= " << a.radius << endl;
  return os;
}
