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
#include "MolecularUnit.h"

#include <algorithm>

#include "Settings.h"

using namespace std;

MolecularUnit::MolecularUnit() {}

MolecularUnit::MolecularUnit(vector<Atom *> atoms) { this->atoms = atoms; }

void MolecularUnit::init_MolecularUnit_A_from_settings(
    PredefinedInteractions *df) {
  Settings *sett = Settings::getInstance();
  readMolecularUnitFromFile(
      sett->MolecularUnitA.file, sett->MolecularUnitA.ignored_rows,
      sett->MolecularUnitA.x_col, sett->MolecularUnitA.y_col,
      sett->MolecularUnitA.z_col, sett->MolecularUnitA.radius_col,
      sett->MolecularUnitA.label_col, sett->MolecularUnitA.atomNo_col);
  this->buildStree();
  if (sett->General
          .candidate_interactions)  // to just use interface atoms for detecting
                                    // boundaries, new contacts etc. from now on
                                    // it will ignore the rest even for
                                    // collision
    this->simplify(*df, true);  // unless sett->Constraint::wholeCollision is
                                // set to true. Then stree which holds all
                                // initial atoms will be used for collison check
}

void MolecularUnit::init_MolecularUnit_B_from_settings(
    PredefinedInteractions *df) {
  Settings *sett = Settings::getInstance();
  readMolecularUnitFromFile(
      sett->MolecularUnitB.file, sett->MolecularUnitB.ignored_rows,
      sett->MolecularUnitB.x_col, sett->MolecularUnitB.y_col,
      sett->MolecularUnitB.z_col, sett->MolecularUnitB.radius_col,
      sett->MolecularUnitB.label_col, sett->MolecularUnitB.atomNo_col);
  this->buildStree();
  if (sett->General.candidate_interactions) this->simplify(*df, true);
}
void MolecularUnit::readMolecularUnitFromFile(
    string filename, vector<int> ignored_rows, unsigned int x_col,
    unsigned int y_col, unsigned int z_col, int radius_col, int label_col,
    int atomNo_col) {
  //    MolecularUnit *this = new MolecularUnit();
  // open file
  ifstream file;
  file.open(filename.c_str());

  // check that the file was good.
  if (!file.is_open()) {
    cout << "File \"" << filename << "\" does not exist." << endl;
    //        delete ret_mu;
    //        return ret_mu;
  }

  // sort the ignored rows vector for simplified read in logic
  sort(ignored_rows.begin(), ignored_rows.end());

  // read file line by line
  for (int row = 0, ignore_idx = 0; !file.eof(); row++) {
    // we've encountered an ignored row
    if (ignore_idx < ignored_rows.size() && row == ignored_rows[ignore_idx]) {
      ignore_idx++;
      continue;
    }

    // read in the line
    string line;
    getline(file, line);

    if (line.empty() || line == " ") continue;

    // read in each column
    vector<string> line_data;
    stringstream line_ss(line);
    string s;
    while (line_ss >> s) line_data.push_back(s);

    // read in column by column
    string label, atomNo, trash;
    double loc[3], radius;
    for (int col = 0; col < line_data.size(); col++) {
      stringstream col_ss(line_data[col]);

      if (col == atomNo_col)
        col_ss >> atomNo;
      else if (col == label_col)
        col_ss >> label;
      else if (col == x_col)
        col_ss >> loc[0];
      else if (col == y_col)
        col_ss >> loc[1];
      else if (col == z_col)
        col_ss >> loc[2];
      else if (col == radius_col)
        col_ss >> radius;
    }

    // default values
    if (radius_col < 0) radius = 1.2;
    if (label_col < 0) label = "";
    if (atomNo_col < 0) atomNo = "";

    // make and add atom to the molecular unit
    Atom *a = new Atom(loc, radius, label, atomNo);
    a->setLine(line_data, x_col, y_col, z_col);
    this->addAtom(a);
  }

  // close file
  file.close();

  // make it calc limits, then return
  this->calcLimits();
  //    return ret_mu;
}

MolecularUnit::MolecularUnit(const MolecularUnit &another) {
  this->atoms = another.atoms;
  for (int i = 0; i < 3; i++) {
    this->hi[i] = another.hi[i];
    this->lo[i] = another.lo[i];
  }
  this->stree = another.stree;
  this->index = another.index;
}

MolecularUnit::~MolecularUnit() {
  for (size_t x = 0; x < this->atoms.size(); x++) {
    delete (this->atoms[x]);
  }
}

void MolecularUnit::addAtom(Atom *a) { this->atoms.push_back(a); }

bool MolecularUnit::contains(Atom *a) {
  return (find(this->atoms.begin(), this->atoms.end(), a) != this->atoms.end());
}

bool MolecularUnit::empty() { return this->atoms.empty(); }

size_t MolecularUnit::size() { return this->atoms.size(); }

Atom *MolecularUnit::getAtomAt(size_t index) {
  // if(index >= 0 && index < atoms.size()) {
  if (index < atoms.size()) {
    return this->atoms.at(index);
  }
  cerr << "getAtomAt: NULL " << index << endl;
  return NULL;
}

Atom *MolecularUnit::getAtomByLabel(string label) {
  for (vector<Atom *>::iterator ait = atoms.begin(); ait != atoms.end();
       ait++) {
    if ((*ait)->getName() == label) return *ait;
  }
  cerr << "can't find the atom with label  " << label << endl;
  return NULL;
}

int MolecularUnit::getIndexOf(Atom *a) {
  int output = -1;
  for (size_t i = 0; (i < this->atoms.size()) && (output < 0); i++) {
    if (this->atoms.at(i) == a) {
      output = (int)i;
    }
  }
  return output;
}

vector<Atom *> MolecularUnit::getAtoms() { return this->atoms; }

vector<Atom *> MolecularUnit::getXFAtoms(double from[][3], double to[][3],
                                         size_t repeat) {
  vector<double> trans_mat = Utils::getTransMatrix(from, to);
  return getXFAtoms(trans_mat, repeat);
}

vector<Atom *> MolecularUnit::getXFAtoms(const vector<double> &trans_mat,
                                         size_t repeat) {
  vector<Atom *> output;
  for (size_t iter = 0; iter < this->atoms.size(); iter++) {
    Atom *current = getXFAtom(iter, trans_mat, repeat);
    output.push_back(current);
  }
  return output;
}

Atom *MolecularUnit::getXFAtom(size_t index, const vector<double> &trans_mat,
                               size_t repeat) {
  Atom *current = new Atom(this->atoms[index]);
  double *l =
      current->getLocation();  // l is the actual location so don't delete
  Vector newL, preL;
  newL = Vector(l).trans(trans_mat);
  for (size_t mer = 0; mer < repeat; mer++) {
    preL = newL;
    newL = preL.trans(trans_mat);
  }
  current->setLocation(newL[0], newL[1], newL[2]);
  return current;

  // double *l = current->getLocation(); // l is the actual location so don't
  // delete double *newL =
  // Utils::matApp(l,from[0],from[1],from[2],to[0],to[1],to[2] ); //this is
  // slower than .trans method used above delete [] newL;
}

double *MolecularUnit::getTransformedCoordinates(
    int index, const std::vector<double> &trans_mat) {
  Atom *current = new Atom(this->atoms[index]);
  double *l =
      current->getLocation();  // l is the actual location so don't delete
  double *ret = new double[3];
  std::memcpy(ret, l, sizeof(double) * 3);
  Vector newL, preL;
  newL = Vector(ret).trans(trans_mat);
  for (int i = 0; i < 3; i++) {
    ret[i] = newL[i];
  }
  delete current;
  return ret;
}

double *MolecularUnit::getCoordinates(int index) {
  double *coordinates = this->getAtomAt(index)->getLocation();

  return coordinates;
}

void MolecularUnit::centerAtoms() {
  size_t iter;
  double sumx = 0, sumy = 0, sumz = 0;
  for (iter = 0; iter < this->atoms.size(); iter++) {
    double *l = (this->atoms[iter])->getLocation();
    sumx += l[0];
    sumy += l[1];
    sumz += l[2];
  }
  sumx = sumx / this->atoms.size();
  sumy = sumy / this->atoms.size();
  sumz = sumz / this->atoms.size();

  for (iter = 0; iter < this->atoms.size(); iter++) {
    double *l = (this->atoms[iter])->getLocation();
    (this->atoms[iter])->setLocation(l[0] - sumx, l[1] - sumy, l[2] - sumz);
  }
}

vector<pair<int, int> > MolecularUnit::getDumbbells() {
  Settings *sett = Settings::getInstance();
  vector<pair<int, int> > output;
  vector<Atom *>::iterator iter1, iter2;

  double minAtom = 0, maxAtom = 0;
  int mainAxis = 2;  // assumes z axis is main axis of helix

  bool middleDumbbells = false;
  if (sett->RootNodeCreation.participatingAtomIndex_low <
      sett->RootNodeCreation.participatingAtomIndex_high)
    middleDumbbells = true;

  if (middleDumbbells) {
    minAtom = this->atoms[sett->RootNodeCreation.participatingAtomIndex_low]
                  ->getLocation()[mainAxis];
    maxAtom = this->atoms[sett->RootNodeCreation.participatingAtomIndex_high]
                  ->getLocation()[mainAxis];
    cout << "minAtom " << minAtom << " maxAtom " << maxAtom << endl;
  }

  int i = 0;
  for (iter1 = this->atoms.begin(); iter1 != this->atoms.end(); iter1++) {
    int j = i;

    for (iter2 = iter1; iter2 != this->atoms.end(); iter2++) {
      Atom *a, *b;
      a = *iter1;
      b = *iter2;

      if (!middleDumbbells ||
          (a->getLocation()[mainAxis] >= minAtom &&
           a->getLocation()[mainAxis] <= maxAtom &&
           b->getLocation()[mainAxis] >= minAtom &&
           b->getLocation()[mainAxis] <=
               maxAtom)) {  // Specifcally looking for  pairs near the middle of
                            // the helices

        double curDist = Utils::dist(a->getLocation(), b->getLocation());

        if (curDist < sett->RootNodeCreation.initial4DContactSeparation_high &&
            curDist > sett->RootNodeCreation.initial4DContactSeparation_low) {
          pair<int, int> result;
          result.first = i;
          result.second = j;
          output.push_back(result);
        }
      }

      j++;
    }
    i++;
  }
  return output;
}

void MolecularUnit::getLimits(double max_xyz[3], double min_xyz[3]) {
  for (int i = 0; i < 3; i++) {
    max_xyz[i] = this->hi[i];
    min_xyz[i] = this->lo[i];
  }
}

void MolecularUnit::calcLimits() {
  for (size_t i = 0; i < this->atoms.size(); i++) {
    double *loc;
    loc = this->atoms[i]->getLocation();
    for (size_t j = 0; j < 3; j++) {
      if (i == 0) {
        this->hi[j] = loc[j];
        this->lo[j] = loc[j];
      }
      if (this->hi[j] < loc[j]) {
        this->hi[j] = loc[j];
      }

      if (this->lo[j] > loc[j]) {
        this->lo[j] = loc[j];
      }
    }
  }
}

void MolecularUnit::buildStree() { this->stree = new Stree(this->atoms); }

void MolecularUnit::simplify(PredefinedInteractions &distfinder,
                             bool first_label) {
  // build an index for fast query

  for (vector<Atom *>::iterator iter = atoms.begin(); iter != atoms.end();
       iter++) {
    index[(*iter)->getName()] = (*iter);
  }

  // clear the old atom array
  this->atoms.clear();

  // for each atom pair in the dist table, find it and
  // insert it into the atom array
  map<string, Atom *>::iterator idx_iter;
  for (PredefinedInteractions::dist_iterator iter = distfinder.dist1begin();
       iter != distfinder.dist1end(); iter++) {
    string label;
    if (first_label)
      label = iter->first.first;
    else
      label = iter->first.second;

    idx_iter = index.find(label);
    if (idx_iter == index.end()) {
      cerr << "can't find atom with label " << label << endl;
      continue;
    } else {
      // check if the atom is alread in the atoms vector
      // void duplication
      bool atom_exist = false;
      for (int i = 0; i < atoms.size(); i++) {
        if (atoms[i]->getName() == label) {
          atom_exist = true;
          break;
        }
      }

      if (!atom_exist) {
        atoms.push_back(idx_iter->second);
        cout << idx_iter->second->getLocation() << endl;
      }
    }
  }

  // re-calculate the bound
  this->calcLimits();
}

Stree *MolecularUnit::getStree() { return this->stree; }

ostream &operator<<(ostream &os, MolecularUnit &h) {
  size_t i;
  for (i = 0; i < h.atoms.size(); i++) {
    os << "[" << i << "]" << *(h.atoms[i]);
  }
  return os;
}
