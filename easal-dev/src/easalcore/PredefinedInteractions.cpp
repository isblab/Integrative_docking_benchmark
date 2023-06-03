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
 * PredefinedInteractions.cpp
 *
 *  Created on: Jun 23, 2010
 *      Author: jrpence
 */

#include "PredefinedInteractions.h"

#include "Atom.h"
#include "Settings.h"
#include "readIn.h"

using namespace std;
/*PredefinedInteractions::PredefinedInteractions(){
                this->dist2Defined = false;
                this->distTableDefined = false;
        }*/

PredefinedInteractions::PredefinedInteractions() {
  Settings *sett = Settings::getInstance();
  vector<vector<string> > data = readData(sett->DistanceData.file);

  if (data.size() == 0) {
    this->dist2Defined = false;
    this->distTableDefined = false;
    return;
  }
  bool matrixInsteadOfColumns = false;
  // Check if the data is in matrix format
  if (data.front().size() == data.size()) {
    for (std::vector<string>::size_type iter = 0; iter != data.size(); iter++) {
      if (data[iter][1] == data[1][iter])
        matrixInsteadOfColumns = true;
      else {
        matrixInsteadOfColumns = false;
        break;
      }
    }
  }

  PredefinedInteractions *output;
  if (matrixInsteadOfColumns) {
    output = new PredefinedInteractions(data);
  } else if (sett->DistanceData.radiusMin_col >= 0 &&
             sett->DistanceData.radiusMax_col >= 0) {
    output = new PredefinedInteractions(
        data, sett->DistanceData.label1_col, sett->DistanceData.label2_col,
        sett->DistanceData.radiusMin_col, sett->DistanceData.radiusMax_col);
  } else if (sett->DistanceData.radius_col >= 0) {
    output = new PredefinedInteractions(data, sett->DistanceData.label1_col,
                                        sett->DistanceData.label2_col,
                                        sett->DistanceData.radius_col);
  } else if (sett->DistanceData.radiusMin_col >= 0) {
    output = new PredefinedInteractions(data, sett->DistanceData.label1_col,
                                        sett->DistanceData.label2_col,
                                        sett->DistanceData.radiusMin_col);
  } else {
    output = new PredefinedInteractions();
  }

  this->assign(output);
  delete output;
}

PredefinedInteractions::PredefinedInteractions(vector<vector<string> > values,
                                               size_t labelAcol,
                                               size_t labelBcol,
                                               size_t radiuscol) {
  this->dist2Defined = false;
  this->distTableDefined = false;
  cout << "distfinder" << endl;
  for (vector<vector<string> >::iterator rowIter = values.begin();
       rowIter != values.end(); rowIter++) {
    for (int i = 0; i < rowIter->size(); i++) {
      cout << rowIter->at(i) << ",";
    }
    cout << endl;
    string labelA = (*rowIter)[labelAcol];
    string labelB = (*rowIter)[labelBcol];
    stringstream mss;
    double min;
    mss << (*rowIter)[radiuscol];
    mss >> min;

    this->dist1[make_pair(labelA, labelB)] = min;
    //		this->minimums[make_pair(labelB,labelA)] = min;
    this->distTableDefined = true;
  }
}

PredefinedInteractions::PredefinedInteractions(vector<vector<string> > values,
                                               size_t labelAcol,
                                               size_t labelBcol,
                                               size_t radiusmincol,
                                               size_t radiusmaxcol) {
  this->dist2Defined = false;
  this->distTableDefined = false;
  for (vector<vector<string> >::iterator rowIter = values.begin();
       rowIter != values.end(); rowIter++) {
    string labelA = (*rowIter)[labelAcol];
    string labelB = (*rowIter)[labelBcol];

    stringstream mss;
    double min, max;
    mss << (*rowIter)[radiusmincol] << " " << (*rowIter)[radiusmaxcol];
    mss >> min >> max;

    this->dist1[make_pair(labelA, labelB)] = min;
    //		this->minimums[make_pair(labelB,labelA)] = min;
    this->distTableDefined = true;
    this->dist2[make_pair(labelA, labelB)] = max;
    //		this->maximums[make_pair(labelB,labelA)] = max;
    this->dist2Defined = true;
  }
}

PredefinedInteractions::PredefinedInteractions(vector<vector<string> > values) {
  this->dist2Defined = false;
  this->distTableDefined = false;
  size_t i, j;
  vector<string> lblA, lblB;  // header

  for (i = 0; i < values.size(); i++)
    for (j = 0; j < values.size(); j++) {
      if (i > 0 && j > 0) {
        stringstream mss;
        double mn;
        mss << values[i][j];
        mss >> mn;

        this->dist1[make_pair(lblA[j], lblB[i])] = mn;
        this->dist1[make_pair(lblB[i], lblA[j])] = mn;
        this->distTableDefined = true;
      }
      if (i == 0) {  // first row
        lblA.push_back(values[i][j]);
      }
      if (j == 0) {  // first column
        lblB.push_back(values[i][j]);
      }
    }
}

PredefinedInteractions::~PredefinedInteractions() {
  // TODO Auto-generated destructor stub
}

double PredefinedInteractions::sumOfRadii(Atom *atomA, Atom *atomB) {
  return atomA->getRadius() + atomB->getRadius();
}

double PredefinedInteractions::collisionLowerBound(Atom *atomA, Atom *atomB) {
  double distance = getDist1(
      atomA, atomB);  // the predetermined pair of atoms that may cause bounds
  if (distance ==
      0)  // 0 if the pair of two atoms are not found in the dist table.
    return collisionLower(sumOfRadii(atomA, atomB));
  else {
    if (this->dist2Defined)
      return distance;
    else
      return collisionLower(distance);
  }
}

double PredefinedInteractions::bondingLowerBound(Atom *atomA, Atom *atomB) {
  Settings *sett = Settings::getInstance();
  if (sett->General.candidate_interactions) {
    double distance =
        getDist1(atomA, atomB);  // only interface atom pairs may cause bounds
    if (distance ==
        0)  // 0 if the pair of two atoms are not found in the dist table.
      return std::numeric_limits<double>::infinity();
    else {
      if (this->dist2Defined)
        return distance;
      else
        return bondLower(distance);
    }

  } else
    return bondLower(sumOfRadii(atomA, atomB));
}

double PredefinedInteractions::bondingUpperBound(Atom *atomA, Atom *atomB) {
  Settings *sett = Settings::getInstance();
  if (sett->General.candidate_interactions) {
    double distance =
        getDist1(atomA, atomB);  // only interface atom pairs may cause bounds
    if (distance ==
        0)  // 0 if the pair of two atoms are not found in the dist table.
      return -std::numeric_limits<double>::infinity();
    else {
      if (this->dist2Defined)
        return getDist2(atomA, atomB);
      else
        return bondUpper(distance);
    }
  } else
    return bondUpper(sumOfRadii(atomA, atomB));
}

double PredefinedInteractions::collisionLower(double distance) {
  Settings *sett = Settings::getInstance();
  return sett->Constraint.collisionLambda * distance +
         sett->Constraint.collisionDelta;
}

double PredefinedInteractions::bondLower(double distance) {
  Settings *sett = Settings::getInstance();
  return sett->Constraint.bondingLowerLambda * distance +
         sett->Constraint.bondingLowerDelta;
}

double PredefinedInteractions::bondUpper(double distance) {
  Settings *sett = Settings::getInstance();
  return sett->Constraint.bondingUpperLambda * distance +
         sett->Constraint.bondingUpperDelta;
}

// todo write a method to check if a pair is bond or not
// and seperate it from mindist. and unify mindist with bondinglowerbound
double PredefinedInteractions::getDist1(Atom *atomA, Atom *atomB) {
  string labelA = atomA->getName();
  string labelB = atomB->getName();
  if (labelA == "" || labelB == "")
    return 0;  // contact distance between arbitrary atoms does not create a
               // child node

  if (!this->distTableDefined) return 0;

  pair<string, string> labels(std::move(labelA), std::move(labelB));
  if (find(this->dist1, labels)) return this->dist1[labels];

  return 0;
}

double PredefinedInteractions::getDist2(Atom *atomA, Atom *atomB) {
  string labelA = atomA->getName();
  string labelB = atomB->getName();
  if (labelA == "" || labelB == "") {
    return 0;  // contact distance between arbitrary atoms does not create a
               // child node
  }

  if (!this->dist2Defined) return 0;

  pair<string, string> labels(std::move(labelA), std::move(labelB));
  if (find(this->dist2, labels)) return this->dist2[labels];

  return 0;
}

void PredefinedInteractions::assign(PredefinedInteractions *other) {
  this->dist1.clear();
  for (map<pair<string, string>, double>::iterator iter = other->dist1.begin();
       iter != other->dist1.end(); iter++) {
    this->dist1[iter->first] = iter->second;
  }
  this->distTableDefined = this->dist1.size() > 0;

  this->dist2.clear();
  for (map<pair<string, string>, double>::iterator iter = other->dist2.begin();
       iter != other->dist2.end(); iter++) {
    this->dist2[iter->first] = iter->second;
  }

  this->dist2Defined = this->dist2.size() > 0;
}

// what about reverse order of val? does not it have the same value?
bool PredefinedInteractions::find(map<pair<string, string>, double> maap,
                                  pair<string, string> val) {
  return (maap.find(val) != maap.end());
}

bool PredefinedInteractions::empty() {
  if (this->distTableDefined) {
    return this->dist1.empty();
  } else {
    return true;
  }
}
