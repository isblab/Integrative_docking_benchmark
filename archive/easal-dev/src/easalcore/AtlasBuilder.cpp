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
 * AtlasBuilder.cpp
 *
 *  Created on: 2008
 *      Author: Aysegul Ozkan
 */

#include "AtlasBuilder.h"

#include "ActiveConstraintGraph.h"
#include "ActiveConstraintRegion.h"
#include "Atlas.h"
#include "CartesianRealizer.h"
#include "CayleyParameterization.h"
#include "CayleyPoint.h"
#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/Geometry"
#include "Eigen/LU"
#include "Eigen/SVD"
#include "MolecularUnit.h"
#include "Settings.h"
#include "Utils.h"
#include "pseudoAtlas.h"
using namespace Eigen;
using Eigen::MatrixXd;

#include <cassert>
#include <fstream>
#include <iterator>
#include <vector>
#include <glog/logging.h>
#include <glog/stl_logging.h>
using namespace std;

#define PI 3.14159265

int num_samples;

string constructSignature(vector<pair<int, int>> part);

vector<pair<int, int>> buildACGFromSignature(string signature) {
  vector<pair<int, int>> part;
  vector<string> parts;
  stringstream signatureStream(signature);
  string bond;

  while (getline(signatureStream, bond, '=')) {
    parts.push_back(bond);
  }

  for (auto bonds : parts) {
    stringstream bondStream(bonds);
    string a, b;
    getline(bondStream, a, '-');
    getline(bondStream, b, '-');
    part.push_back(make_pair(std::stoi(a), std::stoi(b)));
  }
  return part;
}

AtlasBuilder::AtlasBuilder(MolecularUnit* helA, MolecularUnit* helB,
                           SaveLoader* snl, PredefinedInteractions* df,
                           Atlas* atlas)
    : verbose(false) {
  this->snl = snl;
  this->atlas = atlas;

  if (this->verbose) cout << *helA << endl;
  if (this->verbose) cout << *helB << endl;
  if (this->verbose) cout << "AtlasBuilder" << endl;

  this->a = helA;
  this->b = helB;
  this->df = df;
}

AtlasBuilder::AtlasBuilder() : verbose(false) {}

AtlasBuilder::~AtlasBuilder() {
  // TODO Auto-generated destructor stub
}

void AtlasBuilder::setBasin0DNode(ActiveConstraintGraph* acg) {
  this->basin0DNode = acg;
}

Atlas* AtlasBuilder::getAtlas() { return atlas; }

void AtlasBuilder::setData(MolecularUnit* helA, MolecularUnit* helB,
                           SaveLoader* snl, PredefinedInteractions* df,
                           Atlas* atlas) {
  this->snl = snl;
  this->atlas = atlas;

  this->a = helA;
  this->b = helB;
  this->df = df;
}

void AtlasBuilder::create_initial_contactGraphs_for_virusCase() {
  Settings* sett = Settings::getInstance();
  if (this->verbose) cout << "setup" << endl;

  // for each pair of interaction in the distance table
  for (PredefinedInteractions::dist_iterator dit1 = this->df->dist2begin();
       dit1 != this->df->dist2end(); dit1++) {
    for (PredefinedInteractions::dist_iterator dit2 = dit1;
         dit2 != this->df->dist2end(); dit2++) {
      if (dit1 == dit2) continue;

      // get the 4 atoms
      Atom* a1 = a->getAtomByLabel(dit1->first.first);
      Atom* a2 = a->getAtomByLabel(dit2->first.first);

      Atom* b1 = b->getAtomByLabel(dit1->first.second);
      Atom* b2 = b->getAtomByLabel(dit2->first.second);

      if (this->verbose) {
        cout << "checking " << a1->getName() << " " << a2->getName() << endl;
        cout << "\t" << b1->getName() << " " << b2->getName() << endl;
      }

      // get the distance of the atom pairs
      double da = Utils::dist(a1->getLocation(), a2->getLocation());
      double db = Utils::dist(b1->getLocation(), b2->getLocation());

      if (this->verbose) cout << "da db: " << da << " , " << db << endl;

      // if the distance is outside the acceptable range, ignore it
      if (da < sett->RootNodeCreation.initial4DContactSeparation_low ||
          da > sett->RootNodeCreation.initial4DContactSeparation_high ||
          db < sett->RootNodeCreation.initial4DContactSeparation_low ||
          db > sett->RootNodeCreation.initial4DContactSeparation_high)
        continue;

      double ar1 = a1->getRadius();
      double ar2 = a2->getRadius();

      double br1 = b1->getRadius();
      double br2 = b2->getRadius();

      // if ((da + ar1 + ar2 >= db - br1 - br2 -
      // settings::Collision::threshold*2) && (da - ar1 - ar2 <= db + br1 + br2 +
      //settings::Collision::threshold*2) )

      // create the contact graph of pair of atom and add it into the rootGraph
      // list
      vector<pair<int, int>> parts;
      parts.push_back(make_pair(a->getIndexOf(a1), b->getIndexOf(b1)));
      parts.push_back(make_pair(a->getIndexOf(a2), b->getIndexOf(b2)));
      ActiveConstraintGraph* initial = new ActiveConstraintGraph(parts);
      this->rootGraphs.push_back(make_pair(initial, true));
    }
  }
}

void AtlasBuilder::create_initial_contactGraphs(int dimension) {
  Settings* sett = Settings::getInstance();
  dimension = 6 - dimension;
  if (this->verbose) cout << "create_initial_contactGraphs" << endl;

  vector<Atom*> helA = a->getAtoms();
  vector<Atom*> helB = b->getAtoms();
  int m = helA.size();
  int n = helB.size();
  std::vector<pair<int, int>> contacts;

  if (sett->Sampling.sampleAllNodes) {
    for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++) {
        contacts.push_back(make_pair(i, j));
      }
    }

    std::vector<bool> v(contacts.size());
    std::fill(v.begin(), v.begin() + dimension, true);
    do {
      vector<pair<int, int>> parts;  // holds the indices of atoms
      for (size_t i = 0; i < contacts.size(); ++i) {
        if (v[i]) {
          parts.push_back(
              make_pair(std::get<0>(contacts[i]), std::get<1>(contacts[i])));
        }
      }
      ActiveConstraintGraph* initial = new ActiveConstraintGraph(parts);
      this->rootGraphs.push_back(
          make_pair(initial, true));  // contact_graphs with 2 contacts
    } while (std::prev_permutation(v.begin(), v.end()));
  } else {
    if (sett->RootNodeCreation.dimension_of_rootNodes == 5) {
      vector<pair<int, int>> parts;  // holds the indices of atoms
      parts.push_back(make_pair(sett->Sampling.initial_5D_Contact_1,
                                sett->Sampling.initial_5D_Contact_2));
      ActiveConstraintGraph* initial = new ActiveConstraintGraph(parts);
      this->rootGraphs.push_back(
          make_pair(initial, true));  // contact_graphs with 1 contacts
    } else {
      vector<pair<int, int>> parts;  // holds the indices of atoms
      parts.push_back(make_pair(sett->Sampling.initial_5D_Contact_1,
                                sett->Sampling.initial_5D_Contact_2));
      parts.push_back(make_pair(sett->Sampling.initial_5D_Contact_3,
                                sett->Sampling.initial_5D_Contact_4));
      ActiveConstraintGraph* initial = new ActiveConstraintGraph(parts);
      this->rootGraphs.push_back(
          make_pair(initial, true));  // contact_graphs with 1 contacts
    }
  }
}

bool AtlasBuilder::isAncestorOfBasin0D(ActiveConstraintGraph* graph) {
  if (baseAtlasBuilder) return true;
  vector<pair<int, int>> contacts1 = graph->getParticipants();
  vector<pair<int, int>> contacts2 = basin0DNode->getParticipants();

  /*if(contacts1.size() == contacts2.size()) {
    return false;
    }*/

  for (size_t i = 0; i < contacts1.size(); i++) {
    if (find(contacts2.begin(), contacts2.end(), contacts1[i]) ==
        contacts2.end()) {
      return false;
    }
  }
  return true;
}

void AtlasBuilder::findBasins() {
  vector<AtlasNode*> sampledNodes = atlas->getNodes();
  int count = 0;

  Settings* sett = Settings::getInstance();
  mkdir(sett->Basin.BasinDirectory.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
  for (vector<AtlasNode*>::iterator it = sampledNodes.begin();
       it != sampledNodes.end(); it++) {
    if (count > 99) break;

    // We are interested only in the 0D nodes
    if ((*it)->getDim() != 0) continue;

    // Create a New AtlasBuilder

    Atlas* basinAtlas = new Atlas();
#ifdef SERVER
    SaveLoader* basinSNL =
        new SaveLoader(sett->Basin.BasinDirectory + to_string((*it)->getID()),
                       a, b, snl->getMongoDBClient());
#else
    SaveLoader* basinSNL = new SaveLoader(
        sett->Basin.BasinDirectory + to_string((*it)->getID()), a, b);
#endif
    // cout<<"Writing Basin Information to
    // "<<settings::Basin::BasinDirectory+to_string((*it)->getID());

    AtlasBuilder* AB = new AtlasBuilder(a, b, basinSNL, df, basinAtlas);
    AB->setBasin0DNode((*it)->getCG());
    AB->setup();
    AB->startAtlasBuilding();
    BasinAtlases.push_back(AB->getAtlas());
    BasinBottoms.push_back((*it)->getID());
    count++;
  }
}

void AtlasBuilder::generateBasin0DRootNodes() {
  vector<pair<int, int>> contacts = basin0DNode->getParticipants();

  // vector<pair<int, int> > contacts = basin0D->getParticipants();
  int r = 1;

  std::vector<bool> v(contacts.size());
  std::fill(v.begin(), v.begin() + r, true);
  do {
    vector<pair<int, int>> parts;  // holds the indices of atoms
    for (int i = 0; i < (int)contacts.size(); ++i) {
      if (v[i]) {
        parts.push_back(
            make_pair(std::get<0>(contacts[i]), std::get<1>(contacts[i])));
      }
    }
    ActiveConstraintGraph* initial = new ActiveConstraintGraph(parts);
    this->rootGraphs.push_back(
        make_pair(initial, true));  // contact_graphs with 2 contacts
  } while (std::prev_permutation(v.begin(), v.end()));
}

/*void AtlasBuilder::create_initial_4d_contactGraphs_usingDumbbells() {
  if(this->verbose)cout << "create_initial_4d_contactGraphs_usingDumbbells"
<<endl;


  vector<pair<int,int> > vectorA,vectorB; // these are used to store the
dumbbells generated within each helix. vectorA = a->getDumbbells(); // getting
the dumbbell candidates vectorB = b->getDumbbells();

  if(this->verbose)cout << "got dumbbells setA["<<vectorA.size()<<"] setB["<<
vectorB.size()<<"]" << endl; vector<pair<int,int> >::iterator iter,it;

//setup to be bi-incident?
// ar1 ---da---ar2 (helix a)
//  |           |
//  |           |
// br1 ---db---br2 (helix b)
//Each of the variables get set and then the comparison is made.

for(iter = vectorA.begin();iter != vectorA.end();iter++)
for(it = vectorB.begin();it != vectorB.end();it++){

double z_slide_first  = abs( (a->getAtomAt((*iter).first))->getLocation()[2] -
(b->getAtomAt((*it).first))->getLocation()[2]  ) ;   //i can do this because
they are originally aligned to same position in the z axis double z_slide_second
= abs( (a->getAtomAt((*iter).second))->getLocation()[2] -
(b->getAtomAt((*it).second))->getLocation()[2] ) ;   //i can do this because
they are originally aligned to same position in the z axis


double da = Utils::dist((a->getAtomAt((*iter).first))->getLocation(),
(a->getAtomAt((*iter).second))->getLocation() ); double ar1 =
(a->getAtomAt((*iter).first))->getRadius(); double ar2 =
(a->getAtomAt((*iter).second))->getRadius(); double ar12 = ar1 + ar2;
//	double ar12 = a->getAtomAt( (*iter).first )->getMinDist( a->getAtomAt(
(*iter).second) );  //if( ar12 == -1 )  ar12 = ar1 + ar2;

double db = Utils::dist((b->getAtomAt((*it).first))->getLocation(),
(b->getAtomAt((*it).second))->getLocation() ); double br1 =
(b->getAtomAt((*it).first))->getRadius(); double br2 =
(b->getAtomAt((*it).second))->getRadius(); double br12 = br1 + br2;
//	double br12 = b->getAtomAt( (*it).first )->getMinDist( b->getAtomAt(
(*it).second) );

//are they able to be bi-incident?
if(  (da + ar12 >= db - br12) && ( da - ar12 <= db + br12)  )
{

if( !settings::RootNodeCreation::useParticipatingAtomZDistance ||
(z_slide_first<settings::RootNodeCreation::ParticipatingAtomZDistance &&
z_slide_second<settings::RootNodeCreation::ParticipatingAtomZDistance) ) //(
abs((*iter).first - (*it).first)<4  &&   abs((*iter).second - (*it).second)<4  )
//to choose dumbbells proportional to the place in the helix (the other
direction is ignored now)
{
//first paring
vector<pair<int,int> > parts; //holds the indices of atoms
parts.push_back(make_pair((*iter).first,(*it).first) );
parts.push_back(make_pair((*iter).second,(*it).second) );
ActiveConstraintGraph* initial = new
ActiveConstraintGraph(parts,this->a,this->b); this->rootGraphs.push_back(
make_pair(initial, true) );//contact_graphs with 2 contacts

}


if(settings::RootNodeCreation::reversePairDumbbells )	//reverse pairing
{
double z_slide_first_reverse  = abs(
(a->getAtomAt((*iter).first))->getLocation()[2] -
(b->getAtomAt((*it).second))->getLocation()[2]  ) ;   //i can do this because
they are originally aligned to same position in the z axis double
z_slide_second_reverse = abs( (a->getAtomAt((*iter).second))->getLocation()[2] -
(b->getAtomAt((*it).first))->getLocation()[2] ) ;   //i can do this because they
are originally aligned to same position in the z axis

if(!settings::RootNodeCreation::useParticipatingAtomZDistance ||  (
z_slide_first_reverse<settings::RootNodeCreation::ParticipatingAtomZDistance  &&
z_slide_second_reverse<settings::RootNodeCreation::ParticipatingAtomZDistance) )
//settings::RootNodeCreation::closeByDumbbellsAmount=5
{
vector<pair<int,int> > parts;
parts.push_back(make_pair((*iter).first,(*it).second) );
parts.push_back(make_pair((*iter).second,(*it).first) );
ActiveConstraintGraph* initial = new
ActiveConstraintGraph(parts,this->a,this->b); this->rootGraphs.push_back(
make_pair(initial, true) );

}
}
}
}

}*/

/*
void AtlasBuilder::create_initial_5d_contactGraphs() {
        Settings  *sett = Settings::getInstance();
        if(this->verbose)cout << "create_initial_5d_contactGraphs" <<endl;

        if(this->rootGraphs.size() != 0) {
                return;
        }

        vector<Atom*> helA = a->getAtoms();
        vector<Atom*> helB = b->getAtoms();

        for(size_t i=0; i<helA.size(); i++)
        {
                for(size_t j=0; j<helB.size(); j++)
                {

                        double z_slide  = abs( helA[i]->getLocation()[2] -
helB[j]->getLocation()[2]  ) ;   //i can do this because they are originally
aligned to same position in the z axis


                        if(
!sett->RootNodeCreation.useParticipatingAtomZDistance || z_slide <
sett->RootNodeCreation.ParticipatingAtomZDistance) //( abs(i- j)<4 )  //to
choose dumbbells proportional to the place in the helix (the other direction is
ignored now)
                        {
                                vector<pair<int,int> > parts; //holds the
indices of atoms parts.push_back( make_pair(i, j ) ); ActiveConstraintGraph*
initial = new ActiveConstraintGraph(parts); this->rootGraphs.push_back(
make_pair(initial, true) );//contact_graphs with 1 contacts

                        }
                }
        }

}*/

void AtlasBuilder::setup() {
  Settings* sett = Settings::getInstance();
  if (this->rootGraphs.size() != 0) {
    cout << "RootGraph was not empty. Erasing and re-creating." << endl;
    this->rootGraphs.erase(rootGraphs.begin(), rootGraphs.end());
  }
  if (this->verbose) cout << "setup" << endl;
  if (sett->Basin.BasinSampling && !baseAtlasBuilder) {
    generateBasin0DRootNodes();
  } else if (sett->General.candidate_interactions) {
    create_initial_contactGraphs_for_virusCase();
  } else if (sett->RootNodeCreation.dimension_of_rootNodes < 6 &&
             sett->RootNodeCreation.dimension_of_rootNodes > 0) {
    create_initial_contactGraphs(sett->RootNodeCreation.dimension_of_rootNodes);
  } else {
    cout << "Initial root graphs should be at least dimension 1" << endl;
    exit(1);
  }

  cout << "this->rootGraphs.size() " << this->rootGraphs.size() << endl;

  if (this->verbose)
    cout << "Created " << this->rootGraphs.size() << " contactIDs" << endl;
  for (list<pair<ActiveConstraintGraph*, bool>>::iterator iter =
           this->rootGraphs.begin();
       iter != this->rootGraphs.end(); iter++) {
    if (this->verbose) {
      cout << *((*iter).first) << endl;
    }
  }

  // reorder rootGraphs to start from the middle through end and then beginning
  // this is done because the middle rootGraphs are more likely be of interest
  // and we want to see them first.

  int rootGraphSize = this->rootGraphs.size();
  int one_third_size = rootGraphSize/3 - 1;
  list<pair<ActiveConstraintGraph*, bool>>::iterator iterl =
      this->rootGraphs.begin();
  for (int i = 0; i < one_third_size; i++)  // find the 1/3 rootGraph
    iterl++; 

  for (int i = rootGraphSize / 3; i > 0;
       i--)  // from 1/3 through the beginning, push rootGraphs at the end
  {
    this->rootGraphs.push_back(*iterl);
    iterl--;
  }
  for (int i = rootGraphSize / 3; i > 0; i--)  // remove first 1/3 rootGraphs
    this->rootGraphs.pop_front();
}

ActiveConstraintGraph* AtlasBuilder::getNextRootGraph(bool& empty) {
  empty = true;
  ActiveConstraintGraph* nextd;
  for (list<pair<ActiveConstraintGraph*, bool>>::iterator it =
           this->rootGraphs.begin();
       it != this->rootGraphs.end(); it++) {
    if ((*it).second)  // not done
    {
      nextd = (*it).first;
      (*it).second = false;  // set it done
      empty = false;         // able to find one more rootGraph
      break;
    }
  }
  return nextd;
}

void AtlasBuilder::refigureFlip(Orientation* orr,
                                std::vector<std::vector<int>> flipScheme) {
  // cout << "Old flip: " << orr->getFlipNum() << endl;
  // int flipOrig = orr->getFlipNum();

  Settings* sett = Settings::getInstance();

  std::vector<std::vector<int>> tetras = flipScheme;

  double fromBOrr[3][3];
  double toBOrr[3][3];
  orr->getFromTo(fromBOrr, toBOrr);

  vector<double> trans_mat = Utils::getTransMatrix(fromBOrr, toBOrr);

  // mirr list decides the flipNo
  vector<bool> mirr_list;
  vector<double*> tetraCoordinates;

  for (int x = 0; x < tetras.size(); x++) {
    tetraCoordinates.clear();
    for (int y = 0; y < 4; y++) {
      // see CayleyParaterization.h: for index scheme in 'tetrahedra'
      // if index > 5 then atom in helix b, need to do transform
      // if index <=5, then atom is in helix a, just get coordinates

      if (tetras[x][y] > 5) {
        tetraCoordinates.push_back(
            sett->runTimeObjects.muB->getTransformedCoordinates(
                (tetras[x][y] - 6), trans_mat));
      } else {
        tetraCoordinates.push_back(
            sett->runTimeObjects.muA->getCoordinates(tetras[x][y]));
      }
    }

    double* loc_1;
    double* loc_2;
    double* loc_3;
    double* loc_4;

    loc_1 = tetraCoordinates[0];
    loc_2 = tetraCoordinates[1];
    loc_3 = tetraCoordinates[2];
    loc_4 = tetraCoordinates[3];

    double* vec_12 = Utils::vectSub(loc_1, loc_2);
    double* vec_13 = Utils::vectSub(loc_1, loc_3);
    double* vec_14 = Utils::vectSub(loc_1, loc_4);

    double* orthoVector = Utils::crossProd(vec_12, vec_13);

    int tetraFlip = Utils::sign(Utils::dotProd(orthoVector, vec_14));

    if (tetraFlip == -1) {
      tetraFlip = 0;
    }
    mirr_list.push_back(tetraFlip);
    delete[] vec_12;
    delete[] vec_13;
    delete[] vec_14;
    delete[] orthoVector;
  }

  // see mirr_list in CartesianRealizer::computeRealization
  int flip = -1;
  for (int i = 0; i < 8; i++) {
    if (i % 2 == mirr_list[0] && (i / 2) % 2 == mirr_list[1] &&
        (i / 4) % 2 == mirr_list[2]) {
      flip = i;
    }
  }

  if (flip == -1) {
    cout << "ERROR: above should have found flip 0 < i 7" << endl;
  }

  orr->setFlipNum(flip);
  // cout <<"New flip: " << flip << " should match: " << orr->getFlipNum() <<
  // endl;
}

void AtlasBuilder::startAtlasBuilding() {
  Settings* sett = Settings::getInstance();
  int currentrootGraph = 0;
  bool empty;
  ActiveConstraintGraph* nextrootGraph = getNextRootGraph(empty);

  cout << "AtlasBuilder::startAtlasBuilding: this->rootGraphs.size() = "
       << this->rootGraphs.size() << endl;

  while (!empty) {
    if (!sett->AtlasBuilding.stop) {
      int nodenum;
      int success =
          this->atlas->addNode(nextrootGraph, nodenum, ROOT_NODE_PARENT_ID);
      if (success == 1) {
        cout << "AtlasBuilder::startAtlasBuilding: Started sampling rootGraph "
             << currentrootGraph + 1<< " out of " << this->rootGraphs.size()
             << " (node number " << nodenum << ")" << endl;

        AtlasNode* rnode = (*this->atlas)[nodenum];

        if (sett->AtlasBuilding.stop == true) {
          return;
        } else {
          if (sett->AtlasBuilding.breadthFirst) {
            doBreadthFirstSampling(rnode);
          } else {
            this->processAtlasNode(rnode, false, NULL, false, false);
          }
        }

        if (sett->AtlasBuilding.stop == true) {
          return;
        }
      }
      nextrootGraph = getNextRootGraph(empty);
      currentrootGraph++;
    }
  }

  samplingPostProcess();
  // Generate missing parents using 0D nodes through reverse witness and
  // re-sampling.
  if (sett->General.reverseWitness) {
    findAndSampleMissingAncestors();
  }

  cout
      << "FINISHED AtlasBuilder::startAtlasBuilding: this->rootGraphs.size() = "
      << this->rootGraphs.size() << endl;
}

void AtlasBuilder::doBreadthFirstSampling(AtlasNode* rnode) {
  std::pair<std::set<int>::iterator, bool> ret;  // to check duplicates
  std::set<int> myset;
  myset.insert(rnode->getID());

  std::queue<int> myqueue;
  myqueue.push(rnode->getID());

  while (!myqueue.empty()) {
    int nom = myqueue.front();
    rnode = (*atlas)[nom];
    myqueue.pop();
    bool cont = false;
    processAtlasNode(rnode, false, NULL, false, false);
    cont = true;

    vector<int> con = rnode->getConnection();
    size_t dim = rnode->getDim();
    for (vector<int>::iterator it = con.begin(); it != con.end(); it++) {
      AtlasNode* child_node = (*atlas)[*it];
      int child_dim = child_node->getDim();
      if (child_dim < dim) {  // if it is a lesser dimension to make sure it is
                              // child not parent
        ret = myset.insert(*it);
        if (ret.second)  // new element did not exist before
          myqueue.push(*it);
      }
    }
  }
}

void AtlasBuilder::determineStepSizeDynamically(
    ActiveConstraintGraph* acg, ActiveConstraintRegion* region, bool dense,
    CayleyParameterization* cparam) {
  Settings* sett = Settings::getInstance();
  if (this->verbose) cout << "determineStepSizeDynamically" << endl;

  double testStep = 0.4;  // was 0.2 before

  // temporarily disable the steric constraint to allow collision on parameters
  // bool userDefined_stericConstraint = settings::Constraint::stericConstraint;
  // settings::Constraint::stericConstraint = false;

  acg->setStepSize(
      testStep);  // increment this stepsize, if it is causing slowness
  ConvexChart* chart = new ConvexChart(
      acg, dense, cparam, this->df);  // for volume computation, to make blue
                                      // points as minimum as possible

  /**
   * GRID sampling allows some range for contact hence GRID sampling volume per
   * node is proportional with that range of the contact. i.e. Not all root
   * nodes have same volume by GRID sampling. If we are not doing short-range
   * sampling and forcing contact to be exact distance, then in order to have
   * proportional sampling density with GRID sampling, we should set
   * expectedNumberOfSamples per node to be proportional with that ratio
   */
  double contact_lengthUpper =
      df->bondingUpperBound(chart->getAtom(0), chart->getAtom(6));
  double contact_lengthLower =
      df->bondingLowerBound(chart->getAtom(0), chart->getAtom(6));
  double range = pow(contact_lengthUpper, 3) - pow(contact_lengthLower, 3);

  double expectedNumberOfSamples = 500.;  // 10000.; //50000.;
  expectedNumberOfSamples = expectedNumberOfSamples * range;

  /** compute the approximate volume of the region */
  int volume = 0;
  if (chart->initializeChart(false, region) && cparam->is_partial3tree())
    for (; !chart->doneGrid() && !sett->AtlasBuilding.stop; chart->stepGrid())
      volume++;

  delete chart;

  // settings::Constraint::stericConstraint = userDefined_stericConstraint;

  /**
   * tns:total_number_of_steps nspd:number_of_steps_per_dimension s:stepsize
   * nspd = pow(tns, 0.2)
   * nspd1*s1=nspd2*s2
   * s2 = s1* nspd1/nspd2 = s1* pow(tns1, 0.2)/pow(tns2, 0.2) = s1*
   * pow(tns1/tns2, 0.2) = s1 / pow(tns2/tns1, 0.2)
   */

  double voltimes = expectedNumberOfSamples / volume;
  double root = 1. / acg->getDim();
  double pwr = pow(voltimes, root);
  double stepsizeee = testStep / pwr;
  acg->setStepSize(stepsizeee);
}

void AtlasBuilder::createChildAcgFromBoundary(
    list<pair<int, int>>& contactList,
    list<Orientation*>::iterator ori_on_lattice,
    ActiveConstraintGraph* parentACG, bool bret, Orientation* orie_on_boundary,
    AtlasNode* parentNode, bool& noGoodOrientation, int& noPoints,
    CayleyPoint* entryCayleyPoint) {
  Settings* sett = Settings::getInstance();

  bool firstPath = true;
  bool savedOnce = false;
  for (auto contact : contactList) {
    ActiveConstraintGraph* childACG = new ActiveConstraintGraph(parentACG);
    childACG->addContact(contact);

    if (childACG->isDependent()) {
      delete childACG;
      continue;
    }

    Orientation* wit_orr_toSend = new Orientation(orie_on_boundary);

    wit_orr_toSend->setEntryPoint(parentNode->getID(),
                                  entryCayleyPoint->getID(),
                                  orie_on_boundary->getFlipNum());

    std::tuple<int, int, int> test = wit_orr_toSend->getEntryPoint();

    AtlasNode* child_node =
        createChildNode(parentNode, childACG, wit_orr_toSend, savedOnce,
                        noGoodOrientation, noPoints, bret);

    int child_nodeID = child_node->getID();
    orie_on_boundary->addBoundary(child_nodeID);
    (*ori_on_lattice)->addBoundary(child_nodeID);

    delete wit_orr_toSend;
  }
}

void AtlasBuilder::findBoundary(list<Orientation*>::iterator ori_on_lattice,
                                ConvexChart* desc, ActiveConstraintGraph* acg,
                                AtlasNode* rnode,
                                ActiveConstraintRegion* region,
                                ConstraintCheck* constraintChecker, bool bret,
                                bool& noGoodOrientation, int& noPoints,
                                bool& boundary_ori_found_and_saved,
                                CayleyPoint* EntryCayleyPoint) {
  // save the current grid point so we can come back where we were after binary
  // search
  vector<double> pp = desc->getPoint();
  Settings* sett = Settings::getInstance();

  // get the flip on which we want to run binary search
  int flip = (*ori_on_lattice)->getFlipNum();

  // look up and down for one step in each direction around a valid point.
  // If there is a point with collision, then that means there is a contact
  // boundary in between valid point and colliding point.
  // Then you need to do binary search in between to find out that boundary
  // config.
  while (!desc->stepAround() && !sett->AtlasBuilding.stop) {
    bool fail;
    num_samples++;

    vector<vector<int>> flipScheme = rnode->getFlipScheme();
    Orientation* orie_on_boundary =
        desc->computeRealization(flip, fail, flipScheme);
    // if needs maple, then it will always return first root not the one with
    // specified flip ! Orientation* orie_on_boundary =
    // CartesianRealizer::computeRealization(acg, desc, 					flip, fail, flipScheme);

    /// found collision, start binary search
    /// todo ACTUALLY CHECK SHOULD BE DONE AROUND BLUE REGION AS WELL,
    /// IF VOLUME NEGATIVE, SEARCH FOR BOUNDARY IN BETWEEN!!!
    if (sett->Sampling.binarySearch) {
      if (!fail && constraintChecker->stericsViolated(orie_on_boundary)) {
        desc->stepGridBinary(false);  ////walk through valid point

        delete orie_on_boundary;
        num_samples++;
        // orie_on_boundary = CartesianRealizer::computeRealization(acg, desc,
        //		flip, fail, flipScheme);
        orie_on_boundary = desc->computeRealization(flip, fail, flipScheme);

        // check for collision and contacts
        // if there is a collision then it should continue binary stepping
        bool binvalid = !constraintChecker->stericsViolated(orie_on_boundary);
        list<pair<int, int>> contactList2 =
            constraintChecker->getContacts(flip);

        // binary search
        // continue search till you find a valid configuration with new contacts
        // find new contact for new small threshold t // FIXME: what does this
        // line mean?
        int num_bin_iteration = 0;
        while ((contactList2.empty() || !binvalid) && num_bin_iteration < 30) {
          num_bin_iteration++;
          desc->stepGridBinary(binvalid);
          delete orie_on_boundary;
          num_samples++;
          // orie_on_boundary = CartesianRealizer::computeRealization(acg,
          //		desc, flip, fail, flipScheme);
          orie_on_boundary = desc->computeRealization(flip, fail, flipScheme);
          binvalid = !constraintChecker->stericsViolated(orie_on_boundary);
          contactList2 = constraintChecker->getContacts(flip);
        }

        // double check if the valid configuration with new contacts exists.
        // it may not in case it exits the loop because of num_bin_iteration is
        // big
        if (binvalid && contactList2.size() != 0 &&
            (sett->AtlasBuilding.ifBadAngleWitness_createChild ||
             !orie_on_boundary->angleViolated())) {
          // Orientation* orie_on_boundary = sr->getOrienation();
          createChildAcgFromBoundary(contactList2, ori_on_lattice, acg, bret,
                                     orie_on_boundary, rnode, noGoodOrientation,
                                     noPoints, EntryCayleyPoint);

          delete orie_on_boundary;

          break;  // if found a collision in one direction then stop
        }         // if( valid(orie_on_boundary,constraintChecker) &&
                  // contactList2.size() != 0
      }  // if( !fail &&  !valid(orie_on_boundary,constraintChecker)  )else
    }
    if (!boundary_ori_found_and_saved) {
      delete orie_on_boundary;
    }

    desc->setInitialPoint(pp);
  }  // while( !desc->stepGridContact()  && !settings::AtlasBuilding::stop)
  desc->setInitialPoint(pp);
  desc->setDir();
}

/*
 * This is a copy of the initial AtlasBuilder, though the actual AtlasBuilder
 may have changed a bit.
 *
 * Sample(CG(b1,b2,b3...bm) , k) :: m =12 & k=6
 desc = Get description ( CG(b1,b2,b3...bm) , k )
 :: find good parameters (depends on CG)
 & inequalities (depends on identity of  b1,b2,b3...bm)
 FOR each grid point in (convex)desc
 Real=find cartesian realizations w/out helix(CG(b1,b2,b3...bm),k)
 FOR each r in Real
 IF valid w/helix(r)
 getCG(r) = (CG(b1,b2,b3...bm'),k')
 output r with CG
 IF(CG(b1,b2,b3...bm'),k') != ( CG(b1,b2,b3...bm) , k )
 IF k < 6 and notdone(CG(b1,b2,b3...bm'),k')::CG & k used to identify
 Sample(CG(b1,b2,b3...bm'),k')
 FI
 done(CG(b1,b2,b3...bm'),k')
 FI
 FI
 ROF
 ROF
 done(CG(b1,b2,b3...bm),k)
 */
void AtlasBuilder::sampleNonConvexNodes() {
  cout << "Started sampling non-convex nodes" << endl;
  vector<AtlasNode*> nodes = this->atlas->getNodes();
  Settings* sett = Settings::getInstance();
  for (auto it : nodes) {
    bool partial3tree = true;
    if (it->isComplete() && it->isPartial3Tree(partial3tree) == 0) {
      if (!partial3tree) {
        pseudoAtlas* pa = new pseudoAtlas(snl);
        pa->processNodes(it);
        pa->pruneAtlas();
        unordered_map<string, vector<Orientation*>> nodeMap = pa->getNodeMap();
        // Iterate over the roadmap and reassign the points to different nodes;
        //- Delete the points in the original atlas node
        for (auto& region : nodeMap) {
          vector<pair<int, int>> part = buildACGFromSignature(region.first);
          ActiveConstraintGraph* acg = new ActiveConstraintGraph(part);
          CayleyParameterization* desc = new CayleyParameterization(acg, false);
          int child_nodeID;
          // int success = this->atlas->addNode(acg, child_nodeID);
          AtlasNode* child_node = (*this->atlas)[child_nodeID];

          AtlasNode* rnode = atlas->getNodes()[child_nodeID];
          ActiveConstraintRegion* acr = rnode->getACR();

					for(auto& oris:region.second) { 
						CayleyPoint* wp4 = new CayleyPoint(oris, acg->getParamLines(), a, b, acr->nextWitnessPointID());
						refigureFlip((wp4->getOrientations()[0]), rnode->getFlipScheme());
						this->snl->appendWitness(rnode->getID(), wp4);
                        delete wp4;
					}
                    delete child_node;
                    delete rnode;
                    delete acr;
                    delete acg;
                    delete desc;
				}
                delete pa;
      }
    }
  }
  cout << "Finished Sampling Non-Convex nodes" << endl;
}

void AtlasBuilder::samplingPostProcess() {
	//sampleNonConvexNodes();
	writeNonConvexEPTable();
}

void AtlasBuilder::writeNonConvexEPTable () {
  Settings* sett = Settings::getInstance();
  vector<AtlasNode*> nodes = this->atlas->getNodes();
  for (auto node : nodes) {
    bool partial3tree = true;
    if (node->isComplete() && node->isPartial3Tree(partial3tree) == 0) {
      if (!partial3tree) {
	  	if(sett->Paths.implementPathFinding && 
			node->getDim() <= sett->Paths.energyLevelUpperBound) {
			this->snl->writeEntryPointsTable(node);
		}
	  }
	}
  }
}
bool AtlasBuilder::processAtlasNode(AtlasNode* currentNode, bool denseSampling,
                                    Orientation* coming_witness_ori,
                                    bool resumeSampling, bool bret) {
  Settings* sett = Settings::getInstance();
  int currentNodeID = currentNode->getID();

  if (this->verbose) {
    cout << "AtlasBuilder::MySample: Started node " << currentNodeID << endl;
  }

  if (denseSampling) {
    currentNode->setComplete(false);
  }

  ActiveConstraintRegion* region = currentNode->getACR();
  ActiveConstraintGraph* acg = currentNode->getCG();
  CayleyParameterization* cParam = new CayleyParameterization(acg, false);

  if (sett->Sampling.dynamicStepSizeAmong) {
    determineStepSizeDynamically(acg, region, denseSampling, cParam);
  }

  if (sett->Output.writeNodeFiles) snl->saveNodeHeader(currentNode, atlas);

  /*if (!resumeSampling) {
    if (sett->Output.writeNodeFiles) this->snl->appendDimension(currentNode);
  }*/

  if (resumeSampling && !sett->AtlasBuilding.breadthFirst && !denseSampling) {
    vector<int> connectedNodes = currentNode->getConnection();
    size_t dim = currentNode->getDim();

    for (vector<int>::iterator it = connectedNodes.begin();
         it != connectedNodes.end(); it++) {
      AtlasNode* connectedNode = (*this->atlas)[*it];
      if (connectedNode->getDim() < dim && !connectedNode->isComplete()) {
        snl->loadNode(*it, connectedNode->getACR());
        this->processAtlasNode(connectedNode, false, NULL, resumeSampling,
                               false);
      }
    }
  }

  bool noGoodOrientation = true;
  if (resumeSampling) {
    noGoodOrientation = !currentNode->hasAnyGoodOrientation();
  }

  ConvexChart* chart = new ConvexChart(acg, denseSampling, cParam, this->df);
  if (currentNode->getFlipScheme().empty()) {
    if (chart->partial3tree) {
      currentNode->setFlipScheme(chart->getTetras());
    } else {
      std::vector<std::vector<int>> parent_flipScheme =
          this->atlas->getNodes()[currentNode->getFirstParentID()]
              ->getFlipScheme();

      currentNode->setFlipScheme(parent_flipScheme);
    }
  }

  if (coming_witness_ori != NULL) {
    CayleyPoint* wp4 = new CayleyPoint(coming_witness_ori, acg->getParamLines(),
                                       a, b, region->nextWitnessPointID());

    refigureFlip((wp4->getOrientations()[0]), currentNode->getFlipScheme());

    chart->setWitnessPoint(wp4);
    if (!coming_witness_ori->angleViolated()) {
      region->AddWitnessPoint(wp4);
      this->snl->appendWitness(currentNode->getID(), wp4);
      noGoodOrientation = false;
    }
    addReverseWitness(currentNode, acg, chart, region, coming_witness_ori,
                      noGoodOrientation);
  }

  if (bret) {  // bret true means this is child node of 'the node that is
               // sampled
               // by breathFirst'. Hence we need to return back.
    cout << "\nINSIDE bret block" << endl;
    if (this->verbose) cout << "Returning\n" << endl;
    region->trim();
    delete cParam;
    delete chart;
    return true;
  }

  sampleAtlasNode(currentNode, cParam, acg, chart, resumeSampling, region,
                  noGoodOrientation);

  delete cParam;
  delete chart;
#ifdef VERBOSE
  cout << "AtlasBuilder::MySample: Finished node " << currentNodeID << endl;
#endif

  return !noGoodOrientation;
}

void AtlasBuilder::sampleAtlasNode(AtlasNode* currentNode,
                                   CayleyParameterization* cParam,
                                   ActiveConstraintGraph* acg,
                                   ConvexChart* chart, bool resumeSampling,
                                   ActiveConstraintRegion* region,
                                   bool& noGoodOrientation) {
  int noPoints = 0;
  Settings* sett = Settings::getInstance();
  bool bret = sett->AtlasBuilding.breadthFirst;
  ConstraintCheck* constraintChecker = new ConstraintCheck(acg, this->df);
  // resumeSampling : Flag of whether this is a first/new sample or a continued
  // sample.
  if (cParam->is_partial3tree() &&
      chart->initializeChart(resumeSampling, region)) {
    list<Orientation*> realizationList;
    for (; !chart->doneGrid() && !sett->AtlasBuilding.stop; chart->stepGrid()) {
      // find all realizations for the grid point
      realizationList =
          findRealizations(acg, chart, currentNode->getFlipScheme());
      vector<double> currentCayleyParamValues = chart->getPoint();
      CayleyPoint* currentCayleyPoint =
          new CayleyPoint(currentCayleyParamValues);
      currentCayleyPoint->setID(region->nextSamplePointID());
      currentCayleyPoint->setRealizable(!realizationList.empty());
      currentCayleyPoint->SetTetraVolume(chart->TetraVolume());
      currentCayleyPoint->SetTetraEdges(chart->TetraEdges());

      if (!realizationList.empty()) {
        for (list<Orientation*>::iterator ori_on_lattice =
                 realizationList.begin();
             ori_on_lattice != realizationList.end() &&
             !sett->AtlasBuilding.stop;
             ori_on_lattice++) {
          // check for steric constraint violation and for any possible future
          // contacts
          if (!constraintChecker->stericsViolated(*ori_on_lattice)) {
            /// check for bad angle
            if ((*ori_on_lattice)->angleViolated())
              currentCayleyPoint->incrementBadAngleN();
            else
              noGoodOrientation = false;

            bool boundary_ori_found_and_saved = false;

            if (sett->RootNodeCreation.createChildren) {
              findBoundary(ori_on_lattice, chart, acg, currentNode, region,
                           constraintChecker, bret, noGoodOrientation, noPoints,
                           boundary_ori_found_and_saved, currentCayleyPoint);
            }
            if (!boundary_ori_found_and_saved) {
              if (!(*ori_on_lattice)->angleViolated()) {
                currentCayleyPoint->addOrientation((*ori_on_lattice));
              } else {
                delete (*ori_on_lattice);
              }
            } else {
              delete (*ori_on_lattice);
            }
          } else {
            if ((*ori_on_lattice)->angleViolated())
              currentCayleyPoint->incrementBadAngleN();

            delete (*ori_on_lattice);  // collision
            currentCayleyPoint->incrementCollidN();
          }
        }
      }

      region->AddSamplePoint(currentCayleyPoint);
      noPoints++;
      if (noPoints >= sett->Saving.savePointsFrequency) {
        if (sett->Output.writeNodeFiles)
          this->snl->appendSpacePoints(currentNode);
        noPoints = 0;
      }
    }
  }

  if (noPoints != 0) {
    if (sett->Output.writeNodeFiles) this->snl->appendSpacePoints(currentNode);
    noPoints = 0;
  }
  currentNode->setFoundGoodOrientation(!noGoodOrientation);
  if(sett->Paths.implementPathFinding && 
     currentNode->getDim() <= sett->Paths.energyLevelUpperBound) {
  	this->snl->writeEntryPointsTable(currentNode);
  }

  if (!sett->AtlasBuilding.stop) {
    currentNode->setComplete(true);
    region->trim();
  }
  delete constraintChecker;
}

AtlasNode* AtlasBuilder::createChildNode(AtlasNode* rnode,
                                         ActiveConstraintGraph* child_CG,
                                         Orientation* wit_orr_toSend,
                                         bool& savedOnce,
                                         bool& noGoodOrientation, int& noPoints,
                                         bool bret) {
  Settings* sett = Settings::getInstance();
  int from = rnode->getID();

  int child_nodeID;
  int success = this->atlas->addNode(child_CG, child_nodeID, from);
  AtlasNode* child_node = (*this->atlas)[child_nodeID];

	bool resumeSampling = false;
	if(success == 0) {  //existed before : may be done or incomplete
		delete child_CG; child_CG = NULL;
		child_CG = child_node->getCG(); //you need existing CG to get parameter etc. information
		if(!child_node->isComplete() ) //incomplete  //&& !settings::AtlasBuilding::breadthFirst
			resumeSampling = true;
	}

  wit_orr_toSend->addBoundary(from);

  // If resumeSampling==true, it does not mean they are connected, maybe
  // child_node is created by another node than 'rnode' (parent node).
  bool alreadyConnected = this->atlas->isConnected(from, child_nodeID);
  this->atlas->connect(from, child_nodeID);

  /*
  ==================================================
    success resumeSam  bread  convexifyable   action
         0     0       0       0          Add bunch of witnesses
         0     0       0       1          Continue
         0     0       1       0          Add bunch of witnesses
         0     0       1       1          Continue
         0     1       0       0          Not possible
         0     1       0       1          Sample
         0     1       1       0          Not possible
         0     1       1       1          Do not sample; fix noGoodOrientation
         1     0       0       0          Add bunch of witnesses
         1     0       0       1          Sample
         1     0       1       0          Add bunch of witnesses
         1     0       1       1          Do not sample; fix noGoodOrientation
         1     1       0       0          Not possible
         1     1       0       1          Not possible
         1     1       1       0          Not possible
         1     1       1       1          Not possible
  ==================================================
  */
  int case_int = 0;
  case_int = case_int | success;
  case_int = case_int << 1 | resumeSampling;
  case_int = case_int << 1 | sett->AtlasBuilding.breadthFirst;
  bool partial_3_tree = false;
  // At this point child_node should have non-null CG;
  child_node->isPartial3Tree(partial_3_tree);
  case_int = case_int << 1 | partial_3_tree;
  
  // copy of wit_orr_toSend is added to wp4
  CayleyPoint* wp4 = new CayleyPoint(wit_orr_toSend,
                          child_CG->getParamLines(), a, b);
  
  // Set correct flip in orientation wp4->getOrientations()[0].
  // TODO: This code is copied from processAtlasNode(), refactor
  // this code to an appropriate function / class method.
  CayleyParameterization* cParam = new CayleyParameterization(child_CG, false);
  ConvexChart* chart = new ConvexChart(child_CG, /*denseSampling=*/ false, cParam, this->df);
  if (child_node->getFlipScheme().empty()) {
    if (chart->partial3tree) {
      child_node->setFlipScheme(chart->getTetras());
    } else {
      std::vector<std::vector<int>> parent_flipScheme =
          this->atlas->getNodes()[child_node->getFirstParentID()]
              ->getFlipScheme();

      child_node->setFlipScheme(parent_flipScheme);
    }
  }
  

  refigureFlip(wp4->getOrientations()[0], child_node->getFlipScheme());
  const int& child_flip = wp4->getOrientations()[0]->getFlipNum();
  const auto& [parent_id, entry_point_id, entry_point_flip] = 
                                        wit_orr_toSend->getEntryPoint();

//  if (!wit_orr_toSend->angleViolated()) {
  switch(case_int) {
    // Either a convexifiable child_node has just been created or the sampling is being resumed.
    case 5: 
    case 9: if (!savedOnce) {
              // Save parent sample points before going into recursion.
              this->snl->appendSpacePoints(rnode); 
              noPoints = 0;
              savedOnce = true;
            }
            if (sett->Paths.implementPathFinding) {
               rnode->addToEntryPointTable(child_nodeID, child_flip, entry_point_id,
                                         entry_point_flip);
            }

            wp4->setID(child_node->getACR()->nextWitnessPointID());
    		if (sett->Output.writeNodeFiles) this->snl->appendDimension(child_node);
            this->snl->appendWitness(child_node->getID(), wp4);
					  child_node->setFoundGoodOrientation(true);
            
            // Main Recursive Call.
            if (this->processAtlasNode(child_node, false, wit_orr_toSend, resumeSampling,
                                 bret)) {
              // If there is a realization with good angle at the childs, then set
              // noGoodOrientation to false.
              noGoodOrientation = false;
            }
            break;
    // Cases where we add one witness per (child_flip, entry_point_flip) from every parent.
    case 1:
    case 3: if (sett->Paths.implementPathFinding) {
               // If (child_flip, entry_point_flip) is unique, then add the witness/entry points.
               if (!rnode->entryPointExists(child_nodeID, child_flip, entry_point_flip)) { 
					         rnode->addToEntryPointTable(child_nodeID, child_flip, entry_point_id,
                                             entry_point_flip);
                   wp4->setID(child_node->getACR()->nextWitnessPointID());
				   this->snl->appendWitness(child_node->getID(), wp4);
				   child_node->setFoundGoodOrientation(true);
               }
            } else if (sett->Sampling.JacobianSampling) {
              // TODO: Add a lookup of covered child flips. If child_flip is not covered 
              // by any parent then add this witness point.
              LOG(INFO) << "Unimplemented feature";
            }

            break;
    // Cases where we use breadth first search...        
    // TODO: Implement breadth first cases.
    case 7:
    case 11: 
             // rnode->insertChildNodeWitnessSet(child_nodeID
             LOG(INFO) << "Feature unimplemented";
             break;

    // Cases where the node is non-convexifiable, if settings->saveBoundary is
    // set, add all witnesses from the first parent and one witness per flip
    // from every other parent. 
    case 0:
    case 2:
    case 8:
    case 10: { // Ray tracing.
                 bool ray_tracing_add = false;
                 if (sett->Saving.saveBoundary && 
                         child_node->getFirstParentID() == from) { 
                     ray_tracing_add = true;
                 }
                 // TODO: Implement the case for non first-parent, to make sure all
                 // child flips are covered by some parent. 

                 // Find if rnode is first parent.
                 bool path_finding_add = false;
                 if (sett->Paths.implementPathFinding) {
                    path_finding_add = 
                        !rnode->entryPointExists(child_nodeID, child_flip, entry_point_flip);
                 } 
                 if (ray_tracing_add || path_finding_add) {
                    if (sett->Paths.implementPathFinding) {
                        rnode->addToEntryPointTable(child_nodeID, child_flip, entry_point_id,
                                         entry_point_flip);
                    }
                    wp4->setID(child_node->getACR()->nextWitnessPointID());
                    this->snl->appendWitness(child_node->getID(), wp4);
                    child_node->setFoundGoodOrientation(true);
                 }
             }
             break;
    // Cases below should NOT occur.
    case 4:
    case 6:
    case 12:
    case 13:
    case 14:
    case 15: LOG(WARNING) << "Unexpected condition: " 
                          << "Success: " << success
                          << "ResumeSampling: " << resumeSampling
                          << "Breadth First: " << sett->AtlasBuilding.breadthFirst
                          << "Convexifiable: " << partial_3_tree;
  } // End of switch-case.
//  } // End of angleVoilated if.

  wp4->trim_PointMultiD();
  delete wp4;


  if (child_node->hasAnyGoodOrientation()) {  // if the child has any good
                                              // orientation, then parent
                                              // should be displayed.
    noGoodOrientation = false;
  }
  return child_node;
}

void AtlasBuilder::addReverseWitness(AtlasNode* rnode,
                                     ActiveConstraintGraph* acg,
                                     ConvexChart* chart,
                                     ActiveConstraintRegion* region,
                                     Orientation* coming_witness_ori,
                                     bool& noGoodOrientation) {
  Settings* sett = Settings::getInstance();
  if (sett->General.reverseWitness) {
    vector<CgMarker::Edge> edges = this->cgmarker.getEdge(acg);
    cout << "[mark]" << *acg << endl;
    cout << "[mark] getting witness " << edges.size() << endl;
    for (int i = 0; i < edges.size(); i++) {
      Orientation* ornt = edges[i].second;  // DO NOT FORGET TO DELETE ORNT
      CayleyPoint* wp4 = new CayleyPoint(ornt, acg->getParamLines(), a, b);
      chart->setWitnessPoint(wp4);
      if (!coming_witness_ori->angleViolated()) {
        region->AddWitnessPoint(wp4);
        //				acg->AddWitnessPoint(wp4); //  add it to the
        //graph, it will be needed as a starting point in the description class
        // this->snl->appendWitness(rnode->getID(), wp4);// save it too the
        // graph
        noGoodOrientation = false;
      }
    }
    // if acg is a 0-d and currently not in atlas
    // a lot of debugging needed here - Yichi
    if (acg->getDim() == 0 && !rnode->isMarked) {
      this->cgmarker.mark(rnode);  // acg
      rnode->isMarked = true;
      cout << "[mark] mark 0" << endl;
    }
  }
}

void AtlasBuilder::findAndSampleMissingAncestors() {
  vector<AtlasNode*> nodes = this->atlas->getNodes();
  /**
   *
   * load all 0-d nodes with their space points and mark their ancestor
   * the space points are needed to serve as the starting point of the
   * re-sampling
   * TODO: actually only 1 space point is needed, is there a way to avoid
   * loading the whole node?
   */
  /*
  for(vector<AtlasNode*>::iterator iter = nodes.begin();	iter !=
  nodes.end();	iter++){ if((*iter)->getDim() == 0){
                  snl->loadNode((*iter)->getID(), (*iter)->getACR() );
                  // ActiveConstraintGraph* cg = (*iter)->getCG();
                  cgmarker.mark((*iter));
          }
  }
  */

  // for debug - Yichi
  int resampleCount = 0;

  while (!cgmarker.empty()) {
    pair<ActiveConstraintGraph*,
         vector<pair<ActiveConstraintGraph*, Orientation*>>>
        res = cgmarker.pop();
    ActiveConstraintGraph* cg = res.first;
    vector<pair<ActiveConstraintGraph*, Orientation*>>& edges = res.second;

    int nodenum = this->atlas->findNodeNum(cg);
    cout << *(cg) << endl;
    if (nodenum == -1) {
      delete cg;
      continue;
    }

    // for debug - Yichi
    cgmarker.report_size();

    cout << "refining " << nodenum << endl;

    // check each child, if it is not in the roadmap or the connection doesn't
    // exist, redo sample
    for (int i = 0; i < edges.size(); i++) {
      cout << endl << i << "th child" << endl;
      cout << *(res.second[i].first) << endl;
      if (atlas->findNodeNum(res.second[i].first) == -1 ||
          !atlas->isConnected(nodenum,
                              atlas->findNodeNum(res.second[i].first))) {
        AtlasNode* rnod = (*this->atlas)[nodenum];
        snl->loadNode(nodenum, rnod->getACR());

        Orientation* orr = res.second[i].second;

        // debug output, print the Orientation
        if (true) {
          cout << endl;
          cout << "o";

          vector<int> boundary = orr->getBoundary();
          int flip = orr->getFlipNum();

          double fb[3][3], tb[3][3];
          orr->getFromTo(fb, tb);
          for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++) {
              cout << " " << fb[i][j];
            }
          for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++) {
              cout << " " << tb[i][j];
            }
          cout << " " << boundary.size();
          for (int i = 0; i < boundary.size(); i++) cout << " " << boundary[i];
          cout << " " << flip;
          cout << endl;
        }

        // re-sample the node with the new starting Orientation
        processAtlasNode(rnod, false, edges[i].second, false, false);
        resampleCount++;

        // for debug - Yichi
        cgmarker.report_size();

        break;
      }
    }

    // for debug - Yichi
    cout << "Resampling count:" << resampleCount << endl;

    delete cg;
  }
}

list<Orientation*> AtlasBuilder::findRealizations(
    ActiveConstraintGraph* acg, ConvexChart* des,
    vector<vector<int>> flipScheme) {
  list<Orientation*> output;
  num_samples++;

  if (des->partial3tree)
    for (int i = 0; i < 8; i++)  // for each flip
    {
      bool fail;
      // Orientation *relz = CartesianRealizer::computeRealization(acg, des, i,
      // fail, flipScheme);
      Orientation* relz = des->computeRealization(i, fail, flipScheme);
      if (!fail) {
        // We do not do angle check here i.e. if( !relz->angleViolated() )
        // Because after binary search around this orientation with bad angle,
        // there may have orientations with good angle Also child nodes can
        // cause small angles so if this realization cause a new contact, then
        // create the child not and keep it as witness at the child node. After
        // that delete this realization from parent(this node), and make it
        // orange.
        output.push_back(relz);

      } else {  // volume negative

        delete relz;
        relz = NULL;
        break;  // in the next flips it will always get volume negative too
      }
    }

  /********************** commented by BRIAN, seems unused anyways???
#ifdef USE_MATLAB
  //RUIJIN
  if(!des->partial3tree && settings::General::runSolver){
  bool fail;
  num_samples++;
  Orientation *relz = CartesianRealizer::computeRealization(acg, des, 0,
fail);//the first realization created by first constructor, the rest is by
second constructor if( !fail ){	//means maple is able to return in specified
time msol = relz->mapslns; helix_base* solver = relz->solver;
  vector<vector<double> > rts = relz->rts;

  if(!settings::Constraint::checkAngle || relz->angleViolated() ) //small anglee
  output.push_back( relz );
  else{
  //delete relz;
  output.push_back( relz );
  }


  for(int i=1; i<msol  ; i++)
  {
  CartesianRealizer *relz = new CartesianRealizer(des,  solver, rts, i);
  if(!settings::Constraint::checkAngle || relz->angleViolated() ) //small anglee
  output.push_back( relz );
  else{
  //		delete relz;
  output.push_back( relz );
  }
  }

  }
  else
  delete relz;
  }
#endif


   *********************/

  return output;
}
