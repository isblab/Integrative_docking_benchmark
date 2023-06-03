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
 * Thread_BackEnd.cpp
 *
 *  Created on: May 20, 2016
 *      Author: rprabhu
 */

#include "Thread_BackEnd.h"

#include <iostream>
#include <queue>
#include <string>
#include <vector>

#include "ActiveConstraintGraph.h"
#include "ActiveConstraintRegion.h"
#include "AtlasNode.h"
#include "Settings.h"
#include "ThreadShare.h"

using namespace std;

int samplingNodeNumber;

bool Starting() {
  struct stat buffer;
  string node0(ThreadShare::save_loader->relativePath + "Node0.txt");
  //    string
  //    settingsini(ThreadShare::save_loader->relativePath+"settings.ini");
  string settingsini("settings.ini");
  if ((stat(node0.c_str(), &buffer) == 0) &&
      (stat(settingsini.c_str(), &buffer) == 0)) {
    return false;
  } else {
    return true;
  }
}

void setSamplingNodeNumber(int node_num) { samplingNodeNumber = node_num; }

void BackEnd_Thread(int nodeNumber) {
  Settings* sett = Settings::getInstance();
  samplingNodeNumber = nodeNumber;
  while (true) {
    switch (sett->getRunMode()) {
      case ForestSampleAS:
        if (Starting()) {
          ThreadShare::atlas_builder->setup();
          ThreadShare::atlas_builder->startAtlasBuilding();
          ThreadShare::save_loader->saveRoadMap(ThreadShare::atlas_view);
        } else {
          list<size_t> undone = ThreadShare::atlas_view->getUnfinished();
          while (!undone.empty()) {
            AtlasNode* rnod =
                ThreadShare::save_loader->loadAtlasNode(undone.front());
            ActiveConstraintGraph* cont = rnod->getCG();
            ActiveConstraintRegion* region = rnod->getACR();
            ThreadShare::save_loader->loadNode(undone.front(), region);
            if (!rnod->isComplete()) {
              if (region->GetSamplePointsCount() == 0) {
                ThreadShare::atlas_builder->processAtlasNode(rnod, false, NULL,
                                                             false, false);
              } else {
                ThreadShare::atlas_builder->processAtlasNode(rnod, false, NULL,
                                                             true, false);
              }
            }
            undone.pop_front();
            cout << "\nUndone size = " << undone.size() << endl;
          }
          ThreadShare::atlas_builder->startAtlasBuilding();
          ThreadShare::save_loader->saveRoadMap(ThreadShare::atlas_view);
          sett->setRunMode(Stopped);
        }
        break;
      case ForestSample: {
        if (Starting()) continue;
        list<size_t> undone = ThreadShare::atlas_view->getUnfinished();
        while (!undone.empty()) {
          AtlasNode* rnod =
              ThreadShare::save_loader->loadAtlasNode(undone.front());
          ActiveConstraintGraph* cont = rnod->getCG();
          ActiveConstraintRegion* region = rnod->getACR();
          ThreadShare::save_loader->loadNode(undone.front(), region);
          if (!rnod->isComplete()) {
            if (region->GetSamplePointsCount() == 0) {
              ThreadShare::atlas_builder->processAtlasNode(rnod, false, NULL,
                                                           false, false);
            } else {
              ThreadShare::atlas_builder->processAtlasNode(rnod, false, NULL,
                                                           true, false);
            }
          }
          undone.pop_front();
          cout << "\nThread_NodeFinisher: Undone size = " << undone.size()
               << endl;
        }
        ThreadShare::save_loader->saveRoadMap(ThreadShare::atlas_view);
        sett->setRunMode(Stopped);
      } break;
      case TreeSample: {
        if (Starting()) continue;
        ThreadShare::save_loader->saveRoadMap(ThreadShare::atlas_view);
        AtlasNode* atlas_node = (*ThreadShare::atlas_view)[samplingNodeNumber];
        if (atlas_node->isComplete()) {
          sett->setRunMode(Stopped);
          break;
        }

        ThreadShare::save_loader->loadNode(samplingNodeNumber,
                                           atlas_node->getACR());
        ThreadShare::atlas_builder->processAtlasNode(atlas_node, false, NULL,
                                                     false, false);
        ThreadShare::save_loader->saveRoadMap(ThreadShare::atlas_view);
        sett->setRunMode(Stopped);
      } break;
      case BreadthFirst: {
        if (Starting()) {
          continue;
        }
        std::pair<std::set<int>::iterator, bool> ret;  // to check duplicates
        std::set<int> myset;
        myset.insert(samplingNodeNumber);

        std::queue<int> myqueue;
        myqueue.push(samplingNodeNumber);
        sett->AtlasBuilding.breadthFirst = true;

        while (!myqueue.empty() && !sett->AtlasBuilding.stop) {
          int nom = myqueue.front();
          myqueue.pop();
          AtlasNode* rnode = (*ThreadShare::atlas_view)[nom];
          ThreadShare::save_loader->loadNode(
              nom, rnode->getACR());  // for test purpose
          if (!rnode->isComplete()) {
            bool cont = false;
            ThreadShare::atlas_builder->processAtlasNode(rnode, false, NULL,
                                                         cont, false);
            cont = true;

            vector<int> con = rnode->getConnection();  // con = the node numbers
                                                       // connected to this one

            size_t dim = rnode->getDim();
            for (vector<int>::iterator it = con.begin(); it != con.end();
                 it++) {
              AtlasNode* child_node = (*ThreadShare::atlas_view)[*it];
              int child_dim = child_node->getDim();
              if (child_dim < dim) {  // if it is a lesser dimension to make
                                      // sure it is child not parent
                ret = myset.insert(*it);
                if (ret.second)  // new element did not exist before
                  myqueue.push(*it);
              }
            }
          }
        }
        ThreadShare::save_loader->saveRoadMap(ThreadShare::atlas_view);
        cout << "Thread_Breadth: Done sampling." << endl;
        sett->setRunMode(Stopped);
      } break;
      case RefineSampling: {
        if (Starting()) continue;
        ThreadShare::save_loader->saveRoadMap(ThreadShare::atlas_view);
        AtlasNode* atlas_node = (*ThreadShare::atlas_view)[samplingNodeNumber];

        ThreadShare::save_loader->loadNode(samplingNodeNumber,
                                           atlas_node->getACR());
        ThreadShare::atlas_builder->processAtlasNode(atlas_node, true, NULL,
                                                     false, false);
        ThreadShare::save_loader->saveRoadMap(ThreadShare::atlas_view);
        sett->setRunMode(Stopped);
      } break;
      case Stopped:
        continue;
        break;
      default:
        cout << "There is a BUG in the code. New unauthorized mode introduced"
             << endl;
    }
  }
}
