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
 * Atlas.cpp
 *
 *  Created on: Feb 22, 2009
 *      Author: Admin
 */

#include "Atlas.h"

#include <glog/logging.h>
#include <time.h>

#include <algorithm>
#include <fstream>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include "AtlasNode.h"
#include "Settings.h"
#include "ThreadShare.h"

using namespace std;

Atlas::Atlas() { this->nodes.reserve(1000); }

Atlas::Atlas(vector<AtlasNode *> nodes) {
  Settings *sett = Settings::getInstance();
  this->nodes = nodes;
  this->nodes.reserve(1000);

  for (size_t i = 0; i < this->nodes.size(); i++)
    if (this->nodes[i]->getDim() ==
        sett->RootNodeCreation.dimension_of_rootNodes)
      this->rootIndices.push_back(i);
}

Atlas::~Atlas() { cleanAtlas(); }

void Atlas::cleanAtlas() {
  for (vector<AtlasNode *>::iterator iter = this->nodes.begin();
       iter != this->nodes.end(); iter++) {
    (*iter)->trimNode();
    delete *iter;
  }

  this->nodes.clear();
  this->rootIndices.clear();
}

void Atlas::setNodeFlipScheme(vector<vector<int> > parentFlipScheme,
                              int parent) {
  nodes[parent]->setFlipScheme(parentFlipScheme);
}

vector<AtlasNode *> Atlas::getNodes() {
  return this
      ->nodes;  // todo if it is time consuming to pass that vector (if it
                // copies the memory addresses), then pass pointer to the vector
}

AtlasNode *Atlas::getNode(int nodeNum) { return this->nodes[nodeNum]; }

void Atlas::setNodes(std::vector<AtlasNode *> nds) {
  Settings *sett = Settings::getInstance();
  cleanAtlas();

  this->nodes = nds;

  for (size_t i = 0; i < this->nodes.size(); i++)
    if (this->nodes[i]->getDim() ==
        sett->RootNodeCreation.dimension_of_rootNodes)
      this->rootIndices.push_back(i);
}

size_t Atlas::number_of_nodes() { return this->nodes.size(); }

list<size_t> Atlas::getUnfinished() {
  list<size_t> output;
  for (size_t i = 0; i < this->nodes.size(); i++) {
    if (this->nodes[i]->isComplete() == 0) {
      output.push_back(i);
    }
  }
  return output;
}

vector<int> Atlas::getCommonParents(int nodeNum1, int nodeNum2) {
  vector<int> commonParents;
  if (this->nodes[nodeNum1]->getDim() != this->nodes[nodeNum2]->getDim()) {
    cout << "Dimensions not the same." << endl;
  }
  vector<int> parents_1 = getParents(nodeNum1);
  vector<int> parents_2 = getParents(nodeNum2);
  vector<int> minParents;
  vector<int> maxParents;

  if (parents_1.size() < parents_2.size()) {
    minParents = parents_1;
    maxParents = parents_2;
  } else {
    minParents = parents_2;
    maxParents = parents_1;
  }
  /*
   cout << "Parent set 1: " << endl;
   for(int i = 0; i < minParents.size(); i++){
   cout << minParents[i] << endl;
   }

   cout << "Parent set 2: " << endl;
   for(int i = 0; i < maxParents.size(); i++){
   cout << maxParents[i] << endl;
   }
   */
  for (int i = 0; i < minParents.size(); i++) {
    for (int j = 0; j < maxParents.size(); j++) {
      if (minParents[i] == maxParents[j]) {
        commonParents.push_back(minParents[i]);
        break;
      }
    }
  }

  /*
   cout << "Shared parents: " << endl;
   for(int i = 0; i < commonParents.size(); i++){
   cout << commonParents[i] << endl;
   }
   cout << endl;
   */

  return commonParents;
}

list<size_t> Atlas::getChildren(size_t nodeNum, int level) {
  cout << "starting with node " << nodeNum << " looking for level " << level
       << endl;
  list<size_t> output, temp;
  int curdim = this->nodes[nodeNum]->getDim();
  if (curdim != level)  // to prevent it to return itself
    output.push_back(nodeNum);
  for (int x = 0; x < (curdim - level) && !output.empty(); x++) {
    temp.insert(temp.end(), output.begin(), output.end());
    output.clear();

    for (list<size_t>::iterator iter = temp.begin(); iter != temp.end();
         iter++) {
      vector<int> con = this->nodes[*iter]->getConnection();
      int dim = this->nodes[*iter]->getDim();

      for (vector<int>::iterator it = con.begin(); it != con.end(); it++) {
        if (this->nodes[*it]->getDim() < dim &&
            this->nodes[*it]->getDim() >= level) {
          output.push_back(*it);
        }
      }
    }
    output.sort();
    output.unique();

    temp.clear();
    //			if(output.empty()){
    //				break;
    //			}
  }

  output.sort();
  output.unique();
  return output;
}

int Atlas::getSymmetricNode(AtlasNode *node) {  // works for 5d node
  Settings *sett = Settings::getInstance();
  vector<pair<int, int> > particp = node->getCG()->getParticipants();
  if (particp[0].second != particp[0].first) {
    vector<pair<int, int> > parts;
    parts.push_back(make_pair(particp[0].second, particp[0].first));
    ActiveConstraintGraph *reverseNode = new ActiveConstraintGraph(parts);

    // int revnode = getNodeNum(reverseNode);
    int revnode = findNodeNum(reverseNode);
    delete reverseNode;
    return revnode;
  } else
    return -1;
}

bool isEqual(ActiveConstraintGraph *acg1, ActiveConstraintGraph *acg2) {
  vector<pair<int, int> > contacts1 = acg1->getParticipants();
  vector<pair<int, int> > contacts2 = acg2->getParticipants();

  if (contacts1.size() != contacts2.size()) {
    return false;
  }

  for (int i = 0; i < contacts1.size(); i++) {
    if (find(contacts2.begin(), contacts2.end(), contacts1[i]) ==
        contacts2.end()) {
      return false;
    }
  }
  return true;
}

string constructSignature(vector<pair<int, int> > part) {
  sort(part.begin(), part.end());

  std::string graphString;
  for (size_t j = 0; j < part.size(); j++) {
    graphString.append(std::to_string(std::get<0>(part[j])));
    graphString.append("-");
    graphString.append(std::to_string(std::get<1>(part[j])));
    if (j < part.size() - 1) graphString.append("=");
  }
  return graphString;
}

int Atlas::findNodeNum(ActiveConstraintGraph *node) {
  int nodeNum;

  VLOG(3) << "ActiveConstraintGraph: " << *node;

  vector<pair<int, int> > part = node->getParticipants();

  std::string graphString = constructSignature(part);

  try {
    nodeNum = nodeMap.at(graphString);
  } catch (exception e) {
    nodeNum = -1;
  }

  VLOG(3) << "Checking if graphString exists: " << graphString << " "
          << (bool)(nodeNum + 1);

  return nodeNum;
}
/*
int Atlas::findNodeNum(ActiveConstraintGraph *node) {
        int nodeNumber = -1;
        for(int i=0; i<this->nodes.size(); i++) {
                if(this->nodes[i]->getCG()->getDim() != node->getDim()) {
                        continue;
                }
                if(isEqual(this->nodes[i]->getCG(), node)) {
                        return this->nodes[i]->getID();
                }
        }
        return nodeNumber;
}*/

// no need lock for this method, because this method either called if existence
// of input node is certain, or called inside of a lock while adding to roadmap
int Atlas::getNodeNum(ActiveConstraintGraph *node) {
  // this code does not work if you created a node with random specific
  // contacts(if it does not have any parents)
  for (size_t i = 0; i < this->rootIndices.size(); i++) {
    int rootIndex = this->rootIndices[i];
    ActiveConstraintGraph *root_cgk = this->nodes[rootIndex]->getCG();
    if (root_cgk->IsParentOf(node)) {
      if (node->getDim() == root_cgk->getDim()) {
        return rootIndex;
      }

      int found = findNodeAtTheChildsOfParent(rootIndex, node);
      if (found != -1) {
        return found;
      }
    }
  }
  return -1;
}

vector<int> Atlas::getParents(int nodenum) {
  vector<int> parents;
  vector<int> connections = this->nodes[nodenum]->getConnection();

  for (int i = 0; i < connections.size(); i++) {
    if (this->nodes[nodenum]->getDim() ==
        this->nodes[connections[i]]->getDim() - 1) {
      parents.push_back(connections[i]);
    }
  }

  return parents;
}

vector<int> Atlas::getChildren(int nodenum) {
  vector<int> children;
  vector<int> connections = this->nodes[nodenum]->getConnection();

  for (int i = 0; i < connections.size(); i++) {
    if (this->nodes[nodenum]->getDim() ==
        this->nodes[connections[i]]->getDim() + 1) {
      children.push_back(connections[i]);
    }
  }

  return children;
}

vector<AtlasNode *> Atlas::getChildren(int nodenum, char something, bool) {
  vector<AtlasNode *> children;
  vector<int> connections = this->nodes[nodenum]->getConnection();

  for (int i = 0; i < connections.size(); i++) {
    if (this->nodes[nodenum]->getDim() ==
        this->nodes[connections[i]]->getDim() + 1) {
      children.push_back(this->nodes[connections[i]]);
    }
  }

  return children;
}

int Atlas::findNodeAtTheChildsOfParent(int parent,
                                       ActiveConstraintGraph *node) {
  if (node == NULL) return -1;
  vector<int> con =
      this->nodes[parent]
          ->getConnection();  // con = the node numbers connected to this one
  for (vector<int>::iterator it = con.begin(); it != con.end(); it++) {
    ActiveConstraintGraph *cgKk =
        this->nodes[*it]->getCG();  // cgKk is an actual pointer to the graph of
                                    // a connected node.
    if (cgKk->getDim() < this->nodes[parent]->getDim() &&
        cgKk->IsParentOf(node)) {
      if (cgKk->getDim() ==
          node->getDim())  // it is good that you do not compare paramdim,
                           // because parameters of node is not set yet at the
                           // time this method is called :)
        return *it;  // found

      int found = findNodeAtTheChildsOfParent(*it, node);
      if (found != -1) return found;
    }
  }
  return -1;
}

// returns the status. Output contains the node Number.
int Atlas::addNode(ActiveConstraintGraph *nodeA, int &output, int fpID) {
  Settings *sett = Settings::getInstance();
  int added;
  output = this->findNodeNum(nodeA);  // check if it existed before.
  if (output < 0) {
    output = this->nodes.size();
    VLOG(3) << "Creating node with id: " << output;
    vector<int> blankConnect;
    AtlasNode *addNode =
        new AtlasNode(output, fpID, false, false, nodeA->getDim(),
                      blankConnect);  // set empty to be false in order to
                                      // display it initially
    addNode->setCG(nodeA);

    this->nodes.push_back(addNode);

    std::string graphString;
    vector<pair<int, int> > cont = nodeA->getParticipants();
    sort(cont.begin(), cont.end());
    for (size_t j = 0; j < cont.size(); j++) {
      graphString.append(std::to_string(std::get<0>(cont[j])));
      graphString.append("-");
      graphString.append(std::to_string(std::get<1>(cont[j])));
      if (j < cont.size() - 1) {
	  	graphString.append("=");
	  }
    }
    VLOG(3) << "Node: " << output << " NodeMap_key: " << graphString;
    nodeMap.insert(make_pair(graphString, output));

    if (nodeA->getDim() ==
        sett->RootNodeCreation
            .dimension_of_rootNodes) { // todo check if it has a parent in the
                                      // atlas, if not put it to the rootindices
                                      // as well?
      this->rootIndices.push_back(output);
	}
	if (sett->Output.writeNodeFiles) {
	  sett->runTimeObjects.save_loader->appendDimension(addNode);
	}

    added = 1;  // succesfully added
  } else {
    added = 0;  // not added, existed before
  }

  //	this->applyForces();

  return added;
}

void Atlas::connect(int indexA, int indexB) {
  if (indexA > -1 && indexB > -1 && indexA != indexB) {
    bool newConnection = !(this->isConnected(indexA, indexB));
    if (newConnection) {
      this->nodes[indexA]->addConnection(indexB);
      this->nodes[indexB]->addConnection(indexA);
    }
  }
}

bool Atlas::isConnected(int nodeA, int nodeB) {
  bool a = this->nodes[nodeA]->isConnectedTo(nodeB);
  bool b = this->nodes[nodeB]->isConnectedTo(nodeA);
  if (a != b) cerr << "\nasysmetric Roadmap connection\n";
  return a;
}

AtlasNode *Atlas::operator[](size_t id) {
  AtlasNode *output = this->nodes[id];
  return output;
}

void Atlas::BuildTree(vector<AtlasVertex *> &Graph) {
  vector<int> adj;
  Settings *sett = Settings::getInstance();
  for (vector<AtlasNode *>::iterator iter = this->nodes.begin();
       iter != this->nodes.end(); iter++) {
    if ((*iter)->getDim() > 1) {
      continue;
    }
    adj.clear();
    vector<int> push = (*iter)->getConnection();
    for (int i = 0; i < push.size(); i++) {
      if (this->nodes[push[i]]->getDim() <= sett->Paths.energyLevelUpperBound) {
        adj.push_back(push[i]);
      }
    }
    AtlasVertex *v1 = new AtlasVertex((*iter)->getID(), adj);
    Graph.push_back(v1);
  }
}

// Helper method to find the index of a node in the newly built graph
int findindex(vector<AtlasVertex *> Graph, int num) {
  for (int i = 0; i < Graph.size(); i++)
    if (Graph[i]->number == num) return i;
  return -1;
}

void Atlas::BuildTree(vector<vector<int> > &matrix,
                      vector<AtlasVertex *> &Graph) {
  for (int i = 0; i < Graph.size(); i++) {
    vector<int> adj = Graph[i]->adjList;
    for (vector<int>::iterator iter = adj.begin(); iter != adj.end(); iter++) {
      int index = findindex(Graph, *iter);
      if (index != -1) matrix[i][findindex(Graph, *iter)] = 1;
    }
  }
}

std::vector<size_t> Atlas::getRootIndicies() { return this->rootIndices; }

int Atlas::findpath(int src, int dst, std::string relativepath) {
  std::queue<int> Q;
  vector<AtlasVertex *> Graph;
  Q.push(src);
  std::ofstream ofs;
  clock_t begin = clock();
  std::string pathfile = relativepath + "paths.txt";
  ofs.open(pathfile, std::ofstream::out | std::ofstream::app);

  BuildTree(Graph);

  // ofs<<"Finding the shortest path between "<<src<<" and "<<dst<<endl;
  if (findindex(Graph, src) == -1 || findindex(Graph, dst) == -1) {
    ofs << "The Source and Destination have to be 0D or 1D nodes\n" << endl;
    return 0;
  }

  while (!Q.empty()) {
    int current = Q.front();
    Q.pop();
    AtlasVertex *now = Graph[findindex(Graph, current)];
    Graph[findindex(Graph, current)]->visited = true;
    for (int i = 0; i < now->adjList.size(); i++) {
      if (Graph[findindex(Graph, now->adjList[i])]->visited != true) {
        Q.push(now->adjList[i]);
        Graph[findindex(Graph, now->adjList[i])]->visited = true;
        Graph[findindex(Graph, now->adjList[i])]->parent = current;
      }
    }
  }

  int cur = dst;
  bool print = true;
  while (cur != src) {
    if (Graph[findindex(Graph, cur)]->parent == -1) {
      clock_t end = clock();
      double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
      return 0;
    }
    if (print == true) {
      ofs << "The path between " << src << " and " << dst << " is:";
      print = false;
    }
    ofs << cur << "--";

    cur = Graph[findindex(Graph, cur)]->parent;
  }
  clock_t end = clock();
  ofs << cur << endl;
  double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
  ofs << "Found the path in " << elapsed_secs << " seconds" << endl;

  ofs.close();
  for (int i = 0; i < Graph.size(); i++) delete Graph[i];
  return 1;
}

void Atlas::findAllPaths(std::string relativepath) {
  vector<int> zeroDnodes;
  for (vector<AtlasNode *>::iterator iter = this->nodes.begin();
       iter != this->nodes.end(); iter++) {
    if ((*iter)->getDim() == 0) {
      zeroDnodes.push_back((*iter)->getID());
    }
  }

  for (int i = 0; i < zeroDnodes.size(); i++) {
    for (int j = i + 1; j < zeroDnodes.size(); j++) {
      findpath(zeroDnodes[i], zeroDnodes[j], relativepath);
    }
  }
}

void mul(vector<vector<int> > &matrix1, vector<vector<int> > &matrix2) {
  int num_nodes = matrix1.size();
  vector<vector<int> > square(num_nodes);
  for (int i = 0; i < num_nodes; i++) square[i].resize(num_nodes);

  for (int i = 0; i < num_nodes; i++) {
    for (int j = 0; j < num_nodes; j++) {
      for (int k = 0; k < num_nodes; k++) {
        square[i][j] += matrix1[i][k] * matrix2[k][j];
      }
    }
  }

  for (int i = 0; i < num_nodes; i++) {
    for (int j = 0; j < num_nodes; j++) {
      matrix1[i][j] = square[i][j];
    }
  }
}

void Atlas::findNumpaths() {
  Settings *sett = Settings::getInstance();
  std::ofstream ofs;
  clock_t begin = clock();
  clock_t end;
  ofs.open("path_matrix.txt", std::ofstream::out | std::ofstream::app);

  if (sett->Paths.pathLength <= 0) {
    ofs << "Path length must be positive" << endl;
    ofs.close();
    return;
  }

  vector<AtlasVertex *> Graph;
  BuildTree(Graph);

  int num_nodes = Graph.size();
  int n = sett->Paths.pathLength;

  vector<vector<int> > matrix(num_nodes);
  for (int i = 0; i < num_nodes; i++) matrix[i].resize(num_nodes);

  vector<vector<int> > identity(num_nodes);
  for (int i = 0; i < num_nodes; i++) identity[i].resize(num_nodes);

  for (int i = 0; i < num_nodes; i++) {
    for (int j = 0; j < num_nodes; j++) {
      if (i == j) {
        identity[i][j] = 1;
      } else {
        identity[i][j] = 0;
      }
    }
  }

  BuildTree(matrix, Graph);

  if (sett->Paths.pathLength == 1) {
    end = clock();
  } else {
    while (n > 0) {
      if (n % 2 == 0) {
        mul(matrix, matrix);
        n /= 2;
      } else {
        mul(identity, matrix);
        n--;
      }
    }
    end = clock();
  }

  ofs << "The path matrix is" << endl;
  for (int i = 0; i < num_nodes; i++) {
    for (int j = 0; j < num_nodes; j++) {
      ofs << matrix[i][j] << " ";
    }
    ofs << endl;
  }

  double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
  ofs << "Found number of paths in " << elapsed_secs << " seconds" << endl;
  ofs.close();
}

// TODO: Make sure it's not used anywhere and delete.
vector<pair<CayleyPoint *, int> *> Atlas::preProcessEntryPoints(
    int nodeChoice) {
  vector<pair<CayleyPoint *, int> *> allEntryPoints;
  if (nodes[nodeChoice]->getAllEntryPoints().empty()) {
    vector<AtlasNode *> children = this->getChildren(nodeChoice, 'a', true);

    vector<int> flips;

    // TODO:Should be made a pointer
    pair<CayleyPoint *, int> entryP;

    for (int i = 0; i < children.size(); i++) {
      flips = this->findAllFlips(children[i]->getID());

      for (int j = 0; j < flips.size(); j++) {
        bool epStatus =
            nodes[nodeChoice]->findEntryPoint(children[i], flips[j], entryP);

        if (epStatus == true && entryP.second >= 0) {
          /*
                              cout << endl << "0DNode: "
                                   << children[i]->getID() << endl <<"0D flip "
                                   << flips[j] << " to 1D flip " <<
             entryP->second << endl;

                              entryP->first->printData();

                              cout  << endl;

          */

          allEntryPoints.push_back(&entryP);

        } else
          cout << endl
               << "findEntryPoint returned a nullptr w/ nodeID: "
               << children[i]->getID() << " and flip " << flips[j] << endl;
      }
    }
  }
  return allEntryPoints;
}



/*
std::vector<std::pair<CayleyPoint *, int> *> Atlas::testContinuousPath(int src,
                                                                       int dst,
                                                                       int f1,
                                                                       int f2) {
  cout << "src: " << src << endl << "dst: " << dst << endl;

  std::vector<std::pair<CayleyPoint *, int> *> eventpath;

  cout << "Called testContinuousPath" << endl << endl;

  if (!((this->nodes[src]->getDim() == 0) &&
        (this->nodes[dst]->getDim() == 0))) {
    cout << "Source and Destination should be 0D nodes" << endl;
    return eventpath;
  }

  vector<int> parents = getCommonParents(src, dst);
  if (parents.size() < 1) {
    cout << "No common parents." << endl;
    return eventpath;
  }

  int nodeChoice;

  if (parents.size() == 1) {
    cout << "Searching through only common parent " << parents[0] << endl;
    nodeChoice = parents[0];
  }

  else if (parents.size() > 1) {
    cout << "In which common parent do you want to search?" << endl;

    for (int z = 0; z < parents.size(); z++) {
      cout << parents[z] << endl;
    }

    cin >> nodeChoice;
  }

  ActiveConstraintRegion *acr = new ActiveConstraintRegion();
  ThreadShare::save_loader->loadNode(nodeChoice, acr);

  vector<pair<CayleyPoint *, int> *> allEntryPoints =
      preProcessEntryPoints(nodeChoice);
  // this->nodes[nodeChoice]->setAllEntryPoints(allEntryPoints);
  //    if (nodes[nodeChoice]->getAllEntryPoints().empty()) {

  //        vector<AtlasNode*> children = this->getChildren(nodeChoice, 'a',
  //        true);

  //        vector<int> flips;

  //        pair<CayleyPoint*, int>* entryP;

  //        for (int i = 0; i < children.size(); i++) {

  //            flips = this->findAllFlips(children[i]->getID());

  //            for (int j = 0; j < flips.size(); j++) {

  //                entryP = nodes[nodeChoice]->findEntryPoint(children[i],
  //                flips[j]);

  //                if (entryP != nullptr && entryP->second >= 0) {

  //                    cout << endl << "0DNode: "
  //                         << children[i]->getID() << endl <<"0D flip "
  //                         << flips[j] << " to 1D flip " << entryP->second <<
  //                         endl;

  //                    entryP->first->printData();

  //                    cout  << endl;

  //                    allEntryPoints.push_back(entryP);

  //                } else
  //                    cout << endl
  //                         << "findEntryPoint returned a nullptr w/ nodeID: "
  //                         << children[i]->getID() << " and flip " << flips[j]
  //                            << endl;
  //            }

  //        }

  //        this->nodes[nodeChoice]->setAllEntryPoints(allEntryPoints);

  //    }

  //              AHA!!!!!!!!!!!!
  // vector<pair<CayleyPoint*, int>* > entryPoints =
  // nodes[nodeChoice]->getAllEntryPoints();

  vector<pair<CayleyPoint *, int> *> entryPoints = allEntryPoints;

  cout << "Flips of allEntryPoints: " << endl;

  for (int j = 0; j < entryPoints.size(); j++) {
    cout << entryPoints[j]->second << endl;
  }

  std::pair<CayleyPoint *, int> *ep1;
  std::pair<CayleyPoint *, int> *ep2;

  int choice_1;
  int choice_2;

  cout << "# of entryPoints: " << entryPoints.size() << endl;

  for (int i = 0; i < entryPoints.size(); i++) {
    cout << i << ": on flip " << entryPoints[i]->second << ":  ";

    entryPoints[i]->first->printData();

    cout << endl;
  }

  for (int i = 0; i < entryPoints.size() - 1; i++) {
    for (int j = i + 1; j < entryPoints.size(); j++) {
      ep1 = entryPoints[i];
      ep2 = entryPoints[j];
      eventpath = this->nodes[nodeChoice]->findContinuousPath(
          std::get<0>(*ep1), std::get<0>(*ep2), std::get<1>(*ep1),
          std::get<1>(*ep2), this);

      if (!eventpath.empty()) {
        cout << endl
             << "path found between entry points " << i << " and " << j << endl;
      }
    }
  }

  /****
   cin >> choice_1;
   cin >> choice_2;

   ep1 = nodes[nodeChoice]->allEntryPoints[choice_1];
   ep2 = nodes[nodeChoice]->allEntryPoints[choice_2];

   if (ep1 == NULL) {
   cout << "Entry point 1 doesn't exist. DEBUG!!!" << endl;
   return eventpath;
   }

   if (ep2 == NULL) {
   cout << "Entry point 1 doesn't exist. DEBUG!!!" << endl;
   return eventpath;
   }

   eventpath = this->nodes[nodeChoice]->findContinuousPath(
   std::get < 0 > (*ep1), std::get < 0 > (*ep2), std::get < 1 > (*ep1),
   std::get < 1 > (*ep2), this);

   ***************/

  /**************
   for(int i = 0; i < nodes[nodeChoice]->allEntryPoints.size(); i++){
   for(int j = 0; j < nodes[nodeChoice]->allEntryPoints.size(); j++){
   if(i < j){
   ep1 = nodes[nodeChoice]->allEntryPoints[i];
   ep2 = nodes[nodeChoice]->allEntryPoints[j];
   eventpath = this->nodes[nodeChoice]->findContinuousPath(ep1->first,
   ep2->first, ep1->second,  ep2->second, this); if(eventpath.empty()){ cout <<
   "No path between entry points " << i << " and " << j << endl;
   }
   else
   cout << "Path found between entry points "  << i << " and " << j << endl;
   }
   }
   }
   ******************

  nodes[nodeChoice]->printEventForest();
  cout << "End testContinuousPath" << endl << endl;

  /*****
   if(nodeChoice == 4){

   MolecularUnit* MuA = nodes[nodeChoice]->getCG()->getMolecularUnitA();
   MolecularUnit* MuB = nodes[nodeChoice]->getCG()->getMolecularUnitB();
   Orientation* ori_1 =
   nodes[nodeChoice]->getAllEntryPoints()[1]->first->getOrientation(6);
   Orientation* ori_2 =
   nodes[nodeChoice]->getAllEntryPoints()[2]->first->getOrientation(1);

   double dist = nodes[nodeChoice]->findL2Distance(MuA, MuB, ori_1, ori_2);


   }

   *****/

  /*******
  cout << "begin testBifurcation" << endl;
  for (int i = 0; i < 8; i++) {
      nodes[nodeChoice]->testBifurcation(i);
  }
  cout << "end testBifurcation" << endl;

  **************

  // nodes[nodeChoice]->printFlip(1);
  // nodes[nodeChoice]->printFlip(5);

  return eventpath;
}*/

/*
 * For a 0D node, this function finds a flip number that has a valid
 * realization.
 */
vector<int> Atlas::findFlip(int nodenum) {
  vector<int> flips;
  if (nodes[nodenum]->getDim() > 0) {
    return flips;
  }

  ActiveConstraintRegion *acr = new ActiveConstraintRegion();
  ThreadShare::save_loader->loadNode(nodenum, acr);

  if (acr == NULL) {
    delete acr;
    return flips;
  }

  vector<CayleyPoint *> space = acr->getSpace();

  if (space.size() == 0) {
    delete acr;
    return flips;
  }

  vector<Orientation *> orient = space[0]->getOrientations();

  if (orient.size() == 0) {
    delete acr;
    return flips;
  }

  delete acr;

  for (int i = 0; i < orient.size(); i++) {
    flips.push_back(orient[i]->getFlipNum());
  }
  return flips;
}

vector<int> Atlas::findAllFlips(int nodenum) {
  vector<int> flips;
  if (nodes[nodenum]->getDim() > 0) {
    return flips;
  }

  ActiveConstraintRegion *acr = new ActiveConstraintRegion();
  ThreadShare::save_loader->loadNode(nodenum, acr);

  if (acr == NULL) {
    delete acr;
    return flips;
  }

  vector<CayleyPoint *> space = acr->getSpace();

  if (space.size() == 0) {
    delete acr;
    return flips;
  }
  for (int i = 0; i < space.size(); i++) {
    vector<Orientation *> orient = space[i]->getOrientations();
    for (int j = 0; j < orient.size(); j++) {
      flips.push_back(orient[j]->getFlipNum());
    }
  }

  sort(flips.begin(), flips.end());

  auto last = unique(flips.begin(), flips.end());

  flips.erase(last, flips.end());

  delete acr;
  return flips;
}
/*
std::vector<std::pair<CayleyPoint*, int> *>
Atlas::find1DPath(int src, int dst) {
  std::vector<int> srcFlips;
  this->nodes[src]->getWitnessFlips(srcFlips);
  std::vector<int> dstFlips;
  this->nodes[dst]->getWitnessFlips(dstFlips);

  std::vector<std::pair<CayleyPoint*, int> *> ap;
  for(int i=0; i<srcFlips.size(); i++){
    for(int j=0; j<dstFlips.size(); j++) {
       ap = find1DoFPath(src, dst, srcFlips[i], dstFlips[j]);
      if(ap.size() != 0){
        std::cout << "Found a path for flips "<<i<<" and "<<j << std::endl;
      }
    }
  }
  return ap;
}
*/

void Atlas::report_size() {
  cout << "Node num:" << nodes.size() << endl;
  cout << "Nodemap size:" << nodeMap.size() << endl;
}
