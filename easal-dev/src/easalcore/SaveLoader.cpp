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
 * SaveLoader.cpp
 *
 *  Created on: Apr 12, 2009
 *      Author: James Pence
 */

#include "ActiveConstraintGraph.h"
#include "CayleyPoint.h"
#include "MolecularUnit.h"
#include "Orientation.h"
#include "SaveLoader.h"
#include "Settings.h"

#ifdef SQLITE
#include "SQLiteHelper.h"
#endif

#include <time.h>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>



#include <glog/logging.h>
#include <glog/stl_logging.h>

#ifdef SERVER
#include <bsoncxx/builder/basic/array.hpp>
#include <bsoncxx/builder/basic/document.hpp>
#include <bsoncxx/builder/basic/kvp.hpp>
#include <bsoncxx/builder/stream/array.hpp>
#include <bsoncxx/builder/stream/document.hpp>
#include <bsoncxx/builder/stream/helpers.hpp>
#include <bsoncxx/types.hpp>
#include <memory>

#include "json/json.h"
#include "json/value.h"

using namespace bsoncxx;

using bsoncxx::builder::basic::kvp;
Json::Value RoadMap(Json::arrayValue);
int numNodes = 0;

#endif

using namespace std;

time_t SaveLoader::last_saved_time = time(NULL);

SaveLoader::SaveLoader() {}

void SaveLoader::setData(std::string relativePosition, MolecularUnit* a,
                         MolecularUnit* b) {
  this->relativePath = relativePosition;
  this->a = a;
  this->b = b;
}

void SaveLoader::setSessionID(std::string session) {
  this->sessionID = session;
}

#ifdef SERVER
SaveLoader::SaveLoader(std::string relativePosition, MolecularUnit* a,
                       MolecularUnit* b, mongocxx::client* mongoDBC) {
  if (relativePosition.at(relativePosition.length() - 1) != '/') {
    relativePosition.append("/");
  }

  this->relativePath = relativePosition;
#ifdef WIN32
  mkdir(this->relativePath.c_str());
#else
  mkdir(this->relativePath.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
#endif
  this->a = a;
  this->b = b;

  this->mongoDBClient = mongoDBC;
  sessionID = "Run1";
}
#else

SaveLoader::SaveLoader(string relativePosition, MolecularUnit* a,
                       MolecularUnit* b) {
  if (relativePosition.at(relativePosition.length() - 1) != '/') {
    relativePosition.append("/");
  }

  this->relativePath = relativePosition;
#ifdef WIN32
  mkdir(this->relativePath.c_str());
#else
  mkdir(this->relativePath.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
#endif
  this->a = a;
  this->b = b;
#ifdef SERVER
  std::string database_filename = "test.db";
  DeleteDatabase(database_filename);
  CreateDatabase(database_filename);
#endif
}

#endif

SaveLoader::~SaveLoader() {}

void SaveLoader::saveNodeHeader(AtlasNode* node, Atlas* atlas) {
  ofstream outFile;

#ifdef SERVER

  stringstream name;
  stringstream JSONname;

  name << this->relativePath << "RoadMap.txt";

  JSONname << this->relativePath << "RoadMap.json";

  name.flush();
  JSONname.flush();
  outFile.open((name.str()).c_str(), ios_base::app | ios_base::out);

  // writeNodeHeader(node, outFile);
  writeNodetoJSON(node, JSONname.str());

#else
  stringstream name;
  name << this->relativePath << "RoadMap.txt";
  name.flush();
  outFile.open((name.str()).c_str(), ios_base::app | ios_base::out);

  // writeNodeHeader(node, outFile);
#endif

  // mapViewInUse->unlock();

  time_t current_time = time(NULL);
  int time_elapsed = difftime(current_time, last_saved_time);
  if (time_elapsed > 5) {  // 5 seconds
    last_saved_time = current_time;
    saveRoadMap(atlas);
  }
}

#ifdef SERVER
void SaveLoader::writeNodetoJSON(AtlasNode* node, std::string JSONfile) {
  Json::StreamWriterBuilder swb;
  std::unique_ptr<Json::StreamWriter> writer(swb.newStreamWriter());

  std::fstream f;

  f.open(JSONfile.c_str(), ios_base::trunc | ios_base::out);

  ActiveConstraintGraph* cgk = node->getCG();

  // An individual RoadMap Entry - One for every Node in the Atlas
  Json::Value RoadMapEntry;

  // A object to hold pairs of contacts/parameters.
  Json::Value pairs(Json::objectValue);

  // An array to hold contacts
  Json::Value contacts(Json::arrayValue);

  // An array to hold parameters
  Json::Value params(Json::arrayValue);

  // An array to hold locations
  Json::Value loc(Json::arrayValue);

  // An array to hold connected Nodes
  Json::Value connections(Json::arrayValue);

  // Iterate over the contacts and append them to the contacts JSON array
  vector<pair<int, int>> contact = cgk->getParticipants();

  int i = 0;
  for (vector<pair<int, int>>::iterator iter = contact.begin();
       iter != contact.end(); iter++) {
    pairs["one"] = iter->first;
    pairs["two"] = iter->second;
    contacts[i] = pairs;
    i++;
  }

  // Iterate over the parameters and append them to the params JSON array
  vector<pair<int, int>> par = cgk->getParamLines();
  i = 0;
  for (vector<pair<int, int>>::iterator iter = par.begin(); iter != par.end();
       iter++) {
    pairs["one"] = iter->first;
    pairs["two"] = iter->second;
    params[i] = pairs;
    i++;
  }

  // Iterate over the connections and append them to the connections JSON array
  vector<int> con = node->getConnection();
  i = 0;
  for (vector<int>::iterator it = con.begin(); it != con.end(); it++) {
    connections[i] = (*it);
    i++;
  }

  // Populate the RoadMapEntry
  RoadMapEntry["NodeID"] = node->getID();
  RoadMapEntry["FirstParentID"] = node->getID();
  RoadMapEntry["Complete"] = node->isComplete();
  RoadMapEntry["NonEmpty"] = node->hasAnyGoodOrientation();
  RoadMapEntry["Contacts"] = contacts;
  RoadMapEntry["Parameters"] = params;
  RoadMapEntry["StepSize"] = cgk->getStepSize();
  // RoadMapEntry["Location"] = loc;
  RoadMapEntry["ConnectedNodes"] = connections;

  // Append it to the RoadMap
  RoadMap.append(RoadMapEntry);
  numNodes++;

  // f << RoadMapEntry << std::endl;
  writer->write(RoadMap, &f);
  f << endl;

  f.close();

  mongocxx::database db = (*mongoDBClient)["easal"];
  mongocxx::collection session = db[sessionID];

  auto bsonBuilder = bsoncxx::builder::stream::document{};

  Json::FastWriter fastWriter;
  std::string output = fastWriter.write(RoadMap);

  using bsoncxx::builder::stream::close_document;
  using bsoncxx::builder::stream::document;
  using bsoncxx::builder::stream::finalize;
  using bsoncxx::builder::stream::open_document;

  bsoncxx::stdx::optional<bsoncxx::document::value> maybe_result =
      session.find_one(document{} << "RoadMap" << 1 << finalize);
  if (maybe_result) {
    session.update_one(document{} << "RoadMap" << 1 << finalize,
                       document{} << "$set" << open_document << "Nodes"
                                  << output << close_document << finalize);

  } else {
    bsoncxx::document::value samplingInfo =
        bsonBuilder << "RoadMap" << 1 << "Nodes" << output
                    << bsoncxx::builder::stream::finalize;

    bsoncxx::document::view view = samplingInfo.view();

    bsoncxx::stdx::optional<mongocxx::result::insert_one> result =
        session.insert_one(view);
  }
}
#endif

void SaveLoader::writeNodeHeader(AtlasNode* node, ostream& outFile) {
  ActiveConstraintGraph* cgk = node->getCG();

  outFile << "# " << node->getID() << " " << node->getFirstParentID() << " "
          << node->isComplete() << " " << node->hasAnyGoodOrientation() << " ";

  vector<pair<int, int>> contacts = cgk->getParticipants();
  outFile << contacts.size() << " ";
  for (vector<pair<int, int>>::iterator iter = contacts.begin();
       iter != contacts.end(); iter++)
    outFile << "c " << iter->first << " " << iter->second << " ";

  vector<pair<int, int>> par = cgk->getParamLines();
  outFile << par.size() << " ";
  for (vector<pair<int, int>>::iterator iter = par.begin(); iter != par.end();
       iter++)
    outFile << "p " << iter->first << " " << iter->second << " ";

  outFile << cgk->getStepSize() << " ";
  // outFile << node->getLocation()[0] << " " << node->getLocation()[1] << " "
  // << node->getLocation()[2] << " ";
  vector<int> con = node->getConnection();
  outFile << con.size();
  for (vector<int>::iterator it = con.begin(); it != con.end(); it++) {
    outFile << " " << *it;
  }
  outFile << "\n";  // endl;
}

void SaveLoader::saveRecentPointsToFile(AtlasNode* rnode) {
  stringstream name, nameW;
  ofstream outFile, outFileW;
  name << this->relativePath << "node" << rnode->getID() << ".txt";
  nameW << this->relativePath << "node" << rnode->getID() << "_witness.txt";
  name.flush();
  nameW.flush();
#ifdef server
  stringstream jsonname, jsonnamew;
  jsonname << this->relativepath << "node" << rnode->getID() << ".json";
  jsonnamew << this->relativepath << "node" << rnode->getID()
            << "_witness.json";
  jsonname.flush();
  jsonnamew.flush();
  ofstream outjsonfile, outjsonfileW;
  outjsonfile.open((jsonname.str()).c_str(), ios_base::app | ios_base::out);
  outjsonfilew.open((jsonnamew.str()).c_str(), ios_base::app | ios_base::out);
#endif

  outFile.open((name.str()).c_str(), ios_base::app | ios_base::out);
  outFileW.open((nameW.str()).c_str(), ios_base::app | ios_base::out);

  vector<CayleyPoint*> witspc = rnode->getACR()->GetWitnessPoints();
  for (vector<CayleyPoint*>::iterator iter = witspc.begin();
       iter != witspc.end(); iter++) {
    writeCayleyPoint(*iter, outFileW, true);
#ifdef SERVER
    // Write Cayley Points to JSON file here.
    writeCayleyPointtoJSON(*iter, JSONnameW.str(), rnode->getID(), true);
#endif
  }
  vector<CayleyPoint*> spc = rnode->getACR()->GetSamplePoints();
  for (vector<CayleyPoint*>::iterator iter = spc.begin(); iter != spc.end();
       iter++) {
    writeCayleyPoint(*iter, outFile);
#ifdef SERVER
    // Write Cayley Points to JSON file here.
    writeCayleyPointtoJSON(*iter, JSONname.str(), rnode->getID());
#endif
  }

  outFile.flush();
  outFile.close();
  outFileW.flush();
  outFileW.close();

#ifdef SERVER
  outJSONFile.flush();
  outJSONFile.close();
  outJSONFileW.flush();
  outJSONFileW.close();
#endif
  // TODO: Now that all the points have been written to JSON files, just read
  // the JSON file and do one mongodb update query.
  rnode->getACR()->trim();
}

void SaveLoader::appendSpacePoints(AtlasNode* rnode) {
  //	cout << "."; // a progress indicator on the console.
  cout.flush();
  ofstream outFile;
  stringstream name;
  name << this->relativePath << "Node" << rnode->getID() << ".txt";
  name.flush();

#ifdef SERVER
  ofstream outJSONFile;
  stringstream JSONname;
  JSONname << this->relativePath << "Node" << rnode->getID() << ".json";
  JSONname.flush();
  outJSONFile.open((JSONname.str()).c_str(), ios_base::app | ios_base::out);
#endif

  outFile.open((name.str()).c_str(), ios_base::app | ios_base::out);
  vector<CayleyPoint*> space = rnode->getACR()->GetSamplePoints();

  for (size_t i = 0; i < space.size(); i++) {
    if (space[i] == NULL) continue;
    this->writeCayleyPoint(space[i], outFile);
#ifdef SERVER
    this->writeCayleyPointtoJSON(space[i], JSONname.str(), rnode->getID());
#endif
  }
  outFile.flush();
  outFile.close();
#ifdef SERVER
  outJSONFile.flush();
  outJSONFile.close();
#endif
  rnode->getACR()->trim();
}

#ifdef CAF
void SaveLoader::appendSpacePoints(int nodeID,
                                   std::vector<CayleyPointStruct> vecCps) {
  //	cout << "."; // a progress indicator on the console.
  cout.flush();
  ofstream outFile;
  stringstream name;
  name << this->relativePath << "Node" << nodeID << ".txt";
  name.flush();

#ifdef SERVER
  ofstream outJSONFile;
  stringstream JSONname;
  JSONname << this->relativePath << "Node" << nodeID << ".json";
  JSONname.flush();
  outJSONFile.open((JSONname.str()).c_str(), ios_base::app | ios_base::out);
#endif

  outFile.open((name.str()).c_str(), ios_base::app | ios_base::out);

  for (size_t i = 0; i < vecCps.size(); i++) {
    CayleyPoint* cPoint = new CayleyPoint(vecCps[i]);
    this->writeCayleyPoint(cPoint, outFile);
#ifdef SERVER
    this->writeCayleyPointtoJSON(cPoint, JSONname.str(), nodeID);
#endif
    delete cPoint;
  }
  vecCps.erase(vecCps.begin(), vecCps.end());
  outFile.flush();
  outFile.close();
#ifdef SERVER
  outJSONFile.flush();
  outJSONFile.close();
#endif
}

#endif
void SaveLoader::appendWitness(int nodeID, CayleyPoint* wit) {
  ofstream outFile;
  ofstream outJSONFile;

  stringstream name;
  name << this->relativePath << "Node" << nodeID << "_Witness.txt";
  name.flush();
  outFile.open((name.str()).c_str(), ios_base::app | ios_base::out);

  this->writeCayleyPoint(wit, outFile, true);
  outFile.flush();
  outFile.close();

#ifdef SERVER
  stringstream JSONname;
  JSONname << this->relativePath << "Node" << nodeID << "_Witness.json";
  JSONname.flush();
  outJSONFile.open((JSONname.str()).c_str(), ios_base::app | ios_base::out);
  this->writeCayleyPointtoJSON(wit, JSONname.str(), nodeID, true);
  outJSONFile.flush();
  outJSONFile.close();
#endif
}


void SaveLoader::appendDimension(AtlasNode* rnode) {
  Settings* sett = Settings::getInstance();
  if (!rnode->dimWrittenSample) {
    ofstream outFile;
    stringstream name;
    name << this->relativePath << "Node" << rnode->getID() << ".txt";
    name.flush();
    outFile.open((name.str()).c_str(), ios_base::app | ios_base::out);
    outFile << "d " << rnode->getParamDim() << "\n";
    outFile.flush();
    outFile.close();
    rnode->dimWrittenSample = true;
  }

  if (!rnode->dimWrittenWitness) {
    ofstream outFileW;
    stringstream nameW;
    nameW << this->relativePath << "Node" << rnode->getID() << "_Witness.txt";
    nameW.flush();
    outFileW.open((nameW.str()).c_str(), ios_base::app | ios_base::out);
    outFileW << "d " << rnode->getParamDim() << "\n";
    outFileW.flush();
    outFileW.close();
    rnode->dimWrittenWitness = true;
  }

  if(sett->Paths.implementPathFinding && 
  		rnode->getDim() <= sett->Paths.energyLevelUpperBound) {
  	if (!rnode->dimWrittenEPTable) {
      ofstream outFileW;
      stringstream nameW;
      nameW << this->relativePath << "Node" << rnode->getID() << "_EPTable.txt";
      nameW.flush();
      outFileW.open((nameW.str()).c_str(), ios_base::app | ios_base::out);
      outFileW << "d " << rnode->getParamDim() << "\n";
      outFileW << "child_id child_flip parent_flip entry_point_id" << endl;
      outFileW.close();
      rnode->dimWrittenEPTable = true;
    }
  }

}

void SaveLoader::appendDimension(int nodeID, int dim) {
  ofstream outFile;
  stringstream name;
  name << this->relativePath << "Node" << nodeID << ".txt";
  name.flush();
  outFile.open((name.str()).c_str(), ios_base::app | ios_base::out);
  outFile << "d " << dim << "\n";
  outFile.flush();
  outFile.close();

  ofstream outFileW;
  stringstream nameW;
  nameW << this->relativePath << "Node" << nodeID << "_Witness.txt";
  nameW.flush();
  outFileW.open((nameW.str()).c_str(), ios_base::app | ios_base::out);
  outFileW << "d " << dim << "\n";
  outFileW.flush();
  outFileW.close();
}

void SaveLoader::appendVolume(AtlasNode* rnode) {
  double min[6], max[6];
  rnode->getACR()->getSpaceVolume(min, max);
  ofstream outFile;
  stringstream name;
  name << this->relativePath << "Node" << rnode->getID() << ".txt";
  name.flush();
  outFile.open((name.str()).c_str(), ios_base::app | ios_base::out);
  outFile << "v ";
  for (size_t i = 0; i < 6; i++) {
    outFile << min[i] << " ";
  }
  for (size_t i = 0; i < 6; i++) {
    outFile << max[i];
    if (i < 5) {
      outFile << " ";
    }
  }
  outFile << "\n";  // endl;
  outFile.flush();
  outFile.close();
}

void SaveLoader::loadNode(int num, ActiveConstraintRegion* output) {
  ifstream inFile, inFileW;
  stringstream name, nameW;
  name << this->relativePath << "Node" << num << ".txt";
  nameW << this->relativePath << "Node" << num << "_Witness.txt";
  name.flush();
  nameW.flush();

  output->trim();
  // cout << "SaveLoader::loadNode: Reading node " << num << " from "<<
  // name.str() << endl;
  try {
    inFile.open((name.str()).c_str(), ios_base::in);
    inFileW.open((nameW.str()).c_str(), ios_base::in);
    if (inFile.is_open() && inFileW.is_open()) {
      loadActiveConstraintRegion(inFile, name.str(), output);
      loadActiveConstraintRegion(inFileW, nameW.str(), output);
	  output->setSamplePointsCreated(output->GetSamplePointsCount());
      output->setWitnessPointsCreated(output->GetWitnessPointsCount());
	  
    } else {
      cout << "SaveLoader::loadNode: Failed to open " << name.str() << endl;
    }
    inFile.close();
    inFileW.close();
  } catch (std::string e) {
    cerr << e;
    inFile.close();
    inFileW.close();
  } catch (exception& e) {
    cerr << e.what();
    inFile.close();
    inFileW.close();
  }

}

void SaveLoader::saveRoadMap(Atlas* atlas) {
  // int tlock = mapViewInUse->try_lock();  //0 if able to lock
  // if( tlock == 0 ){
  cout.flush();
  ofstream outFile;
  stringstream name;
  name << this->relativePath << "RoadMap.txt";
  name.flush();
  for (size_t i = 3; i > 0; i--) {
    stringstream oldername, newername;
    oldername << this->relativePath << "Roadmap" << i << ".txt";
    oldername.flush();
    newername << this->relativePath << "Roadmap" << (i - 1) << ".txt";
    newername.flush();
    if (i > 1) {
      remove(oldername.str().c_str());
      rename(newername.str().c_str(), oldername.str().c_str());
    } else {
      remove(oldername.str().c_str());
      rename(name.str().c_str(), oldername.str().c_str());
    }
  }

  outFile.open(name.str().c_str());
  if (!outFile.is_open())
    cout << "SaveLoader::saveMapView: Failed to create " << name.str() << endl;
  outFile
      << "/NodeID firstParentID Complete nonEmpty ContactSize \"contacts # #\" "
         " ParamDimension \"parameters # #\"  StepSize \"Location x y z\" "
         "NumberOfConnections \"Nodes this node is connected to # # # ...\" "
      << "\n";  // endl;

  size_t size = atlas->number_of_nodes();
  for (size_t iter = 0; iter < size;
       iter++)  // since roadmap is extending during saving process by caller
                // thread,  it is possible that there will be connection nodes,
                // where the node itself is not saved yet.
  {
    AtlasNode* node = (*atlas)[iter];
    writeNodeHeader(node, outFile);
  }

  outFile.flush();
  outFile.close();

  // mapViewInUse->unlock();
  //}
  // else{
  //   cout << "saveRoadMap lock not achieved" << endl;
  //    }
}

void SaveLoader::saveAtlas(Atlas* atlas) {
  size_t size = atlas->number_of_nodes();
  for (size_t num = 0; num < size; num++) {
    AtlasNode* node = (*atlas)[num];
    appendDimension(node);  // on 21 Feb 13

    if (node->getACR()->GetSpaceSize() > 0) saveRecentPointsToFile(node);
  }
}

void SaveLoader::loadMapView(Atlas* output) {
  bool debug = false;
  if (debug) cout << "loadMapView: Begin." << endl;

  // mapViewInUse->lock();

  if (output != NULL) {
    if (output->number_of_nodes() > 0) {
      output->cleanAtlas();
      // cerr << "deleted nodes ";
    }
  } else
    output = new Atlas();

  ifstream inFile;
  stringstream name;
  vector<AtlasNode*>* nodes = new vector<AtlasNode*>();
  std::string trash;

  name << this->relativePath << "RoadMap.txt";
  name.flush();
  inFile.open(name.str().c_str());
  inFile >> trash;

  if (debug) cout << "loadMapView: Opening file " << name.str() << endl;

  while ((trash[0] == '#' || trash[0] == '/') && !inFile.eof()) {
    if (debug) cout << "Trash: " << trash << endl;

    if (trash[0] == '#') {
      // cout << 1 << endl;
      AtlasNode* tmp = loadNextAtlasNode(inFile);
      nodes->push_back(tmp);
      // cout << 2 << endl;
    } else {
      std::string trsh;
      getline(inFile, trsh);
    }
    inFile >> trash;  // read #
  }
  inFile.close();

  output->setNodes(*nodes);

  if (debug) cout << "loadMapView: Unlocking." << endl;

  // mapViewInUse->unlock();
  //    return output;
}

AtlasNode* SaveLoader::loadAtlasNode(int nodenum) {
  AtlasNode* tmp;
  // mapViewInUse->lock();
  ifstream inFile;
  stringstream name;
  std::string trash;
  int id;
  bool found = false;
  name << this->relativePath << "RoadMap.txt";
  name.flush();
  inFile.open((name.str()).c_str());
  if (inFile.is_open()) {
    inFile >> trash;

    while ((trash[0] == '#' || trash[0] == '/') && !inFile.eof()) {
      if (trash[0] == '#') {
        inFile >> id;
        if (id == nodenum) {
          inFile.unget();  // to make loadAtlasNode able to read the id again
          tmp = loadNextAtlasNode(inFile);
          found = true;
          break;
        } else
          getline(inFile, trash);
      } else
        getline(inFile, trash);

      inFile >> trash;  // read #
    }
    inFile.close();
  } else {
    cout << "SaveLoader::loadAtlasNode: Could not open roadmap file "
         << name.str() << endl;
  }

  // mapViewInUse->unlock();

  if (!found) tmp = new AtlasNode();
  return tmp;
}

AtlasNode* SaveLoader::loadNextAtlasNode(ifstream& inFile) {
  std::string trash;
  int id;
  int pid;
  bool complete;
  bool empty;
  int contactSize;
  vector<pair<int, int>> partic;
  int paramDim;
  vector<pair<int, int>> params;

  double stepsize;
  // double* location = new double[3];
  int numberOfConnections;
  vector<int> con;

  inFile >> id >> pid >> complete >> empty >> contactSize;

  for (int i = 0; i < contactSize; i++) {
    inFile >> trash;  // read 'c'
    pair<int, int> contact;
    inFile >> contact.first >> contact.second;
    partic.push_back(contact);
  }

  inFile >> paramDim;
  for (int i = 0; i < paramDim; i++) {
    inFile >> trash;  // read 'p'
    pair<int, int> par;
    inFile >> par.first >> par.second;
    params.push_back(par);
  }

  inFile >> stepsize;
  // inFile >> location[0] >> location[1] >> location[2];

  inFile >> numberOfConnections;
  for (int it = 0; it < numberOfConnections; it++) {
    int val;
    inFile >> val;
    con.push_back(val);
  }

  AtlasNode* tmp =
      new AtlasNode(id, pid, complete, empty, 6 - contactSize, con);

  ActiveConstraintGraph* output =
      new ActiveConstraintGraph(partic, params);  // 16SEP13

  if (empty)
    tmp->setFoundGoodOrientation(true);
  else
    tmp->setFoundGoodOrientation(false);

  output->setStepSize(stepsize);

  tmp->setCG(
      output);  // necessary since i need contact information to see if a new
                // cgk is sampled before or not (need to go through whole cgks)

  Settings* sett = Settings::getInstance();
  if(sett->Paths.implementPathFinding && paramDim <= sett->Paths.energyLevelUpperBound ) {
	 loadEntryPointsTable(tmp);
  }

  return tmp;
}

void SaveLoader::loadActiveConstraintRegion(istream& file,
                                            std::string filename,
                                            ActiveConstraintRegion* output) {
  time_t Start_t = time(NULL);

  int oris_in_sparseregion = 0, totaloris = 0;
  bool debug = false;
  std::string element;
  CayleyPoint* currentPoint = NULL;
  CayleyPoint* lastValidPoint;
  CayleyPoint* prevPoint = NULL;
  vector<CayleyPoint*> witnesses;
  int paramDim;
  size_t line = 0;
  //int samplePointsCreated = 0, witnessPointsCreated = 0;
  try {
    while (!file.eof()) {
      line++;

      if (file.fail()) {  // can happen if the previous file reading operation
                          // "file >> something" cause error such as mismatch of
                          // types
        throw std::string("previous");
      }

      file >> element;
      if (file.eof()) break;

      if (file.fail()) {
        throw std::string("read element");
      }

      switch (element[0]) {
        case 'd':  // param dimension
        {
          if (debug) cout << "d";
          file >> paramDim;
        } break;

        case '*':  // Point
        {
          if (debug) cout << "*";

          currentPoint = readCayleyPoint(file, paramDim);
          if (currentPoint != NULL) {
            output->AddSamplePoint(currentPoint);
            lastValidPoint = currentPoint;
            //samplePointsCreated++;
          }

        } break;

        case 'w':  // witness point
        {
          currentPoint =
              NULL;  // set it null otherwise when there file read error, and it
                     // crashes, you are gonna delete previous valid point
          if (debug) cout << "w";
          try {
            currentPoint = readCayleyPoint(file, paramDim);
          } catch (std::string err) {
            file.clear();
            std::string trash;
            getline(file, trash);
            cerr << "\n[" << line << "]Bad witness " << err << "\n";
            delete currentPoint;
            currentPoint = NULL;
            break;
          }
          if (currentPoint != NULL) {
            //							output->AddWitnessPoint(currentPoint);
            witnesses.push_back(currentPoint);
            lastValidPoint = currentPoint;
            //witnessPointsCreated++;
          }

        } break;

        case 'o':  // orientation
        {
          totaloris++;
          if (debug) cout << "o";
          Orientation* temp = NULL;
          try {
            temp = readOrientation(file);
          } catch (std::string err) {
            file.clear();
            cerr << "\n[" << line << "]Bad orientation " << err << "\n";
            file.ignore(numeric_limits<streamsize>::max(), '\n');
            if (temp != NULL) {
              delete temp;
              temp = NULL;
            }
            break;
          }
          if ((temp != NULL) && (currentPoint != NULL)) {
            currentPoint->addOrientation(temp);
          } else {
            delete temp;
            temp = NULL;
            cerr << "\n[" << line
                 << "] orientation deleted because of bad point \n";
          }
        } break;
        case '/': {
          if (debug) cout << "/";
          file.ignore(numeric_limits<streamsize>::max(), '\n');
        } break;

        // for old type node files
        case 'f':
        case 'r':
        case 'b':
        case 'c':
        case 'p':
        case 't':
        case 's':
        case 'v':
        case '-':
          file.ignore(numeric_limits<streamsize>::max(), '\n');
          break;

        default:
          cout << "\nline = \"" << line << " element " << element << "\""
               << "\n";  // endl;
          cout.flush();
      }
      if (debug) cout.flush();
    }

    // added by Brian
    /*if (currentPoint != NULL) {
      output->setSamplePointsCreated(samplePointsCreated);
      output->setWitnessPointsCreated(witnessPointsCreated);
    } else {
      // output->setPointsCreated(lastValidPoint->getID());
      // cout << endl << "ERRRORRRR line 705 SaveLoader.cpp" << endl ;
    }*/

    for (size_t i = 0; i < witnesses.size(); i++) {
      output->AddWitnessPoint(witnesses[i]);
    }
  } catch (std::string err) {
    delete currentPoint;
    cerr << "\n" << filename << ":" << line << "file error - " << err << endl;
  }

  // cout << "totaloris " <<  totaloris << endl;

  time_t t1 = time(NULL);
  int time_1 = difftime(t1, Start_t);
  // cout << "time_1 " << time_1 << endl;
}

CayleyPoint* SaveLoader::readCayleyPoint(istream& file, size_t dim) {
  CayleyPoint* output = NULL;
  vector<double> pos;
  bool t;
  int badAngleN, collideN, ID;

  double tmp;
  for (size_t i = 0; i < dim; i++) {
    file >> tmp;
    if (file.fail()) throw std::string("read pointpos");
    pos.push_back(tmp);
  }

  file >> t;

  int axis = 0;  // for jacobian sampling
  file >> axis;

  std::string temp;
  file >> temp;
  if (temp[0] == 'o')  // somehow dimension is 1 but that parameter is not
                       // written to the file !!!
  {
    badAngleN = t;
    t = tmp;
    file.unget();
  } else
    badAngleN = atoi(temp.c_str());

  file >> collideN;

  for (size_t i = 0; i < pos.size(); i++) {
    if (std::isnan(pos[i])) {
      throw(std::string("nan value in Pos"));
    }
  }
  file >> ID;

  // Read tetra info.
  std::string str_value;
  std::array<double, 3> tetra_volume;
  std::array<std::array<double, 6>, 3> tetra_edges;
  for (int i=0; i<3; i++) {
    // Read edge lenths.
    for (int j=0; j<6; j++) {
      file >> str_value;
      tetra_edges[i][j] = atof(str_value.c_str());
    }
    // Read volume.
    file >> str_value;
    tetra_volume[i] = atof(str_value.c_str());
  }



  output = new CayleyPoint(pos);

  //	if(file.fail()){
  //		//throw string("read point-tet");   //I COMMENDED BECAUSE COULDNOT
  //FIND THE ERROR
  //	}
  output->setRealizable(t);
  output->setBadAngleN(badAngleN);
  output->setCollidN(collideN);
  output->zIndex = axis;
  output->setID(ID);
  //output->SetIsTetraVolumeZero(is_tetra_volume_zero);
  output->SetTetraVolume(tetra_volume);
  output->SetTetraEdges(tetra_edges);

  return output;
}

Orientation* SaveLoader::readOrientation(istream& file) {
  double fb[3][3], tb[3][3];
  std::string strValue;
  bool good = true;
  //	try{

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) {
      file >> strValue;
      if (strValue == "nan") {
        //				throw string("read orient fromB");
        good = false;
      } else {
        fb[i][j] = atof(strValue.c_str());
      }
    }

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) {
      file >> strValue;
      if (strValue == "nan") {
        //				throw string("read orient toB");
        good = false;
      } else {
        tb[i][j] = atof(strValue.c_str());
      }
    }

  //	}
  //	catch(ifstream::failure &e){
  if (file.fail()) {
    throw std::string("read orient fromto");
  }
  int boundsize, flip;
  vector<int> boundary;
  //	try{
  file >> boundsize;
  int bndry;
  for (int i = 0; i < boundsize; i++) {
    file >> bndry;
    boundary.push_back(bndry);
  }
  file >> flip;

  bool hasEntryPoint;
  file >> hasEntryPoint;

  int nodeID, cayleyPointID, flipNo;

  if (hasEntryPoint) {
    file >> nodeID;
    file >> cayleyPointID;
    file >> flipNo;
  }

  //	}
  //	catch(ifstream::failure &e){
  if (file.fail()) {
    throw std::string("read orient boundflip");
  }
  //	for(int i = 0;i<3;i++)for(int j = 0;j<3;j++){
  //		if( isnan(fa[i][j]) || isnan(fb[i][j]) || isnan(ta[i][j]) ||
  //isnan(tb[i][j]) ){ 			throw(string("nan value in fromto"));
  //		}
  //	}

  Orientation* output = NULL;
  if (good) {
    output = new Orientation(fb, tb);
    output->setBoundary(boundary);
    output->setFlipNum(flip);
    if (hasEntryPoint) output->setEntryPoint(nodeID, cayleyPointID, flipNo);
  }

  return output;
}

#ifdef SERVER
/*
void SaveLoader(std::string file){
    ifstream ifile(file);
    Json::Value CayleyPoints;
    Json::CharReaderBuilder builder;
    std::string errs;

    if(ifile.good()) {
        bool ok = Json::parseFromStream(builder, ifile, &CayleyPoints, &errs);
    }

    ifile.close();

mongocxx::database db = (*mongoDBClient)["easal"];
  mongocxx::collection session = db[sessionID];

  auto bsonBuilder = bsoncxx::builder::stream::document{};

  using bsoncxx::builder::stream::document;
  using bsoncxx::builder::stream::finalize;
  using bsoncxx::builder::stream::open_document;
  using bsoncxx::builder::stream::close_document;


  bsoncxx::stdx::optional<bsoncxx::document::value> maybe_result =
          session.find_one(document{} << "Node" << nodeID << finalize);
  if(maybe_result) {
          session.update_one(document{} << "Node" << nodeID << finalize,
                          document{} << "$set" << open_document <<
                          "Sampling" << output << close_document << finalize);
  } else {
          bsoncxx::document::value samplingInfo= bsonBuilder
                  << "Node" << nodeID
                  << "Sampling" << output
                  << bsoncxx::builder::stream::finalize;

          bsoncxx::document::view view = samplingInfo.view();

          bsoncxx::stdx::optional<mongocxx::result::insert_one> result =
session.insert_one(view);
}
*/

void SaveLoader::writeCayleyPointtoJSON(CayleyPoint* pnt, std::string file,
                                        int nodeID, bool witness) {
  ifstream ifile(file);
  Json::Value CayleyPoints;
  Json::CharReaderBuilder builder;
  std::string errs;

  if (ifile.good()) {
    bool ok = Json::parseFromStream(builder, ifile, &CayleyPoints, &errs);
  }

  ifile.close();

  Json::StreamWriterBuilder swb;
  std::unique_ptr<Json::StreamWriter> writer(swb.newStreamWriter());

  // A Cayley Point
  Json::Value CayleyPointJSON;

  // An array to hold Cayley Point
  Json::Value CayleyParameterValues(Json::arrayValue);

  // An array to hold Orientations
  Json::Value Orientations(Json::arrayValue);

  // A object to hold individual orientation
  Json::Value orientation(Json::objectValue);

  // An array to hold all boundaries of a Cayley point from all orientations.
  Json::Value cayleyPointBoundaryJSON(Json::arrayValue);

  // An array to hold EntryPoints
  Json::Value entryPointArrayJSON(Json::arrayValue);

  vector<int> cayleyPointBoundary;

  for (int i = 0; i < pnt->dim(); i++) {
    CayleyParameterValues[i] = ((*pnt)[i]);
  }

  vector<Orientation*> ornt = pnt->getOrientations();
  int orients = 0;
  for (vector<Orientation*>::iterator iter = ornt.begin(); iter != ornt.end();
       iter++) {
    // An array to hold from
    Json::Value fbJSON(Json::arrayValue);

    // An array to hold to
    Json::Value tbJSON(Json::arrayValue);

    // An array to hold boundary
    Json::Value boundaryJSON(Json::arrayValue);

    double fb[3][3], tb[3][3];
    (*iter)->getFromTo(fb, tb);
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        Json::Value val = fb[i][j];
        fbJSON.append(val);
      }
    }

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        Json::Value val = tb[i][j];
        tbJSON.append(val);
      }
    }

    vector<tuple<int, CayleyPoint*, Orientation*>> ep;
    //	= (*iter)->getEntryPoints();

    for (int i = 0; i < ep.size(); i++) {
      Json::Value entryPointJSON;
      Json::Value parentid = std::get<0>(ep[i]);
      Json::Value epid = std::get<1>(ep[i])->getID();
      Json::Value flipNumber = (int)std::get<2>(ep[i])->getFlipNum();
      entryPointJSON["ParentID"] = parentid;
      entryPointJSON["EntryPointID"] = epid;
      entryPointJSON["FlipNumber"] = flipNumber;

      entryPointArrayJSON.append(entryPointJSON);
    }

    vector<int> boundary = (*iter)->getBoundary();
    for (int i = 0; i < boundary.size(); i++) {
      boundaryJSON[i] = (boundary[i]);
      cayleyPointBoundary.push_back(boundary[i]);
    }

    orientation["fb"] = fbJSON;
    orientation["tb"] = tbJSON;
    orientation["boundaries"] = boundaryJSON;
    orientation["flip"] = (int)(*iter)->getFlipNum();
  }

  std::sort(cayleyPointBoundary.begin(), cayleyPointBoundary.end());
  std::vector<int>::iterator it;
  it = std::unique(cayleyPointBoundary.begin(), cayleyPointBoundary.end());
  cayleyPointBoundary.resize(std::distance(cayleyPointBoundary.begin(), it));

  for (std::vector<int>::iterator iter = cayleyPointBoundary.begin();
       iter != cayleyPointBoundary.end(); iter++) {
    Json::Value bound = *iter;
    cayleyPointBoundaryJSON.append(bound);
  }

  // Populate the CayleyPoint
  CayleyPointJSON["NodeID"] = nodeID;
  CayleyPointJSON["Witness"] = witness;
  CayleyPointJSON["CayleyPoint"] = CayleyParameterValues;
  CayleyPointJSON["Realizable"] = pnt->isRealizable();
  CayleyPointJSON["Z-index"] = pnt->zIndex;
  CayleyPointJSON["Orientations"] = Orientations;
  CayleyPointJSON["BadAngle"] = pnt->getBadAngleN();
  CayleyPointJSON["Collision"] = pnt->getCollidN();
  CayleyPointJSON["Boundaries"] = cayleyPointBoundaryJSON;
  CayleyPointJSON["EntryPoints"] = entryPointArrayJSON;

  CayleyPoints.append(CayleyPointJSON);

  ofstream ofile(file, std::ofstream::out | std::ofstream::trunc);

  Json::FastWriter fastWriter;
  std::string output = fastWriter.write(CayleyPoints);

  writer->write(CayleyPoints, &ofile);
  ofile.close();

  mongocxx::database db = (*mongoDBClient)["easal"];
  mongocxx::collection session = db[sessionID];

  auto bsonBuilder = bsoncxx::builder::stream::document{};

  using bsoncxx::builder::stream::close_document;
  using bsoncxx::builder::stream::document;
  using bsoncxx::builder::stream::finalize;
  using bsoncxx::builder::stream::open_document;

  bsoncxx::stdx::optional<bsoncxx::document::value> maybe_result =
      session.find_one(document{} << "Node" << nodeID << finalize);
  if (maybe_result) {
    session.update_one(document{} << "Node" << nodeID << finalize,
                       document{} << "$set" << open_document << "Sampling"
                                  << output << close_document << finalize);
  } else {
    bsoncxx::document::value samplingInfo =
        bsonBuilder << "Node" << nodeID << "Sampling" << output
                    << bsoncxx::builder::stream::finalize;

    bsoncxx::document::view view = samplingInfo.view();

    bsoncxx::stdx::optional<mongocxx::result::insert_one> result =
        session.insert_one(view);
  }
}

#endif

void SaveLoader::writeCayleyPoint(CayleyPoint* pnt, ostream& file,
                                  bool witness) {
  if (witness) {
    file << "w";
  } else {
    file << "*";
  }
  size_t dim;

  for (size_t i = 0; i < pnt->dim(); i++) {
    file << " " << (*pnt)[i];
  }

  file << " " << pnt->isRealizable();

  file << " " << pnt->zIndex;  // for jacobian sampling

  // file << " " << pnt->getBadAngle(); //if its volume is negative, then it
  // will not look forward and check for grid constraint and so keep angle to be
  // true!!!
  file << " " << pnt->getBadAngleN();
  file << " " << pnt->getCollidN();
  file << " " << pnt->getID();
  
  // Write tetra info.
  const auto& tetra_volume = pnt->TetraVolume();
  const auto& tetra_edges = pnt->TetraEdges();
  for (int i=0; i<3; i++) {
    // Write edge lenths.
    for (int j=0; j<6; j++) {
      file << " " << tetra_edges[i][j];
    }
    file << " " << tetra_volume[i];
  }

/*
  // Helpful in debugging GetTetrahedralBoundaryFlips 
  vector<vector<int>> boundary_flips = pnt->GetTetrahedralBoundaryFlips();
  
  file << " | ";
  for (const auto& partition : boundary_flips) {
    for (const auto& flip : partition) {
      file << " " << flip;
    }
    file << " | ";
  }
*/

  file << "\n";  // endl;
  vector<Orientation*> ornt = pnt->getOrientations();
  for (vector<Orientation*>::iterator iter = ornt.begin(); iter != ornt.end();
       iter++) {
    writeOrientation(*iter, file);

    //		delete *iter;
  }
  //	delete pnt;
}

void SaveLoader::writeOrientation(Orientation* ornt, ostream& file) {
  int flip = ornt->getFlipNum();

  file << "o";

  vector<int> boundary = ornt->getBoundary();

  double fb[3][3], tb[3][3];
  ornt->getFromTo(fb, tb);

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) {
      file << " " << fb[i][j];
    }

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) {
      file << " " << tb[i][j];
    }
  file << " " << boundary.size();
  for (int i = 0; i < boundary.size(); i++) file << " " << boundary[i];
  file << " " << flip;  // the value is incorrect if it witness orientation.

  if (ornt->getHasEntryPoint() == true) {
    file << " " << ornt->getHasEntryPoint();

    std::tuple<int, int, int> entryPoints;
    entryPoints = ornt->getEntryPoint();

    file << " " << std::get<0>(entryPoints) << " " << std::get<1>(entryPoints)
         << " " << std::get<2>(entryPoints);
  } else {
    file << " " << 0;
  }

  file << "\n";  // endl;
}

void SaveLoader::writeOrientationInPDBformat(MolecularUnit* hA,
                                             MolecularUnit* hB,
                                             Orientation* orien,
                                             std::string filename) {
  double fb[3][3], tb[3][3];
  orien->getFromTo(fb, tb);
  vector<Atom*> atms = hA->getAtoms();
  vector<Atom*> btms = hB->getXFAtoms(fb, tb);

  ofstream file;
  file.open(filename.c_str());
  for (size_t i = 0; i < atms.size(); i++) {
    file << atms[i]->getLine() << "\n";  // endl;
    cout << atms[i]->getLine() << "\n";  // endl;
  }
  for (size_t i = 0; i < btms.size(); i++) {
    file << btms[i]->getLine() << "\n";  // endl;
    cout << btms[i]->getLine() << "\n";  // endl;
  }
  file.flush();
  file.close();

  for (size_t i = 0; i < btms.size(); i++) {
    delete btms[i];
  }
}

#ifdef SERVER
mongocxx::client* SaveLoader::getMongoDBClient() { return mongoDBClient; }
#endif


void SaveLoader::writeEntryPointsTable(AtlasNode* rnode) {
  stringstream fileName;
  ofstream outFile;
  fileName << this->relativePath << "Node" << rnode->getID() << "_EPTable.txt";
  fileName.flush();

  outFile.open((fileName.str()).c_str(), ios_base::app | ios_base::out);
  

  const std::unordered_map<std::pair<int, int>, unordered_map<int, int>, AtlasNode::hash_pair>& entryPointTable = rnode->getEntryPointTable();

  for(auto& entryPoints: entryPointTable) {
 	auto& [child_id, child_flip] = entryPoints.first;
	for(auto& entryPoint: entryPoints.second) {
		auto& [parent_flip, entry_point_id] = entryPoint;
   			outFile<<child_id<<" "<<child_flip<<" " 
   					<<parent_flip<<" "<<entry_point_id<<endl;
      }
  }

  outFile.close();

}

void SaveLoader::loadEntryPointsTable(AtlasNode* rnode) {
  ifstream inFile;
  stringstream fileName;
  fileName << this->relativePath << "Node" << rnode->getID()<< "_EPTable.txt";
  fileName.flush();

  try {
    inFile.open((fileName.str()).c_str(), ios_base::in);
    if (inFile.is_open()) {
      std::string trsh;
	  getline(inFile, trsh);
	  getline(inFile, trsh);
	  int child_id, child_flip, entry_point_id, parent_flip;
      while(getline(inFile, trsh)) {
	  	if(trsh == "")
			break;
		std::stringstream ss(trsh);
		ss >> child_id >> child_flip >> parent_flip >> entry_point_id ;
		rnode->addToEntryPointTable(child_id, child_flip, 
			entry_point_id, parent_flip);
	  }
    } else {
      LOG(WARNING) << "SaveLoader::loadEPTable: Failed to open " << fileName.str();
    }
    inFile.close();
  } catch (exception& e) {
    cerr << e.what();
    inFile.close();
  }

}

