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

#include "SamplerActor.h"

#include <glog/logging.h>

#include <chrono>

#include "CartesianRealizer.h"
#include "CayleyParameterization.h"
#include "ConstraintCheck.h"
#include "ConvexChart.h"
#include "Settings.h"

using namespace std;
using namespace std::chrono;
#define cout aout(self)

string constructSignature(vector<pair<int, int> > part);
void determineStepSizeDynamically(Sampler::stateful_pointer<Sampler_state> self,
                                  ActiveConstraintGraph* cgK,
                                  ActiveConstraintRegion* region, bool dense,
                                  CayleyParameterization* cparam) {
  Settings* sett = Settings::getInstance();

  VLOG(3) << "determineStepSizeDynamically";
#ifdef VERBOSE
  cout << "determineStepSizeDynamically" << endl;
#endif
  double testStep = 0.4;  // was 0.2 before

  cgK->setStepSize(
      testStep);  // increment this stepsize, if it is causing slowness
  ConvexChart* chart = new ConvexChart(
      cgK, dense, cparam,
      sett->runTimeObjects.df);  // for volume computation, to make blue points
                                 // as minimum as possible

  /**
   * GRID sampling allows some range for contact hence GRID sampling volume per
   * node is proportional with that range of the contact. i.e. Not all root
   * nodes have same volume by GRID sampling. If we are not doing short-range
   * sampling and forcing contact to be exact distance, then in order to have
   * proportional sampling density with GRID sampling, we should set
   * expectedNumberOfSamples per node to be proportional with that ratio
   */
  double contact_lengthUpper = sett->runTimeObjects.df->bondingUpperBound(
      chart->getAtom(0), chart->getAtom(6));
  double contact_lengthLower = sett->runTimeObjects.df->bondingLowerBound(
      chart->getAtom(0), chart->getAtom(6));
  double range = pow(contact_lengthUpper, 3) - pow(contact_lengthLower, 3);

  double expectedNumberOfSamples = 500.;  // 10000.; //50000.;
  expectedNumberOfSamples = expectedNumberOfSamples * range;

  /** compute the approximate volume of the region */
  int volume = 0;
  if (chart->initializeChart(false, region) && cparam->is_partial3tree())
    for (; !chart->doneGrid(); chart->stepGrid()) volume++;

  delete chart;

  // settings::Constraint::stericConstraint = userDefined_stericConstraint;

  double voltimes = expectedNumberOfSamples / volume;
  double root = 1. / cgK->getDim();
  double pwr = pow(voltimes, root);
  double stepsizeee = testStep / pwr;
  cgK->setStepSize(stepsizeee);
}

void refigureFlip(Sampler::stateful_pointer<Sampler_state> self,
                  Orientation* orr, std::vector<std::vector<int> > flipScheme) {
  Settings* sett = Settings::getInstance();

  std::vector<std::vector<int> > tetras = flipScheme;

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

  LOG_IF(WARNING, (flip == -1))
      << "ERROR: above should have found flip 0 < i 7";
#ifdef VERBOSE
  if (flip == -1) {
    cout << "ERROR: above should have found flip 0 < i 7" << endl;
  }
#endif

  orr->setFlipNum(flip);
  tetraCoordinates.clear();
}

list<Orientation*> findRealizations(ConvexChart* chart,
                                    vector<vector<int> > flipScheme) {
  list<Orientation*> output;

  if (chart->partial3tree)
    for (int i = 0; i < 8; i++)  // for each flip
    {
      bool fail;
      // Orientation *relz = CartesianRealizer::computeRealization(cgK, des, i,
      // fail, flipScheme);
      Orientation* relz = chart->computeRealization(i, fail, flipScheme);
      if (!fail) {
        output.push_back(relz);
      } else {  // volume negative

        delete relz;
        relz = NULL;
        break;  // in the next flips it will always get volume negative too
      }
    }

  return output;
}

void createChildContactGraphs_fromTheBoundary(
    Sampler::stateful_pointer<Sampler_state> self,
    list<pair<int, int> >& contactList,
    list<Orientation*>::iterator ori_on_lattice, ActiveConstraintGraph* cgK,
    Orientation* orie_on_boundary, AtlasNode* rnode, bool& noGoodOrientation,
    int& noPoints, CayleyPoint* EntryCayleyPoint) {
  ActiveConstraintGraph* cgKprime[contactList.size()];
  bool firstPath = true;
  bool savedOnce = false;
  for (size_t c = 0; c < contactList.size(); c++) {  // TRY ALL CONTACT LIST
    cgKprime[c] = new ActiveConstraintGraph(cgK);
    pair<int, int> contact = contactList.front();
    cgKprime[c]->addContact(contact);
    contactList.pop_front();

    if (cgKprime[c]->constraintSize() > 6 || cgKprime[c]->isDependent() ||
        !cgKprime[c]->mayHaveThirdAtom()) {
      delete cgKprime[c];
      cgKprime[c] = NULL;
      contactList.push_back(contact);
      continue;
    }

    //--------------
    Orientation* wit_orr_toSend = new Orientation(
        orie_on_boundary);  // orie_on_boundary->getOrienation() ;

    // Insert Entry point
    wit_orr_toSend->setEntryPoint(rnode->getID(), EntryCayleyPoint->getID(),
                                  orie_on_boundary->getFlipNum());

    // wit_orr_toSend->addBoundary(rnode->getID());
    vector<pair<int, int> > part = cgKprime[c]->getParticipants();
    string graphString= constructSignature(part);
    Settings* sett = Settings::getInstance();
    self->send(sett->ActorSystem.AB, self->state.rnode->getID(), *(cgKprime[c]),
               *wit_orr_toSend, EntryCayleyPoint->getID(),
               int((*ori_on_lattice)->getFlipNum()),
               self->state.rnode->getFlipScheme());
    self->state.numNewRegionsFound++;
    self->state.numMessagesInABQ++;

    // TODO: Optimize by checking if we have sent this region before. This will
    // reduce the number of messages sent.

    VLOG(2) << "SamplerActor sent a new region to AB Actor from node: "
            << self->state.rnode->getID()
            << " ContactGraph: " << graphString 
            << " CayleyPoint: " << EntryCayleyPoint->getID()
            << " Flip: " << int((*ori_on_lattice)->getFlipNum());
#ifdef VERBOSE
    cout << "SamplerActor sent a new region to AB Actor from node: "
         << self->state.rnode->getID()
         << " ContactGraph: " << graphString 
         << " CayleyPoint: " << EntryCayleyPoint->getID()
         << " Flip: " << int((*ori_on_lattice)->getFlipNum()) << endl;
#endif

    // AtlasNode* child_node = createChildNode(self, rnode, cgKprime[c],
    //		wit_orr_toSend, savedOnce, leafWitness, noGoodOrientation,
    //		noPoints);
    // int child_nodeID = child_node->getID();
    // orie_on_boundary->addBoundary(child_nodeID);
    //(*ori_on_lattice)->addBoundary(child_nodeID);

    contactList.push_back(contact);  // to let 'for loop' to stop
    delete wit_orr_toSend;
    delete cgKprime[c];
  }
}

void findBoundary(Sampler::stateful_pointer<Sampler_state> self,
                  list<Orientation*>::iterator ori_on_lattice,
                  ConvexChart* desc, ActiveConstraintGraph* cgK,
                  AtlasNode* rnode, ActiveConstraintRegion* region,
                  ConstraintCheck* detector, bool& noGoodOrientation,
                  int& noPoints, bool& boundary_ori_found_and_saved,
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
  while (!desc->stepAround()) {
    bool fail;
    vector<vector<int> > flipScheme = rnode->getFlipScheme();
    Orientation* orie_on_boundary =
        desc->computeRealization(flip, fail, flipScheme);
    // Orientation* orie_on_boundary =
    //	CartesianRealizer::computeRealization(cgK, desc, flip, fail,
    //flipScheme); //if needs maple, then it will always return first root not
    //the one with specified flip !

    /// found collision, start binary search
    /// todo ACTUALLY CHECK SHOULD BE DONE AROUND BLUE REGION AS WELL,
    /// IF VOLUME NEGATIVE, SEARCH FOR BOUNDARY IN BETWEEN!!!
    if (sett->Sampling.binarySearch) {
      if (!fail && detector->stericsViolated(orie_on_boundary)) {
        desc->stepGridBinary(false);  ////walk through valid point

        delete orie_on_boundary;
        orie_on_boundary = desc->computeRealization(flip, fail, flipScheme);

        // check for collision and contacts
        // if there is a collision then it should continue binary stepping
        bool binvalid = !detector->stericsViolated(orie_on_boundary);
        list<pair<int, int> > contactList2 = detector->getContacts(flip);

        // binary search
        // continue search till you find a valid configuration with new contacts
        // find new contact for new small threshold t // FIXME: what does this
        // line mean?
        int num_bin_iteration = 0;
        while ((contactList2.empty() || !binvalid) && num_bin_iteration < 30) {
          num_bin_iteration++;
          desc->stepGridBinary(binvalid);
          delete orie_on_boundary;
          orie_on_boundary = desc->computeRealization(flip, fail, flipScheme);
          binvalid = !detector->stericsViolated(orie_on_boundary);
          contactList2 = detector->getContacts(flip);
        }

        // double check if the valid configuration with new contacts exists.
        // it may not in case it exits the loop because of num_bin_iteration is
        // big
        if (binvalid && contactList2.size() != 0 &&
            (sett->AtlasBuilding.ifBadAngleWitness_createChild ||
             !orie_on_boundary->angleViolated())) {
          // Orientation* orie_on_boundary = sr->getOrienation();
          createChildContactGraphs_fromTheBoundary(
              self, contactList2, ori_on_lattice, cgK, orie_on_boundary, rnode,
              noGoodOrientation, noPoints, EntryCayleyPoint);

          delete orie_on_boundary;  // TODO: Check back after valgrind run

          break;  // if found a collision in one direction then stop
        }  // if( valid(orie_on_boundary,detector) && contactList2.size() != 0
      }    // if( !fail &&  !valid(orie_on_boundary,detector)  )else
    }
    if (!boundary_ori_found_and_saved) {
      delete orie_on_boundary;
    }

    desc->setInitialPoint(pp);
  }  // while( !desc->stepGridContact()  && !settings::AtlasBuilding::stop)
  desc->setInitialPoint(pp);
  desc->setDir();
}

bool sampleAtlasNode(Sampler::stateful_pointer<Sampler_state> self, bool dense,
                     bool continu,
                     std::vector<std::vector<int> >& parentFlipScheme) {
  Settings* sett = Settings::getInstance();

  int from = self->state.rnode->getID();

  VLOG(1) << "AtlasBuilder::MySample: Started node " << from;
#ifdef VERBOSE
  cout << "AtlasBuilder::MySample: Started node " << from << endl;
#endif
  if (dense) self->state.rnode->setComplete(false);

  ActiveConstraintRegion* region = self->state.rnode->getACR();
  ActiveConstraintGraph* cgK = self->state.rnode->getCG();
  ConstraintCheck* detector = new ConstraintCheck(
      cgK, sett->runTimeObjects.df);  // Create a fresh ConstraintCheck detector
  CayleyParameterization* desc =
      new CayleyParameterization(cgK, false);  // Create a fresh description

  if (sett->Sampling.dynamicStepSizeAmong)
    determineStepSizeDynamically(self, cgK, region, dense, desc);

  bool noGoodOrientation = true;
  if (continu) noGoodOrientation = !self->state.rnode->hasAnyGoodOrientation();

  ConvexChart* chart =
      new ConvexChart(cgK, dense, desc, sett->runTimeObjects.df);

  // Do flip reconciliation.

  if (self->state.rnode->getFlipScheme().empty()) {
    if (chart->partial3tree) {
      VLOG(2) << "chart->partial3tree: " << chart->partial3tree;
#ifdef VERBOSE
      cout << "chart->partial3tree: " << chart->partial3tree << endl;
#endif
      self->state.rnode->setFlipScheme(chart->getTetras());
    } else {
      self->state.rnode->setFlipScheme(parentFlipScheme);
    }
  }

  // if the coming_witness_ori is not NULL, calculate the parameter value
  if (self->state.ori != NULL) {
    VLOG(2) << "We got a non-NULL orientation";
#ifdef VERBOSE
    cout << "We got a non-NULL orientation" << endl;
#endif
    CayleyPoint* wp4 = new CayleyPoint(
        self->state.ori, cgK->getParamLines(), sett->runTimeObjects.muA,
        sett->runTimeObjects.muB); // deletion of coming_witness_ori
                                   // is handled WHERE IT IS CREATED.

    refigureFlip(self, (wp4->getOrientations()[0]),
                 self->state.rnode->getFlipScheme());
    vector<Orientation> sendOriVecToWitnessWriter{*(self->state.ori)};
    chart->setWitnessPoint(wp4);
    if (!self->state.ori->angleViolated()) {
      region->AddWitnessPoint(wp4);
      // TODO:Call witness writer to write this
      self->send(self->state.writer, from, self->state.rnode->getParamDim(),
                 sendOriVecToWitnessWriter, cgK->getParamLines());
      self->send(sett->ActorSystem.AB, Write_v);
      noGoodOrientation = false;
    }
    // addReverseWitness(self->state.rnode, cgK, chart, region, self->state.ori,
    // noGoodOrientation);
  }

  int noPoints = 0;

  // continu : Flag of whether this is a first/new sample or a continued sample.
  if (desc->is_partial3tree() && chart->initializeChart(continu, region)) {
    VLOG(2) << "Started setting things up for sampling.";
#ifdef VERBOSE
    cout << "Started setting things up for sampling." << endl;
#endif
    list<Orientation*> real;
    for (; !chart->doneGrid(); chart->stepGrid()) {
      // find all realization for the grid point
      // real = findRealizations(cgK, chart,
      // self->state.rnode->getFlipScheme());
      real = findRealizations(chart, self->state.rnode->getFlipScheme());
      vector<double> out = chart->getPoint();
      CayleyPoint* p4d = new CayleyPoint(out);
      p4d->setID(region->nextSamplePointID());
      p4d->setRealizable(!real.empty());

      if (!real.empty()) {
        for (list<Orientation*>::iterator ori_on_lattice = real.begin();
             ori_on_lattice != real.end(); ori_on_lattice++) {
          // check for steric constraint violation and for any possible future
          // contacts
          if (!detector->stericsViolated(*ori_on_lattice)) {
            // check for bad angle
            if ((*ori_on_lattice)->angleViolated())
              p4d->incrementBadAngleN();
            else
              noGoodOrientation = false;

            bool boundary_ori_found_and_saved = false;

            if (sett->RootNodeCreation.createChildren) {
              // find the boundary by binary search
              findBoundary(self, ori_on_lattice, chart, cgK, self->state.rnode,
                           region, detector, noGoodOrientation, noPoints,
                           boundary_ori_found_and_saved, p4d);
            }

            // save both 'ori_on_lattice' and 'orie_on_boundary'
            if (!boundary_ori_found_and_saved) {
              if (!(*ori_on_lattice)->angleViolated()) {
                p4d->addOrientation((*ori_on_lattice));
              } else {
                delete (*ori_on_lattice);  // orange
              }
            } else {
              delete (*ori_on_lattice);
            }
          } else {
            if ((*ori_on_lattice)->angleViolated()) p4d->incrementBadAngleN();

            delete (*ori_on_lattice);  // collision
            p4d->incrementCollidN();
          }

        }  // endif real iterator
      }
      region->AddSamplePoint(p4d);
      noPoints++;
    }
    VLOG(2) << "Done with stepping through the grid.";
#ifdef VERBOSE
    cout << "Done with stepping through the grid." << endl;
#endif
  } else {
    VLOG(2) << "Didn't sample anything because:";
    VLOG(2) << "Partial 3-tree: " << desc->is_partial3tree()
            << " ChartInitialization: "
            << chart->initializeChart(continu, region);
#ifdef VERBOSE
    cout << "Didn't sample anything because:" << endl;
    cout << "Partial 3-tree: " << desc->is_partial3tree()
         << " ChartInitialization: " << chart->initializeChart(continu, region)
         << endl;
#endif
  }

  self->state.rnode->setFoundGoodOrientation(!noGoodOrientation);

  // TODO:This probably needs to go inside pause
  /*if(!sett->AtlasBuilding.stop ) {
          self->state.rnode->setComplete(true);
          region->trim();
  }*/

  delete detector;
  delete desc;
  delete chart;

  VLOG(1) << "AtlasBuilder::MySample: Finished node " << from;
#ifdef VERBOSE
  cout << "AtlasBuilder::MySample: Finished node " << from << endl;
#endif

  return !noGoodOrientation;
}

void sendCayleyPointsToWriter(
    const Sampler::stateful_pointer<Sampler_state>& self) {
  std::vector<CayleyPointStruct> vecCps;
  Settings* sett = Settings::getInstance();
  for (CayleyPoint* cp : self->state.rnode->getACR()->GetSamplePoints()) {
    if (cp == NULL) continue;
    CayleyPointStruct cps = cp->getCayleyPointStruct();
    vecCps.push_back(std::move(cps));
  }
  self->send(self->state.writer, self->state.rnode->getID(),
             self->state.rnode->getParamDim(), vecCps);
  self->send(sett->ActorSystem.AB, Write_v);
  vecCps.erase(vecCps.begin(), vecCps.end());
}

Sampler::behavior_type typedSampler(
    Sampler::stateful_pointer<Sampler_state> self, Writer writer) {
  Settings* sett = Settings::getInstance();

  self->state.writer = writer;
  return {
      [=](AtlasNodeStruct ans, Orientation ori,
          std::vector<std::vector<int> > parentFlipScheme) {
        VLOG(2) << "SamplerActor received Sampling message for node: "
                << ans.numID;
#ifdef VERBOSE
        cout << "SamplerActor received Sampling message for node: " << ans.numID
             << " at TIME :" << high_resolution_clock::now() << endl;
#endif
        self->state.rnode = new AtlasNode(ans);
        self->state.ori = new Orientation(ori);
        bool dense = false;
        bool continu = false;
        sampleAtlasNode(self, dense, continu, parentFlipScheme);
        VLOG(2) << "SamplerActor sent sampling results for node: " << ans.numID;
#ifdef VERBOSE
        cout << "SamplerActor sent sampling results for node: " << ans.numID
             << " at TIME :" << high_resolution_clock::now() << endl;
#endif
        if (self->state.rnode->getParamDim() == 0 ||
            self->state.numNewRegionsFound == 0) {
          sendCayleyPointsToWriter(self);
          VLOG(2) << "SamplerActor about to send DONE message for node: "
                  << self->state.rnode->getID();
#ifdef VERBOSE
          cout << "SamplerActor about to send DONE message for node: "
               << self->state.rnode->getID()
               << " at TIME :" << high_resolution_clock::now() << endl;
#endif

          self->send(sett->ActorSystem.AB, Done_v,
                     self->state.rnode->getID(),
                     self->state.rnode->hasAnyGoodOrientation());
        }
        // return make_tuple(ans.numID, self->state.childNodes, true,
        // self->state.rnode->hasAnyGoodOrientation(),
        // self->state.rnode->getFlipScheme());
      },

      [=](AtlasNodeStruct ans,
          std::vector<std::vector<int> > parentFlipScheme) {
        VLOG(2) << "SamplerActor received Sampling message for node "
                << ans.numID;
#ifdef VERBOSE
        // cout <<"Started Sampling Node number "<<ans.numID<<endl;
        cout << "SamplerActor received Sampling message for node " << ans.numID
             << " at TIME :" << high_resolution_clock::now() << endl;
#endif
        self->state.rnode = new AtlasNode(ans);
        self->state.ori = NULL;
        bool dense = false;
        bool continu = false;
        sampleAtlasNode(self, dense, continu, parentFlipScheme);
        VLOG(2) << "SamplerActor sent sampling results for node: " << ans.numID;
#ifdef VERBOSE
        cout << "SamplerActor sent sampling results for node: " << ans.numID
             << " at TIME :" << high_resolution_clock::now() << endl;
#endif
        if (self->state.rnode->getParamDim() == 0 ||
            self->state.numNewRegionsFound == 0) {
          sendCayleyPointsToWriter(self);
          VLOG(2) << "SamplerActor about to send DONE message for node: "
                  << self->state.rnode->getID();
#ifdef VERBOSE
          cout << "SamplerActor about to send DONE message for node: "
               << self->state.rnode->getID()
               << " at TIME :" << high_resolution_clock::now() << endl;
#endif

          self->send(sett->ActorSystem.AB, Done_v,
                     self->state.rnode->getID(),
                     self->state.rnode->hasAnyGoodOrientation());
        }
        // return make_tuple(ans.numID, self->state.childNodes, true,
        // self->state.rnode->hasAnyGoodOrientation(),
        // self->state.rnode->getFlipScheme());
      },

      [=](vector<tuple<int, int, int>> boundaryList) -> result<Done, int, bool> {
#ifdef VERBOSE
        cout << "SamplerActor received boundaryList message for node "
             << self->state.rnode->getID()
             << " at TIME :" << high_resolution_clock::now() << endl;
#endif
        vector<CayleyPoint*> region = self->state.rnode->getACR()->getSpace();
        for (auto it = boundaryList.begin(); it != boundaryList.end(); it++) {
          // CayleyPoint *cp =
          // self->state.rnode->getACR()->getEntryPoint(std::get<0>(*it));
          CayleyPoint* cp = region[std::get<0>(*it)];
          Orientation* ori = cp->getOrientation(std::get<1>(*it));
          if (ori != NULL) {  // Angle wasn't violated.
            ori->addBoundary(std::get<2>(*it));
          }
        }
        sendCayleyPointsToWriter(self);
        VLOG(2) << "Saving points to the disk.";
#ifdef VERBOSE
        cout << "Saving points to the disk." << endl;
#endif
        VLOG(1) << "SamplerActor sent DONE message for node: "
                << self->state.rnode->getID();
#ifdef VERBOSE
        cout << "SamplerActor sent DONE message for node: "
             << self->state.rnode->getID()
             << " at TIME :" << high_resolution_clock::now() << endl;
#endif
        return {Done_v, self->state.rnode->getID(),
                               self->state.rnode->hasAnyGoodOrientation()};
      },

      [=](int cpID, int flipNum, int nodeNum) {
        self->state.numMessagesInABQ--;
        VLOG(2) << "SamplerActor received boundary message for node: "
                << self->state.rnode->getID() << " CayleyPoint: " << cpID
                << " Flip: " << flipNum
                << " BoundaryNode: " << nodeNum
                << " NumMessagesInABQ: " << self->state.numMessagesInABQ;
#ifdef VERBOSE
        cout << "SamplerActor received boundary message for node: "
             << self->state.rnode->getID() << " for CayleyPoint: " << cpID
             << " flipNum:" << flipNum
             << " boundary: " << nodeNum
             << " NumMessagesInABQ: " << self->state.numMessagesInABQ
             << " at TIME :" << high_resolution_clock::now() << endl;
#endif
        // vector<CayleyPoint*> region =
        self->state.rnode->getACR()->setBoundaryPoint(cpID, flipNum, nodeNum);
        /*CayleyPoint *cp = region[cpID];
        Orientation* ori = cp->getOrientation(flipNum);
        if (ori!=NULL) { //Angle wasn't violated.
            ori->addBoundary(nodeNum);
        }*/

        if (self->state.numMessagesInABQ == 0) {
          VLOG(2) << "Saving points to the disk.";
#ifdef VERBOSE
          cout << "Saving points to the disk." << endl;
#endif
          sendCayleyPointsToWriter(self);
          VLOG(1) << "SamplerActor about to send DONE message for node: "
                  << self->state.rnode->getID();
#ifdef VERBOSE
          cout << "SamplerActor about to send DONE message for node: "
               << self->state.rnode->getID()
               << " at TIME :" << high_resolution_clock::now() << endl;
#endif

          self->send(sett->ActorSystem.AB, Done_v,
                     self->state.rnode->getID(),
                     self->state.rnode->hasAnyGoodOrientation());
          // TODO: Take care of HAS ANY GOOD ORIENTATION
        }
      },

      [=](Pause) -> result<int, std::unordered_map<std::string, 
                                   std::pair<ActiveConstraintGraph,
                                   std::vector<std::tuple<Orientation, int, int>>>>,
                                   bool, bool, std::vector<std::vector<int>>> {
        // TODO: Fix the situation of immidiate sending
        VLOG(1) << "SamplerActor received PAUSE message for node: "
                << self->state.rnode->getID();
#ifdef VERBOSE
        cout << "SamplerActor received PAUSE message for node: "
             << self->state.rnode->getID()
             << " at TIME :" << high_resolution_clock::now() << endl;
#endif
        LOG(INFO) << "Pausing";
        cout << "Pausing" << endl;
#ifdef VERBOSE
        cout << "SamplerActor sent sampling results for node: "
             << self->state.rnode->getID()
             << " at TIME :" << high_resolution_clock::now() << endl;
#endif
        return {self->state.rnode->getID(), self->state.childNodes,
                          false, self->state.rnode->hasAnyGoodOrientation(),
                          self->state.rnode->getFlipScheme()};
      },
      [=](Kill) {
        VLOG(1) << "SamplerActor received KILL message for node: "
                << self->state.rnode->getID();
#ifdef VERBOSE
        cout << "SamplerActor received KILL message for node: "
             << self->state.rnode->getID()
             << " at TIME :" << high_resolution_clock::now() << endl;
        cout << "SamplerActor sent EXIT message for node: "
             << self->state.rnode->getID()
             << " at TIME :" << high_resolution_clock::now() << endl;
#endif
        delete self->state.rnode;
        anon_send_exit(self, exit_reason::normal);
      }};
}
