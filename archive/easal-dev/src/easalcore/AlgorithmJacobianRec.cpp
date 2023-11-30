/*
 * AlgorithmJacobianRec.cpp
 *
 *  Created on: May 27, 2014
 *      Author: Aysegul Ozkan
 */

#include "AlgorithmJacobianRec.h"

#include "ActiveConstraintGraph.h"
#include "ActiveConstraintRegion.h"
#include "Atlas.h"
#include "CartesianRealizer.h"
#include "CayleyParameterization.h"
#include "CayleyPoint.h"
#include "CgMarker.h"
#include "ConstraintCheck.h"
#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/Geometry"
#include "Eigen/LU"
#include "Eigen/SVD"
#include "MolecularUnit.h"
#include "SaveLoader.h"
#include "Settings.h"
#include "Utils.h"
//#include <Eigen/JacobiSVD.h>
using namespace Eigen;
using Eigen::MatrixXd;
using Eigen::Quaterniond;

//#include<stdlib.h>
#include <cassert>
#include <cmath>
#include <iterator>
#include <limits>
#include <list>
#include <set>
#include <vector>
using namespace std;

#define PI 3.14159265

double AlgorithmJacobianRec::cayleyStep[6];
MatrixXd AlgorithmJacobianRec::cayleyDirections;
// MatrixXd AlgorithmJacobianRec::cayleyDirections = MatrixXd::Identity(5, 5)*
// Settings::Sampling::stepSize;
bool AlgorithmJacobianRec::aDebug = true;

static double infinit;
AlgorithmJacobianRec::AlgorithmJacobianRec(MolecularUnit* a, MolecularUnit* b,
                                           SaveLoader* snl,
                                           PredefinedInteractions* df,
                                           Atlas* atlas)
    : AtlasBuilder(a, b, snl, df, atlas) {
  Settings* set = Settings::getInstance();
  for (int i = 0; i < 6; i++) {
    cayleyStep[i] = set->Sampling.stepSize;
  }
  infinit = std::numeric_limits<float>::infinity();
}

AlgorithmJacobianRec::AlgorithmJacobianRec() {
  Settings* set = Settings::getInstance();
  for (int i = 0; i < 6; i++) {
    cayleyStep[i] = set->Sampling.stepSize;
  }
  infinit = std::numeric_limits<float>::infinity();
}

void AlgorithmJacobianRec::startAtlasBuilding() {
  Settings* set = Settings::getInstance();

  cout << "AlgorithmJacobianRec this->rootGraphs.size()"
       << this->rootGraphs.size() << endl;
  int currentrootGraph = 0;

  bool empty;
  ActiveConstraintGraph* nextrootGraph = getNextRootGraph(empty);

  while (
      !empty &&
      !set->AtlasBuilding.stop) {  // rootGraphs is the list of the pregenerated
                                   // rootGraphs which will have many or one.

    // nextrootGraph->hybrid = true;  ////commented on 4jan13 to make it exact
    // 5d node
    int nodenum;
    int success = this->atlas->addNode(
        nextrootGraph, nodenum,
        -1);  // ATTENTION: firstParentID has been hardcoded -1 for root nodes.
    if (success == 1) {
      cout << "we start sampling rootGraph " << currentrootGraph << " out of "
           << this->rootGraphs.size() << endl;
      AtlasNode* rnode = (*this->atlas)[nodenum];

      int symmNodeNum = this->atlas->getSymmetricNode(rnode);
      if (symmNodeNum != -1) {
        AtlasNode* symmNode = (*this->atlas)[symmNodeNum];

        vector<pair<int, int> > particp = symmNode->getCG()->getParticipants();
        vector<pair<int, int> > sparticp;
        sparticp.push_back(make_pair(particp[0].second, particp[0].first));

        vector<pair<int, int> > plines = symmNode->getCG()->getParamLines();
        vector<pair<int, int> > splines;
        for (int i = 0; i < plines.size(); i++)
          splines.push_back(make_pair(plines[i].second, plines[i].first));

        ActiveConstraintGraph* reverseNode =
            new ActiveConstraintGraph(sparticp, splines);
        //				reverseNode->use_existing_parametrization =
        //true;
        rnode->setCG(reverseNode);

        //				string symm_param =
        //symmNode->getCG()->getParams(); 				string reverse_param = ""; 				for(int
        //i=0; i<symm_param.size();  i=i+5){ //to start from reverse if you put
        //contacts at the end 					reverse_param += symm_param.substr(i+2,2) ;
        //					reverse_param += symm_param.substr(i,2)
        //; 					reverse_param += " ";
        //				}
        //				rnode->getCG()->setParams(reverse_param);
      }

      // this->MySample(rnode, false, NULL, false, false);
      MyJacobianSampleRec(rnode, false, NULL, false,
                          set->AtlasBuilding.breadthFirst, false);
    }
    nextrootGraph = getNextRootGraph(empty);

    currentrootGraph++;

    //		if( this->atlas->getNodeNum( this->rootGraphs.front() ) == -1  ){
    ////if( !this->rootGraphs.front()->isComplete()  ){ 			int nodenum;
    //			this->atlas->addNode( this->rootGraphs.front(),  nodenum
    //); 			snl->saveNode( this->atlas->getNode(nodenum) );
    //			this->MySample(this->rootGraphs.front(), false, NULL, false, false,
    //false);
    //		}

    //		this->rootGraphs.pop_front(); // that rootGraph is removed from the
    //list and the list size is smaller.
  }
}

map<int, double> AlgorithmJacobianRec::stepRatios(list<Orientation*> oldreal,
                                                  list<Orientation*> real) {
  map<int, double> output;
  bool found = false;
  for (list<Orientation*>::iterator ori = real.begin(); ori != real.end();
       ori++) {
    Vector3d eaB = (*ori)->eaB;
    Vector3d TB = (*ori)->TB;

    for (list<Orientation*>::iterator old_ori = oldreal.begin();
         old_ori != real.end(); old_ori++) {
      if ((*ori)->getFlipNum() == (*old_ori)->getFlipNum()) {
        Vector3d old_eaB = (*old_ori)->eaB;
        Vector3d old_TB = (*old_ori)->TB;

        Vector3d diff_TB = TB - old_TB;
        Vector3d diff_eaB = eaB - old_eaB;
        diff_eaB = diff_eaB / 10;  /// normalize

        double dist =
            sqrt(diff_TB[0] * diff_TB[0] + diff_TB[1] * diff_TB[1] +
                 diff_TB[2] * diff_TB[2] + diff_eaB[0] * diff_eaB[0] +
                 diff_eaB[1] * diff_eaB[1] + diff_eaB[2] * diff_eaB[2]);
        // output.push_back( dist );
        output[(*ori)->getFlipNum()] = dist;
        found = true;
        break;
      }
    }
    //		if(!found)
    //			output.push_back( -1);
  }

  return output;
}

/*
//newest
list<CartesianRealizer*> AlgorithmJacobianRec::findRealizations(ConvexChart*
des, list<CartesianRealizer*> oldreal, list<CartesianRealizer*> &real)
{
        list<CartesianRealizer*> output;
        //output = real; //include initial ones as well
        bool map = des->mapl;
        int msol=0;
        double a[6][4];
        string ss;
        int dim = des->getParams().size()/5;

        if( !map )
        {
                int orisize = real.size();
                int cori = 0;
                for (list<CartesianRealizer*>::iterator ori = real.begin(); ori
!= real.end(); ori++)
                {
                        cori++;
                        orisize = real.size();

                        Vector3d eaB = (*ori)->eaB;
                        Vector3d TB = (*ori)->TB;

                        double pconnections[12][12];
                        for(int i = 0; i<12; i++)	for(int j = 0; j<12;
j++)	pconnections[i][j] = (*ori)->pconnections[i][j];


                        for (list<CartesianRealizer*>::iterator old_ori =
oldreal.begin(); old_ori != real.end(); old_ori++)
                        {
                                if( (*ori)->getFlip() == (*old_ori)->getFlip() )
                                {
                                        Vector3d old_eaB = (*old_ori)->eaB;
                                        Vector3d old_TB = (*old_ori)->TB;

                                        Vector3d diff_TB = TB - old_TB;
                                        Vector3d diff_eaB = eaB - old_eaB;
                                        for(int h=0; h<3; h++) if(
abs(diff_eaB[h]) > 180 )  diff_eaB[h] = -1 * Utils::sign(diff_eaB[h]) * (360 -
abs(diff_eaB[h]) );

                                        diff_eaB = diff_eaB/10;  ///normalize

                                        double dist = sqrt(
diff_TB[0]*diff_TB[0] + diff_TB[1]*diff_TB[1] + diff_TB[2]*diff_TB[2] +
diff_eaB[0]*diff_eaB[0] + diff_eaB[1]*diff_eaB[1] + diff_eaB[2]*diff_eaB[2] );


                                        int ratio = round(dist);

                                        if(dist < 0.5) //if consecutive 2
orientation is too close, you need to delete the new one.
                                        {
                                                delete *ori;
                                                real.erase(ori);
                                                ori--;
                                                cout << "deleted " << cori << "
out of " << orisize << endl;
                                        }
                                        for(int s=1; s<ratio; s++) //we need
ratio times intermediate steps
                                        {

                                                double old_pconnections[12][12];
                                                for(int i = 0; i<12; i++)
for(int j = 0; j<12; j++)	old_pconnections[i][j] =
(*old_ori)->pconnections[i][j];

                                                for(int i=0; i<12; i++)
                                                        for(int j=0; j<12; j++)
                                                                old_pconnections[i][j]
+= s*(pconnections[i][j] - old_pconnections[i][j])/ratio;

                                                CartesianRealizer *relz = new
CartesianRealizer(des, (*ori)->getFlip(), old_pconnections); if(
relz->getVolumeTest() ){ if( !Settings::Collision::checkAngle ||
!relz->angleViolated() ) //small anglee  //relz->angleBetweenTwoMolecularUnits()
                                                                output.push_back(
relz ); else{
                                                                //delete relz;
                                                                output.push_back(
relz ); //child nodes can cause small angles so if this realization cause a new
contact, then create the child not and keep it as witness at the child node.
After that delete this realization from parent(this node), and make it orange.

                                                        }
                                                }
                                                else{ //volume negative
                                                        delete relz;
                                                        break; //in the next
flips it will always get volume negative too
                                                }
                                        }



                                        break;
                                }
                        }

                }
        }

        return output;
}
*/

CayleyPoint* AlgorithmJacobianRec::findStartPoint(
    ActiveConstraintGraph* cgK, ConvexChart* desc, int flip,
    ActiveConstraintRegion* region,
    AtlasNode* rnode)  // AtlasNode *rnode
{
  Settings* set = Settings::getInstance();
  cout << "trying to find a start point" << endl;

  Orientation* r = NULL;
  Orientation* temp_r = NULL;
  CayleyPoint* cayleyPoint = NULL;
  CayleyPoint* cayleyPoint_t = NULL;
  // ActiveConstraintGraph *cgK = rnode->getCG();

  ConstraintCheck* detector = new ConstraintCheck(
      desc->currentGraph, df);  // Create a fresh Collision detector

  // ConvexChart *desc = new ConvexChart( cgK, false ); // Create a fresh
  // description

  int dim = desc->parameters.size();
  MatrixXd gridSteps;
  gridSteps.setIdentity(dim, dim);
  for (int g = 0; g < dim; g++) {
    int d = desc->independent_directions[g];
    gridSteps(g, g) = set->Sampling.gridSteps[d];
  }

  desc->stepSize = 0.4;
  bool found1feasible = false;
  if (desc->initializeChart(false, region)) {
    for (; !desc->doneGrid() && !set->AtlasBuilding.stop; desc->stepGrid()) {
      // r = new CartesianRealizer(desc, flip);
      bool fail;
      r = CartesianRealizer::computeRealization(cgK, desc, flip, fail,
                                                rnode->getFlipScheme());
      vector<double> cayleyData = desc->getParamPoint();
      cayleyPoint = new CayleyPoint(cayleyData);
      cayleyPoint->setRealizable(true);
      cayleyPoint->addOrientation(r);

      if (!fail) {
        // cout << "1";
        bool myvalid = !detector->stericsViolated(r);
        if (myvalid) {  // checks against rest of the helix
          if (!r->angleViolated()) {
            if (!found1feasible) {
              // temp_r = new CartesianRealizer(desc, flip);
              temp_r = CartesianRealizer::computeRealization(
                  cgK, desc, flip, fail, rnode->getFlipScheme());
              found1feasible = true;

              vector<double> cayleyData_t = desc->getParamPoint();
              cayleyPoint_t = new CayleyPoint(cayleyData_t);
              cayleyPoint_t->setRealizable(true);
              cayleyPoint_t->addOrientation(temp_r);
            }
            bool succeed;
            MatrixXd JacInv =
                computeJacobian(rnode, cgK, desc, cayleyPoint, flip, succeed,
                                -1, gridSteps, infinit);
            if (succeed) {
              // end the for loop
              break;  // FOUND A VALID POINT INSIDE THE GRID BOUNDARIES and can
                      // compute Jacobian at that point
            }
          }
        }
      }

      // delete r;
      // r = NULL;
      cayleyPoint->trim_PointMultiD();
      delete cayleyPoint;
      cayleyPoint = NULL;
    }
  }
  delete detector;
  // delete desc;

  if (cayleyPoint != NULL) {
    cout << "FOUND A VALID POINT INSIDE THE GRID BOUNDARIES" << endl;
    cayleyPoint_t->trim_PointMultiD();
    delete cayleyPoint_t;
    cayleyPoint_t = NULL;
  } else {
    cout << "FAILED to find a start point" << endl;
    if (cayleyPoint_t != NULL) {
      cout << " but will still return a feasible realization with bad jacobian"
           << endl;
      cayleyPoint = cayleyPoint_t;
    }
  }
  return cayleyPoint;
}

static MatrixXd besttCayleyDirections;
static MatrixXd bestGridSteps;
bool AlgorithmJacobianRec::MyJacobianSampleRec(AtlasNode* rnode, bool dense,
                                               Orientation* orrw, bool continu,
                                               bool breath, bool bret) {
  Settings* set = Settings::getInstance();

  this->cayleyDirections =
      MatrixXd::Identity(rnode->getDim(), rnode->getDim()) *
      set->Sampling.stepSize;
  cout << "Dimension is " << rnode->getDim() << endl;
  if (AlgorithmJacobianRec::aDebug) cout << "Sample Started" << endl;

  ActiveConstraintGraph* cgK = rnode->getCG();
  ActiveConstraintRegion* region = rnode->getACR();

  cout << *cgK;

  // Settings::Sampling::stepSize = 0.2; // sqrt(
  // (Settings::Sampling::gridSteps[0]*Settings::Sampling::gridSteps[0]) /
  // cgK->getDim() );  //  dim * s^2 = g^2
  cout << "Settings::Sampling::stepSize set to be " << set->Sampling.stepSize
       << endl;

  int from = rnode->getID();  // atlas->getNodeNum(cgK); // from is the node
                              // number of the incoming contact graph as stored
                              // in the Roadmap object.
  if (AlgorithmJacobianRec::aDebug) cout << "from  " << from << endl;

  cout << "MyJacobianSampleRec " << from << endl;
  bool mynewstart = continu;  // Flag of whether this is a first/new sample or a
                              // continued sample.

  ConstraintCheck* detector =
      new ConstraintCheck(cgK, df);  // Create a fresh Collision detector
  if (AlgorithmJacobianRec::aDebug) cout << "Created detector" << endl;

  CayleyParameterization* desc =
      new CayleyParameterization(cgK, dense);  // Create a fresh description
  if (AlgorithmJacobianRec::aDebug) cout << "Created description" << endl;
  int dim = desc->getParameters().size();

  cout << set->Sampling.gridSteps[0] << " " << endl;

  snl->saveNodeHeader(rnode, this->atlas);  // save node after parameters is
                                            // set.

  if (!continu) this->snl->appendDimension(rnode);

  ConvexChart* chart = new ConvexChart(cgK, dense, desc, df);

  // commented out to sample 3D nodes. however, this particular 3D node (0-0 1-0
  // 2-1) IS a partial 3 tree. meaning computeRealization does not use the
  // flipScheme vector. so it doesn't matter if this is commented out or not.
  if (rnode->getFlipScheme().empty()) {
    // cout << "chart -> partial3tree: " << chart -> partial3tree << endl;

    if (chart->partial3tree) {
      rnode->setFlipScheme(chart->getTetras());

    } else {
      int rnum = rnode->getID();

      vector<int> parents = this->atlas->getParents(rnum);

      // probably a faster way to do this w/ .sort or similar
      // choose the parent w/ the lowest node ID, as that will be the first
      // parent (atlas is built w/dfs) is it the case here that rnode always
      // only has one parent?
      int firstParent = parents[0];
      for (int h = 1; h < parents.size(); h++) {
        if (parents[h] < firstParent) {
          firstParent = parents[h];
        }
      }

      std::vector<std::vector<int> > parent_flipScheme =
          this->atlas->getNodes()[firstParent]->getFlipScheme();
      rnode->setFlipScheme(parent_flipScheme);
    }
  }

  if (orrw != NULL) {
    CayleyPoint* wp4 = new CayleyPoint(
        orrw, cgK->getParamLines(), a, b,
        region->nextWitnessPointID());  // deletion of orrw is handled
                                             // WHERE IT IS CREATED? idk

    refigureFlip((wp4->getOrientations()[0]), rnode->getFlipScheme());

    chart->setWitnessPoint(wp4);
    if (!orrw->angleViolated()) {
      region->AddWitnessPoint(wp4);
      // cgK -> AddWitnessPoint(wp4); // add it to the graph // it may be needed
      // as a starting point to sample in the description class note: it will
      // 'not' be saved to file 'twice' in next steps again. only space points
      // are saved at the next steps, there is no more witness and witness
      // saving job later.
      this->snl->appendWitness(rnode->getID(),
                               wp4);  // save it to the nodeX.txt
                                      // noGoodOrientation = false;
    }
  }

  // don't add reverse witness for now (addReverseWitness). the setting is false
  // for this.

  int d1 = chart->pr1 / 2;
  int d2 = chart->pr2 / 2;
  int d3 = chart->pr3 / 2;
  int d4 = chart->pr4 / 2;
  int d5 = chart->pr5 / 2;
  int d6 = chart->pr6 / 2;

  int end_coordinats[6] = {chart->pr1, chart->pr2, chart->pr3,
                           chart->pr4, chart->pr5, chart->pr6};

  for (int i = 0; i < 6; i++) chart->end_coordinats[i] = end_coordinats[i];

  int noPoints = 0;
  int loop = 0;
  // int flip = 4;
  for (int flip = 0; flip < 8; flip++) {
    cout << "FLIPPPPPPPP " << flip << endl;
    chart->gridDone = false;
    for (size_t i = 0; i < chart->end_coordinats[0]; i++)
      for (size_t j = 0; j < chart->end_coordinats[1]; j++)
        for (size_t k = 0; k < chart->end_coordinats[2]; k++)
          for (size_t l = 0; l < chart->end_coordinats[3]; l++)
            for (size_t m = 0; m < chart->end_coordinats[4]; m++)
              for (size_t n = 0; n < chart->end_coordinats[5]; n++)
                (chart->visited[i][j][k][l][m][n]) = false;

    chart->initializeChart(false, region);
    bool fail;
    // Orientation *sreal = CartesianRealizer::computeRealization(cgK, chart,
    // flip, fail, rnode -> getFlipScheme());

    /// step grid until you find a valid realization.
    Orientation* sreal;
    for (; !chart->doneGrid() && !set->AtlasBuilding.stop; chart->stepGrid()) {
      sreal = CartesianRealizer::computeRealization(cgK, chart, flip, fail,
                                                    rnode->getFlipScheme());
      /// break as soon as you find a valid realization.
      if (sreal != NULL) break;
    }
    /// what if you absolutely can't find one? continue with the current loop.
    if (sreal == NULL) continue;

    // return false;
    // desc -> independent_directions.clear(); // then you will have the issue
    // to which axis to project among a bunch of sets.

    /// bool fail;
    /// Orientation* sreal = CartesianRealizer::computeRealization(cgK, chart,
    /// flip, fail, rnode -> getFlipScheme());
    vector<double> cayleyData;
    cayleyData = chart->getParamPoint();
    CayleyPoint* cayleyPoint_s = new CayleyPoint(cayleyData);
    cayleyPoint_s->setRealizable(true);
    cayleyPoint_s->addOrientation(sreal);

    /**
    // when did this happen? and why? because it causes an error later when the
    method jacobian is called if(cayleyPoint_s -> getOrientations().size() == 0)
            continue;
    */

    // ----------------------------------
    bool enteredNegVal;
    bool gridStepsChange[6];

    MatrixXd J = jacobian(cgK, chart, cayleyPoint_s, flip, -1, gridStepsChange,
                          cayleyDirections, enteredNegVal, rnode);

    vector<int> indep = computeIndependentColumns(
        J, dim);  // to allow directions to be initialized

    chart->independent_directions.clear();
    // desc -> independent_directions = indep; // somehow these directions are
    // not good!! so I commented it

    if (chart->independent_directions.size() ==
        0)  // if still not initialized because of neg vol
      for (int i = 0; i < dim; i++) chart->independent_directions.push_back(i);
    /**
    if(dim == 2) { // somehow for dim2 for the specific node i am using as a
    test case, this directions work best! chart ->
    independent_directions.clear(); chart ->
    independent_directions.push_back(0); chart ->
    independent_directions.push_back(2);
    }
    */
    cgK->independent_directions.clear();
    cgK->independent_directions.assign(chart->independent_directions.begin(),
                                       chart->independent_directions.end());
    // ----------------------------------

    /**
    if(dim >= 5) {
            delete cayleyPoint_s;
            cayleyPoint_s = NULL;
            cayleyPoint_s = findStartPoint(cgK, chart, flip, region, rnode);
            if(cayleyPoint_s != NULL)
                    sreal = cayleyPoint_s -> getOrientations().front();
    }

    if(cayleyPoint_s == NULL)
            continue;
    */

    cout << "sreal TB_EAB " << sreal->TB[0] << " " << sreal->TB[1] << " "
         << sreal->TB[2] << " " << sreal->eaB[0] << " " << sreal->eaB[1] << " "
         << sreal->eaB[2] << " " << endl;

    MatrixXd gridSteps;
    gridSteps.setIdentity(dim, dim);
    for (int g = 0; g < dim; g++) {
      int d = chart->independent_directions[g];
      gridSteps(g, g) = set->Sampling.gridSteps[d];
    }

    besttCayleyDirections = cayleyDirections;
    //		bestGridSteps = gridSteps;

    bool succeed;
    MatrixXd JacInv = computeJacobian(
        rnode, cgK, chart, cayleyPoint_s, flip, succeed, -1, gridSteps,
        infinit);  // to make sure JacInv will really lead to gridSteps at least
                   // in terms grid directions. i.e. coming_gridDirections will
                   // have same direction as gridSteps.
    // since specificDirection is -1, it will NOT reverse grid directions, hence
    // no need to find 'gridStepsLater' and update 'coming_gridDirections'

    //    	MatrixXd gridStepsLater = bestGridSteps;
    JacInv = besttCayleyDirections;

    vector<bool> coming_gridDirections;
    for (int i = 0; i < dim; i++) coming_gridDirections.push_back(true);
    //		for(int h=0; h<dim; h++) //make coming_gridDirections consistent with
    //bestGridSteps since gridSteps may have changed in computeJacobian. SO that
    //coming_gridDirections will be consistent with JacInv 			if(
    //gridStepsLater(h,h)*gridSteps(h,h) < 0 ) 				coming_gridDirections[h] =
    //!coming_gridDirections[h];

    //---------------------
    //		int
    //coordinats[6]={desc->pr/2,desc->pr/2,desc->pr/2,desc->pr/2,desc->pr/2,desc->pr/2};
    //		int coordinats[6] = {desc->pr1/2, desc->pr2/2, desc->pr3/2,
    //desc->pr4/2, desc->pr5/2, desc->pr6/2};

    int temp_coor[6];
    temp_coor[0] = round((sreal->TB[0] + set->Sampling.GridXY) /
                         set->Sampling.gridSteps[0]);  // 19.5
    temp_coor[1] = round((sreal->TB[1] + set->Sampling.GridXY) /
                         set->Sampling.gridSteps[1]);  // 19.5
    temp_coor[2] =
        round(sreal->TB[2] + set->Sampling.GridZ / set->Sampling.gridSteps[2]);
    temp_coor[3] =
        round((sreal->eaB[0] + 180) /
              set->Sampling
                  .gridSteps[3]);  // added +1 to prevent floating point errors.
    temp_coor[4] =
        round((sreal->eaB[1] /**- 0.86602540378*/) /
              set->Sampling
                  .gridSteps[4]);  // allows you 4 steps, 0.866, 0.91, 0.954, 1
    temp_coor[5] = round((sreal->eaB[2] + 180) / set->Sampling.gridSteps[5]);
    //			temp_coor[4] =floor(
    //(sreal->eaB[1]-0.883333)/Settings::Sampling::gridSteps[4] );

    //			double cos_eaB1 = cos(-sreal->eaB[1]*PI/180);
    //			if(cos_eaB1<.91)  			//0.883333 costheta,    27.952927333
    //degrees 				d5=0; 			else if(cos_eaB1<.955)  	//0.93 costheta,  21.565185015
    //degrees 				d5=1;
    //			else					//0.976667 costheta,  12.401408289
    //degrees 				d5=2;

    // 2*Settings::Sampling::GridXY/Settings::Sampling::gridSteps[0] => the
    // number of steps within grid range
    temp_coor[0] +=
        (chart->pr1 - (2 * set->Sampling.GridXY / set->Sampling.gridSteps[0])) /
        2;  // since the grid size increased +2 to prevent rounding issues.
    temp_coor[1] +=
        (chart->pr2 - (2 * set->Sampling.GridXY / set->Sampling.gridSteps[0])) /
        2;
    temp_coor[2] +=
        (chart->pr3 - (2 * set->Sampling.GridZ / set->Sampling.gridSteps[2])) /
        2;
    temp_coor[3] += (chart->pr4 - 36) / 2;
    temp_coor[4] += (chart->pr5 - 5) / 2;
    temp_coor[5] += (chart->pr6 - 36) / 2;

    //		temp_coor[0]++;  //since the grid size increased +2 to prevent
    //rounding issues. 		temp_coor[1]++; 		temp_coor[2]++; 		temp_coor[3]++;
    //		temp_coor[4]++;
    //		temp_coor[5]++;

    int coordinats[6] = {d1, d2, d3, d4, d5, d6};
    for (int i = 0; i < dim; i++) {
      int gI = chart->independent_directions.at(i);
      coordinats[gI] = temp_coor[gI];
    }

    /***/
    cout << "independent directions ";
    for (int i = 0; i < dim; i++) {
      cout << chart->independent_directions.at(i) << " ";
    }
    cout << endl;

    cout << "temp coor ";
    for (int i = 0; i < (sizeof(temp_coor) / sizeof(*temp_coor)); i++) {
      cout << temp_coor[i] << " ";
    }
    cout << endl;
    /***/

    cout << "start coordinats ";
    for (int i = 0; i < 6; i++) {
      cout << coordinats[i] << " ";
    }
    cout << endl;

    for (int i = 0; i < 6; i++) {
      chart->start_coordinats[i] = coordinats[i];
    }

    vector<double> dum;
    for (int i = 0; i < dim; i++) dum.push_back(set->Sampling.stepSize);
    findRealizationsJacobianRec(sreal, chart, cayleyPoint_s, flip, dum, 0,
                                chart->start_coordinats, detector, cgK,
                                noPoints, rnode, JacInv, coming_gridDirections);
    // findBoundary(); /// TODO: finish this method's input parameter, and all
    // should be done. 2018/9/14

    if (noPoints != 0) {
      this->snl->appendSpacePoints(rnode);
      noPoints = 0;
    }
  }

  /***/
  // delete visited after using it
  /**
  for (size_t i = 0; i < pr1; i++) {
          this -> visited[i] = new bool ****[pr2];
  }

  for (size_t i = 0; i < pr1; i++) {
          for (size_t j = 0; j < pr2; j++) {
                  this -> visited[i][j] = new bool ***[pr3];
          }
  }

  for (size_t i = 0; i < pr1; i++) {
          for (size_t j = 0; j < pr2; j++) {
                  for (size_t k = 0; k < pr3; k++) {
                          this -> visited[i][j][k] = new bool **[pr4];
                  }
          }
  }

  for (size_t i = 0; i < pr1; i++) {
          for (size_t j = 0; j < pr2; j++) {
                  for (size_t k = 0; k < pr3; k++) {
                          for (size_t l = 0; l < pr4; l++) {
                                  this -> visited[i][j][k][l] = new bool *[pr5];
                          }
                  }
          }
  }
  */
  for (size_t i = 0; i < end_coordinats[0]; i++) {
    for (size_t j = 0; j < end_coordinats[1]; j++) {
      for (size_t k = 0; k < end_coordinats[2]; k++) {
        for (size_t l = 0; l < end_coordinats[3]; l++) {
          for (size_t m = 0; m < end_coordinats[4]; m++) {
            delete[](chart->visited[i][j][k][l][m]);
          }
        }
      }
    }
  }

  for (size_t i = 0; i < end_coordinats[0]; i++) {
    for (size_t j = 0; j < end_coordinats[1]; j++) {
      for (size_t k = 0; k < end_coordinats[2]; k++) {
        for (size_t l = 0; l < end_coordinats[3]; l++) {
          delete[](chart->visited[i][j][k][l]);
        }
      }
    }
  }

  for (size_t i = 0; i < end_coordinats[0]; i++) {
    for (size_t j = 0; j < end_coordinats[1]; j++) {
      for (size_t k = 0; k < end_coordinats[2]; k++) {
        delete[](chart->visited[i][j][k]);
      }
    }
  }

  for (size_t i = 0; i < end_coordinats[0]; i++) {
    for (size_t j = 0; j < end_coordinats[1]; j++) {
      delete[](chart->visited[i][j]);
    }
  }

  for (size_t i = 0; i < end_coordinats[0]; i++) {
    delete[](chart->visited[i]);
  }

  delete[](chart->visited);

  /// int currentSamplingDimension = tempListOfChildNodes.begin() -> first ->
  /// getDim();
  /**
  if(previousSamplingDimension == currentSamplingDimension) {
          return true;
  }
  */

  /***/
  vector<int> list_of_sampled_child_nodes;
  if (set->RootNodeCreation.createChildren) {
    list<pair<AtlasNode*, Orientation*> > tempListOfChildNodes =
        rnode->getListOfChildNodes();
    for (list<pair<AtlasNode*, Orientation*> >::iterator i =
             tempListOfChildNodes.begin();
         i != tempListOfChildNodes.end(); i++) {
      bool outer_continue = false;
      /**
      if(i -> first -> getDim() == 0) {
              continue;
      }
      if(i -> first -> getDim() == 1) {
              // sampleTheNode(i -> first, false, i -> second, continu, bret);
              MyJacobianSampleRec(i -> first, false, i -> second, false, false,
      false); continue;
      }
      */
      /// we don't want to sample the same nodes over and over again
      /// so we include a test condition here
      /// which allows us to continue if the IDs match
      if (list_of_sampled_child_nodes.size() == 0)
        list_of_sampled_child_nodes.push_back(i->first->getID());
      else {
        for (int j = 0; j < list_of_sampled_child_nodes.size(); j++) {
          if (i->first->getID() == list_of_sampled_child_nodes[j]) {
            outer_continue = true;
            break;
          }
        }
        if (!outer_continue)
          list_of_sampled_child_nodes.push_back(i->first->getID());
      }
      if (outer_continue) continue;
      /// include base case here
      if (i->first->getDim() == 0) continue;
      MyJacobianSampleRec(i->first, false, i->second, false, false, false);
    }
    rnode->clearListOfChildNodes();
  }
  /***/

  if (noPoints != 0) {
    this->snl->appendSpacePoints(rnode);
    noPoints = 0;
  }

  rnode->setFoundGoodOrientation(true);

  if (!set->AtlasBuilding.stop && !breath) {
    rnode->setComplete(true);

    region->trim();
  }

  delete detector;
  delete desc;
  delete chart;

  return true;
}
/**
void AlgorithmJacobianRec::calculate_jacobian_in_0D
(
                AtlasNode*   rnode,
                CayleyPoint* cpoint,
                bool         dense
) {
        ActiveConstraintGraph*  cgK      = rnode -> getCG();
        ActiveConstraintRegion* region   = rnode -> getACR();
        CayleyParameterization* desc     = new CayleyParameterization(cgK,
dense);	// create a fresh description ConstraintCheck*        detector = new
ConstraintCheck(cgK, df);			// create a fresh collision
detector ConvexChart*            chart    = new ConvexChart(cgK, dense, desc,
df); Orientation*            ori      = NULL;

        chart -> initializeChart(false, region);

        snl -> saveNodeHeader (rnode, this -> atlas);
        snl -> appendDimension(rnode);

        /// retrieve the position of the parameters in "param_length" variable.
        /// use a nested for loop, traverse "param_length". end for loop as soon
as the first instance is found. double         cayley_param = (chart ->
getPoint())[0]; pair<int, int> cayley_param_pos; for(size_t i = 0; i < (a ->
size() + b -> size()); i++) { bool exit = false; for(size_t j = 0; j < (a ->
size() + b -> size()); j++) { if((chart -> param_length)[i][j] == cayley_param)
{ cayley_param_pos = std::make_pair(i, j); exit = true; break;
                        }
                }
                if(exit)
                        break;
        }

        /// retrieve the upper and lower bound of the cayley parameter.
        /// then, update cayley parameter length at cayley_param_pos. we'll
start with the lower bound.
        /// "param_length" is a symmetric matrix, so update cayley parameter at
2 positions, (x, y) and (y, x). double upper = chart -> get_maxOfMax_first();
        double lower = chart -> get_minOfMin_first();
        double current_value = lower;
        (chart -> param_length)[cayley_param_pos.first][cayley_param_pos.second]
= current_value; (chart ->
param_length)[cayley_param_pos.second][cayley_param_pos.first] = current_value;

        /// TODO: create an array similar to the "visited" array, but of smaller
size
        /// vector of two-tuples; first element is the cayley parameter, second
element is an integer boolean
        /// 0 means not visited, 1 means visited.
        //  is this actually needed? we are going in one direction anyway.
        //  vector<pair<double, int> >
visited.push_back(std::make_pair(current_value, 1));

        /// variables x, y, z, phi, cos(theta), psi for transformation and
rotation vector<double> step_size; for(int i = 0; i < 6; i++) {
                step_size.push_back(Settings::Sampling::gridSteps[i]);
        }
        vector<double> old_trans_rot;
        vector<double> new_trans_rot;
        vector<double> ratio_k;
        vector<double> cayley_data;
        CayleyPoint*   cayley_point;

        bool fail;
        for(int flip = 0; flip < 8; flip++) {
                ori = CartesianRealizer::computeRealization(cgK, chart, flip,
fail, rnode -> getFlipScheme());
                /// do until a feasible realization is found, because some
cayley parameters don't give feasible realization lower =
get_first_parameter_with_feasible_orientation(ori, lower, upper, 64,
cayley_param_pos, cgK, chart, flip, fail, rnode); upper =
get_last_parameter_with_feasible_orientation(ori, lower, upper, 64,
cayley_param_pos, cgK, chart, flip, fail, rnode);
                /// found a feasible realization
                if(current_value != NULL) {
                        chart -> getParamPoint(cayley_data);
                        cayley_point =  new CayleyPoint(cayley_data);
                        cayley_point -> setRealizable(true);
                        cayley_point -> addOrientation(ori);
                } else {
                        cout << "There are no more feasible realizations in the
current node." << endl; return;
                }
                old_trans_rot = translation_rotation_getter(cayley_point, flip,
cgK);	// use chart to obtain the 6 variables, here chart has "lower"
                /// update chart to use the new current_value, then find the
"new_trans_rot" as mid point of lower and upper lower = current_value;
                current_value = (upper + lower) / 2;
                (chart ->
param_length)[cayley_param_pos.first][cayley_param_pos.second] = current_value;
                (chart ->
param_length)[cayley_param_pos.second][cayley_param_pos.first] = current_value;
                ori = CartesianRealizer::computeRealization(cgK, chart, flip,
fail, rnode -> getFlipScheme());
                /// check if it fails. if it does, return. if not, create new
point. if(!fail) { cayley_data.clear(); chart -> getParamPoint(cayley_data);
                        cayley_point = new CayleyPoint(cayley_data);
                        cayley_point -> setRealizable(true);
                        cayley_point -> addOrientation(ori);
                } else
                        return;
                new_trans_rot = translation_rotation_getter(cayley_point, flip,
cgK);	// use chart to obtain the 6 variables, here chart has "current_value"
                k.clear();
// reset k for(int i = 0; i < 6; i++) {
// calculate ks using the old and new translation/rotation
                        ratio_k.push_back((new_trans_rot[i] - old_trans_rot[i])
/ step_size[i]);
                }
                while(true) {
                        while(norm(k) < 1. - 0.1 || norm(k) > 1. + 0.1) {
                                k.clear();
                                if(norm(k) < 1. - 0.1) {
// too close
                                        lower         = current_value;
// updata lower bound
                                        current_value = (lower + upper) / 2;
// binary search (chart ->
param_length)[cayley_param_pos.first][cayley_param_pos.second] = current_value;
                                        (chart ->
param_length)[cayley_param_pos.second][cayley_param_pos.first] = current_value;
                                        ori =
CartesianRealizer::computeRealization(cgK, chart, flip, fail, rnode ->
getFlipScheme()); if(ori != NULL) { cayley_data.clear(); chart ->
getParamPoint(cayley_data); cayley_point = new CayleyPoint(cayley_data);
                                                cayley_point ->
setRealizable(true); cayley_point -> addOrientation(ori);
                                        }
                                        new_trans_rot =
translation_rotation_getter(cayley_point, flip, cgK);		// get new
position for(int i = 0; i < 6; i++) { ratio_k.push_back((new_trans_rot[i] -
old_trans_rot[i]) / step_size[i]);
                                        }
                                } else if(norm(k) > 1. + 0.1) {
// too far upper         = current_value;
// update upper bound
                                        current_value = (lower + upper) / 2;
// binary search (chart ->
param_length)[cayley_param_pos.first][cayley_param_pos.second] = current_value;
                                        (chart ->
param_length)[cayley_param_pos.second][cayley_param_pos.first] = current_value;
                                        ori =
CartesianRealizer::computeRealization(cgK, chart, flip, fail, rnode ->
getFlipScheme()); if(ori != NULL) { cayley_data.clear(); chart ->
getParamPoint(cayley_data); cayley_point = new CayleyPoint(cayley_data);
                                                cayley_point ->
setRealizable(true); cayley_point -> addOrientation(ori);
                                        }
                                        new_trans_rot =
translation_rotation_getter(cayley_point, flip, cgK);		// get new
position for(int i = 0; i < 6; i++) { ratio_k.push_back((new_trans_rot[i] -
old_trans_rot[i]) / step_size[i]);
                                        }
                                }
                        }
                        save(new_trans_rot);
                        for(int i = 0; i < 6; i++) {
                                old_trans_rot[i] = new_trans_rot[i];
                        }
                        lower = current_value;
// re-initialize lower, this is where we ended upper = chart ->
get_maxOfMax_first();				// re-initialize upper (chart ->
param_length)[cayley_param_pos.first][cayley_param_pos.second] = upper; (chart
-> param_length)[cayley_param_pos.second][cayley_param_pos.first] = upper; ori =
CartesianRealizer::computeRealization(cgK, chart, flip, fail, rnode ->
getFlipScheme()); if(ori != NULL) { cayley_data.clear(); chart ->
getParamPoint(cayley_data); cayley_point = new CayleyPoint(cayley_data);
                                cayley_point -> setRealizable(true);
                                cayley_point -> addOrientation(ori);
                        }
                        new_trans_rot =
translation_rotation_getter(cayley_point, flip, cgK); k.clear(); for(int i = 0;
i < 6; i++) { ratio_k.push_back((new_trans_rot[i] - old_trans_rot[i]) /
step_size[i]);
                        }
                        /// when there is not enough distance length left
between current_value and upper,
                        /// terminate the search.
                        if(norm(k) < 1 - 0.1)
                                break;

                        /// else, there's still enough distance
                        current_value = (lower + upper) / 2;
                        (chart ->
param_length)[cayley_param_pos.first][cayley_param_pos.second] = current_value;
                        (chart ->
param_length)[cayley_param_pos.second][cayley_param_pos.first] = current_value;
                        ori = CartesianRealizer::computeRealization(cgK, chart,
flip, fail, rnode -> getFlipScheme()); if(ori != NULL) { cayley_data.clear();
                                chart -> getParamPoint(cayley_data);
                                cayley_point = new CayleyPoint(cayley_data);
                                cayley_point -> setRealizable(true);
                                cayley_point -> addOrientation(ori);
                        }
                        new_trans_rot =
translation_rotation_getter(cayley_point, flip, cgK);
                        /// prepare for the loop to start again.
                        k.clear();
                        for(int i = 0; i < 6; i++) {
                                ratio_k.push_back((new_trans_rot[i] -
old_trans_rot[i]) / step_size[i]);
                        }
                }
        }

        return;
}

double AlgorithmJacobianRec::norm(vector<double> k) {
        double output = 0;
        for(int i = 0; i < k.size(); i++) {
                output += k[i] * k[i];
        }
        return sqrt(output);
}

double AlgorithmJacobianRec::get_first_parameter_with_feasible_orientation(
        Orientation &ori, double &lower, double &upper, int iteration,
pair<double> cayley_param_pos, ActiveConstraintGraph* cgK, ConvexChart* chart,
int flip, bool fail, AtlasNode* rnode) { int counter = 0;
        /// repeat until ori is not null or max iteration is met
        while(ori == NULL && counter < iteration) {
                /// sterics or angle violation, current_value not feasible
                cout << "Cayley parameter length not realizable, proceed to next
point." << endl; lower = lower + (lower + upper) / 64;		// partition the
range into 64 equal distance pieces.
                /// update "param_length"
                (chart ->
param_length)[cayley_param_pos.first][cayley_param_pos.second] = lower; (chart
-> param_length)[cayley_param_pos.second][cayley_param_pos.first] = lower;
                /// get orientation again
                ori = CartesianRealizer::computeRealization(cgK, chart, flip,
fail, rnode -> getFlipScheme()); counter++;
        }
        if(ori == NULL && counter >= iteration)
                return NULL;
        else
                return current_value;
}

double AlgorithmJacobianRec::get_last_parameter_with_feasible_orientation(
        Orientation &ori, double &lower, double &upper, int iteration,
pair<double> cayley_param_pos, ActiveConstraintGraph* cgK, ConvexChart* chart,
int flip, bool fail, AtlasNode* rnode) { int counter = 0;
        /// repeat until ori is not null or max iteration is met
        while(ori == NULL && counter < iteration) {
                /// sterics or angle violation, current_value not feasible
                cout << "Cayley parameter length not realizable, proceed to next
point." << endl; upper = upper - (lower + upper) / 64;		// partition the
range into 64 equal distance pieces.
                /// update "param_length"
                (chart ->
param_length)[cayley_param_pos.first][cayley_param_pos.second] = upper; (chart
-> param_length)[cayley_param_pos.second][cayley_param_pos.first] = upper;
                /// get orientation again
                ori = CartesianRealizer::computeRealization(cgK, chart, flip,
fail, rnode -> getFlipScheme()); counter++;
        }
        if(ori == NULL && counter >= iteration)
                return NULL;
        else
                return current_value;
}

vector<double> AlgorithmJacobianRec::translation_rotation_getter(CayleyPoint*
cp, int flip, ActiveConstraintGraph* cgK) { vector<Atom*> helB = b ->
getAtoms(); vector<double> pos;

        vector<Orientation*> sol = cp -> getOrientations();
        for(size_t a = 0; a < sol.size(); a++) {
                if(sol[a] -> getFlipNum() == flip) {
                        double fb[3][3], tb[3][3];
                        sol[a] -> getFromTo(fb, tb);

                        Vector3d p1(fb[0][0], fb[0][1], fb[0][2]);
                        Vector3d p2(fb[1][0], fb[1][1], fb[1][2]);
                        Vector3d p3(fb[2][0], fb[2][1], fb[2][2]);

                        Vector3d P1(tb[0][0], tb[0][1], tb[0][2]);
                        Vector3d P2(tb[1][0], tb[1][1], tb[1][2]);
                        Vector3d P3(tb[2][0], tb[2][1], tb[2][2]);

                        Vector3d v1 = p2 - p1;
                        Vector3d v2 = p3 - p1;
                        Vector3d v3 = v1.cross(v2);

                        Vector3d V1 = P2 - P1;
                        Vector3d V2 = P3 - P1;
                        Vector3d V3 = V1.cross(V2);

                        Matrix3d m, M, R;
                        m << v1(0), v2(0), v3(0),	v1(1), v2(1), v3(1),
v1(2), v2(2), v3(2); M << V1(0), V2(0), V3(0),	V1(1), V2(1), V3(1),	V1(2),
V2(2), V3(2);

                        R = M * m.inverse(); 		// rotation matrix
                        Vector3d t = P1 - R * p1;	// translation

                        Quaterniond q(R);

                        // compute mean
                        Vector3d mean(0, 0, 0);
                        for(size_t iter = 0; iter < helB.size(); iter++) {
                                double* l = helB[iter] -> getLocation();
                                Vector3d p(l[0], l[1], l[2]);
                                mean += R * p;
                        }
                        mean = mean / helB.size();

                        Vector3d TB = t + mean;
                        Matrix3d RB = quaterniontoRotation(q.w(), q.x(), q.y(),
q.z());

                        Vector3d eaB =  Utils::RotMatToEuler(R);
                        eaB[0] = eaB[0] * 180 / PI;
                        eaB[2] = eaB[2] * 180 / PI;
                        if(cgK -> independent_directions.size() > 0)
                                for(vector<int >::iterator iter = cgK ->
independent_directions.begin(); iter != cgK -> independent_directions.end();
iter++) { if((*iter) < 3) pos.push_back(TB(*iter)); else pos.push_back(eaB(*iter
- 3)); } else { pos.push_back(TB(0));	pos.push_back(TB(1));
pos.push_back(TB(2)); pos.push_back(eaB(0));	pos.push_back(eaB(1));
pos.push_back(eaB(2));
                        }
                }
        }
        return pos;

}

Matrix3d AlgorithmJacobianRec::quaterniontoRotation(double q0, double q1, double
q2, double q3) { Matrix3d R; R <<
        1 - 2*q2*q2 - 2*q3*q3,	2*q1*q2 + 2*q0*q3,		2*q1*q3 -
2*q0*q2, 2*q1*q2 - 2*q0*q3,		1 - 2*q1*q1 - 2*q3*q3,	2*q2*q3 +
2*q0*q1, 2*q1*q3 + 2*q0*q2,		2*q2*q3 - 2*q0*q1,		1 -
2*q1*q1 - 2*q2*q2; return R;
}
*/
/*
void AlgorithmJacobianRec::findRealizationsJacobianRec(CartesianRealizer *
orij_real, ConvexChart* des, CartesianRealizer * real, int flip, vector<double>
diff_point, int coming_direction, int coordinats[6], ConstraintCheck *detector,
ActiveConstraintGraph *cgK, int &noPoints, AtlasNode *rnode)
{
        bool debug = false;
        if(debug) cout << "findRealizationsJacobianRec" << endl;
        if(debug) {for(int i=0; i<6; i++) cout << " " << coordinats[i];    cout
<< endl; }


        bool out = false;
        for(int i=0; i<6; i++)
                if(coordinats[i]<0 || coordinats[i] >= des->end_coordinats[i] )
                        out = true;
        if(out) return;


         if( noPoints >= 1000 ) //100
        {
                this->snl->appendSpacePoints( rnode );
                noPoints=0;
        }

         if( Settings::AtlasBuilding::stop )
                return;



                int dim = des->paramsOrdered.size()/5;

                Vector3d old_eaB = real->eaB;
                Vector3d old_TB = real->TB;

//		MatrixXd gridSteps;
//		gridSteps.setIdentity(dim,dim);
//		for(int g=0; g<dim; g++)
//		{
//			int d = des->independent_directions[g];
//			gridSteps(g,g) = Settings::Sampling::gridSteps[d];
//		}
//		bool succeed;
//		MatrixXd JacInv = computeJacobian( des, real, flip, succeed,
gridSteps);


        bool jac_computed_once = false;
        MatrixXd aValidJacInv;
        MatrixXd JacInv;
        for(int dir=0; dir<2; dir++)
        {
                for(int g=0; g<dim; g++)  //walk along the grid, x,y,z,angles
//no only first 'dim' grid parameters
                {
                        if( Settings::AtlasBuilding::stop )
                                        return;


//			if(g==dim-1 && dir==1)
//					cout << endl;

                        int move=1;
                        if(dir==1) move = -1;
                        int curr_direction = move * (g+1); //added 1 to overcome
nullish after multiplication

                        int gI = des->independent_directions.at(g);
                        double change=0;

                        int new_coordinats[6];
                        for(int i=0; i<6; i++)
                        {
                                new_coordinats[i] = coordinats[i];
                                if(gI==i){ //g==i
                                        new_coordinats[i] += move;
                                        change =
Settings::Sampling::gridSteps[gI] *
(new_coordinats[i]-des->start_coordinats[i]); //des->pr/2
                                }
                        }

//			if( new_coordinats[3] > 28 && new_coordinats[1] < 27 &&
gI==3 && dir==0 )
//				cout << endl;


                        bool out = false;
                        for(int i=0; i<6; i++)
                                if(new_coordinats[i]<0 || new_coordinats[i] >=
des->end_coordinats[i] ) out = true; if(out) continue;



                        if( des->visited[ new_coordinats[0] ][ new_coordinats[1]
][ new_coordinats[2] ][ new_coordinats[3] ][ new_coordinats[4] ] [
new_coordinats[5] ] ) //already visited continue;


                        des->visited[ new_coordinats[0] ][ new_coordinats[1] ][
new_coordinats[2] ][ new_coordinats[3] ][ new_coordinats[4] ] [
new_coordinats[5] ] = true;



                        double orij_step;
                        if(gI<3)
                                orij_step =  orij_real->TB[gI] + change
-real->TB[gI]; else orij_step =  orij_real->eaB[gI-3] + change -real->eaB[gI-3]
;

//			orij_step = orij_step/
Settings::Sampling::gridSteps[gI];


                        MatrixXd gridSteps;
                        gridSteps.setIdentity(dim,dim);
                        for(int gg=0; gg<dim; gg++)
                        {
                                int dd = des->independent_directions[gg];
                                gridSteps(gg,gg) =
move*Settings::Sampling::gridSteps[dd];
                        }
                        gridSteps(g,g) = orij_step;

                        bool succeed;
                        if(!jac_computed_once)
                                JacInv = computeJacobian( des, real, flip,
succeed, g, gridSteps); else JacInv = computeJacobian( des, real, flip, succeed,
g, gridSteps, JacInv); //NOT aValidJacInv, if previous JacInv come closer but
could not succeed, it should continue from previous JacInv. But if previous
JacInv was too bad, it may cause the new Jac to go bad too.

                        if(succeed){
//				AlgorithmJacobianRec::cayleyDirections = JacInv;
//				jac_computed_once = true;
                        }

                        bool hit;  bool negVol;
                        vector<double> prev_point = real->getPoint();
                        CartesianRealizer* relg = binaryJump(JacInv, des, real,
flip, g, hit, orij_step, negVol); if(relg == NULL &&
curr_direction==coming_direction) relg = jumpByPreviousDirection(diff_point,
des, real, flip, g, hit, orij_step, negVol);

                        relg = adaptiveStepSize2JumpRec(orij_real, relg, des, g,
g, flip, new_coordinats, old_TB, old_eaB, orij_step); //false hence=>keep new
adapted version even if g becomes worse    //hit=true => cancel_if_gBecameWorser



//			cout << "5" << endl;
                        if( relg!= NULL  ){ //&& !hit
                                jac_computed_once = true;
                                AlgorithmJacobianRec::cayleyDirections = JacInv;
                                //aValidJacInv = JacInv;

                                        vector<double> out = relg->getPoint();
//should change the point for each real PointMultiD* p4d = new PointMultiD(out);
                                        p4d->tetra( true );

                                        bool myvalid = valid(relg, detector);

                                        if(myvalid)//checks against rest of the
helix
                                        {
                                                Orientation* ori_on_lattice =
relg->getOrienation();//19april2013

                                                if(
!Settings::Collision::checkAngle || !ori_on_lattice->angleViolated()){ // for
test reason p4d->addOrientation(ori_on_lattice);

                                                        vector<double>
diff_point; for(int i=0; i<out.size(); i++)
diff_point.push_back(out[i]-prev_point[i]);
                                                        findRealizationsJacobianRec(orij_real,
des, relg, flip, diff_point, curr_direction, new_coordinats, detector, cgK,
noPoints, rnode); p4d->colorCode = relg->jacColor; //colorcode is set as a
result of recursive call
                                                }
                                                else{
                                                        delete ori_on_lattice;
//orange p4d->angleViolated() = true; p4d->colorCode = "6"; //bad angle
                                                }
                                        }
                                        else
                                                p4d->colorCode = "5"; //REAL
COLLISION


                                        cgK->AddSamplePoint(p4d);
                                        noPoints++;
                                        delete relg;

                        }
                        else{ //volume negative
//				des->isBoundary[g] = true;
                                delete relg;
                                des->visited[ new_coordinats[0] ][
new_coordinats[1] ][ new_coordinats[2] ][ new_coordinats[3] ][ new_coordinats[4]
][ new_coordinats[5] ] = false;

                                if(hit && !negVol)
                                        real->jacColor = "1";
                                else if(negVol && !hit)
                                        real->jacColor = "2";
                                else if(negVol && hit)
                                        real->jacColor = "3";
                                else
                                        real->jacColor = "4"; //couldnot pass
adaptiveSampling. hence bad jacobian sampling.
                        }


                }
        }


        MatrixXd gridSteps;
        gridSteps.setIdentity(dim,dim);
        for(int gg=0; gg<dim; gg++)
        {
                int dd = des->independent_directions[gg];
                gridSteps(gg,gg) = Settings::Sampling::gridSteps[dd];
        }
        bool succeed;
        if(!jac_computed_once)
                JacInv = computeJacobian( des, real, flip, succeed,  -1,
gridSteps); else JacInv = computeJacobian( des, real, flip, succeed,  -1,
gridSteps, JacInv); if(succeed) AlgorithmJacobianRec::cayleyDirections = JacInv;

        findRealizationsJacobianRec_atCorner( orij_real,  des,  real,  flip,
coordinats, detector, cgK, noPoints, 0, JacInv, rnode);
        findRealizationsJacobianRec_atCorner( orij_real,  des,  real,  flip,
coordinats, detector, cgK, noPoints, 1, JacInv, rnode);
        findRealizationsJacobianRec_atCorner( orij_real,  des,  real,  flip,
coordinats, detector, cgK, noPoints, 2, JacInv, rnode);
}
*/

// void AlgorithmJacobianRec::findRealizationsJacobianRec(CartesianRealizer *
// orij_real, ConvexChart* des, CartesianRealizer * real, int flip,
// vector<double> diff_point, int coming_direction, int coordinats[6],
// ConstraintCheck *detector, ActiveConstraintGraph *cgK, int &noPoints,
// AtlasNode *rnode, MatrixXd coming_cayleyDirections)
//{
//	bool debug = false;
//	if(debug) cout << "findRealizationsJacobianRec" << endl;
//	if(debug) {for(int i=0; i<6; i++) cout << " " << coordinats[i];    cout
//<< endl; }
//
//
//	bool out = false;
//	for(int i=0; i<6; i++)
//		if(coordinats[i]<0 || coordinats[i] >= des->end_coordinats[i] )
//			out = true;
//
////	if(coordinats[4]<2  ) //just for test reason to prevent deep recursion
////		out = true;
//
//	if(out) return;
//
//
//	 if( noPoints >= 1000 ) //100
//	{
//		this->snl->appendSpacePoints( rnode );
//		noPoints=0;
//	}
//
//	 if( Settings::AtlasBuilding::stop )
//		return;
//
//
//
//		int dim = des->paramsOrdered.size()/5;
//
//		Vector3d old_eaB = real->eaB;
//		Vector3d old_TB = real->TB;
//
//
//
//	MatrixXd JacInv;
//	for(int dir=0; dir<2; dir++)
//	{
//
//		MatrixXd gridSteps;
//		gridSteps.setIdentity(dim,dim);
//		int move=1;  if(dir==1) move = -1;
////		move = move * Utils::sign(coming_direction); //to walk circular
///with hope less deep recursive
//
//		double orij_step[dim];
//		for(int g=0; g<dim; g++)  //walk along the grid, x,y,z,angles //no
//only first 'dim' grid parameters
//		{
//			int gI = des->independent_directions.at(g);
//			double change = Settings::Sampling::gridSteps[gI] *
//(coordinats[gI] + move -des->start_coordinats[gI]);
//
//			if(gI<3)
//				orij_step[g] =  orij_real->TB[gI] + change
//-real->TB[gI]; 			else 				orij_step[g] =  orij_real->eaB[gI-3] + change
//-real->eaB[gI-3] ;
//
//			gridSteps(g,g) = orij_step[g];
//		}
//
//
//		bool succeed;
//		JacInv = computeJacobian( des, real, flip, succeed, -1, gridSteps,
//coming_cayleyDirections, 0, 1000, coming_cayleyDirections); 		if(succeed)
//AlgorithmJacobianRec::cayleyDirections = JacInv;
////		else JacInv = coming_cayleyDirections;
//
//
//
//		for(int g=0; g<dim; g++)  //walk along the grid, x,y,z,angles //no
//only first 'dim' grid parameters
//		{
//			if( Settings::AtlasBuilding::stop )
//					return;
//
//			int curr_direction = move * (g+1); //added 1 to overcome
//nullish after multiplication
//
//			int gI = des->independent_directions.at(g);
//
//			int new_coordinats[6];
//			for(int i=0; i<6; i++)	new_coordinats[i] =
//coordinats[i]; 			new_coordinats[gI] += move;
//
//
//			bool out = false;
//			for(int i=0; i<6; i++)
//				if(new_coordinats[i]<0 || new_coordinats[i] >=
//des->end_coordinats[i] ) 					out = true;
//
////			if(coordinats[4]<2  ) //just for test reason to prevent deep
///recursion /					out = true;
//
//			if(out) continue;
//
//
//
//			if( des->visited[ new_coordinats[0] ][ new_coordinats[1] ][
//new_coordinats[2] ][ new_coordinats[3] ][ new_coordinats[4] ] [
//new_coordinats[5] ] ) //already visited 				continue;
//
//
//			des->visited[ new_coordinats[0] ][ new_coordinats[1] ][
//new_coordinats[2] ][ new_coordinats[3] ][ new_coordinats[4] ] [
//new_coordinats[5] ] = true;
//
//
//
//			bool hit;  bool negVol;
//			vector<double> prev_point = real->getPoint();
//			CartesianRealizer* relg = binaryJump(JacInv, des, real,
//flip, g, hit, orij_step[g], negVol); 			if(relg == NULL &&
//curr_direction==coming_direction) 				relg = jumpByPreviousDirection(diff_point,
//des, real, flip, g, hit, orij_step[g], negVol);
//
//			if( relg == NULL  )
//			{
//				VectorXd cayleyySteps = VectorXd::Ones(dim) *
//Settings::Sampling::stepSize; 				VectorXd gridStep(dim);
//gridStep(g)=gridSteps(g,g);
//
//				VectorXd JacInv1Dir = computeJacobian1Dir( des, real,
//flip, succeed, -1, gridStep, cayleyySteps, 0, 1000, cayleyySteps);
//				vector<double> jac_diff_point;
//				for(int i=0; i<dim; i++) jac_diff_point.push_back(
//JacInv1Dir(i) ); 				relg = jumpByPreviousDirection(jac_diff_point, des, real,
//flip, g, hit, orij_step[g], negVol);
//			}
//
//
//			if( relg == NULL  )
//			{
//				JacInv = computeJacobian( des, real, flip, succeed, g,
//gridSteps, coming_cayleyDirections, 0, 1000, coming_cayleyDirections); 				relg =
//binaryJump(JacInv, des, real, flip, g, hit, orij_step[g], negVol);
//			}
//
//
//
//			relg = adaptiveStepSize2JumpRec(orij_real, relg, des, g, g,
//flip, new_coordinats, old_TB, old_eaB, orij_step[g]); //false hence=>keep new
//adapted version even if g becomes worse    //hit=true =>
//cancel_if_gBecameWorser
//
//
//
//			if( relg!= NULL  ){
//
//				AlgorithmJacobianRec::cayleyDirections = JacInv;
//
//
//				vector<double> out = relg->getPoint(); //should change
//the point for each real 				PointMultiD* p4d = new PointMultiD(out); 				p4d->tetra(
//true );
//
//				bool myvalid = valid(relg, detector);
//
//				if(myvalid)//checks against rest of the helix
//				{
//					Orientation* ori_on_lattice =
//relg->getOrienation();//19april2013
//
//					if(  !Settings::Collision::checkAngle ||
//!ori_on_lattice->angleViolated()){ // for test reason
//						p4d->addOrientation(ori_on_lattice);
//
//						vector<double> diff_point;
//						for(int i=0; i<out.size(); i++)
//diff_point.push_back(out[i]-prev_point[i]);
//						findRealizationsJacobianRec(orij_real,
//des, relg, flip, diff_point, curr_direction, new_coordinats, detector, cgK,
//noPoints, rnode, JacInv); 						p4d->colorCode = relg->jacColor; //colorcode is set
//as a result of recursive call
//					}
//					else{
//						delete ori_on_lattice; //orange
//						p4d->angleViolated() = true;
//						p4d->colorCode = "6"; //bad
//angle
//					}
//				}
//				else
//					p4d->colorCode = "5"; //REAL COLLISION
//
//
//				cgK->AddSamplePoint(p4d);
//				noPoints++;
//				delete relg;
//
//			}
//			else{ //volume negative
////				des->isBoundary[g] = true;
//				delete relg;
//				des->visited[ new_coordinats[0] ][ new_coordinats[1]
//][ new_coordinats[2] ][ new_coordinats[3] ][ new_coordinats[4] ][
//new_coordinats[5] ] = false;
//
//				if(hit && !negVol)
//					real->jacColor = "1";
//				else if(negVol && !hit)
//					real->jacColor = "2";
//				else if(negVol && hit)
//					real->jacColor = "3";
//				else
//					real->jacColor = "4"; //couldnot pass
//adaptiveSampling. hence bad jacobian sampling.
//			}
//
//
//		}
//	}
//
//
////	MatrixXd gridSteps;
////	gridSteps.setIdentity(dim,dim);
////	for(int gg=0; gg<dim; gg++)
////	{
////		int dd = des->independent_directions[gg];
////		gridSteps(gg,gg) = Settings::Sampling::gridSteps[dd];
////	}
//
//	MatrixXd gridSteps;
//	gridSteps.setIdentity(dim,dim);
//	for(int g=0; g<dim; g++)  //walk along the grid, x,y,z,angles //no only
//first 'dim' grid parameters
//	{
//		int gI = des->independent_directions.at(g);
//		double change = Settings::Sampling::gridSteps[gI] *
//(coordinats[gI] + 1 -des->start_coordinats[gI]);
//
//		if(gI<3)
//			gridSteps(g,g) =  orij_real->TB[gI] + change
//-real->TB[gI]; 		else 			gridSteps(g,g) =  orij_real->eaB[gI-3] + change
//-real->eaB[gI-3] ;
//	}
//
//
//
//	bool succeed;
////	if(!jac_computed_once)
//		JacInv = computeJacobian( des, real, flip, succeed,  -1,
//gridSteps);
////	else
////		JacInv = computeJacobian(rnode, des, real, flip, succeed,  -1,
///gridSteps, JacInv);
//	if(succeed) AlgorithmJacobianRec::cayleyDirections = JacInv;
////	else JacInv = coming_cayleyDirections;
//
//	findRealizationsJacobianRec_atCorner( orij_real,  des,  real,  flip,
//coordinats, detector, cgK, noPoints, 0, JacInv, rnode);
//	findRealizationsJacobianRec_atCorner( orij_real,  des,  real,  flip,
//coordinats, detector, cgK, noPoints, 1, JacInv, rnode);
//	findRealizationsJacobianRec_atCorner( orij_real,  des,  real,  flip,
//coordinats, detector, cgK, noPoints, 2, JacInv, rnode);
//
////	findRealizationsJacobianRec_atCorner_merged( orij_real,  des,  real,
///flip,  coordinats, detector, cgK, noPoints, 0, JacInv, rnode);
//
//}

void AlgorithmJacobianRec::findRealizationsJacobianRec(
    Orientation* orij_real, ConvexChart* des, CayleyPoint* cp_r, int flip,
    vector<double> diff_point, int coming_direction, int coordinats[6],
    ConstraintCheck* detector, ActiveConstraintGraph* cgK, int& noPoints,
    AtlasNode* rnode, MatrixXd coming_cayleyDirections,
    vector<bool> coming_gridDirections) {
  Settings* set = Settings::getInstance();

  bool debug = false;
  if (debug) cout << "findRealizationsJacobianRec" << endl;
  if (debug) {
    for (int i = 0; i < 6; i++) cout << " " << coordinats[i];
    cout << endl;
  }

  ActiveConstraintRegion* region = rnode->getACR();

  Orientation* real = cp_r->getOrientations().front();

  bool out = false;
  for (int i = 0; i < 6; i++)
    if (coordinats[i] < 0 || coordinats[i] >= des->end_coordinats[i]) {
      cout << "index for coordinats is " << i << " and i is " << coordinats[i]
           << endl;
      out = true;
    }
  if (out) return;

  /// this can be changed to the save point frequency variable in settings.ini
  if (noPoints >= 1000) {  // 100
    this->snl->appendSpacePoints(rnode);
    noPoints = 0;
  }

  if (set->AtlasBuilding.stop) return;

  int dim = des->parameters.size();

  Vector3d old_eaB = real->eaB;
  Vector3d old_TB = real->TB;

  MatrixXd initial_coming_cayleyDirections = coming_cayleyDirections;

  MatrixXd JacInv;
  // for(int dir = 0; dir < 2; dir++)
  for (int circle_move = 0; circle_move < 2 * dim; circle_move++) {
    // int circle_direction = circle_move;
    int circle_direction = (coming_direction + circle_move) % (2 * dim);
    int g = circle_direction % dim;
    int move = 1;
    if (circle_direction >= dim) move = -1;
    // int move = 1; if(dir == 1) move = -1;

    // for(int g = 0; g < dim; g++) // walk along the grid, x,y,z,angles // no
    // only first 'dim' grid parameters
    {
      if (set->AtlasBuilding.stop) return;  /// added 2018/9/14

      int gI = des->independent_directions.at(g);
      int curr_direction =
          circle_direction;  // move * (g + 1); // added 1 to overcome nullish
                             // after multiplication

      int new_coordinats[6];
      for (int i = 0; i < 6; i++) new_coordinats[i] = coordinats[i];
      new_coordinats[gI] += move;

      bool out = false;
      for (int i = 0; i < 6; i++)
        if (new_coordinats[i] < 0 ||
            new_coordinats[i] >= des->end_coordinats[i]) {
          cout << "index for coordinats is " << i << " and i is "
               << coordinats[i] << endl;
          out = true;
        }
      if (out) continue;

      /**
      /// new_coordinats can contain negative entries. therefore we must create
      a bijection between
      /// the set of natural numbers and integers
      /// 0	1	2	3	4 (set of natural numbers including 0)
      /// 0  -1	1  -2	2 (set of integers)
      new_coordinats[0] = new_coordinats[0] < 0 ? new_coordinats[0] * -2 - 1 :
      new_coordinats[0] * 2; new_coordinats[1] = new_coordinats[1] < 0 ?
      new_coordinats[1] * -2 - 1 : new_coordinats[1] * 2; new_coordinats[2] =
      new_coordinats[2] < 0 ? new_coordinats[2] * -2 - 1 : new_coordinats[2] *
      2; new_coordinats[3] = new_coordinats[3] < 0 ? new_coordinats[3] * -2 - 1
      : new_coordinats[3] * 2; new_coordinats[4] = new_coordinats[4] < 0 ?
      new_coordinats[4] * -2 - 1 : new_coordinats[4] * 2; new_coordinats[5] =
      new_coordinats[5] < 0 ? new_coordinats[5] * -2 - 1 : new_coordinats[5] *
      2;

      /// use another check
      int pr[6] = {des -> pr1, des -> pr2, des -> pr3, des -> pr4, des -> pr5,
      des -> pr6}; for(int i = 0; i < 6; i++) { if(new_coordinats[i] > pr[i]) {
                      return;
              }
      }
      */
      if (des->visited[new_coordinats[0]][new_coordinats[1]][new_coordinats[2]]
                      [new_coordinats[3]][new_coordinats[4]]
                      [new_coordinats[5]])  // already visited
        continue;

      des->visited[new_coordinats[0]][new_coordinats[1]][new_coordinats[2]]
                  [new_coordinats[3]][new_coordinats[4]][new_coordinats[5]] =
          true;

      double orij_step;
      double change = set->Sampling.gridSteps[gI] *
                      (new_coordinats[gI] - des->start_coordinats[gI]);
      if (gI < 3)
        orij_step = orij_real->TB[gI] + change - real->TB[gI];
      else
        orij_step = orij_real->eaB[gI - 3] + change - real->eaB[gI - 3];

      // orij_step = move * Settings::Sampling::gridSteps[gI]; // to show
      // incliment

      MatrixXd gridSteps;
      gridSteps.setIdentity(dim, dim);
      // for(int gg = 0; gg < dim; gg++)
      // {
      // int dd = des -> independent_directions[gg];
      // gridSteps(gg, gg) = move * Settings::Sampling::gridSteps[dd];
      // }
      // gridSteps(g, g) = orij_step;

      for (int gg = 0; gg < dim;
           gg++) {  // make grid consistent with coming_gridDirections
        int comindir = 1;
        if (!coming_gridDirections[gg]) comindir = -1;
        int dd = des->independent_directions[gg];
        gridSteps(gg, gg) = comindir * set->Sampling.gridSteps[dd];
      }
      gridSteps(g, g) = orij_step;

      coming_cayleyDirections = initial_coming_cayleyDirections;
      if ((coming_gridDirections[g] && move == -1) ||
          (!coming_gridDirections[g] &&
           move == 1))  // if coming and current direction are opposite to each
                        // other than reverse the cayley direction
        for (int hh = 0; hh < dim; hh++) {
          coming_cayleyDirections(hh, g) = -coming_cayleyDirections(hh, g);
        }

      besttCayleyDirections = coming_cayleyDirections;
      bestGridSteps = gridSteps;
      bool succeed;
      JacInv = computeJacobian(rnode, cgK, des, cp_r, flip, succeed, g,
                               gridSteps, infinit, coming_cayleyDirections, 0,
                               coming_cayleyDirections);
      MatrixXd gridStepsLater = bestGridSteps;
      JacInv = besttCayleyDirections;

      vector<bool> temp_coming_gridDirections = coming_gridDirections;
      if (move == 1)
        temp_coming_gridDirections[g] = true;
      else
        temp_coming_gridDirections[g] = false;

      for (int h = 0; h < dim;
           h++)  // make coming_gridDirections consistent with bestGridSteps
                 // since gridSteps may have changed in computeJacobian. SO that
                 // coming_gridDirections will be consistent with JacInv
        if (gridStepsLater(h, h) * gridSteps(h, h) <
            0)  // gridStepsLater(h, h) != gridSteps(h, h)
          temp_coming_gridDirections[h] = !temp_coming_gridDirections[h];

      // if(succeed) AlgorithmJacobianRec::cayleyDirections = JacInv; //
      // commented 23 may, because samples increased else JacInv =
      // coming_cayleyDirections;

      bool hit;
      bool negVol;
      // vector<double> prev_point = real -> getPoint();
      double real_point[6];
      cp_r->getPoint(real_point);

      // CartesianRealizer* relg = binaryJump(JacInv, des, real, flip, g, hit,
      // orij_step, negVol);
      CayleyPoint* cayleyPoint = binaryJump(cgK, JacInv, des, cp_r, flip, g,
                                            hit, orij_step, negVol, rnode);

      // CartesianRealizer* relg = walkByJacobian(des, real -> pconnections, g,
      // flip, JacInv, 1, hit); if(hit) { delete relg; relg = NULL; }

      // cayley step is not the final step!, so you need this.
      if (cayleyPoint == NULL && curr_direction == coming_direction)
        cayleyPoint = jumpByPreviousDirection(cgK, diff_point, des, cp_r, flip,
                                              g, hit, orij_step, negVol, rnode);

      if (cayleyPoint == NULL) {
        VectorXd cayleyySteps = VectorXd::Ones(dim) * set->Sampling.stepSize;
        VectorXd gridStep(dim);
        gridStep(g) = gridSteps(g, g);

        VectorXd JacInv1Dir =
            computeJacobian1Dir(cgK, des, cp_r, flip, succeed, -1, gridStep,
                                cayleyySteps, 0, infinit, cayleyySteps, rnode);
        vector<double> jac_diff_point;
        for (int i = 0; i < dim; i++) jac_diff_point.push_back(JacInv1Dir(i));
        cayleyPoint =
            jumpByPreviousDirection(cgK, jac_diff_point, des, cp_r, flip, g,
                                    hit, orij_step, negVol, rnode);
      }

      // no need to specifically jump by previous CayleySteps because
      // computeJacobian already considers that. but what if current point
      // deteriorates new CayleySteps!
      if (cayleyPoint == NULL)
        cayleyPoint = binaryJump(cgK, coming_cayleyDirections, des, cp_r, flip,
                                 g, hit, orij_step, negVol, rnode);

      bool wasNull = false;
      if (cayleyPoint == NULL) wasNull = true;

      cayleyPoint = adaptiveStepSize2JumpRec(
          rnode, cgK, orij_real, cayleyPoint, des, g, g, flip, new_coordinats,
          old_TB, old_eaB, orij_step, JacInv,
          temp_coming_gridDirections);  // false hence => keep new adapted
                                        // version even if g becomes worse //
                                        // hit = true => cancel_if_gBecameWorser

      if (cayleyPoint == NULL) {
        // des -> isBoundary[g] = true;
        des->visited[new_coordinats[0]][new_coordinats[1]][new_coordinats[2]]
                    [new_coordinats[3]][new_coordinats[4]][new_coordinats[5]] =
            false;

        // commented to clean jaccolor from pointmultid.h
        // if(hit && !negVol)
        // cp_r -> jacColor = "1";
        // else if(negVol && !hit)
        // cp_r -> jacColor = "2";
        // else if(negVol && hit)
        // cp_r -> jacColor = "3";
        // else
        // cp_r -> jacColor = "4"; // could not pass adaptiveSampling. hence bad
        // jacobian sampling.
      }

      bool enteredFarAway = false;
      // not sure to put this "|| cayleyPoint == NULL" part, i believe it should
      // not jump far away when jacobian is bad, then you may not have a
      // complete coverage. Maybe jumpfarAway only if it is initial point?
      // Actually not, it will not hurt completeness, since it is recursive
      // sampling, and it was already stuck if you would not jump far away.
      // Maybe one concern is it will distort uniform sampling a little bit,
      // because of 'jump away' puts the point NOT on exact grid point.
      if (wasNull || cayleyPoint == NULL) {
        cayleyPoint =
            jumpFarAway(cgK, des, real_point, g, flip, coming_cayleyDirections,
                        hit, coordinats, new_coordinats, rnode);
        enteredFarAway = true;
      }

      /**
      /// 2018/10/2, findBoundary algorithm added here
      /// after finding a cayley point, walk around the cayley point. if there
      is a collision, do binary search in between to find the configuration.
      /// if cayley point is not null, then we can find boundary.
      if(Settings::RootNodeCreation::createChildren) {
              /// create a list of orientation pointers, fill it
              list<Orientation*> real;
              real = AtlasBuilder::findRealizations(cgK, des, rnode ->
      getFlipScheme());
              /// create a list iterator
              list<Orientation*>::iterator ori_on_lattice = real.begin();


              /// chart is given
              /// cgK is given
              /// rnode is given
              /// region is given
              /// detector is given
              /// bret is false, because we are not doing breathfirst search
              /// noGoodOriention is passed to
      createChildContactGraps_fromTheBoundary. set this to true for now? bool
      noGoodOrientation = true;
              /// noPoints uhhhh give it a zero
              int noPoints = 0;
              /// create a boolean variable
              bool boundary_ori_found_and_saved = false;
              /// the cayley point entry should be the cayley point found above
              cout << "You are about to enter findBounary" << endl;
              findBoundary(ori_on_lattice, des, cgK, rnode, region, detector,
      false, noGoodOrientation, noPoints, boundary_ori_found_and_saved,
      cayleyPoint);

              /// if i set ifBadAngleWitness_createChild to false, does that
      mean the boundary points won't be created either?

      }
      */

      if (cayleyPoint != NULL) {
        // cout << des ->
        // visited[new_coordinats[0]][new_coordinats[1]][new_coordinats[2]][new_coordinats[3]][new_coordinats[4]][new_coordinats[5]]
        // << " " << new_coordinats[0] << " " << new_coordinats[1] << " " <<
        // new_coordinats[2] << " " << new_coordinats[3] << " " <<
        // new_coordinats[4] << " " << new_coordinats[5] << endl;

        // if(!enteredFarAway)
        // AlgorithmJacobianRec::cayleyDirections = JacInv; // commented 23 may,
        // because samples increased

        // vector<double> out = relg -> getPoint(); // should change the point
        // for each real PointMultiD* p4d = new PointMultiD(out); p4d ->
        // tetra(true);

        double relg_point[6];
        cayleyPoint->getPoint(relg_point);

        Orientation* ori_on_lattice = cayleyPoint->getOrientations().front();
        bool myvalid = !detector->stericsViolated(ori_on_lattice);

        bool outerGrid = detector->outerGrid;

        if (outerGrid) {
          des->visited[new_coordinats[0]][new_coordinats[1]][new_coordinats[2]]
                      [new_coordinats[3]][new_coordinats[4]]
                      [new_coordinats[5]] = false;
          // cp_r -> jacColor = "2"; // commented to clean pointmultid from
          // jaccolor
          cayleyPoint->trim_PointMultiD();
        } else if (myvalid)  // checks against rest of the helix
        {
          // Orientation* ori_on_lattice = relg -> getOrienation(); // 19 april
          // 2013

          if (!ori_on_lattice->angleViolated()) {  // for test reason
            // p4d -> addOrientation(ori_on_lattice);

            vector<double> diff_point;
            int dim = des->parameters.size();
            for (int i = 0; i < dim; i++)
              diff_point.push_back(relg_point[i] - real_point[i]);

            if (enteredFarAway) {
              vector<bool> temp_coming_gridDirections = coming_gridDirections;
              if (move == 1)
                temp_coming_gridDirections[g] = true;
              else
                temp_coming_gridDirections[g] = false;

              findRealizationsJacobianRec(
                  orij_real, des, cayleyPoint, flip, diff_point, curr_direction,
                  new_coordinats, detector, cgK, noPoints, rnode,
                  coming_cayleyDirections, temp_coming_gridDirections);
            } else {
              findRealizationsJacobianRec(
                  orij_real, des, cayleyPoint, flip, diff_point, curr_direction,
                  new_coordinats, detector, cgK, noPoints, rnode, JacInv,
                  temp_coming_gridDirections);
            }
            // cayleyPoint -> colorCode = relg -> jacColor; // colorcode is set
            // as a result of recursive call
          } else {
            // delete ori_on_lattice; // orange
            cayleyPoint->trim_PointMultiD();
            cayleyPoint->incrementBadAngleN();
          }
        } else {
          // Orientation* ori_on_lattice = relg -> getOrienation();
          // p4d -> addOrientation(ori_on_lattice); // to see collision regions
          if (ori_on_lattice->angleViolated())
            cayleyPoint->incrementBadAngleN();
          cayleyPoint->trim_PointMultiD();
          cayleyPoint->incrementCollidN();
        }

        /***/
        bool no_good_o = true;
        list<Orientation*> o;
        o = findRealizations(cgK, des, rnode->getFlipScheme());
        cayleyPoint->setID(region->nextSamplePointID());
        cayleyPoint->setRealizable(!o.empty());
        if (!o.empty()) {
          for (list<Orientation*>::iterator lattice_o = o.begin();
               lattice_o != o.end() && !set->AtlasBuilding.stop; lattice_o++) {
            if (!detector->stericsViolated(*lattice_o)) {
              if ((*lattice_o)->angleViolated()) {
                cayleyPoint->incrementBadAngleN();
              } else {
                no_good_o = false;
              }
              bool boundary_o_found_and_saved = false;
              if (set->RootNodeCreation.createChildren) {
                findBoundary(lattice_o, des, cgK, rnode, region, detector,
                             false, no_good_o, noPoints,
                             boundary_o_found_and_saved, cayleyPoint);
              }
              if (!boundary_o_found_and_saved) {
                if (!(*lattice_o)->angleViolated()) {
                  cayleyPoint->addOrientation((*lattice_o));
                } else {
                  delete (*lattice_o);
                }
              } else {
                delete (*lattice_o);
              }
            } else {
              if ((*lattice_o)->angleViolated()) {
                cayleyPoint->incrementBadAngleN();
              }
              delete (*lattice_o);
              cayleyPoint->incrementCollidN();
            }
          }
        }
        /***/
        region->AddSamplePoint(cayleyPoint);
        noPoints++;
        // delete relg;
      }
      // else { // volume negative
      // des -> isBoundary[g] = true;
      // delete relg;
      // des ->
      // visited[new_coordinats[0]][new_coordinats[1]][new_coordinats[2]][new_coordinats[3]][new_coordinats[4]][new_coordinats[5]]
      // = false;
      //
      // if(hit && !negVol)
      // real -> jacColor = "1";
      // else if(negVol && !hit)
      // real -> jacColor = "2";
      // else if(negVol && hit)
      // real -> jacColor = "3";
      // else
      // real -> jacColor = "4"; // couldnot pass adaptiveSampling. hence bad
      // jacobian sampling.
      // }
    }
  }

  /// added 2018/9/14
  return;  // since jumpFarAway good enough! but still keep corner sampling in
           // case bad jacobian.

  // MatrixXd gridSteps;
  // gridSteps.setIdentity(dim, dim);
  // for(int gg = 0; gg < dim; gg++)
  // {
  // int dd = des -> independent_directions[gg];
  // gridSteps(gg, gg) = Settings::Sampling::gridSteps[dd];
  // }

  MatrixXd gridSteps;
  gridSteps.setIdentity(dim, dim);
  for (int g = 0; g < dim; g++)  // walk along the grid, x,y,z,angles // no only
                                 // first 'dim' grid parameters
  {
    int gI = des->independent_directions.at(g);
    double change = set->Sampling.gridSteps[gI] *
                    (coordinats[gI] + 1 - des->start_coordinats[gI]);

    if (gI < 3)
      gridSteps(g, g) = orij_real->TB[gI] + change - real->TB[gI];
    else
      gridSteps(g, g) = orij_real->eaB[gI - 3] + change - real->eaB[gI - 3];
  }

  bool succeed;
  // if(!jac_computed_once)
  // JacInv = computeJacobian(des, real, flip, succeed, -1, gridSteps);
  JacInv = computeJacobian(rnode, cgK, des, cp_r, flip, succeed, -1, gridSteps,
                           infinit, coming_cayleyDirections, 0,
                           coming_cayleyDirections);

  // else
  // JacInv = computeJacobian(des, real, flip, succeed, -1, gridSteps, JacInv);
  // if(succeed) AlgorithmJacobianRec::cayleyDirections = JacInv; // commented
  // 23 may, because samples increased else JacInv = coming_cayleyDirections;

  findRealizationsJacobianRec_atCorner(orij_real, des, cp_r, flip, coordinats,
                                       detector, cgK, noPoints, 0, JacInv,
                                       rnode);
  findRealizationsJacobianRec_atCorner(orij_real, des, cp_r, flip, coordinats,
                                       detector, cgK, noPoints, 1, JacInv,
                                       rnode);
  findRealizationsJacobianRec_atCorner(orij_real, des, cp_r, flip, coordinats,
                                       detector, cgK, noPoints, 2, JacInv,
                                       rnode);

  // findRealizationsJacobianRec_atCorner_merged(orij_real, des, real, flip,
  // coordinats, detector, cgK, noPoints, 0, JacInv, rnode);
}

// this is for jacobian stepping.
// 1: before computing for the next magnitude and direction, step in the
// previous direction and magnitude first. 2: after computing for the next
// magnitude and direction, step in the negation of that magnitude and direction
// first. return what stericsViolated returns?
/**
bool AlgorithmJacobianRec::step_around_jacobian(int choice, CayleyPoint
cayley_point, MatrixXd magnitude_and_direction) { if(choice == 0) { // walk in
previous magnitude and direction CayleyPoint cayley_point =
walkByJacobian(previous_direction); findRealizations(cayley_point);
                // however, each cayley point can correspond to as many as 8
realizations. so what to do here, if some of them violate sterics constraints?
                for(all the realizations of cayley_point) {
                        if(sterics_violated(one_of_the_realizations) {
                                binary_search;
                                create_child_node;
                        }
                }
        } else if(choice == 1) { // walk in negation of next magnitude and
direction new_direction = compute_next_direction; negation_direction =
-new_direction; CayleyPoint cayley_point = walkByJacobian(negation_direction);
                findRealizations(cayley_point);
                for(all the realizations of cayley_point) {
                        if(sterics_violated(one_of_the_realizations) {
                                binary_search;
                                create_child_node;
                        }
                }
        }
        return false;
}
*/

CayleyPoint* AlgorithmJacobianRec::jumpFarAway(ActiveConstraintGraph* cgK,
                                               ConvexChart* des, double data[6],
                                               int g, int flip, MatrixXd JacInv,
                                               bool& hit, int coordinats[6],
                                               int new_coordinats[6],
                                               AtlasNode* rnode) {
  // return NULL;
  Settings* set = Settings::getInstance();

  //	CartesianRealizer* rel = walkByJacobian( des, paramconnections, g, flip,
  //JacInv, 1, hit);
  CayleyPoint* cayleyPoint =
      walkByJacobian(cgK, des, data, g, flip, JacInv, 1, hit, rnode);

  if (cayleyPoint == NULL) return NULL;

  Orientation* rel = cayleyPoint->getOrientations().front();

  int dim = des->parameters.size();

  int coor[6];

  coor[0] = round((rel->TB[0] + set->Sampling.GridXY) /
                  set->Sampling.gridSteps[0]);  // 19.5
  coor[1] = round((rel->TB[1] + set->Sampling.GridXY) /
                  set->Sampling.gridSteps[1]);  // 19.5
  coor[2] =
      round(rel->TB[2] + set->Sampling.GridZ) / set->Sampling.gridSteps[2];
  coor[3] =
      round((rel->eaB[0] + 180) /
            set->Sampling
                .gridSteps[3]);  // added +1 to prevent floating point errors.
  coor[4] = round(
      (rel->eaB[1] - 0.86602540378) /
      set->Sampling.gridSteps[4]);  // allows you 4 steps, 0.866, 0.91, 0.954, 1
  coor[5] = round((rel->eaB[2] + 180) / set->Sampling.gridSteps[5]);

  coor[0] +=
      (des->pr1 - (2 * set->Sampling.GridXY) / (set->Sampling.gridSteps[0])) /
      2;  // since the grid size increased +2 to prevent rounding issues.
  coor[1] +=
      (des->pr2 - (2 * set->Sampling.GridXY) / (set->Sampling.gridSteps[1])) /
      2;
  coor[2] +=
      (des->pr3 - (2 * set->Sampling.GridZ) / (set->Sampling.gridSteps[2])) / 2;
  coor[3] += (des->pr4 - 36) / 2;
  coor[4] += (des->pr5 - 5) / 2;
  coor[5] += (des->pr6 - 36) / 2;

  //	coor[0] += (des->pr1 - 40)/2;  //since the grid size increased +2 to
  //prevent rounding issues. 	coor[1] += (des->pr2 - 40)/2; 	coor[2] += (des->pr3
  //- 7)/2; 	coor[3] += (des->pr4 - 36)/2; 	coor[4] += (des->pr5 - 3)/2; 	coor[5]
  //+= (des->pr6 - 36)/2;

  //	coor[0]++;
  //	coor[1]++;
  //	coor[2]++;
  //	coor[3]++;
  //	coor[4]++;
  //	coor[5]++;

  //	for(int i=0; i<dim; i++)
  //	{
  //		int gI = des->independent_directions.at(i);
  //		coor[gI] = round( des->start_coordinats[gI] + (coor[gI] -
  //des->start_coordinats[gI])/Settings::Sampling::gridSteps[gI] );
  //	}

  for (int i = 0; i < dim; i++)  // do not touch rest of coordinates, they are
                                 // dependent on the first dim coordinates.
  {
    int gI = des->independent_directions.at(i);
    new_coordinats[gI] = coor[gI];
  }

  /**
                  /// new_coordinats can contain negative entries. therefore we
     must create a bijection between
                  /// the set of natural numbers and integers
                  /// 0	1	2	3	4 (set of natural numbers
     including 0)
                  /// 0  -1	1  -2	2 (set of integers)
                  new_coordinats[0] = new_coordinats[0] < 0 ? new_coordinats[0]
     * -2 - 1 : new_coordinats[0] * 2; new_coordinats[1] = new_coordinats[1] < 0
     ? new_coordinats[1] * -2 - 1 : new_coordinats[1] * 2; new_coordinats[2] =
     new_coordinats[2] < 0 ? new_coordinats[2] * -2 - 1 : new_coordinats[2] * 2;
                  new_coordinats[3] = new_coordinats[3] < 0 ? new_coordinats[3]
     * -2 - 1 : new_coordinats[3] * 2; new_coordinats[4] = new_coordinats[4] < 0
     ? new_coordinats[4] * -2 - 1 : new_coordinats[4] * 2; new_coordinats[5] =
     new_coordinats[5] < 0 ? new_coordinats[5] * -2 - 1 : new_coordinats[5] * 2;

                  /// use another check
                  int pr[6] = {des -> pr1, des -> pr2, des -> pr3, des -> pr4,
     des -> pr5, des -> pr6}; for(int i = 0; i < 6; i++) { if(new_coordinats[i]
     > pr[i]) { return NULL;
                          }
                  }
  */

  bool out = false;
  for (int i = 0; i < 6; i++)
    if (new_coordinats[i] < 0 || new_coordinats[i] >= des->end_coordinats[i]) {
      cout << "index for coordinats is " << i << " and i is " << coordinats[i]
           << endl;
      out = true;
    }
  if (out) {
    cayleyPoint->trim_PointMultiD();
    delete cayleyPoint;
    cayleyPoint = NULL;
    return NULL;
  }

  if (des->visited[new_coordinats[0]][new_coordinats[1]][new_coordinats[2]]
                  [new_coordinats[3]][new_coordinats[4]]
                  [new_coordinats[5]])  // already visited
  {
    cayleyPoint->trim_PointMultiD();
    delete cayleyPoint;
    cayleyPoint = NULL;
    return NULL;
  }

  int coor_diff = 0;
  for (int i = 0; i < dim; i++) {
    int gI = des->independent_directions.at(i);
    coor_diff += abs(new_coordinats[gI] - coordinats[gI]);
  }

  if (coor_diff > 1) {
    des->visited[new_coordinats[0]][new_coordinats[1]][new_coordinats[2]]
                [new_coordinats[3]][new_coordinats[4]][new_coordinats[5]] =
        true;
    return cayleyPoint;
  } else {
    cayleyPoint->trim_PointMultiD();
    delete cayleyPoint;
    cayleyPoint = NULL;
    return NULL;
  }
}

/*
void
AlgorithmJacobianRec::findRealizationsJacobianRec_atCorner_merged(CartesianRealizer
* orij_real, ConvexChart* des, CartesianRealizer * real, int flip, int
coordinats[6], ConstraintCheck *detector, ActiveConstraintGraph *cgK, int
&noPoints, int corner, MatrixXd JacInvb, AtlasNode *rnode)
{
        bool debug = false;
        if(debug) cout << "findRealizationsJacobianRec" << endl;
        if(debug) {for(int i=0; i<6; i++) cout << " " << coordinats[i];    cout
<< endl; }


        bool out = false;
        for(int i=0; i<6; i++)
                if(coordinats[i]<0 || coordinats[i] >= des->end_coordinats[i] )
                        out = true;

//	if(coordinats[4]<2  ) //just for test reason to prevent deep recursion
//			out = true;

        if(out) return;


         if( noPoints >= 1000 ) //100
        {
                this->snl->appendSpacePoints( rnode );
                noPoints=0;
        }

         if( Settings::AtlasBuilding::stop )
                return;



                int dim = des->paramsOrdered.size()/5;

                Vector3d old_eaB = real->eaB;
                Vector3d old_TB = real->TB;



        MatrixXd JacInv;
        for(int dir=0; dir<2; dir++)
        {

                MatrixXd gridSteps;
                gridSteps.setIdentity(dim,dim);
                int move=1;  if(dir==1) move = -1;
                double orij_step[dim];
                for(int g=0; g<dim; g++)  //walk along the grid, x,y,z,angles
//no only first 'dim' grid parameters
                {
                        int gI = des->independent_directions.at(g);
                        double change = Settings::Sampling::gridSteps[gI] *
(coordinats[gI] + move -des->start_coordinats[gI]);

                        if(gI<3)
                                orij_step[g] =  orij_real->TB[gI] + change
-real->TB[gI]; else orij_step[g] =  orij_real->eaB[gI-3] + change
-real->eaB[gI-3] ;

                        gridSteps(g,g) = orij_step[g];

                        int next = (g+dim-1) % dim;
                        gridSteps(g,next) = orij_step[g];
                }


                bool succeed;
                JacInv = computeJacobian( des, real, flip, succeed, -1,
gridSteps); if(succeed) AlgorithmJacobianRec::cayleyDirections = JacInv;



                for(int g=0; g<dim; g++)  //walk along the grid, x,y,z,angles
//no only first 'dim' grid parameters
                {
                        if( Settings::AtlasBuilding::stop )
                                        return;

                        int gI = des->independent_directions.at(g);

                        int new_coordinats[6];
                        for(int i=0; i<6; i++)	new_coordinats[i] =
coordinats[i]; new_coordinats[gI] += move;

                        int gnext = (g+1) % dim;
                        int gInext = des->independent_directions.at(gnext);
                        new_coordinats[gInext] += move;


                        bool out = false;
                        for(int i=0; i<6; i++)
                                if(new_coordinats[i]<0 || new_coordinats[i] >=
des->end_coordinats[i] ) out = true;

//			if(coordinats[4]<2  ) //just for test reason to prevent
deep recursion
//					out = true;

                        if(out) continue;



                        if( des->visited[ new_coordinats[0] ][ new_coordinats[1]
][ new_coordinats[2] ][ new_coordinats[3] ][ new_coordinats[4] ] [
new_coordinats[5] ] ) //already visited continue;





//			bool hit;  bool negVol;
//			CartesianRealizer* relg = binaryJump(JacInv, des, real,
flip, g, hit, orij_step[g], negVol);
//			relg = adaptiveStepSize2JumpRec(orij_real, relg, des, g,
g, flip, new_coordinats, old_TB, old_eaB, orij_step[g]); //false hence=>keep new
adapted version even if g becomes worse    //hit=true => cancel_if_gBecameWorser


                        bool hit;  bool negVol;
                        CartesianRealizer* relg =
binaryJump_Corner_merged(JacInv, des, real, flip, g, gnext, hit, orij_step[g],
orij_step[gnext], negVol); double curr_prev_step =
sqrt(orij_step[g]*orij_step[g] + orij_step[gnext]*orij_step[gnext]); relg =
adaptiveStepSize2JumpRec(orij_real, relg, des, g, gnext, flip, new_coordinats,
old_TB, old_eaB, curr_prev_step, JacInv); //false hence=>keep new adapted
version even if g becomes worse    //hit=true => cancel_if_gBecameWorser



                        if( relg!= NULL  ){

                                des->visited[ new_coordinats[0] ][
new_coordinats[1] ][ new_coordinats[2] ][ new_coordinats[3] ][ new_coordinats[4]
] [ new_coordinats[5] ] = true;

                                AlgorithmJacobianRec::cayleyDirections = JacInv;


                                vector<double> out = relg->getPoint(); //should
change the point for each real PointMultiD* p4d = new PointMultiD(out);
                                p4d->tetra( true );

                                bool myvalid = valid(relg, detector);

                                if(myvalid)//checks against rest of the helix
                                {
                                        Orientation* ori_on_lattice =
relg->getOrienation();//19april2013

                                        if(  !Settings::Collision::checkAngle ||
!ori_on_lattice->angleViolated()){ // for test reason
                                                p4d->addOrientation(ori_on_lattice);
                                                findRealizationsJacobianRec(orij_real,
des, relg, flip, out, 0, new_coordinats, detector, cgK, noPoints, rnode,
JacInv); p4d->colorCode = relg->jacColor; //colorcode is set as a result of
recursive call
                                        }
                                        else{
                                                delete ori_on_lattice; //orange
                                                p4d->angleViolated() = true;
                                                p4d->colorCode = "6"; //bad
angle
                                        }
                                }
                                else
                                        p4d->colorCode = "5"; //REAL COLLISION


                                cgK->AddSamplePoint(p4d);
                                noPoints++;
                                delete relg;


                        }
                        else{ //volume negative
//				des->isBoundary[g] = true;
                                delete relg;

                                if(hit && !negVol)
                                        real->jacColor = "1";
                                else if(negVol && !hit)
                                        real->jacColor = "2";
                                else if(negVol && hit)
                                        real->jacColor = "3";
                                else
                                        real->jacColor = "4"; //couldnot pass
adaptiveSampling. hence bad jacobian sampling.
                        }


                }
        }


}
*/

/*
void AlgorithmJacobianRec::findRealizationsJacobianRec_atCorner(Orientation *
orij_real, ConvexChart* des,  PointMultiD * cp_r , int flip, int coordinats[6],
ConstraintCheck *detector, ActiveConstraintGraph *cgK, int &noPoints, int
corner, MatrixXd JacInv , AtlasNode *rnode){

        bool out = false;
        for(int i=0; i<6; i++)
                if(coordinats[i]<0 || coordinats[i] >= des->end_coordinats[i] )
                        out = true;
        if(out) return;

        int dim = des->paramsOrdered.size()/5;


        Orientation * real = cp_r->getOrientations().front();

        Vector3d old_eaB = real->eaB;
        Vector3d old_TB = real->TB;


        bool was_stuck;
        double previous_step;
        int prev_move=1;
        int prev_g=0;
        for(int dir=0; dir<2; dir++)
        {
                for(int g=0; g<dim; g++)  //walk along the grid, x,y,z,angles
//no only first 'dim' grid parameters
                {
                        if( Settings::AtlasBuilding::stop )
                                        return;

                        int move=1;
                        if(dir==1) move = -1;

                        int gI = des->independent_directions.at(g);
                        double change;

                        int new_coordinats[6];
                        for(int i=0; i<6; i++)
                        {
                                new_coordinats[i] = coordinats[i];
                                if(gI==i){
                                        new_coordinats[i] += move;
                                        change =
Settings::Sampling::gridSteps[gI] *
(new_coordinats[i]-des->start_coordinats[i]);
                                }
                        }

                        out = false;
                        for(int i=0; i<6; i++)
                                if(new_coordinats[i]<0 || new_coordinats[i] >=
des->end_coordinats[i] ) out = true; if(out) {was_stuck=false; continue;}


                        bool stuck_now = true;
                        if( des->visited[ new_coordinats[0] ][ new_coordinats[1]
][ new_coordinats[2] ][ new_coordinats[3] ][ new_coordinats[4] ][
new_coordinats[5] ]  ) //already visited stuck_now = false;


                        double curr_step;
                        if(gI<3)
                                curr_step =  orij_real->TB[gI] + change
-real->TB[gI]; else curr_step =  orij_real->eaB[gI-3] + change -real->eaB[gI-3]
;

//			curr_step = move*Settings::Sampling::gridSteps[gI]; //to
show incliment

                        if(dir==0 && g==0)
                        {
                                was_stuck = stuck_now;
                                previous_step = curr_step;
                                prev_move = move;
                                prev_g = g;
                                continue;  //skip the first direction
                        }


                        int middle_coordinats[6];
                        for(int i=0; i<6; i++)
                        {
                                middle_coordinats[i] = coordinats[i];
                                int prev_gI =
des->independent_directions.at(prev_g); if(gI==i) middle_coordinats[i] += move;
                                if( prev_gI==i)
                                        middle_coordinats[i] += prev_move;
                        }


                        out = false;
                        for(int i=0; i<6; i++)
                                if(middle_coordinats[i]<0 ||
middle_coordinats[i] >= des->end_coordinats[i] ) out = true; if(out) continue;


                        if(!out && was_stuck && stuck_now &&  !des->visited[
middle_coordinats[0] ][ middle_coordinats[1] ][ middle_coordinats[2] ][
middle_coordinats[3] ][ middle_coordinats[4] ][ middle_coordinats[5] ] )
                        {

//				if(true) cout << "real TB_EAB " << real->TB[0]
<< " " << real->TB[1] << " " << real->TB[2] << " " << real->eaB[0] << " " <<
real->eaB[1] << " " << real->eaB[2] << " " << endl;

                                double curr_step_corner=curr_step,
previous_step_corner=previous_step; if(corner == 1) curr_step_corner = curr_step
- move*Settings::Sampling::gridSteps[gI]/2; else if(corner == 2){ int prev_gI =
des->independent_directions.at(prev_g); previous_step_corner = previous_step -
prev_move*Settings::Sampling::gridSteps[prev_gI]/2;
                                }



//				MatrixXd gridSteps;
//				gridSteps.setIdentity(dim,dim);
//				for(int gg=0; gg<dim; gg++)
//				{
//					int dd =
des->independent_directions[gg];
//					gridSteps(gg,gg) =
Settings::Sampling::gridSteps[dd];
//				}
//				gridSteps(prev_g,prev_g) = previous_step_corner;
//				gridSteps(g,g) = curr_step_corner;
//				bool succeed;
//				MatrixXd JacInv = computeJacobian( des, real,
flip, succeed, gridSteps);




                                bool hit;  bool negVol;
                                PointMultiD* p_relg = binaryJump_Corner(JacInv,
des, cp_r, flip, prev_g, g, hit, previous_step_corner, curr_step_corner,
negVol); double curr_prev_step = sqrt(curr_step_corner*curr_step_corner +
previous_step_corner*previous_step_corner); p_relg =
adaptiveStepSize2JumpRec(orij_real, p_relg, des, g, prev_g, flip,
middle_coordinats, old_TB, old_eaB, curr_prev_step, JacInv); //false hence=>keep
new adapted version even if g becomes worse    //hit=true =>
cancel_if_gBecameWorser



                                if( p_relg!= NULL  )
                                {
                                        des->visited[ middle_coordinats[0] ][
middle_coordinats[1] ][ middle_coordinats[2] ][ middle_coordinats[3] ][
middle_coordinats[4] ][ middle_coordinats[5] ] = true;

                                        double out[6];
                                        p_relg->getPoint(out);
                                        vector<double> jac_diff_point;
                                        for(int i=0; i<dim; i++)
jac_diff_point.push_back( out[i] );

//					vector<double> out = relg->getPoint();
//should change the point for each real
//					PointMultiD* p4d = new PointMultiD(out);
//					p4d->tetra( true );

                                        Orientation * ori_on_lattice =
p_relg->getOrientations().front(); bool myvalid = valid(ori_on_lattice,
detector); bool outerGrid = detector->outerGrid;

                                        if(outerGrid)
                                        {
                                                des->visited[
middle_coordinats[0] ][ middle_coordinats[1] ][ middle_coordinats[2] ][
middle_coordinats[3] ][ middle_coordinats[4] ][ middle_coordinats[5] ] = false;
                                                cp_r->jacColor = "2";
                                                p_relg->trim_PointMultiD();
                                        }
                                        else
                                        if(myvalid)//checks against rest of the
helix
                                        {
//						Orientation* ori_on_lattice =
relg->getOrienation();//19april2013

                                                if(
!Settings::Collision::checkAngle || !ori_on_lattice->angleViolated()){ // for
test reason
//
p4d->addOrientation(ori_on_lattice); findRealizationsJacobianRec(orij_real, des,
p_relg, flip, jac_diff_point, 0, middle_coordinats, detector, cgK, noPoints,
rnode, JacInv);
//							p4d->colorCode =
relg->jacColor; //colorcode is set as a result of recursive call
                                                }
                                                else{
//							delete ori_on_lattice;
//orange p_relg->trim_PointMultiD();

                                                        p_relg->angleViolated()
= true; p_relg->colorCode = "6"; //bad angle
                                                }
                                        }
                                        else{
                                                p_relg->colorCode = "5"; //REAL
COLLISION p_relg->trim_PointMultiD();
                                        }

                                        cgK->AddSamplePoint(p_relg);
                                        noPoints++;


//					delete relg;

                                }
                                else{ //volume negative
                                        delete p_relg;
                                }


                        }

                        was_stuck = stuck_now;
                        previous_step = curr_step;
                        prev_move = move;
                        prev_g = g;
                }
        }


}
*/

// make sure that there is no change in directions other than g.
CayleyPoint* AlgorithmJacobianRec::adaptiveStepSize2JumpRec(
    AtlasNode* rnode, ActiveConstraintGraph* cgK, Orientation* orij_real,
    CayleyPoint* cp, ConvexChart* des, int g, int g2, int flip,
    int coordinats[6], Vector3d old_TB, Vector3d old_eaB,
    double grid_step_expected, MatrixXd comingJacInv,
    vector<bool> coming_gridDirections, bool cancel_if_gBecameWorser) {
  Settings* set = Settings::getInstance();
  // return relg;
  bool debug = false;

  if (cp == NULL) return NULL;

  Orientation* relg;  // = cp->getOrientations().front();

  double change[6] = {1, 1, 1, 1, 1, 1};
  int gI = des->independent_directions.at(g);
  int g2I = des->independent_directions.at(g2);

  int dim = des->parameters.size();

  for (int r = 0; r < dim; r++) {
    relg = cp->getOrientations().front();

    if (r == g || r == g2)  // fix stepping in other directions (caused by
                            // sampling along direction g)
      continue;

    int rI = des->independent_directions.at(r);
    //			if(rI<3){
    //				change[r] =  (relg->TB[rI] - old_TB[rI])  ;
    //				if(debug && AlgorithmJacobianRec::aDebug)  cout << change[r] << "
    //= " << relg->TB[rI] << " - " << old_TB[rI] << endl;
    //			}
    //			else
    //				change[r] =  (relg->eaB[rI-3] - old_eaB[rI-3]) ; //SHOULD DIVIDE
    //TO GRID STEPSIZE?

    double expected_change = set->Sampling.gridSteps[rI] *
                             (coordinats[rI] - des->start_coordinats[rI]);
    if (rI < 3)
      change[r] = relg->TB[rI] - (orij_real->TB[rI] + expected_change);
    else
      change[r] =
          relg->eaB[rI - 3] - (orij_real->eaB[rI - 3] + expected_change);

    double ratio =
        change[r] / set->Sampling.gridSteps[rI];  // grid_step_expected;
    ratio = abs(ratio);
    if (abs(ratio) > 0.2)  // 0.2
    // if(ratio >= 1.5 || (ratio <= 0.5  ) ) //&& ratio != 0  // ratio nun dusuk
    // olmasi yanlis jac hesabindan deil, contak vb boundar carpmasi sonucu
    // stepin degistirilmesi olabilir
    //			if( abs(change[r]) > 0.1 )
    // if( change[r] != 0 )
    {
      // if(true)  cout << "adaptive 2 change " << r << " " << change[r] << " "
      // ;
      //				CartesianRealizer* reld =
      //walkByJacobian( des,  r,  flip, JacInv, ratio[r]); 				reld =
      //adaptiveStepSize(relg, des, r, flip, JacInv, old_TB, old_eaB, ratio[r]);
      //
      //				delete relg; relg = NULL;
      //				double change;
      //				if(r<3)
      //					change =  (reld->TB[r] -
      //old_TB[r]) ; 				else 					change =  (reld->eaB[r-3] - old_eaB[r-3]) ; 				if(reld !=
      //NULL) 				{delete reld; reld = NULL; }

      change[r] = -change[r];

      MatrixXd gridSteps;
      gridSteps.setIdentity(dim, dim);

      //			for(int g=0; g<dim; g++)
      //			{
      //				int d = des->independent_directions[g];
      //				gridSteps(g,g) =
      //Settings::Sampling::gridSteps[d]  * ratio;
      //			}
      //			gridSteps(r,r) = change[r];

      for (int gg = 0; gg < dim;
           gg++)  // make grid consistent with coming_gridDirections
      {
        int comindir = 1;
        if (!coming_gridDirections[gg]) comindir = -1;
        int dd = des->independent_directions[gg];
        gridSteps(gg, gg) = comindir * set->Sampling.gridSteps[dd] * ratio;
      }
      gridSteps(r, r) = change[r];
      int move = 1;
      if (change[r] < 0) move = -1;

      ////				MatrixXd gridSteps;
      ////				gridSteps.setIdentity(dim,dim);
      ////				for(int g=0; g<dim; g++)  //walk along the
      ///grid, x,y,z,angles //no only first 'dim' grid parameters /
      ///{ /					int d =
      ///des->independent_directions.at(g);
      ////					double changed =
      ///Settings::Sampling::gridSteps[d] * (coordinats[d] + 1
      ///-des->start_coordinats[d]);
      ////
      ////					if(d<3)
      ////						gridSteps(g,g) =
      ///orij_real->TB[d] + changed    -relg->TB[d]; /
      ///else
      ////						gridSteps(g,g) =
      ///orij_real->eaB[d-3] + changed -relg->eaB[d-3] ; /
      ///} /				gridSteps(r,r) = change[r];

      MatrixXd coming_cayleyDirections = comingJacInv;
      if ((coming_gridDirections[r] && move == -1) ||
          (!coming_gridDirections[r] &&
           move == 1))  // if coming and current direction are opposite to each
                        // other than reverse the cayley direction
        for (int hh = 0; hh < dim; hh++) {
          coming_cayleyDirections(hh, r) = -coming_cayleyDirections(hh, r);
        }

      bool tempBound;
      bool negVol;
      bool succeed;
      bool hit;

      besttCayleyDirections = coming_cayleyDirections;
      MatrixXd JacInv = computeJacobian(
          rnode, cgK, des, cp, flip, succeed, r, gridSteps, infinit,
          coming_cayleyDirections, 0, coming_cayleyDirections);
      JacInv = besttCayleyDirections;

      //			MatrixXd JacInv = comingJacInv;

      // MatrixXd JacInv = computeJacobian( des, relg, flip, succeed, r,
      // gridSteps);

      //				CartesianRealizer* reld =
      //binaryJump(JacInv, des, relg, flip, r, tempBound, change[r], negVol);
      CayleyPoint* cayleyPoint = binaryJump(
          cgK, JacInv, des, cp, flip, r, tempBound, change[r], negVol, rnode);

      //				if( reld == NULL  )
      ////computeJacobian1Dir is not that accurate. it may make it worser than
      //adapting to its correct place.
      //				{
      //					VectorXd cayleyySteps =
      //VectorXd::Ones(dim) * Settings::Sampling::stepSize; 					VectorXd
      //gridStep(dim); gridStep(r)=gridSteps(r,r);
      //
      //					VectorXd JacInv1Dir =
      //computeJacobian1Dir( des, relg, flip, succeed, -1, gridStep,
      //cayleyySteps, 0, 1000, cayleyySteps); 					vector<double> jac_diff_point;
      //					for(int i=0; i<dim; i++)
      //jac_diff_point.push_back( JacInv1Dir(i) ); 					reld =
      //jumpByPreviousDirection(jac_diff_point, des, relg, flip, r, hit,
      //change[r], negVol);
      //				}

      if (cayleyPoint != NULL)  // WE SHOULD NOT RETURN NOW, SHOULD RETURN AFTER
                                // THE FOR LOOP DONE.
      {
        cp->trim_PointMultiD();
        delete cp;
        cp = NULL;
        // return reld;
        cp = cayleyPoint;
      }
    }
  }

  //		bool gBecameWorser = false; //TRUE
  //
  //		double chan;
  //		if(gI<3)
  //			chan =  (relg->TB[gI] - (old_TB[gI]+grid_step_expected) )
  //; 		else 			chan =  (relg->eaB[gI-3] - (old_eaB[gI-3]+grid_step_expected) ) ;
  //
  //		if(abs(chan) > 0.1  ) //while making that there is no change in
  //directions other than g, WE SHOULD keep g as it is, otherwise cancel the job
  //you done. 			gBecameWorser = true;     //added 0.1, to allow adaptive2jump to
  //be executed, since it is beneficial in other directions than g.

  //		bool other_dims_still_changed = false; //after trying to fix later
  //dimensions, the first ones become corrupted again 		for(int r=0; r<dim; r++)
  //		{
  //			if(r == g) //fix stepping in other directions (caused by
  //sampling along direction g) 				continue;
  //
  //			double chan;
  //			int rI = des->independent_directions.at(r);
  //			if(rI<3)
  //				chan =  (relg->TB[rI] - old_TB[rI])  ;
  //			else
  //				chan =  (relg->eaB[rI-3] - old_eaB[rI-3]) ; //SHOULD DIVIDE
  //TO GRID STEPSIZE?
  //
  //			double ratio = chan / grid_step_expected;
  //			if( abs(ratio) > 0.1  )
  //				other_dims_still_changed = true;
  //		}
  //
  //		if(other_dims_still_changed ){
  //			delete relg;
  //			return NULL;
  //		}

  bool other_dims_still_changed =
      false;  // after trying to fix later dimensions, the first ones become
              // corrupted again
  for (int r = 0; r < dim; r++) {
    if (r == g || r == g2)  // fix stepping in other directions (caused by
                            // sampling along direction g)
      continue;

    int rI = des->independent_directions.at(r);
    double expected_change = set->Sampling.gridSteps[rI] *
                             (coordinats[rI] - des->start_coordinats[rI]);
    double chan;
    if (rI < 3)
      chan = orij_real->TB[rI] + expected_change - relg->TB[rI];
    else
      chan = orij_real->eaB[rI - 3] + expected_change - relg->eaB[rI - 3];

    double ratio = chan / set->Sampling.gridSteps[rI];
    //			if( abs(ratio) > 0.9  )
    //			if( abs(ratio) > 0.3  )
    if (abs(ratio) > 0.2) other_dims_still_changed = true;
  }

  if (other_dims_still_changed) {
    cp->trim_PointMultiD();
    delete cp;
    return NULL;
  }

  double ratio;
  double chan;
  if (gI < 3)
    chan = relg->TB[gI] - old_TB[gI];
  else
    chan = relg->eaB[gI - 3] - old_eaB[gI - 3];

  if (g == g2)
    ratio = chan / grid_step_expected;
  else {
    double change0 = 1;
    if (g2I < 3)
      change0 = (relg->TB[g2I] - old_TB[g2I]);  // abs
    else
      change0 = (relg->eaB[g2I - 3] - old_eaB[g2I - 3]);  // abs

    double change = sqrt(change0 * change0 + chan * chan);
    ratio = change / grid_step_expected;
  }

  bool gBecameWorser = true;
  //		if( ratio>0.1 )
  if (ratio < 1.2 && ratio > 0.8)
    //		if( ratio<1.3 && ratio>0.7 )
    //		if( ratio<1.5 && ratio>0.5 )
    gBecameWorser = false;

  if (gBecameWorser) {
    cp->trim_PointMultiD();
    delete cp;
    return NULL;
  }

  return cp;
}

/*
//make sure that there is no change in directions other than g.
CartesianRealizer*
AlgorithmJacobianRec::adaptiveStepSize2JumpRecCombined(CartesianRealizer *
orij_real, CartesianRealizer* relg, ConvexChart* des, int g, int g2, int flip,
int coordinats[6], Vector3d old_TB, Vector3d old_eaB, double grid_step_expected,
MatrixXd JacInv, bool cancel_if_gBecameWorser)
{
        //return relg;
        bool debug = false;

        if( relg!= NULL  )
        {

                double change[6] = {0,0,0,0,0,0};
                double ratio[6] = {0,0,0,0,0,0};
                double sumRatio = 0;
                int gI = des->independent_directions.at(g);
                int g2I = des->independent_directions.at(g2);

                int dim = des->paramsOrdered.size()/5;
                for(int r=0; r<dim; r++)
                {
                        if(r == g || r == g2) //fix stepping in other directions
(caused by sampling along direction g) continue;

                        int rI = des->independent_directions.at(r);


                        double expected_change =
Settings::Sampling::gridSteps[rI] * (coordinats[rI]-des->start_coordinats[rI]);
                        if(rI<3)
                                change[r] =  relg->TB[rI]  -  (orij_real->TB[rI]
+ expected_change)    ; else change[r] =  relg->eaB[rI-3]- (orij_real->eaB[rI-3]
+ expected_change)  ;


                        ratio[r] = change[r] /
Settings::Sampling::gridSteps[rI]; //grid_step_expected; sumRatio +=
abs(ratio[r]);
                }

                if( sumRatio > 0.2*(dim-1)  )
                {

                        MatrixXd gridSteps;
                        gridSteps.setIdentity(dim,dim);
                        for(int g=0; g<dim; g++)
                        {
                                int d = des->independent_directions[g];
                                gridSteps(g,g) = Utils::sign(-change[g]) *
Settings::Sampling::gridSteps[d];
                        }

//			for(int g=0; g<dim; g++)
//				gridSteps(g,0) = -change[g];  //make the first
column to be combined changes that you want to reset


                        VectorXd vvv(dim);
                        for(int i=0; i<dim; i++) vvv(i) = -ratio[i];
//-change[i];

                        bool tempBound; bool negVol; bool succeed; bool hit;
                        //MatrixXd JacInv = computeJacobian( des, relg, flip,
succeed, g, gridSteps);

                        VectorXd rrr = JacInv * vvv;
                        vector<double> diff_point;
                        for(int i=0; i<dim; i++) diff_point.push_back( rrr(i) );
                        double multiplier_of_diffpoint=1;
                        CartesianRealizer* reld = walkByPreviousDirection( des,
relg->pconnections, g,  flip, diff_point, multiplier_of_diffpoint, hit);


//			MatrixXd JacInv = computeJacobian( des, relg, flip,
succeed, 0, gridSteps);
//			CartesianRealizer* reld = walkByJacobian( des,
relg->pconnections, 0,  flip, JacInv, 1, hit);
//			CartesianRealizer* reld =  binaryJump(JacInv, des, relg,
flip, 0, tempBound, change[r], negVol);

                        if( reld == NULL  ){
                                delete reld;
                        }
                        else //WE SHOULD NOT RETURN NOW, SHOULD RETURN AFTER THE
FOR LOOP DONE.
                        {
                                delete relg; relg = NULL;
                                //return reld;
                                relg = reld;
                        }

                }



                bool other_dims_still_changed = false; //after trying to fix
later dimensions, the first ones become corrupted again for(int r=0; r<dim; r++)
                {
                        if(r == g || r == g2) //fix stepping in other directions
(caused by sampling along direction g) continue;

                        int rI = des->independent_directions.at(r);
                        double expected_change =
Settings::Sampling::gridSteps[rI] * (coordinats[rI]-des->start_coordinats[rI]);
                        double chan;
                        if(rI<3)
                                chan =  orij_real->TB[rI] + expected_change
-relg->TB[rI]; else chan =  orij_real->eaB[rI-3] + expected_change
-relg->eaB[rI-3] ;

                        double ratio = chan / Settings::Sampling::gridSteps[rI];
//			if( abs(ratio) > 0.5  )
//			if( abs(ratio) > 0.3  )
                        if( abs(ratio) > 0.2  )
                                other_dims_still_changed = true;
                }

                if(other_dims_still_changed ){
                        delete relg;
                        return NULL;
                }


                double ratioo;
                double chan;
                if(gI<3)
                        chan =  relg->TB[gI] - old_TB[gI];
                else
                        chan =  relg->eaB[gI-3] - old_eaB[gI-3] ;


                if(g == g2)
                        ratioo = chan / grid_step_expected;
                else
                {
                        double change0 = 1;
                        if(g2I<3)
                                change0 =  (relg->TB[g2I]-old_TB[g2I])  ;  //abs
                        else
                                change0 =  (relg->eaB[g2I-3]-old_eaB[g2I-3]) ;
//abs

                        double change = sqrt(change0*change0 + chan*chan);
                        ratioo = change / grid_step_expected;
                }


                bool gBecameWorser = true;
                if( ratioo<1.2 && ratioo>0.8 )
//		if( ratio<1.3 && ratio>0.7 )
//		if( ratio<1.5 && ratio>0.5 )
                        gBecameWorser = false;


                if( gBecameWorser ){
                        delete relg;
                        return NULL;
                }
        }
        return relg;

}
*/

CayleyPoint* AlgorithmJacobianRec::walkByPreviousDirection(
    ActiveConstraintGraph* cgK, ConvexChart* des, double out[6], int g,
    int flip, vector<double> diff_point, double& gridstep, bool& hit,
    AtlasNode* rnode) {
  Settings* set = Settings::getInstance();
  bool debug = false;
  if (debug && AlgorithmJacobianRec::aDebug)
    cout << "_____walkByByPreviousDirection on grid " << g << endl;
  double pconnections[12][12];

  for (int i = 0; i < 12; i++)
    for (int j = 0; j < 12; j++)
      pconnections[i][j] =
          des->param_length[i][j];  // des->paramconnections[i][j];

  //	double minGridFactor = 1;
  vector<double> cayleyData;
  hit = false;
  for (size_t p = 0; p < des->parameters.size(); p++) {
    int v1 = des->parameters[p].first;
    int v2 = des->parameters[p].second;

    if (debug && AlgorithmJacobianRec::aDebug)
      cout << v1 << "-" << v2 << " " << pconnections[v1][v2];
    // pconnections[v1][v2] += des->direction[g] * Settings::Sampling::stepSize
    // * (diff_TB_eaB_list[p/5][g] /total_TB_eaB[g]) * (gridSteps[g]
    // /total_TB_eaB[g]);  //(gridSteps[g] /total_TB_eaB[g]) is for having same
    // step size with grid steps

    pconnections[v1][v2] = out[p] + diff_point[p] * gridstep;

    int isContact = des->find_pair_in_vector(des->contacts, des->parameters[p]);

    Atom* atomA = des->getAtom(v1);
    Atom* atomB = des->getAtom(v2);
    double bondingLowerBound = df->bondingLowerBound(atomA, atomB);
    double bondingUpperBound = df->bondingUpperBound(atomA, atomB);
    double collisionLowerBound = df->collisionLowerBound(atomA, atomB);

    // free cayley params
    if (isContact < 0 &&
        pconnections[v1][v2] <
            collisionLowerBound) {  // tetrahedral boundary is handled by volume
                                    // test, but collison should be handled
                                    // seperately
      pconnections[v1][v2] = collisionLowerBound;
      //			isBoundary[g] = true;
      if (debug && AlgorithmJacobianRec::aDebug) cout << "f";
      hit = true;

      //			double changeFactor =
      //abs(paramconnections[v1][v2] - pconnections[v1][v2]) / abs(JacInv(p/5,g)
      //* gridstep); 			if(changeFactor < minGridFactor) minGridFactor =
      //changeFactor;
    }

    if (isContact >= 0 && set->Sampling.short_range_sampling) {
      if (pconnections[v1][v2] <
          bondingLowerBound) {  // tetrahedral boundary is handled by volume
                                // test, but collison should be handled
                                // seperately
        pconnections[v1][v2] = bondingLowerBound;
        //			isBoundary[g] = true;
        if (debug && AlgorithmJacobianRec::aDebug) cout << "c";
        hit = true;
      } else if (pconnections[v1][v2] > bondingUpperBound) {
        pconnections[v1][v2] = bondingUpperBound;
        //			isBoundary[g] = true;
        if (debug && AlgorithmJacobianRec::aDebug) cout << "b";
        hit = true;
      }
    }

    pconnections[v2][v1] = pconnections[v1][v2];
    if (debug && AlgorithmJacobianRec::aDebug)
      cout << " " << pconnections[v1][v2] << ", ";

    cayleyData.push_back(pconnections[v1][v2]);
  }

  //	CartesianRealizer *relg = new CartesianRealizer(des, flip,
  //pconnections);
  bool fail;
  Orientation* relg = CartesianRealizer::computeRealization(
      cgK, des, flip, pconnections, fail, rnode->getFlipScheme());

  if (debug && AlgorithmJacobianRec::aDebug && !fail)
    cout << " TB_EAB " << relg->TB[0] << " " << relg->TB[1] << " "
         << relg->TB[2] << " " << relg->eaB[0] << " " << relg->eaB[1] << " "
         << relg->eaB[2] << " ";

  if (!fail) {
    if (debug && AlgorithmJacobianRec::aDebug) cout << "vol pos" << endl;

    /**
    ActiveConstraintRegion* region = rnode -> getACR();
    CayleyParameterization* desc = new CayleyParameterization(cgK, false);
    ConstraintCheck* detector = new ConstraintCheck(cgK, this -> df);
    int noPoints = 0;

    bool noGoodOrientation = false;
    list<Orientation*> real;
    // real = findRealizations(cgK, des, rnode -> getFlipScheme());
    real.push_back(relg);
    CayleyPoint* cayleyPoint = new CayleyPoint(cayleyData);
    cayleyPoint -> setID(region -> getPointsCreated());
    region -> incrementPointsCreated();
    cayleyPoint -> setRealizable(!real.empty());
    cayleyPoint -> addOrientation(relg);

    if(!real.empty()) {
            for(list<Orientation*>::iterator ori_on_lattice = real.begin();
    ori_on_lattice != real.end()&& !Settings::AtlasBuilding::stop;
    ori_on_lattice++) {
                    // check for steric constraint violation and for any
    possible future contacts if(!detector -> stericsViolated(*ori_on_lattice)) {
                            /// check for bad angle
                            if((*ori_on_lattice) -> angleViolated())
                                    cayleyPoint -> incrementBadAngleN();
                            else
                                    noGoodOrientation = false;

                            bool boundary_ori_found_and_saved = false;
                            if(Settings::RootNodeCreation::createChildren) {
                                    findBoundary(ori_on_lattice, des, cgK,
    rnode, region, detector, false, noGoodOrientation, noPoints,
    boundary_ori_found_and_saved, cayleyPoint);
                            }

                            if(!boundary_ori_found_and_saved) {
                                    if(!(*ori_on_lattice) -> angleViolated()) {
                                            cayleyPoint ->
    addOrientation((*ori_on_lattice));
                                    }
                                    else {
                                            delete(*ori_on_lattice); // orange
                                            // p4d -> badAngle = true; //
    commented since incrementBadAngleN done above
                                    }
                            }
                            else {
                                    delete(*ori_on_lattice);
                            }
                    }
                    else {
                            if((*ori_on_lattice) -> angleViolated())
                                    cayleyPoint -> incrementBadAngleN();

                            delete(*ori_on_lattice); // collision
                            cayleyPoint -> incrementCollidN();
                    }
            }
    }
    */
    CayleyPoint* cayleyPoint = new CayleyPoint(cayleyData);
    cayleyPoint->setRealizable(true);
    cayleyPoint->addOrientation(relg);
    return cayleyPoint;
  } else {
    // hit = false;?
    if (debug && AlgorithmJacobianRec::aDebug) cout << "vol neg" << endl;
    /**
    ActiveConstraintRegion* region = rnode -> getACR();
    bool noGoodOrientation = true;
    CayleyParameterization* desc = new CayleyParameterization(cgK, false);
    ConstraintCheck* detector = new ConstraintCheck(cgK, this -> df);
    int noPoints = 0;
    vector<vector<int> > flipScheme = rnode -> getFlipScheme();

    Orientation* realization = CartesianRealizer::computeRealization(cgK, des,
    flip, fail, flipScheme); list<Orientation*> ori; ori.push_back(realization);
    list<Orientation*>::iterator ori_on_lattice = ori.begin();

    Orientation* orie_on_boundary = CartesianRealizer::computeRealization(cgK,
    des, flip, fail, flipScheme);

    bool binvalid = !detector -> stericsViolated(orie_on_boundary);
    list<pair<int, int> > contactList2 = detector -> getContacts(flip);
                    // binary search
                    // continue search till you find a valid configuration with
    new contacts
                    // find new contact for new small threshold t // FIXME: what
    does this line mean? int num_bin_iteration = 0; while ((contactList2.empty()
    || !binvalid) && num_bin_iteration < 30) { num_bin_iteration++; des ->
    stepGridBinary(binvalid); delete orie_on_boundary;
                            // num_samples++;
                            orie_on_boundary =
    CartesianRealizer::computeRealization(cgK, des, flip, fail, flipScheme);
                            binvalid = !detector ->
    stericsViolated(orie_on_boundary); contactList2 = detector ->
    getContacts(flip);
                    }

                    cout << "____________________AlgorithmJacobianRec
    walkByJacobian" << endl;

                    // cout << binvalid << endl;
                    // cout << contactList2.size() << endl;
                    // cout <<
    Settings::AtlasBuilding::ifBadAngleWitness_createChild << endl;
                    // cout << orie_on_boundary -> angleViolated() << endl;

                    // double check if the valid configuration with new contacts
    exists.
                    // it may not in case it exits the loop because of
    num_bin_iteration is big
                    // --------------------
                    // IF STATEMENT
                    // --------------------
                    cout << "binvalid is " << binvalid << " contactList2 size is
    " << contactList2.size() << endl; if (binvalid && contactList2.size() != 0
    && (Settings::AtlasBuilding::ifBadAngleWitness_createChild ||
    !orie_on_boundary -> angleViolated())) { cout << "binvalid &&
    contactList2.size() != 0 &&
    (Settings::AtlasBuilding::ifBadAngleWitness_createChild || !orie_on_boundary
    -> angleViolated())" << endl;
                            // Orientation* orie_on_boundary = sr ->
    getOrienation(); CayleyPoint* cayleyPoint = new CayleyPoint(cayleyData);
                            createChildContactGraps_fromTheBoundary(contactList2,
    ori_on_lattice, cgK, false, orie_on_boundary, rnode, noGoodOrientation,
    noPoints, cayleyPoint);

            if(Settings::Saving::saveBoundary) { // (since we do not save it to
    child) // Do not add witness point to parent node! cout <<
    "Settings::Saving::saveBoundary" << endl; vector<double> outt = des ->
    getPoint(); //to get exact boundary position instead of grid position
                                    CayleyPoint* p4 = new CayleyPoint(outt);
    //use p4 INSTEAD OF p4d  to save exact position of point after binary search
                p4 -> setID(region -> getPointsCreated());
                region -> incrementPointsCreated();

                                    ///
                                    bool cayleyPointEqual = true;
                                    bool orientationEqual = false;
                                    list<pair<CayleyPoint*, Orientation*> > temp
    = rnode -> getListOfVisitedPoints(); for(list<pair<CayleyPoint*,
    Orientation*> >::iterator candidate = temp.begin(); candidate != temp.end();
    candidate++) { cayleyPointEqual = true; orientationEqual = false; for(int i
    = 0; i < p4 -> getData().size(); i++) { if(p4 -> getData()[i] != candidate
    -> first -> getData()[i]) { cayleyPointEqual = false; break;
                                                    }
                                            }
                                            if(orie_on_boundary ->
    isEqual(candidate -> second)) { orientationEqual = true;
                                            }
                                            if(cayleyPointEqual &&
    orientationEqual) { break;
                                            }
                                    }
                                    ///
                                    if(!cayleyPointEqual || !orientationEqual) {
                                            rnode -> pushBackVisitedPoint(p4,
    orie_on_boundary); if(!orie_on_boundary -> angleViolated()) { p4 ->
    addOrientation(orie_on_boundary); cout << "walkByPreviousDirection has found
    a suitable orientation on boundary......................................."
    << endl;
                                                    //
    boundary_ori_found_and_saved = true; } else { delete orie_on_boundary; //
    orange
                                                    // p4 -> badAngle = true;
                                                    p4 -> incrementBadAngleN();
                                            }
                                            // pthread_mutex_lock(&space_sync);
                                            region -> AddSamplePoint(p4); // should
    it be insertwitness?!!! no, witness is the one saved at the child node. we
    are here saving the boundary config to current node.
                                            //
    pthread_mutex_unlock(&space_sync);

                                            noPoints++;
                                    }
                            } else
                                    delete orie_on_boundary;

                            // break; // if found a collision in one direction
    then stop } // if( valid(orie_on_boundary,detector) && contactList2.size()
    != 0
    */
    delete relg;
    return NULL;
  }
}

// CartesianRealizer* AlgorithmJacobianRec::walkByJacobian(ConvexChart* des,
// double paramconnections[][12], int g, int flip, MatrixXd JacInv, double &
// gridstep, bool & hit)
CayleyPoint* AlgorithmJacobianRec::walkByJacobian(
    ActiveConstraintGraph* cgK, ConvexChart* des, double out[6], int g,
    int flip, MatrixXd JacInv, double ratio, bool& hit, AtlasNode* rnode) {
  Settings* set = Settings::getInstance();
  bool debug = false;
  debug = true;

  if (debug && AlgorithmJacobianRec::aDebug)
    cout << "_____walkByJacobian on grid " << g << endl;
  double pconnections[12][12];

  for (int i = 0; i < 12; i++)
    for (int j = 0; j < 12; j++)
      pconnections[i][j] =
          des->param_length[i][j];  // des->paramconnections[i][j];

  //	double minGridFactor = 1;
  vector<double> cayleyData;
  hit = false;
  for (size_t p = 0; p < des->parameters.size(); p++) {
    int v1 = des->parameters[p].first;
    int v2 = des->parameters[p].second;

    if (debug && AlgorithmJacobianRec::aDebug)
      cout << v1 << "-" << v2 << " " << pconnections[v1][v2];
    // pconnections[v1][v2] += des->direction[g] * Settings::Sampling::stepSize
    // * (diff_TB_eaB_list[p/5][g] /total_TB_eaB[g]) * (gridSteps[g]
    // /total_TB_eaB[g]);  //(gridSteps[g] /total_TB_eaB[g]) is for having same
    // step size with grid steps

    //		pconnections[v1][v2] +=  JacInv(p/5,g) * gridstep ;

    int gI = des->independent_directions.at(g);
    pconnections[v1][v2] =
        out[p] +
        JacInv(p, g) * ratio;  //* gridstep / Settings::Sampling::gridSteps[gI];

    int isContact = des->find_pair_in_vector(des->contacts, des->parameters[p]);

    Atom* atomA = des->getAtom(v1);
    Atom* atomB = des->getAtom(v2);
    double bondingLowerBound = df->bondingLowerBound(atomA, atomB);
    double bondingUpperBound = df->bondingUpperBound(atomA, atomB);
    double collisionLowerBound = df->collisionLowerBound(atomA, atomB);

    // free cayley params
    if (isContact < 0 &&
        pconnections[v1][v2] <
            collisionLowerBound) {  // tetrahedral boundary is handled by volume
                                    // test, but collison should be handled
                                    // seperately
      pconnections[v1][v2] = collisionLowerBound;
      //			isBoundary[g] = true;
      if (debug && AlgorithmJacobianRec::aDebug) cout << "f";
      hit = true;

      //			double changeFactor =
      //abs(paramconnections[v1][v2] - pconnections[v1][v2]) / abs(JacInv(p/5,g)
      //* gridstep); 			if(changeFactor < minGridFactor) minGridFactor =
      //changeFactor;
    }

    if (isContact >= 0 && set->Sampling.short_range_sampling) {
      if (pconnections[v1][v2] <
          bondingLowerBound) {  // tetrahedral boundary is handled by volume
                                // test, but collison should be handled
                                // seperately
        pconnections[v1][v2] = bondingLowerBound;
        //			isBoundary[g] = true;
        if (debug && AlgorithmJacobianRec::aDebug) cout << "c";
        hit = true;
      } else if (pconnections[v1][v2] > bondingUpperBound) {
        pconnections[v1][v2] = bondingUpperBound;
        //			isBoundary[g] = true;
        if (debug && AlgorithmJacobianRec::aDebug) cout << "b";
        hit = true;
      }
    }

    //		double min = des->minconnections[v1][v2]; // =0; temporarily set it 0
    //to allow collision on free cayley params 		double max =
    //des->maxconnections[v1][v2]; 		if( max < pconnections[v1][v2]){
    //			pconnections[v1][v2]  = max ;
    //			if(debug && AlgorithmJacobianRec::aDebug)  cout << " x";
    //			hit = true;
    //
    //			double changeFactor = abs(paramconnections[v1][v2] -
    //pconnections[v1][v2]) / abs(JacInv(p/5,g) * gridstep); 			if(changeFactor <
    //minGridFactor) minGridFactor = changeFactor;
    //
    //		}
    //		else if( min > pconnections[v1][v2]){
    //			pconnections[v1][v2]  = min ;
    //			if(debug && AlgorithmJacobianRec::aDebug)  cout << " n"
    //; 			hit = true;
    //
    //			double changeFactor =  abs(paramconnections[v1][v2] -
    //pconnections[v1][v2]) / abs(JacInv(p/5,g) * gridstep); 			if(changeFactor <
    //minGridFactor) minGridFactor = changeFactor;
    //
    //		}

    pconnections[v2][v1] = pconnections[v1][v2];
    if (debug && AlgorithmJacobianRec::aDebug)
      cout << " " << pconnections[v1][v2] << ", ";

    //		if(p==15 && pconnections[v1][v2]<9.8)
    //			{
    //				cout << endl;
    //			}

    cayleyData.push_back(pconnections[v1][v2]);
  }

  //	if(hit)
  //		des->isBoundary[g] = true;

  //-----------allows you have point which is projection on cayley space when
  //you pass it. But is not whole grid step. so it cause dense sampling
  //	CartesianRealizer *relg;
  //	if(hit && minGridFactor !=0 && minGridFactor !=1 && false)  //SET IT
  //FALSE TO PREVENT INFINITE loop
  //	{
  //		if(debug && AlgorithmJacobianRec::aDebug)  cout << "with
  //minGridFactor " << minGridFactor << " "; 		gridstep = gridstep *
  //minGridFactor; 		relg = walkByJacobian(des,  paramconnections, g, flip,
  //JacInv, gridstep);
  //	}
  //	else
  //		relg = new CartesianRealizer(des, flip, pconnections);

  // CartesianRealizer *relg = new CartesianRealizer(des, flip, pconnections);
  bool fail;
  Orientation* relg = CartesianRealizer::computeRealization(
      cgK, des, flip, pconnections, fail, rnode->getFlipScheme());

  if (debug && AlgorithmJacobianRec::aDebug && !fail)
    cout << " TB_EAB " << relg->TB[0] << " " << relg->TB[1] << " "
         << relg->TB[2] << " " << relg->eaB[0] << " " << relg->eaB[1] << " "
         << relg->eaB[2] << " ";

  //	relg->axis = g;

  if (!fail) {
    if (debug && AlgorithmJacobianRec::aDebug) cout << "vol pos" << endl;

    /**
    ActiveConstraintRegion* region = rnode -> getACR();
    CayleyParameterization* desc = new CayleyParameterization(cgK, false);
    ConstraintCheck* detector = new ConstraintCheck(cgK, this -> df);
    int noPoints = 0;

    bool noGoodOrientation = false;
    list<Orientation*> real;
    // real = findRealizations(cgK, des, rnode -> getFlipScheme());
    real.push_back(relg);
    CayleyPoint* cayleyPoint = new CayleyPoint(cayleyData);
    cayleyPoint -> setID(region -> getPointsCreated());
    region -> incrementPointsCreated();
    cayleyPoint -> setRealizable(!real.empty());
    cayleyPoint -> addOrientation(relg);

    if(!real.empty()) {
            for(list<Orientation*>::iterator ori_on_lattice = real.begin();
    ori_on_lattice != real.end()&& !Settings::AtlasBuilding::stop;
    ori_on_lattice++) {
                    // check for steric constraint violation and for any
    possible future contacts if(!detector -> stericsViolated(*ori_on_lattice)) {
                            /// check for bad angle
                            if((*ori_on_lattice) -> angleViolated())
                                    cayleyPoint -> incrementBadAngleN();
                            else
                                    noGoodOrientation = false;

                            bool boundary_ori_found_and_saved = false;
                            if(Settings::RootNodeCreation::createChildren) {
                                    findBoundary(ori_on_lattice, des, cgK,
    rnode, region, detector, false, noGoodOrientation, noPoints,
    boundary_ori_found_and_saved, cayleyPoint);
                            }

                            if(!boundary_ori_found_and_saved) {
                                    if(!(*ori_on_lattice) -> angleViolated()) {
                                            cayleyPoint ->
    addOrientation((*ori_on_lattice));
                                    }
                                    else {
                                            delete(*ori_on_lattice); // orange
                                            // p4d -> badAngle = true; //
    commented since incrementBadAngleN done above
                                    }
                            }
                            else {
                                    delete(*ori_on_lattice);
                            }
                    }
                    else {
                            if((*ori_on_lattice) -> angleViolated())
                                    cayleyPoint -> incrementBadAngleN();

                            delete(*ori_on_lattice); // collision
                            cayleyPoint -> incrementCollidN();
                    }
            }
    }
    */
    CayleyPoint* cayleyPoint = new CayleyPoint(cayleyData);
    cayleyPoint->setRealizable(true);
    cayleyPoint->addOrientation(relg);
    return cayleyPoint;
  } else {
    // hit = false;?
    /// TODO add code to recursively sample. ignore the repeating cayley points
    /// first.
    if (debug && AlgorithmJacobianRec::aDebug) cout << "vol neg" << endl;
    /**
    ActiveConstraintRegion* region = rnode -> getACR();
    bool noGoodOrientation = true;
    CayleyParameterization* desc = new CayleyParameterization(cgK, false);
    ConstraintCheck* detector = new ConstraintCheck(cgK, this -> df);
    int noPoints = 0;
    vector<vector<int> > flipScheme = rnode -> getFlipScheme();

    Orientation* realization = CartesianRealizer::computeRealization(cgK, des,
    flip, fail, flipScheme); list<Orientation*> ori; ori.push_back(realization);
    list<Orientation*>::iterator ori_on_lattice = ori.begin();

    Orientation* orie_on_boundary = CartesianRealizer::computeRealization(cgK,
    des, flip, fail, flipScheme);

    bool binvalid = !detector -> stericsViolated(orie_on_boundary);
    list<pair<int, int> > contactList2 = detector -> getContacts(flip);
                    // binary search
                    // continue search till you find a valid configuration with
    new contacts
                    // find new contact for new small threshold t // FIXME: what
    does this line mean? int num_bin_iteration = 0; while ((contactList2.empty()
    || !binvalid) && num_bin_iteration < 30) { num_bin_iteration++; des ->
    stepGridBinary(binvalid); delete orie_on_boundary;
                            // num_samples++;
                            orie_on_boundary =
    CartesianRealizer::computeRealization(cgK, des, flip, fail, flipScheme);
                            binvalid = !detector ->
    stericsViolated(orie_on_boundary); contactList2 = detector ->
    getContacts(flip);
                    }

                    cout << "____________________AlgorithmJacobianRec
    walkByJacobian" << endl;

                    // cout << binvalid << endl;
                    // cout << contactList2.size() << endl;
                    // cout <<
    Settings::AtlasBuilding::ifBadAngleWitness_createChild << endl;
                    // cout << orie_on_boundary -> angleViolated() << endl;

                    // double check if the valid configuration with new contacts
    exists.
                    // it may not in case it exits the loop because of
    num_bin_iteration is big
                    // --------------------
                    // IF STATEMENT
                    // --------------------
                    cout << "binvalid is " << binvalid << " contactList2 size is
    " << contactList2.size() << endl; if (binvalid && contactList2.size() != 0
    && (Settings::AtlasBuilding::ifBadAngleWitness_createChild ||
    !orie_on_boundary -> angleViolated())) { cout << "binvalid &&
    contactList2.size() != 0 &&
    (Settings::AtlasBuilding::ifBadAngleWitness_createChild || !orie_on_boundary
    -> angleViolated())" << endl;
                            // Orientation* orie_on_boundary = sr ->
    getOrienation(); CayleyPoint* cayleyPoint = new CayleyPoint(cayleyData);
                            if(rnode -> getDim() == 1) {
                                    cout << "yes" << endl;
                            }
                            createChildContactGraps_fromTheBoundary(contactList2,
    ori_on_lattice, cgK, false, orie_on_boundary, rnode, noGoodOrientation,
    noPoints, cayleyPoint);

            if(Settings::Saving::saveBoundary) { // (since we do not save it to
    child) // Do not add witness point to parent node! cout <<
    "Settings::Saving::saveBoundary" << endl; vector<double> outt = des ->
    getPoint(); //to get exact boundary position instead of grid position
                                    CayleyPoint* p4 = new CayleyPoint(outt);
    //use p4 INSTEAD OF p4d  to save exact position of point after binary search
                p4 -> setID(region -> getPointsCreated());
                region -> incrementPointsCreated();

                                    ///
                                    bool cayleyPointEqual = true;
                                    bool orientationEqual = false;
                                    list<pair<CayleyPoint*, Orientation*> > temp
    = rnode -> getListOfVisitedPoints(); for(list<pair<CayleyPoint*,
    Orientation*> >::iterator candidate = temp.begin(); candidate != temp.end();
    candidate++) { cayleyPointEqual = true; orientationEqual = false; for(int i
    = 0; i < p4 -> getData().size(); i++) { if(p4 -> getData()[i] != candidate
    -> first -> getData()[i]) { cayleyPointEqual = false; break;
                                                    }
                                            }
                                            if(orie_on_boundary ->
    isEqual(candidate -> second)) { orientationEqual = true;
                                            }
                                            if(cayleyPointEqual &&
    orientationEqual) { break;
                                            }
                                    }
                                    ///
                                    if(!cayleyPointEqual || !orientationEqual) {
                                            rnode -> pushBackVisitedPoint(p4,
    orie_on_boundary); if(!orie_on_boundary -> angleViolated()) { p4 ->
    addOrientation(orie_on_boundary); cout << "walkByJacobian has found a
    suitable orientation on boundary......................................." <<
    endl;
                                                    //
    boundary_ori_found_and_saved = true; } else { delete orie_on_boundary; //
    orange
                                                    // p4 -> badAngle = true;
                                                    p4 -> incrementBadAngleN();
                                            }
                                            // pthread_mutex_lock(&space_sync);
                                            region -> AddSamplePoint(p4); // should
    it be insertwitness?!!! no, witness is the one saved at the child node. we
    are here saving the boundary config to current node.
                                            //
    pthread_mutex_unlock(&space_sync);

                                            noPoints++;
                                    }
                            } else
                                    delete orie_on_boundary;

                            // break; // if found a collision in one direction
    then stop } // if( valid(orie_on_boundary,detector) && contactList2.size()
    != 0
    */

    delete relg;
    return NULL;
  }
}

/*
//CartesianRealizer* AlgorithmJacobianRec::walkByJacobian(ConvexChart* des,
double paramconnections[][12], int g, int flip, MatrixXd JacInv, double &
gridstep, bool & hit) CartesianRealizer*
AlgorithmJacobianRec::walkByJacobian(ConvexChart* des,  double
paramconnections[][12], int g, int flip, MatrixXd JacInv, double ratio, bool &
hit)
{
        bool debug = false;
        debug = true;

        if(debug && AlgorithmJacobianRec::aDebug)  cout << "_____walkByJacobian
on grid " << g << endl; double pconnections[12][12];



        for(int i = 0; i<12; i++)	for(int j = 0; j<12; j++)
pconnections[i][j] = paramconnections[i][j]; //des->paramconnections[i][j];

//	double minGridFactor = 1;
        hit = false;
        for(int p = des->paramsOrdered.size()-5; p>=0; p=p-5) //start change
from reverse direction, because of the order of dependence
        {
                //des->updateBoundaries();
                des->updateBoundariesSpecific(p/5);

                int v1 = atoi( des->paramsOrdered.substr(p,2).c_str() ) ;
                int v2 = atoi( des->paramsOrdered.substr(p+2,2).c_str() ) ;

                double min = des->minconnections[v1][v2];
                double max = des->maxconnections[v1][v2];


                double prev = pconnections[v1][v2];

                if(debug && AlgorithmJacobianRec::aDebug)  cout << v1 << "-" <<
v2 << " " << pconnections[v1][v2] ;
                //pconnections[v1][v2] += des->direction[g] *
Settings::Sampling::stepSize * (diff_TB_eaB_list[p/5][g] /total_TB_eaB[g]) *
(gridSteps[g] /total_TB_eaB[g]);  //(gridSteps[g] /total_TB_eaB[g]) is for
having same step size with grid steps


//		pconnections[v1][v2] +=  JacInv(p/5,g) * gridstep ;

                int gI = des->independent_directions.at(g);
                pconnections[v1][v2] +=  JacInv(p/5,g) * ratio ; gridstep /
Settings::Sampling::gridSteps[gI];

                double curr = pconnections[v1][v2];


                int isContact = des->getCurContacts().find(
des->paramsOrdered.substr(p,4) , 0);

                //free cayley params
                if( isContact < 0 && pconnections[v1][v2] <
des->mygetRadiusSum(v1, v2) * Settings::Collision::contactLowerThreshold ){
//tetrahedral boundary is handled by volume test, but collison should be handled
seperately pconnections[v1][v2] = des->mygetRadiusSum(v1, v2) *
Settings::Collision::contactLowerThreshold;
//			isBoundary[g] = true;
                        if(debug && AlgorithmJacobianRec::aDebug)  cout << "f";
                        hit = true;

//			double changeFactor = abs(paramconnections[v1][v2] -
pconnections[v1][v2]) / abs(JacInv(p/5,g) * gridstep);
//			if(changeFactor < minGridFactor) minGridFactor =
changeFactor;
                }




                if( curr > max ){
                        hit = true;
                        pconnections[v1][v2]  = max ;
                        pconnections[v2][v1]  = max ;
                        if(debug && AlgorithmJacobianRec::aDebug)  cout <<
"max";
//			if(max == prev) //make it different than prev point
//			{
//				pconnections[v1][v2]  = max -
Settings::Sampling::stepSize;
//				pconnections[v2][v1]  = max -
Settings::Sampling::stepSize;
//			}
                }
                else if( curr < min ){
                        hit = true;
                        pconnections[v1][v2]  = min ;
                        pconnections[v2][v1]  = min ;
                        if(debug && AlgorithmJacobianRec::aDebug)  cout <<
"min";
//			if(min == prev) //make it different than prev point
//			{
//				pconnections[v1][v2]  = min +
Settings::Sampling::stepSize;
//				pconnections[v2][v1]  = min +
Settings::Sampling::stepSize;
//			}
                }











//		//ACTUALLT THIS ARE NOT BOU
//		//contact for 6 d
//		if( isContact>=0 && pconnections[v1][v2] <
des->mygetRadiusSum(v1, v2) * Settings::Collision::contactLowerThreshold ){
//tetrahedral boundary is handled by volume test, but collison should be handled
seperately
//			pconnections[v1][v2] = des->mygetRadiusSum(v1, v2) *
Settings::Collision::contactLowerThreshold;
////			isBoundary[g] = true;
//			if(debug && AlgorithmJacobianRec::aDebug)  cout << "c";
//			hit = true;
//		}
//		else if( isContact>=0 && pconnections[v1][v2] >
des->mygetRadiusSum(v1, v2) + Settings::Collision::contactUpperThreshold )
//		{
//			pconnections[v1][v2] =  des->mygetRadiusSum(v1, v2) +
Settings::Collision::contactUpperThreshold;
////			isBoundary[g] = true;
//			if(debug && AlgorithmJacobianRec::aDebug)  cout << "b";
//			hit = true;
//		}



//		double min = des->minconnections[v1][v2]; // =0; temporarily set
it 0 to allow collision on free cayley params
//		double max = des->maxconnections[v1][v2];
//		if( max < pconnections[v1][v2]){
//			pconnections[v1][v2]  = max ;
//			if(debug && AlgorithmJacobianRec::aDebug)  cout << " x";
//			hit = true;
//
//			double changeFactor = abs(paramconnections[v1][v2] -
pconnections[v1][v2]) / abs(JacInv(p/5,g) * gridstep);
//			if(changeFactor < minGridFactor) minGridFactor =
changeFactor;
//
//		}
//		else if( min > pconnections[v1][v2]){
//			pconnections[v1][v2]  = min ;
//			if(debug && AlgorithmJacobianRec::aDebug)  cout << " n"
;
//			hit = true;
//
//			double changeFactor =  abs(paramconnections[v1][v2] -
pconnections[v1][v2]) / abs(JacInv(p/5,g) * gridstep);
//			if(changeFactor < minGridFactor) minGridFactor =
changeFactor;
//
//		}


                pconnections[v2][v1] = pconnections[v1][v2];
                if(debug && AlgorithmJacobianRec::aDebug)    cout << " " <<
pconnections[v1][v2] << ", ";

//		if(p==15 && pconnections[v1][v2]<9.8)
//			{
//				cout << endl;
//			}
        }



//	if(hit)
//		des->isBoundary[g] = true;


//-----------allows you have point which is projection on cayley space when you
pass it. But is not whole grid step. so it cause dense sampling
//	CartesianRealizer *relg;
//	if(hit && minGridFactor !=0 && minGridFactor !=1 && false)  //SET IT
FALSE TO PREVENT INFINITE loop
//	{
//		if(debug && AlgorithmJacobianRec::aDebug)  cout << "with
minGridFactor " << minGridFactor << " ";
//		gridstep = gridstep * minGridFactor;
//		relg = walkByJacobian(des,  paramconnections, g, flip, JacInv,
gridstep);
//	}
//	else
//		relg = new CartesianRealizer(des, flip, pconnections);


        CartesianRealizer *relg = new CartesianRealizer(des, flip,
pconnections);

        if(debug && AlgorithmJacobianRec::aDebug)    cout << " TB_EAB " <<
relg->TB[0] << " " << relg->TB[1] << " " << relg->TB[2] << " " << relg->eaB[0]
<< " " << relg->eaB[1] << " " << relg->eaB[2] << " "  ;

//	relg->axis = g;

        if( relg->getVolumeTest() ){
                if(debug && AlgorithmJacobianRec::aDebug)  cout << "vol pos" <<
endl; return relg;
        }
        else{
                //hit = false;?
                if(debug && AlgorithmJacobianRec::aDebug)  cout << "vol neg" <<
endl; delete relg; return NULL;
        }
}
*/

/*
//CartesianRealizer* AlgorithmJacobianRec::walkByJacobian2(ConvexChart* des,
double paramconnections[][12], int g1, int g2, int flip, MatrixXd JacInv, double
gridstep1, double gridstep2, bool & hit) CayleyPoint*
AlgorithmJacobianRec::walkByJacobian2(ActiveConstraintGraph* cgK, ConvexChart*
des,  double out[6], int g1, int g2, int flip, MatrixXd JacInv, double ratio1,
double ratio2, bool & hit, AtlasNode* rnode)
{
        bool debug = false;
        if(debug && AlgorithmJacobianRec::aDebug)  cout << "_____walkByJacobian
on grid " << g1 << " " << g2 << endl; double pconnections[12][12];

        if(debug && AlgorithmJacobianRec::aDebug)  cout <<  "JacInv \n" <<
JacInv << endl;

        for(int i = 0; i<12; i++)	for(int j = 0; j<12; j++)
pconnections[i][j] = des->param_length[i][j]; //des->paramconnections[i][j];

        vector<double> cayleyData;
        hit = false;
        for(size_t p = 0; p<des->parameters.size(); p++)
        {
                int v1 = des->parameters[p].first;
                int v2 = des->parameters[p].second;

                if(debug && AlgorithmJacobianRec::aDebug)  cout << v1 << "-" <<
v2 << " " << pconnections[v1][v2] ;

//		pconnections[v1][v2] +=  JacInv(p/5,g1) * gridstep1 +
JacInv(p/5,g2) * gridstep2;

                int gI1 = des->independent_directions.at(g1);
                int gI2 = des->independent_directions.at(g2);
//		pconnections[v1][v2] +=  JacInv(p/5,g1) * gridstep1 /
Settings::Sampling::gridSteps[gI1] + JacInv(p/5,g2) * gridstep2 /
Settings::Sampling::gridSteps[gI2]; pconnections[v1][v2] =out[p] +  JacInv(p,g1)
* ratio1 + JacInv(p,g2) * ratio2;

                int isContact = des->find_pair_in_vector( des->contacts,
des->parameters[p]);

                Atom * atomA = des->getAtom(v1);
                Atom * atomB = des->getAtom(v2);
                double bondingLowerBound = df->bondingLowerBound(atomA, atomB);
                double bondingUpperBound = df->bondingUpperBound(atomA, atomB);
                double collisionLowerBound = df->collisionLowerBound(atomA,
atomB);

                //free cayley params
                if( isContact<0 && pconnections[v1][v2] <  collisionLowerBound
){ //tetrahedral boundary is handled by volume test, but collison should be
handled seperately pconnections[v1][v2] = collisionLowerBound; if(debug &&
AlgorithmJacobianRec::aDebug)  cout << "f"; hit = true;
                }

                if(isContact>=0 && Settings::Sampling::short_range_sampling)
                {
                        if(  pconnections[v1][v2] <  bondingLowerBound ){
//tetrahedral boundary is handled by volume test, but collison should be handled
seperately pconnections[v1][v2] = bondingLowerBound;
        //			isBoundary[g] = true;
                                if(debug && AlgorithmJacobianRec::aDebug)  cout
<< "c"; hit = true;
                        }
                        else if(  pconnections[v1][v2] >  bondingUpperBound )
                        {
                                pconnections[v1][v2] =  bondingUpperBound;
        //			isBoundary[g] = true;
                                if(debug && AlgorithmJacobianRec::aDebug)  cout
<< "b"; hit = true;
                        }
                }

                pconnections[v2][v1] = pconnections[v1][v2];
                if(debug && AlgorithmJacobianRec::aDebug)  cout << " " <<
pconnections[v1][v2] << ", ";

                cayleyData.push_back( pconnections[v1][v2] );
        }

//	CartesianRealizer *relg = new CartesianRealizer(des, flip,
pconnections); bool fail; Orientation * relg =
CartesianRealizer::computeRealization(cgK, des, flip, pconnections, fail);


        if(debug && AlgorithmJacobianRec::aDebug && !fail)  cout << " TB_EAB "
<< relg->TB[0] << " " << relg->TB[1] << " " << relg->TB[2] << " " <<
relg->eaB[0] << " " << relg->eaB[1] << " " << relg->eaB[2] << " " << endl;

        if(  !fail )
        {
                if(debug && AlgorithmJacobianRec::aDebug)  cout << "vol pos" <<
endl;

                CayleyPoint* cayleyPoint = new CayleyPoint(cayleyData);
                cayleyPoint->setRealizable(true);
                cayleyPoint->addOrientation(relg);

                return cayleyPoint;
        }
        else{
                delete relg;
                return NULL;
        }
}

*/

vector<int> AlgorithmJacobianRec::computeIndependentColumns(
    MatrixXd A, int dim) {  // A is 6 by dim
  cout << "computeIndependentColumns" << endl;
  vector<int> indep;

  cout << "A \n" << A << endl;

  double maxDet = -1000000;
  MatrixXd B(dim, dim);

  for (int i = 0; i < 6; i++) {
    for (int n = 0; n < dim; n++) B(0, n) = A(i, n);
    for (int j = i + 1; j < 6; j++) {
      if (dim > 1)
        for (int n = 0; n < dim; n++) B(1, n) = A(j, n);
      for (int k = j + 1; k < 6; k++) {
        if (dim > 2)
          for (int n = 0; n < dim; n++) B(2, n) = A(k, n);
        for (int l = k + 1; l < 6; l++) {
          if (dim > 3)
            for (int n = 0; n < dim; n++) B(3, n) = A(l, n);
          for (int m = l + 1; m < 6; m++) {
            if (dim > 4)
              for (int n = 0; n < dim; n++) B(4, n) = A(m, n);

            double determ = B.determinant();
            if (maxDet < determ) {
              maxDet = determ;
              indep.clear();
              indep.push_back(i);
              indep.push_back(j);
              indep.push_back(k);
              indep.push_back(l);
              indep.push_back(m);
              for (int n = 0; n < 5 - dim; n++)
                indep.pop_back();  // delete last 5-dim indices, since indep is
                                   // size dim
            }
            if (dim <= 4) break;
          }
          if (dim <= 3) break;
        }
        if (dim <= 2) break;
      }
      if (dim <= 1) break;
    }
  }

  for (int n = 0; n < indep.size(); n++) cout << indep[n] << " " << endl;

  return indep;

  //	MatrixXd U = A.svd().matrixU();
  //	MatrixXd V = A.svd().matrixV();
  //	VectorXd svl =  A.svd().singularValues();
  //	cout <<  "U \n" << U << endl;
  //	cout <<  "V \n" << V << endl;
  //	cout <<  "svl \n" << svl << endl;

  //	MatrixXd R = A.qr().matrixR();
  //	MatrixXd Q = A.qr().matrixQ();
  //	cout <<  "Q \n" << Q << endl;
  //	cout <<  "R \n" << R << endl;
  //	int size = R.size();
  //	int i = 1;
  //	for(int j =0; j<size; j++)
  //	{
  //		if( R(i,j) != 0  )
  //		{
  //			des->independent_directions.push_back(j);
  //			i++;
  //			cout << j << " ";
  //		}
  //		if(i>size)
  //			break;
  //	}
  //	cout << endl;
}

/*
MatrixXd AlgorithmJacobianRec::computeJacobianfromCartesian(ConvexChart* des,
CartesianRealizer * real, int flip)
{
        bool debug = false;
        debug = true;

        if(debug && AlgorithmJacobianRec::aDebug)  cout <<
"________computeJacobianfromCartesian " << endl; int dim =
des->paramsOrdered.size()/5; MatrixXd JacInv(dim,dim); for(int i=0; i < dim;
i++) for(int j=0; j < dim; j++) JacInv(i,j) = 0;


if(des->independent_directions.size() > 0 ){
        vector<Atom*> helA = a->getAtoms(), helB = b->getAtoms();
        vector<Atom*> outputB;
        Atom *current;
        vector<pair<int,int> > PL = des->currentGraph->getParamLines();
//ConvexChart::mysetParamLine()

        for(int i=0; i<dim; i++)
        {
                int d = des->independent_directions[i];

                Vector3d eaB = real->eaB;
                Vector3d TB = real->TB;

                if(d<3)
                        TB[d] = TB[d] + Settings::Sampling::gridSteps[d];
                else
                        eaB[d-3] = eaB[d-3] + Settings::Sampling::gridSteps[d];

                eaB = eaB * PI / 180;
                Matrix3d RB;
                RB = AngleAxisd(eaB[0], Vector3d::UnitZ())* AngleAxisd(eaB[1],
Vector3d::UnitX())* AngleAxisd(eaB[2], Vector3d::UnitZ());


                Vector3d meanB(0,0,0);
                for(size_t iter = 0;   iter < helB.size();  iter++)
                {
                        current = new Atom(helB[iter]);
                        double *l = current->getLocation();
                        Vector3d p(l[0], l[1], l[2]);
                        Vector3d r = RB * p;
                        current->setLocation(r(0), r(1), r(2));
                        outputB.push_back(current);
                        meanB += r;
                }
                meanB = meanB /  helB.size();
                for(size_t iter = 0;   iter < outputB.size();  iter++)
                {
                        double *l = outputB[iter]->getLocation();
                        Vector3d r(l[0], l[1], l[2]);
                        Vector3d h = r + TB - meanB;
                        outputB[iter]->setLocation(h(0), h(1), h(2));
                //	cout << h(0) << " " <<  h(1) << " " << h(2) << endl;
                }




//	   for(size_t i = 0;i < this->paramlines.size();i++){ //Paramlines needs
to be defined for this. (in description class)
//		}

                for(size_t p = 0; p<des->paramsOrdered.size(); p=p+5)
                {
                        int v1 = atoi( des->paramsOrdered.substr(p,2).c_str() )
; int v2 = atoi( des->paramsOrdered.substr(p+2,2).c_str() ) ;

                        Atom* atomA = helA[PL[p/5].first];
                        Atom* atomB = outputB[PL[p/5].second];
                        double distance =
Utils::dist(atomA->getLocation(),atomB->getLocation()); cout << v1<<"-" << v2 <<
" PL:" << PL[p/5].first << "-" << PL[p/5].second << ", dist:" << distance  << ",
real_dist:" << real->pconnections[v1][v2] << endl;

                        JacInv(p/5,i) = (distance - real->pconnections[v1][v2]
)/ Settings::Sampling::gridSteps[d];
//			cout << (distance - real->pconnections[v1][v2] )/
Settings::Sampling::gridSteps[d] << " " ;

                        cout << "5-5: " << real->pconnections[0][6] << " " <<
Utils::dist(helA[5]->getLocation(),outputB[5]->getLocation()) << endl;
                }
//		cout << endl;



                   vector<Atom*>::iterator itb;
                   for(itb = outputB.begin();itb != outputB.end(); itb++){
                                   delete (*itb);
                   }
                   outputB.clear();
        }
}

if(debug && AlgorithmJacobianRec::aDebug)  cout <<  "JacInv car \n" << JacInv <<
endl; return JacInv;

}
*/

MatrixXd AlgorithmJacobianRec::jacobian(
    ActiveConstraintGraph* cgK, ConvexChart* des, CayleyPoint* cayleyPoint,
    int flip, int specificGridDiretion, bool gridStepsChange[6],
    MatrixXd cayleyDirections, bool& enteredNegVal, AtlasNode* rnode) {
  Settings* set = Settings::getInstance();
  bool debug = false;
  debug = true;

  for (int y = 0; y < 6; y++) gridStepsChange[y] = false;

  Orientation* real = cayleyPoint->getOrientations().front();
  int dim = des->parameters.size();

  // cayleyDirections.conservativeResize(dim,dim); //=
  // cayleyDirections.topLeftCorner<dim, dim>() ; MatrixXd temp(dim,dim);
  // for(int i=0; i<dim; i++) for(int j=0; j<dim; j++) temp(i,j) =
  // cayleyDirections(i,j); cayleyDirections = temp;

  Vector3d old_eaB = real->eaB;
  Vector3d old_TB = real->TB;

  vector<vector<double> > diff_TB_eaB_list;

  double out[6];
  cayleyPoint->getPoint(out);

  double pconnections[12][12];
  for (size_t i = 0; i < des->parameters.size(); i++) {
    int v1 = des->parameters[i].first;
    int v2 = des->parameters[i].second;
    // copy only parameters from real, not the rest of the connections because
    // they should be updated and changed later according to updated params.
    // if you copy the rest, they wont be updated in realization class since
    // there already exist a value for them
    des->param_length[v1][v2] = out[i];  // real->pconnections[v1][v2];
    des->param_length[v2][v1] = out[i];  // real->pconnections[v1][v2];
  }

  for (size_t i = 0; i < des->parameters.size(); i++) {
    int v1 = des->parameters[i].first;
    int v2 = des->parameters[i].second;
    if (debug && AlgorithmJacobianRec::aDebug)
      cout << i << ": " << v1 << "-" << v2 << " " << des->param_length[v1][v2]
           << ", ";
  }

  if (debug && AlgorithmJacobianRec::aDebug)
    cout << "TB_EAB " << old_TB[0] << " " << old_TB[1] << " " << old_TB[2]
         << " " << old_eaB[0] << " " << old_eaB[1] << " " << old_eaB[2] << " "
         << endl;

  int negVolDirection = -1;
  for (size_t direction = 0; direction < dim; direction++) {
    for (int k = 0; k < 12; k++) {
      for (int j = 0; j < 12; j++) {
        pconnections[k][j] = des->param_length[k][j];
      }
    }

    for (size_t i = 0; i < des->parameters.size(); i++) {
      int v1 = des->parameters[i].first;
      int v2 = des->parameters[i].second;
      pconnections[v1][v2] += cayleyDirections(
          i, direction);  // cayleySteps[i/5]; // Settings::Sampling::stepSize;
      pconnections[v2][v1] = pconnections[v1][v2];
    }

    // CartesianRealizer *relc = new CartesianRealizer(des, flip, pconnections);
    bool fail;
    Orientation* relc = CartesianRealizer::computeRealization(
        cgK, des, flip, pconnections, fail, rnode->getFlipScheme());
    if (fail) {  // volume neg
      // negVolDirection = direction;
      // enteredNegVal = true;
      // delete relc; relc = NULL;
      // break;

      if (specificGridDiretion != direction && specificGridDiretion != -1) {
        for (size_t i = 0; i < des->parameters.size(); i++) {
          int v1 = des->parameters[i].first;
          int v2 = des->parameters[i].second;
          pconnections[v1][v2] -=
              2 * cayleyDirections(i, direction);  // try reverse direction
          pconnections[v2][v1] = pconnections[v1][v2];
        }

        delete relc;
        relc = NULL;
        // relc = new CartesianRealizer(des, flip, pconnections);
        relc = CartesianRealizer::computeRealization(
            cgK, des, flip, pconnections, fail, rnode->getFlipScheme());
        if (!fail) {
          for (size_t i = 0; i < des->parameters.size(); i++)
            cayleyDirections(i, direction) =
                -cayleyDirections(i, direction);  // try reverse direction

          gridStepsChange[direction] = true;
        }
      }

      if (fail) {  // volume neg again
        negVolDirection = direction;
        enteredNegVal = true;
        delete relc;
        relc = NULL;
        break;

        if (debug && AlgorithmJacobianRec::aDebug)
          cout << "J volume negative when increased by stepsize" << endl;

        // for(int k = 0; k < 12; k++) for(int j = 0; j < 12; j++)
        // pconnections[k][j] = des -> paramconnections[k][j]; // des ->
        // paramconnections[i][j];

        // des -> updateBoundaries();
        // for(size_t i = 0; i < des -> paramsOrdered.size(); i = i + 5)
        for (int i = des->parameters.size() - 1; i >= 0;
             i--) {  // start change from reverse direction, because of the
                     // order of dependence
          // des -> updateBoundaries();
          des->updateBoundariesSpecific(i);
          int v1 = des->parameters[i].first;
          int v2 = des->parameters[i].second;

          double min = des->param_lengthLower[v1][v2];
          double max = des->param_lengthUpper[v1][v2];

          double prev = des->param_length[v1][v2];
          double curr = pconnections[v1][v2];
          if (debug && AlgorithmJacobianRec::aDebug)
            cout << "J " << i << ": " << v1 << "-" << v2 << " max:" << max
                 << " min:" << min << " prev:" << prev << " curr:" << curr;

          if (curr > max) {
            pconnections[v1][v2] = max;
            pconnections[v2][v1] = max;
            if (max == prev) {  // make it different than prev point
              pconnections[v1][v2] = max - set->Sampling.stepSize;
              pconnections[v2][v1] = max - set->Sampling.stepSize;
            }
          } else if (curr < min) {
            pconnections[v1][v2] = min;
            pconnections[v2][v1] = min;
            if (min == prev) {  // make it different than prev point
              pconnections[v1][v2] = min + set->Sampling.stepSize;
              pconnections[v2][v1] = min + set->Sampling.stepSize;
            }
          }
          if (debug && AlgorithmJacobianRec::aDebug)
            cout << " now: " << pconnections[v1][v2] << endl;

          des->param_length[v1][v2] = pconnections[v1][v2];
          des->param_length[v2][v1] = pconnections[v1][v2];
        }

        // take it back the changes you have just done above.
        for (size_t i = 0; i < des->parameters.size(); i++) {
          int v1 = des->parameters[i].first;
          int v2 = des->parameters[i].second;
          // copy only parameters from real, not the rest of the connections
          // because they should be updated and changed later according to
          // updated params. if you copy the rest, they wont be updated in
          // realization class since there already exist a value for them
          des->param_length[v1][v2] = out[i];  // real -> pconnections[v1][v2];
          des->param_length[v2][v1] = out[i];  // real -> pconnections[v1][v2];
        }

        bool fail;
        delete relc;
        relc = NULL;
        // relc = new CartesianRealizer(des, flip, pconnections);
        relc = CartesianRealizer::computeRealization(
            cgK, des, flip, pconnections, fail, rnode->getFlipScheme());

        if (!fail) {
          for (size_t i = 0; i < des->parameters.size(); i++) {
            int v1 = des->parameters[i].first;
            int v2 = des->parameters[i].second;
            cayleyDirections(i, direction) =
                pconnections[v1][v2] - des->param_length[v1][v2];
          }

          if (debug && AlgorithmJacobianRec::aDebug)
            cout << "cayleyDirections after neg vol fix: \n"
                 << cayleyDirections << endl;

        } else {
          negVolDirection = direction;
          enteredNegVal = true;
          delete relc;
          relc = NULL;
          break;
        }
      }
    }

    if (debug && AlgorithmJacobianRec::aDebug)
      cout << "TB_EAB after walking on cayley combined direction" << direction
           << " is: " << relc->TB[0] << " " << relc->TB[1] << " " << relc->TB[2]
           << " " << relc->eaB[0] << " " << relc->eaB[1] << " " << relc->eaB[2]
           << " " << endl;

    vector<double> diff_TB_eaB(6);  /// 6 zero-initialized elements
    if (!fail) {  // do not check angle now, maybe when you change the direction
                  // of sampling using jacobian, the new realization will have
                  // good angle
      Vector3d eaB = relc->eaB;
      Vector3d TB = relc->TB;

      Vector3d diff_TB = TB - old_TB;
      Vector3d diff_eaB = eaB - old_eaB;
      // if(debug && AlgorithmJacobianRec::aDebug) cout << "J " << eaB[0] << " "
      // << old_eaB[0] << " " << diff_eaB[0] << " " ;
      for (int h = 0; h < 3; h++) {
        if (abs(diff_eaB[h]) > 180) {
          diff_eaB[h] =
              -1 * Utils::sign(diff_eaB[h]) *
              (360 - abs(diff_eaB[h]));  // to keep it between (-180, 180)
        }
      }
      // for(int h = 0; h < 3; h++) diff_eaB[h] = ((int) diff_eaB[h] + 180) %
      // 360 - 180; // to keep it between (-180, 180) if(debug &&
      // AlgorithmJacobianRec::aDebug) cout << diff_eaB[0] << endl;

      for (int h = 0; h < 3; h++) diff_TB_eaB[h] = diff_TB[h];
      for (int h = 0; h < 3; h++) diff_TB_eaB[h + 3] = diff_eaB[h];

      // double min = diff_TB_eaB[3];
      // if(abs(diff_TB_eaB[3] + 180) < abs(min)) min = diff_TB_eaB[3] + 180;
      // if(abs(diff_TB_eaB[3] - 180) < abs(min)) min = diff_TB_eaB[3] - 180;
      // diff_TB_eaB[3] = min;

      // min = diff_eaB[2];
      // if(diff_eaB[2] + 180 < min) min = diff_eaB[2] + 180;
      // if(diff_eaB[2] - 180 < min) min = diff_eaB[2] - 180;
      // diff_eaB[2] = min;

      delete relc;
      relc = NULL;
    } else {
      if (debug && AlgorithmJacobianRec::aDebug)
        cout << direction
             << " J OOOPS vol neg still.. then it is CayleyBoundary" << endl;
      delete relc;
      relc = NULL;  // delete in either case
    }

    diff_TB_eaB_list.push_back(
        diff_TB_eaB);  // if volume neg, then (0,,,0) is added to diff_TB_eaB
                       // for current params
  }

  MatrixXd A(6, dim);

  if (!enteredNegVal)
    for (int i = 0; i < 6; i++)
      for (int j = 0; j < dim; j++)
        A(i, j) = diff_TB_eaB_list[j][i];  // cayleySteps[j];

  return A;
}

MatrixXd AlgorithmJacobianRec::computeJacobian(
    AtlasNode* rnode, ActiveConstraintGraph* cgK, ConvexChart* des,
    CayleyPoint* cayleyPoint, int flip, bool& succeed, int specificGridDiretion,
    MatrixXd gridSteps, double bestNorm, MatrixXd cayleyDirections, int loop,
    MatrixXd bestCayleyDirections) {
  bool debug = false;
  debug = true;
  int loopConst = 15;

  if (loop >= loopConst)  // base case
  {
    if (debug && AlgorithmJacobianRec::aDebug)
      cout << "SUCCEED FALSE!!!!!!!!!!!" << endl;

    succeed = false;
    return bestCayleyDirections;
  }

  succeed = true;

  bool enteredNegVal = false;
  if (debug && AlgorithmJacobianRec::aDebug)
    cout << "________computeJacobian " << loop << endl;
  if (debug && AlgorithmJacobianRec::aDebug)
    cout << "cayleyDirections: \n" << cayleyDirections << endl;
  int dim = des->parameters.size();

  // // cayleyDirections.conservativeResize(dim,dim); // =
  // cayleyDirections.topLeftCorner<dim, dim>(); MatrixXd temp(dim,dim); for(int
  // i=0; i<dim; i++) for(int j=0; j<dim; j++) temp(i,j) =
  // cayleyDirections(i,j); cayleyDirections = temp;

  MatrixXd prev_gridSteps = gridSteps;
  MatrixXd prev_cayleyDirections = cayleyDirections;

  bool gridStepsChange[6];
  MatrixXd J =
      jacobian(cgK, des, cayleyPoint, flip, specificGridDiretion,
               gridStepsChange, cayleyDirections, enteredNegVal, rnode);

  for (int y = 0; y < 6; y++) {
    if (gridStepsChange[y]) gridSteps(y, y) = -gridSteps(y, y);
  }

  if (enteredNegVal) {
    if (debug && AlgorithmJacobianRec::aDebug)
      cout << "enteredNegVal during compute Jacobian" << endl;

    for (size_t i = 0; i < dim; i++)
      for (size_t j = 0; j < dim; j++)
        // if(j == negVolDirection)
        cayleyDirections(i, j) = cayleyDirections(i, j) / 2;  // try halfing

    // for(size_t i = 0; i < dim; i++) // it did not work because it diminishes
    // that dimension from cayleyDirections as well.
    // {
    // cayleyDirections(i, negVolDirection) = cayleyDirections(i,
    // negVolDirection) / 2; // try halfing gridSteps(i, negVolDirection) =
    // gridSteps(i, negVolDirection) / 2;
    // }

    MatrixXd tempCayley = computeJacobian(
        rnode, cgK, des, cayleyPoint, flip, succeed, specificGridDiretion,
        gridSteps, bestNorm, cayleyDirections, ++loop, bestCayleyDirections);
    if (succeed) bestCayleyDirections = tempCayley;

  } else {
    MatrixXd Jac(dim, dim);  // remove dependent columns by setting dim x dim

    for (int i = 0; i < dim; i++) {
      for (int j = 0; j < dim; j++) {
        int d = des->independent_directions[i];
        Jac(i, j) = J(d, j);  // cayleySteps[j];
      }
    }
    if (debug && AlgorithmJacobianRec::aDebug) cout << "Jac\n" << Jac << endl;

    MatrixXd JacInv = Jac.inverse();
    if (debug && AlgorithmJacobianRec::aDebug)
      cout << "JacInv\n" << JacInv << endl;

    MatrixXd JacInv_gridsteps =
        JacInv *
        gridSteps;  // i assume jac and jacinv have same signs(or directions)
    // jac is the cartesian changes. if jac has same direction(or sign) with
    // gridSteps, then their multiplication will be positive and hence after
    // cayleyDirections * JacInv_gridsteps,  the directions of cayleyDirections
    // will be kept as it is. that is what we want. if we were getting the
    // correct directions, we want to keep them as they are.

    MatrixXd identity;
    identity.setIdentity(dim, dim);
    // MatrixXd diff = JacInv - identity;
    MatrixXd diff = JacInv_gridsteps - identity;
    double norm_cayleyMultiplier = diff.norm();

    cayleyDirections = cayleyDirections * JacInv_gridsteps;  // JacInv;
    // even if after neg vol arrangenments, the direction is reversed and hence
    // corresponding cayley direction reversed, it will have neg effect on jac
    // and hence jacInv. After multiplication cayley with jacinv, that neg
    // effect will be notralized and cayley will turn to be pos again. actually
    // you should not go neg direction of cayleyDirections, it will not give you
    // correct directions on cartesian. then how will you make sure at the
    // beginning you got correct signed cayley directions, actually you do not
    // need to know. If it is wrong cayley direction at the beginning, then jac
    // will be negative(or opposite direction than it is supposed to be), hence
    // After multiplication cayley with jacinv, that neg effect will be
    // notralized and cayley will turn to be correct direction again. and it is
    // good that we multiply JacInv with gridSteps, so that we know which
    // directions that we want to go.

    double norm_val_cay = cayleyDirections.norm();
    if (norm_val_cay < 0.1)  // if cayleyDirections goes nullish // base case
      succeed = false;
    else {
      if (norm_cayleyMultiplier <
          bestNorm) {  // since we are returning besttCayleyDirections, even if
                       // it oscillates and doesnot converge, we will return the
                       // best one on the way.
        bestNorm = norm_cayleyMultiplier;
        bestCayleyDirections = prev_cayleyDirections;  // cayleyDirections;

        besttCayleyDirections = prev_cayleyDirections;
        bestGridSteps = prev_gridSteps;  // gridSteps;
      }

      // MatrixXd diff_cayleyDirections = cayleyDirections -
      // prev_cayleyDirections;

      // double norm_diff_cayleyDirections = diff_cayleyDirections.norm();
      // if(debug && AlgorithmJacobianRec::aDebug) cout <<
      // "norm_diff_cayleyDirections " << norm_diff_cayleyDirections << endl;
      if (debug && AlgorithmJacobianRec::aDebug)
        cout << "norm_cayleyMultiplier " << norm_cayleyMultiplier << endl;

      if (norm_cayleyMultiplier >
          0.05 * dim * dim) {  // && norm_diff_cayleyDirections>0.1 // &&
                               // norm_diff_cayleyDirections>dim/2.)
        if (debug && AlgorithmJacobianRec::aDebug)
          cout << "cayleyDirections " << cayleyDirections << endl;
        if (debug && AlgorithmJacobianRec::aDebug)
          cout << "prev_cayleyDirections " << prev_cayleyDirections << endl;
        // if(debug && AlgorithmJacobianRec::aDebug)cout <<
        // "diff_cayleyDirections " << diff_cayleyDirections << endl;
        MatrixXd tempCayley =
            computeJacobian(rnode, cgK, des, cayleyPoint, flip, succeed,
                            specificGridDiretion, gridSteps, bestNorm,
                            cayleyDirections, ++loop, bestCayleyDirections);
        if (succeed) bestCayleyDirections = tempCayley;

      } else {
        if (debug && AlgorithmJacobianRec::aDebug)
          cout << "SUCCEED :))" << endl;
      }
    }
  }
  return bestCayleyDirections;
}

/*   july2014
MatrixXd AlgorithmJacobianRec::computeJacobian(ConvexChart* des, PointMultiD *
cayleyPoint, int flip, bool & succeed, int specificGridDiretion, MatrixXd
gridSteps, double bestNorm,  MatrixXd cayleyDirections, int loop, MatrixXd
bestCayleyDirections)
{
//	if(loop>1)
//		return bestCayleyDirections;


        Orientation * real = cayleyPoint->getOrientations().front();

        bool debug = false;
        debug = true;
        int loopConst = 15;

//	succeed = false;
        succeed = true;

        bool enteredNegVal = false;
        if(debug && AlgorithmJacobianRec::aDebug)  cout <<
"________computeJacobian " << loop << endl; if(debug &&
AlgorithmJacobianRec::aDebug) cout << "cayleyDirections: \n" << cayleyDirections
<< endl; int dim = des->paramsOrdered.size()/5;

//	//cayleyDirections.conservativeResize(dim,dim); //=
cayleyDirections.topLeftCorner<dim, dim>() ;
//	MatrixXd temp(dim,dim);
//	for(int i=0; i<dim; i++) for(int j=0; j<dim; j++) temp(i,j) =
cayleyDirections(i,j);
//	cayleyDirections = temp;

        MatrixXd prev_gridSteps = gridSteps;
        MatrixXd prev_cayleyDirections = cayleyDirections;

        Vector3d old_eaB = real->eaB;
        Vector3d old_TB = real->TB;

        vector<vector<double> > diff_TB_eaB_list;

        double out[6];
        cayleyPoint->getPoint(out);


        double pconnections[12][12];
        for(size_t i = 0; i<des->paramsOrdered.size(); i=i+5){
                int v1 = atoi( des->paramsOrdered.substr(i,2).c_str() ) ;
                int v2 = atoi( des->paramsOrdered.substr(i+2,2).c_str() ) ;
                des->paramconnections[v1][v2] =  out[i/5];
//real->pconnections[v1][v2]; //copy only parameters from real, not the rest of
the connections because they should be updated and changed later according to
updated params. if you copy the rest, they wont be updated in realization class
since there already exist a value for them des->paramconnections[v2][v1] =
out[i/5]; //real->pconnections[v1][v2];
        }

        for(size_t i = 0; i<des->paramsOrdered.size(); i=i+5){
                int v1 = atoi( des->paramsOrdered.substr(i,2).c_str() ) ;
                int v2 = atoi( des->paramsOrdered.substr(i+2,2).c_str() ) ;
                if(debug && AlgorithmJacobianRec::aDebug)  cout << i/5 << ": "
<< v1 << "-" << v2 << " " << des->paramconnections[v1][v2] << ",  ";
        }
        if(debug && AlgorithmJacobianRec::aDebug)  cout << "TB_EAB " <<
old_TB[0] << " " << old_TB[1] << " " << old_TB[2] << " " << old_eaB[0] << " " <<
old_eaB[1] << " " << old_eaB[2] << " " << endl;


        int negVolDirection = -1;
        for(size_t direction = 0; direction<des->paramsOrdered.size()/5;
direction++)
        {
                for(int k = 0; k<12; k++)	for(int j = 0; j<12; j++)
pconnections[k][j] = des->paramconnections[k][j]; for(size_t i = 0;
i<des->paramsOrdered.size(); i=i+5)
                {
                        int v1 = atoi( des->paramsOrdered.substr(i,2).c_str() )
; int v2 = atoi( des->paramsOrdered.substr(i+2,2).c_str() ) ;
                        pconnections[v1][v2] += cayleyDirections(i/5,
direction);//  cayleySteps[i/5]; //Settings::Sampling::stepSize;
                        pconnections[v2][v1] = pconnections[v1][v2];
                }

                //CartesianRealizer *relc = new CartesianRealizer(des, flip,
pconnections); bool fail; Orientation *relc =
CartesianRealizer::computeRealization(des, flip, pconnections, fail); if( fail )
//volume neg
                {
//			negVolDirection = direction;
//			enteredNegVal = true;
//			delete relc; relc = NULL;
//			break;

                        if(specificGridDiretion != direction &&
specificGridDiretion != -1)
                        {
                                for(size_t i = 0; i<des->paramsOrdered.size();
i=i+5)
                                {
                                        int v1 = atoi(
des->paramsOrdered.substr(i,2).c_str() ) ; int v2 = atoi(
des->paramsOrdered.substr(i+2,2).c_str() ) ; pconnections[v1][v2] -= 2*
cayleyDirections(i/5, direction); //try reverse direction pconnections[v2][v1] =
pconnections[v1][v2];
                                }

                                delete relc; relc = NULL;
//				relc = new CartesianRealizer(des, flip,
pconnections); relc = CartesianRealizer::computeRealization(des, flip,
pconnections, fail); if( !fail ){ for(size_t i = 0; i<des->paramsOrdered.size();
i=i+5) cayleyDirections(i/5, direction) = -cayleyDirections(i/5, direction);
//try reverse direction

                                        gridSteps(direction, direction) =
-gridSteps(direction, direction);
                                }
                        }

                        if( fail ) //volume neg again
                        {
                                negVolDirection = direction;
                                enteredNegVal = true;
                                delete relc; relc = NULL;
                                break;

                                if(debug && AlgorithmJacobianRec::aDebug)  cout
<< "J volume negative when increased by stepsize" << endl;

//				for(int k = 0; k<12; k++)	for(int j = 0;
j<12; j++)	pconnections[k][j] = des->paramconnections[k][j];
//des->paramconnections[i][j];

//				des->updateBoundaries();
//				for(size_t i = 0; i<des->paramsOrdered.size();
i=i+5) for(int i = des->paramsOrdered.size()-5; i>=0; i=i-5) //start change from
reverse direction, because of the order of dependence
                                {
                                        //des->updateBoundaries();
                                        des->updateBoundariesSpecific(i/5);
                                        int v1 = atoi(
des->paramsOrdered.substr(i,2).c_str() ) ; int v2 = atoi(
des->paramsOrdered.substr(i+2,2).c_str() ) ;

                                        double min =
des->minconnections[v1][v2]; double max = des->maxconnections[v1][v2];


                                        double prev =
des->paramconnections[v1][v2]; double curr = pconnections[v1][v2]; if(debug &&
AlgorithmJacobianRec::aDebug)  cout << "J " << i/5 << ": " << v1 << "-" << v2 <<
" max:" << max << " min:" << min  << " prev:"  << prev << " curr:" << curr ;


                                        if( curr > max ){
                                                pconnections[v1][v2]  = max ;
                                                pconnections[v2][v1]  = max ;
                                                if(max == prev) //make it
different than prev point
                                                {
                                                        pconnections[v1][v2]  =
max - Settings::Sampling::stepSize; pconnections[v2][v1]  = max -
Settings::Sampling::stepSize;
                                                }
                                        }
                                        else if( curr < min ){
                                                pconnections[v1][v2]  = min ;
                                                pconnections[v2][v1]  = min ;
                                                if(min == prev) //make it
different than prev point
                                                {
                                                        pconnections[v1][v2]  =
min + Settings::Sampling::stepSize; pconnections[v2][v1]  = min +
Settings::Sampling::stepSize;
                                                }
                                        }
                                        if(debug &&
AlgorithmJacobianRec::aDebug)  cout << " now:" << pconnections[v1][v2] << endl;

                                        des->paramconnections[v1][v2] =
pconnections[v1][v2]; des->paramconnections[v2][v1] = pconnections[v1][v2];


                                }



                                //take it back the changes you have just done
above. for(size_t i = 0; i<des->paramsOrdered.size(); i=i+5){ int v1 = atoi(
des->paramsOrdered.substr(i,2).c_str() ) ; int v2 = atoi(
des->paramsOrdered.substr(i+2,2).c_str() ) ; des->paramconnections[v1][v2] =
out[i/5]; //real->pconnections[v1][v2]; //copy only parameters from real, not
the rest of the connections because they should be updated and changed later
according to updated params. if you copy the rest, they wont be updated in
realization class since there already exist a value for them
                                        des->paramconnections[v2][v1] =
out[i/5]; //real->pconnections[v1][v2];
                                }



                                bool fail;
                                delete relc; relc = NULL;
//				relc = new CartesianRealizer(des, flip,
pconnections); relc = CartesianRealizer::computeRealization(des, flip,
pconnections, fail);


                                if( !fail )
                                {
                                        for(size_t i = 0;
i<des->paramsOrdered.size(); i=i+5){ int v1 = atoi(
des->paramsOrdered.substr(i,2).c_str() ) ; int v2 = atoi(
des->paramsOrdered.substr(i+2,2).c_str() ) ; cayleyDirections(i/5, direction) =
pconnections[v1][v2] - des->paramconnections[v1][v2];
                                        }

                                        if(debug &&
AlgorithmJacobianRec::aDebug) cout << "cayleyDirections after neg vol fix: \n"
<< cayleyDirections << endl;

                                }
                                else
                                {
                                        negVolDirection = direction;
                                        enteredNegVal = true;
                                        delete relc; relc = NULL;
                                        break;
                                }

                        }
                }

                if(debug && AlgorithmJacobianRec::aDebug)  cout << "TB_EAB after
walking on cayley combined direction" << direction << " is: " << relc->TB[0] <<
" " << relc->TB[1] << " " << relc->TB[2] << " " << relc->eaB[0] << " " <<
relc->eaB[1] << " " << relc->eaB[2] << " " << endl;


                vector<double> diff_TB_eaB(6); /// 6 zero-initialized elements
                if( !fail ) //do not check angle now, maybe when you change the
direction of sampling using jacobian, the new realization will have good angle
                {
                        Vector3d eaB = relc->eaB;
                        Vector3d TB = relc->TB;

                        Vector3d diff_TB = TB - old_TB;
                        Vector3d diff_eaB = eaB - old_eaB;
                        //if(debug && AlgorithmJacobianRec::aDebug)  cout << "J
" << eaB[0] << " " << old_eaB[0] << " " << diff_eaB[0] << " " ; for(int h=0;
h<3; h++) if( abs(diff_eaB[h]) > 180 )  diff_eaB[h] = -1 *
Utils::sign(diff_eaB[h]) * (360 - abs(diff_eaB[h]) );//to keep it between (-180,
180)
                        //for(int h=0; h<3; h++)   diff_eaB[h] =  ( (int)
diff_eaB[h]+180)%360 - 180 ;//to keep it between (-180, 180)
                        //if(debug && AlgorithmJacobianRec::aDebug)  cout <<
diff_eaB[0] << endl;

                        for(int h=0; h<3; h++) diff_TB_eaB[h] = diff_TB[h];
                        for(int h=0; h<3; h++) diff_TB_eaB[h+3] = diff_eaB[h];



//			double min = diff_TB_eaB[3];
//			if( abs(diff_TB_eaB[3]+180) < abs(min) ) min =
diff_TB_eaB[3]+180;
//			if( abs(diff_TB_eaB[3]-180) < abs(min) ) min =
diff_TB_eaB[3]-180;
//			diff_TB_eaB[3] = min;

//			min = diff_eaB[2];
//			if(diff_eaB[2]+180 < min) min = diff_eaB[2]+180;
//			if(diff_eaB[2]-180 < min) min = diff_eaB[2]-180;
//			diff_eaB[2] = min;



                        delete relc; relc = NULL;
                }
                else
                {
                        if(debug && AlgorithmJacobianRec::aDebug)  cout <<
direction << " J OOOPS vol neg still.. then it is CayleyBoundary" << endl;
                        delete relc; relc = NULL; // delete in either case
                }

                diff_TB_eaB_list.push_back(diff_TB_eaB); //if volume neg, then
(0,,,0) is added to diff_TB_eaB for current params


        }



        if(enteredNegVal && loop<loopConst )
        {

                if(debug && AlgorithmJacobianRec::aDebug)  cout <<
"enteredNegVal during compute Jacobian" << endl;

                for(size_t i = 0; i<dim; i++)
                        for(size_t j = 0; j<dim; j++)
        //			if(j==negVolDirection)
                                cayleyDirections(i, j) = cayleyDirections(i,
j)/2; //try halfing


//		for(size_t i = 0; i<dim; i++) //it did not work because it
diminishes that dimension from cayleyDirections as well.
//		{
//			cayleyDirections(i, negVolDirection) =
cayleyDirections(i, negVolDirection)/2; //try halfing
//			gridSteps(i, negVolDirection) = gridSteps(i,
negVolDirection)/2;
//		}


                MatrixXd tempCayley = computeJacobian(des, cayleyPoint,  flip,
succeed, specificGridDiretion, gridSteps, bestNorm, cayleyDirections, ++loop,
bestCayleyDirections); if(succeed){ cayleyDirections = tempCayley;
                        bestCayleyDirections = tempCayley;
                }

        }
        else if(!enteredNegVal)  //it may enter here when enteredNegVal (if
recursive computeJacobian keep entering NegVal), when this is the case
        {


                if(des->independent_directions.size()==0)
                {
                        MatrixXd A(6,dim);
                        for(int i=0; i<6; i++)
                                for(int j=0; j<dim; j++)
                                        A(i,j) = diff_TB_eaB_list[j][i]; //
/cayleySteps[j];

                        vector<int> indep = computeIndependentColumns(A, dim) ;
                        des->independent_directions.clear();
                        des->independent_directions = indep;

                        if(dim == 3)
                        {
                                des->independent_directions.clear();
                                des->independent_directions.push_back(0);
                                des->independent_directions.push_back(1);
                                des->independent_directions.push_back(2);
                        }
                        if(dim == 5)
                        {
                                des->independent_directions.clear();
                                des->independent_directions.push_back(0);
                                des->independent_directions.push_back(1);
                                des->independent_directions.push_back(2);
                                des->independent_directions.push_back(3);
                                des->independent_directions.push_back(4);
                        }
                        if(dim == 6)
                        {
                                des->independent_directions.clear();
                                des->independent_directions.push_back(0);
                                des->independent_directions.push_back(1);
                                des->independent_directions.push_back(2);
                                des->independent_directions.push_back(3);
                                des->independent_directions.push_back(4);
                                des->independent_directions.push_back(5);
                        }
                }



                MatrixXd Jac(dim,dim); //remove dependent columns by setting dim
x dim

                for(int i=0; i < dim; i++)
                {
                        for(int j=0; j < dim; j++){
                                int d = des->independent_directions[i];
                                Jac(i,j) = diff_TB_eaB_list[j][d];// /
cayleySteps[j];
                        }
                }
                if(debug && AlgorithmJacobianRec::aDebug)  cout <<  "Jac \n" <<
Jac << endl;

                MatrixXd JacInv = Jac.inverse();
                if(debug && AlgorithmJacobianRec::aDebug)  cout <<  "JacInv \n"
<< JacInv << endl;


                MatrixXd JacInv_gridsteps = JacInv * gridSteps; //i assume jac
and jacinv have same signs(or directions)
                //jac is the cartesian changes. if jac has same direction(or
sign) with gridSteps, then their multiplication will be positive
                //and hence after cayleyDirections * JacInv_gridsteps,  the
directions of cayleyDirections will be kept as it is. that is what we want. if
we were getting the correct directions, we want to keep them as they are.


                MatrixXd identity;
                identity.setIdentity(dim,dim);
//		MatrixXd diff = JacInv - identity;
                MatrixXd diff = JacInv_gridsteps - identity;
                double norm_val = diff.norm();

                cayleyDirections = cayleyDirections * JacInv_gridsteps;
//JacInv; //even if after neg vol arrangenments, the direction is reversed and
hence corresponding cayley direction reversed, it will have neg effect on jac
and hence jacInv. After multiplication cayley with jacinv, that neg effect will
be notralized and cayley will turn to be pos again.
                //actually you should not go neg direction of cayleyDirections,
it will not give you correct directions on cartesian.
                //then how will you make sure at the beginning you got correct
signed cayley directions, actually you do not need to know. If it is wrong
cayley direction at the beginning, then jac will be negative(or opposite
direction than it is supposed to be), hence  After multiplication cayley with
jacinv, that neg effect will be notralized and cayley will turn to be correct
direction again.
                //and it is good that we multiply JacInv with gridSteps, so that
we know which directions that we want to go.


                if(norm_val < bestNorm)
                {
                        bestNorm = norm_val;
                        bestCayleyDirections =  prev_cayleyDirections;
//cayleyDirections; //


                        double norm_val_cay = prev_cayleyDirections.norm();
                        if(norm_val_cay >= 0.1) //if cayleyDirections does not
go nullish
                        {
                                besttCayleyDirections =  prev_cayleyDirections;
                                bestGridSteps = prev_gridSteps; //gridSteps;  //
                        }
                }


                MatrixXd diff_cayleyDirections = cayleyDirections -
prev_cayleyDirections;

                double norm_val_cayleyDirections = diff_cayleyDirections.norm();
                if(debug && AlgorithmJacobianRec::aDebug) cout <<
"norm_val_cayleyDirections " << norm_val_cayleyDirections << endl; if(debug &&
AlgorithmJacobianRec::aDebug) cout << "norm_val " << norm_val << endl;
                if(norm_val > 0.05*dim*dim  && loop<loopConst  &&
norm_val_cayleyDirections>0.1)  // && norm_val_cayleyDirections>dim/2.)
                {
                        if(debug && AlgorithmJacobianRec::aDebug)cout <<
"cayleyDirections " << cayleyDirections << endl; if(debug &&
AlgorithmJacobianRec::aDebug)cout << "prev_cayleyDirections " <<
prev_cayleyDirections << endl; if(debug && AlgorithmJacobianRec::aDebug)cout <<
"diff_cayleyDirections " << diff_cayleyDirections << endl; MatrixXd tempCayley =
computeJacobian(des, cayleyPoint,  flip, succeed, specificGridDiretion,
gridSteps, bestNorm, cayleyDirections, ++loop, bestCayleyDirections);
                        if(succeed){
                                cayleyDirections = tempCayley;
                                bestCayleyDirections = tempCayley;
                                //bestGridSteps = gridSteps;
                        }
                }
                else{
                        cayleyDirections = prev_cayleyDirections;
                        if(norm_val < 0.05*dim*dim){
                                if(debug && AlgorithmJacobianRec::aDebug) cout
<< "SUCCEED :))" << endl;
                        }
                        else{
//				succeed = false;
                                if(debug && AlgorithmJacobianRec::aDebug) cout
<< "not SUCCEED :((" << endl;
                        }
                }

//		if(norm_val < 0.05*dim*dim)
//			succeed = true;


        }
        else{
                succeed = false;
                if(debug && AlgorithmJacobianRec::aDebug) cout << "SUCCEED
FALSE!!!!!!!!!!!" << endl;

//		MatrixXd JacInv(dim,dim);
//		for(int i=0; i < dim; i++)
//			for(int j=0; j < dim; j++)
//				JacInv(i,j) = 0;
//		return JacInv;
                return bestCayleyDirections;
        }



//	AlgorithmJacobianRec::cayleyDirections=cayleyDirections;

        double norm_val_cay = cayleyDirections.norm();
        if(norm_val_cay< 0.1) //if cayleyDirections goes nullish
                succeed = false;

//	cout <<  "final cayleyDirections \n" << cayleyDirections << endl <<
endl;


//	return cayleyDirections;
        return bestCayleyDirections;

}
*/

VectorXd AlgorithmJacobianRec::computeJacobian1Dir(
    ActiveConstraintGraph* cgK, ConvexChart* des, CayleyPoint* cayleyPoint,
    int flip, bool& succeed, int specificGridDiretion, VectorXd gridSteps,
    VectorXd cayleyDirections, int loop, double bestNorm,
    VectorXd bestCayleyDirections, AtlasNode* rnode) {
  Settings* set = Settings::getInstance();
  bool debug = false;
  debug = true;

  //	succeed = false;
  succeed = true;

  Orientation* real = cayleyPoint->getOrientations().front();

  bool enteredNegVal = false;
  if (debug && AlgorithmJacobianRec::aDebug)
    cout << "________computeJacobian 1 direction " << loop << endl;
  if (debug && AlgorithmJacobianRec::aDebug)
    cout << "cayleyDirections: \n" << cayleyDirections << endl;
  int dim = des->parameters.size();

  //	bool zerofound = false;
  //	for(size_t i = 0; i<dim; i++)
  //			for(size_t j = 0; j<dim; j++)
  //				if(AlgorithmJacobianRec::cayleyDirections(i, j) == 0
  //) 					zerofound = true;

  //	//cayleyDirections.conservativeResize(dim,dim); //=
  //cayleyDirections.topLeftCorner<dim, dim>() ; 	MatrixXd temp(dim,dim); 	for(int
  //i=0; i<dim; i++) for(int j=0; j<dim; j++) temp(i,j) = cayleyDirections(i,j);
  //	cayleyDirections = temp;

  VectorXd prev_cayleyDirections = cayleyDirections;

  Vector3d old_eaB = real->eaB;
  Vector3d old_TB = real->TB;

  vector<vector<double> > diff_TB_eaB_list;

  double out[6];
  cayleyPoint->getPoint(out);

  double pconnections[12][12];
  for (size_t i = 0; i < des->parameters.size(); i++) {
    int v1 = des->parameters[i].first;
    int v2 = des->parameters[i].second;
    des->param_length[v1][v2] =
        out[i];  // real->pconnections[v1][v2]; //copy only parameters from
                 // real, not the rest of the connections because they should be
                 // updated and changed later according to updated params. if you
                 // copy the rest, they wont be updated in realization class
                 // since there already exist a value for them
    des->param_length[v2][v1] = out[i];  // real->pconnections[v1][v2];
  }

  for (size_t i = 0; i < des->parameters.size(); i++) {
    int v1 = des->parameters[i].first;
    int v2 = des->parameters[i].second;
    if (debug && AlgorithmJacobianRec::aDebug)
      cout << i << ": " << v1 << "-" << v2 << " " << des->param_length[v1][v2]
           << ",  ";
  }
  if (debug && AlgorithmJacobianRec::aDebug)
    cout << "TB_EAB " << old_TB[0] << " " << old_TB[1] << " " << old_TB[2]
         << " " << old_eaB[0] << " " << old_eaB[1] << " " << old_eaB[2] << " "
         << endl;

  int negVolDirection = -1;
  for (size_t direction = 0; direction < des->parameters.size(); direction++) {
    for (int k = 0; k < 12; k++)
      for (int j = 0; j < 12; j++) pconnections[k][j] = des->param_length[k][j];

    int v1 = des->parameters[direction].first;
    int v2 = des->parameters[direction].second;
    pconnections[v1][v2] += cayleyDirections(
        direction);  //  cayleySteps[i/5]; //Settings::Sampling::stepSize;
    pconnections[v2][v1] = pconnections[v1][v2];

    // CartesianRealizer *relc = new CartesianRealizer(des, flip, pconnections);
    bool fail;
    Orientation* relc = CartesianRealizer::computeRealization(
        cgK, des, flip, pconnections, fail, rnode->getFlipScheme());
    if (fail)  // volume neg
    {
      negVolDirection = direction;
      enteredNegVal = true;
      delete relc;
      relc = NULL;
      break;

      //			if(specificGridDiretion != direction)
      //			{
      //				for(size_t i = 0;
      //i<des->paramsOrdered.size(); i=i+5)
      //				{
      //					int v1 = atoi(
      //des->paramsOrdered.substr(i,2).c_str() ) ; 					int v2 = atoi(
      //des->paramsOrdered.substr(i+2,2).c_str() ) ; 					pconnections[v1][v2] -= 2*
      //cayleyDirections(i/5, direction); //try reverse direction
      //					pconnections[v2][v1] =
      //pconnections[v1][v2];
      //				}
      //
      //				delete relc; relc = NULL;
      //				relc = new CartesianRealizer(des, flip,
      //pconnections); 				if( relc->getVolumeTest() ){ 					for(size_t i = 0;
      //i<des->paramsOrdered.size(); i=i+5) 						cayleyDirections(i/5, direction) =
      //-cayleyDirections(i/5, direction); //try reverse direction
      //				}
      //			}

      if (fail)  // volume neg again
      {
        //				enteredNegVal = true;
        //				delete relc; relc = NULL;
        //				break;

        if (debug && AlgorithmJacobianRec::aDebug)
          cout << "J volume negative when increased by stepsize" << endl;

        //				for(int k = 0; k<12; k++)	for(int j = 0;
        //j<12; j++)	pconnections[k][j] = des->paramconnections[k][j];
        ////des->paramconnections[i][j];

        //				des->updateBoundaries();
        //				for(size_t i = 0; i<des->paramsOrdered.size();
        //i=i+5)
        for (int i = des->parameters.size() - 1; i >= 0;
             i--)  // start change from reverse direction, because of the order
                   // of dependence
        {
          // des->updateBoundaries();
          des->updateBoundariesSpecific(i);
          int v1 = des->parameters[i].first;
          int v2 = des->parameters[i].second;

          double min = des->param_lengthLower[v1][v2];
          double max = des->param_lengthUpper[v1][v2];

          double prev = des->param_length[v1][v2];
          double curr = pconnections[v1][v2];
          if (debug && AlgorithmJacobianRec::aDebug)
            cout << "J " << i << ": " << v1 << "-" << v2 << " max:" << max
                 << " min:" << min << " prev:" << prev << " curr:" << curr;

          if (curr > max) {
            pconnections[v1][v2] = max;
            pconnections[v2][v1] = max;
            if (max == prev)  // make it different than prev point
            {
              pconnections[v1][v2] = max - set->Sampling.stepSize;
              pconnections[v2][v1] = max - set->Sampling.stepSize;
            }
          } else if (curr < min) {
            pconnections[v1][v2] = min;
            pconnections[v2][v1] = min;
            if (min == prev)  // make it different than prev point
            {
              pconnections[v1][v2] = min + set->Sampling.stepSize;
              pconnections[v2][v1] = min + set->Sampling.stepSize;
            }
          }
          if (debug && AlgorithmJacobianRec::aDebug)
            cout << " now:" << pconnections[v1][v2] << endl;

          des->param_length[v1][v2] = pconnections[v1][v2];
          des->param_length[v2][v1] = pconnections[v1][v2];
        }

        // take it back the changes you have just done above.
        for (size_t i = 0; i < des->parameters.size(); i++) {
          int v1 = des->parameters[i].first;
          int v2 = des->parameters[i].second;
          des->param_length[v1][v2] =
              out[i];  // real->pconnections[v1][v2]; //copy only parameters
                       // from real, not the rest of the connections because they
                       // should be updated and changed later according to
                       // updated params. if you copy the rest, they wont be
                       // updated in realization class since there already exist
                       // a value for them
          des->param_length[v2][v1] = out[i];  // real->pconnections[v1][v2];
        }

        bool fail;
        delete relc;
        relc = NULL;
        //				relc = new CartesianRealizer(des, flip,
        //pconnections);
        relc = CartesianRealizer::computeRealization(
            cgK, des, flip, pconnections, fail, rnode->getFlipScheme());

        if (!fail) {
          for (size_t i = 0; i < des->parameters.size(); i++) {
            int v1 = des->parameters[i].first;
            int v2 = des->parameters[i].second;
            cayleyDirections(i, direction) =
                pconnections[v1][v2] - des->param_length[v1][v2];
          }

          if (debug && AlgorithmJacobianRec::aDebug)
            cout << "cayleyDirections after neg vol fix: \n"
                 << cayleyDirections << endl;

        } else {
          negVolDirection = direction;
          enteredNegVal = true;
          delete relc;
          relc = NULL;
          break;
        }
      }
    }

    if (debug && AlgorithmJacobianRec::aDebug)
      cout << "TB_EAB after walking on cayley combined direction" << direction
           << " is: " << relc->TB[0] << " " << relc->TB[1] << " " << relc->TB[2]
           << " " << relc->eaB[0] << " " << relc->eaB[1] << " " << relc->eaB[2]
           << " " << endl;

    vector<double> diff_TB_eaB(6);  /// 6 zero-initialized elements
    if (!fail)  // do not check angle now, maybe when you change the direction
                // of sampling using jacobian, the new realization will have good
                // angle
    {
      Vector3d eaB = relc->eaB;
      Vector3d TB = relc->TB;

      Vector3d diff_TB = TB - old_TB;
      Vector3d diff_eaB = eaB - old_eaB;
      // if(debug && AlgorithmJacobianRec::aDebug)  cout << "J " << eaB[0] << "
      // " << old_eaB[0] << " " << diff_eaB[0] << " " ;
      for (int h = 0; h < 3; h++)
        if (abs(diff_eaB[h]) > 180)
          diff_eaB[h] =
              -1 * Utils::sign(diff_eaB[h]) *
              (360 - abs(diff_eaB[h]));  // to keep it between (-180, 180)
      // for(int h=0; h<3; h++)   diff_eaB[h] =  ( (int) diff_eaB[h]+180)%360 -
      // 180 ;//to keep it between (-180, 180) if(debug &&
      // AlgorithmJacobianRec::aDebug)  cout << diff_eaB[0] << endl;

      for (int h = 0; h < 3; h++) diff_TB_eaB[h] = diff_TB[h];
      for (int h = 0; h < 3; h++) diff_TB_eaB[h + 3] = diff_eaB[h];

      //			double min = diff_TB_eaB[3];
      //			if( abs(diff_TB_eaB[3]+180) < abs(min) ) min =
      //diff_TB_eaB[3]+180; 			if( abs(diff_TB_eaB[3]-180) < abs(min) ) min =
      //diff_TB_eaB[3]-180; 			diff_TB_eaB[3] = min;

      //			min = diff_eaB[2];
      //			if(diff_eaB[2]+180 < min) min = diff_eaB[2]+180;
      //			if(diff_eaB[2]-180 < min) min = diff_eaB[2]-180;
      //			diff_eaB[2] = min;

      delete relc;
      relc = NULL;
    } else {
      if (debug && AlgorithmJacobianRec::aDebug)
        cout << direction
             << " J OOOPS vol neg still.. then it is CayleyBoundary" << endl;
      delete relc;
      relc = NULL;  // delete in either case
    }

    diff_TB_eaB_list.push_back(
        diff_TB_eaB);  // if volume neg, then (0,,,0) is added to diff_TB_eaB
                       // for current params
  }

  if (enteredNegVal && loop < 15) {
    if (debug && AlgorithmJacobianRec::aDebug)
      cout << "enteredNegVal during compute Jacobian" << endl;

    cayleyDirections(negVolDirection) =
        cayleyDirections(negVolDirection) / 2;  // try halfing

    VectorXd tempCayley = computeJacobian1Dir(
        cgK, des, cayleyPoint, flip, succeed, specificGridDiretion, gridSteps,
        cayleyDirections, ++loop, bestNorm, bestCayleyDirections, rnode);
    if (succeed) {
      cayleyDirections = tempCayley;
      bestCayleyDirections = tempCayley;
    }

  } else if (!enteredNegVal)  // it may enter here when enteredNegVal (if
                              // recursive computeJacobian keep entering NegVal),
                              // when this is the case
  {
    MatrixXd Jac(dim, dim);  // remove dependent columns by setting dim x dim

    for (int i = 0; i < dim; i++) {
      for (int j = 0; j < dim; j++) {
        int d = des->independent_directions[i];
        Jac(i, j) = diff_TB_eaB_list[j][d];  // / cayleySteps[j];
      }
    }
    if (debug && AlgorithmJacobianRec::aDebug) cout << "Jac \n" << Jac << endl;

    MatrixXd JacInv = Jac.inverse();
    if (debug && AlgorithmJacobianRec::aDebug)
      cout << "JacInv \n" << JacInv << endl;

    VectorXd JacInv_gridsteps =
        JacInv *
        gridSteps;  // i assume jac and jacinv have same signs(or directions)
    // jac is the cartesian changes. if jac has same direction(or sign) with
    // gridSteps, then their multiplication will be positive and hence after
    // cayleyDirections * JacInv_gridsteps,  the directions of cayleyDirections
    // will be kept as it is. that is what we want. if we were getting the
    // correct directions, we want to keep them as they are.

    VectorXd identity;
    identity.setOnes(dim);
    //		MatrixXd diff = JacInv - identity;
    VectorXd diff = JacInv_gridsteps - identity;
    double norm_val = diff.norm();

    if (norm_val < bestNorm) {
      bestNorm = norm_val;
      bestCayleyDirections = prev_cayleyDirections;
    }

    for (int s = 0; s < dim; s++)
      cayleyDirections(s) =
          cayleyDirections(s) *
          JacInv_gridsteps(
              s);  // JacInv; //even if after neg vol arrangenments, the
                   // direction is reversed and hence corresponding cayley
                   // direction reversed, it will have neg effect on jac and
                   // hence jacInv. After multiplication cayley with jacinv, that
                   // neg effect will be notralized and cayley will turn to be
                   // pos again.
    // actually you should not go neg direction of cayleyDirections, it will not
    // give you correct directions on cartesian. then how will you make sure at
    // the beginning you got correct signed cayley directions, actually you do
    // not need to know. If it is wrong cayley direction at the beginning, then
    // jac will be negative(or opposite direction than it is supposed to be),
    // hence  After multiplication cayley with jacinv, that neg effect will be
    // notralized and cayley will turn to be correct direction again. and it is
    // good that we multiply JacInv with gridSteps, so that we know which
    // directions that we want to go.

    VectorXd diff_cayleyDirections = cayleyDirections - prev_cayleyDirections;

    double norm_val_cayleyDirections = diff_cayleyDirections.norm();
    if (debug && AlgorithmJacobianRec::aDebug)
      cout << "norm_val_cayleyDirections " << norm_val_cayleyDirections << endl;
    if (debug && AlgorithmJacobianRec::aDebug)
      cout << "norm_val " << norm_val << endl;
    if (norm_val > 0.05 * dim * dim && loop < 15 &&
        norm_val_cayleyDirections >
            0.1)  // && norm_val_cayleyDirections>dim/2.)
    {
      if (debug && AlgorithmJacobianRec::aDebug)
        cout << "cayleyDirections " << cayleyDirections << endl;
      if (debug && AlgorithmJacobianRec::aDebug)
        cout << "prev_cayleyDirections " << prev_cayleyDirections << endl;
      if (debug && AlgorithmJacobianRec::aDebug)
        cout << "diff_cayleyDirections " << diff_cayleyDirections << endl;
      VectorXd tempCayley = computeJacobian1Dir(
          cgK, des, cayleyPoint, flip, succeed, specificGridDiretion, gridSteps,
          cayleyDirections, ++loop, bestNorm, bestCayleyDirections, rnode);
      if (succeed) {
        cayleyDirections = tempCayley;
        bestCayleyDirections = tempCayley;
      }
    } else {
      cayleyDirections = prev_cayleyDirections;
      if (norm_val < 0.05 * dim * dim) {
        if (debug && AlgorithmJacobianRec::aDebug)
          cout << "SUCCEED :))" << endl;
      } else {
        //				succeed = false;
        if (debug && AlgorithmJacobianRec::aDebug)
          cout << "not SUCCEED :((" << endl;
      }
    }

    //		if(norm_val < 0.05*dim*dim)
    //			succeed = true;

  } else {
    succeed = false;
    if (debug && AlgorithmJacobianRec::aDebug)
      cout << "SUCCEED FALSE!!!!!!!!!!!" << endl;

    //		MatrixXd JacInv(dim,dim);
    //		for(int i=0; i < dim; i++)
    //			for(int j=0; j < dim; j++)
    //				JacInv(i,j) = 0;
    //		return JacInv;
    return bestCayleyDirections;
  }

  //	AlgorithmJacobianRec::cayleyDirections=cayleyDirections;

  double norm_val_cay = cayleyDirections.norm();
  if (norm_val_cay < 0.1)  // if cayleyDirections goes nullish
    succeed = false;

  //	cout <<  "final cayleyDirections \n" << cayleyDirections << endl <<
  //endl;

  //	return cayleyDirections;
  return bestCayleyDirections;
}

/*
MatrixXd AlgorithmJacobianRec::computeJacobian3(ConvexChart* des,
CartesianRealizer * real, int flip, bool & succeed, double cayleySteps[6])
{
        bool debug = false;
        debug = true;

        succeed = true;
        bool enteredNegVal = false;
        if(debug && AlgorithmJacobianRec::aDebug)  cout <<
"________computeJacobian " << endl; int dim = des->paramsOrdered.size()/5;
        //double cayleySteps[6] =
{Settings::Sampling::stepSize,Settings::Sampling::stepSize,Settings::Sampling::stepSize,
Settings::Sampling::stepSize, Settings::Sampling::stepSize,
Settings::Sampling::stepSize};

        Vector3d old_eaB = real->eaB;
        Vector3d old_TB = real->TB;

        vector<vector<double> > diff_TB_eaB_list;

        double pconnections[12][12];
        for(size_t i = 0; i<des->paramsOrdered.size(); i=i+5){
                int v1 = atoi( des->paramsOrdered.substr(i,2).c_str() ) ;
                int v2 = atoi( des->paramsOrdered.substr(i+2,2).c_str() ) ;
                des->paramconnections[v1][v2] = real->pconnections[v1][v2];
//copy only parameters from real, not the rest of the connections because they
should be updated and changed later according to updated params. if you copy the
rest, they wont be updated in realization class since there already exist a
value for them des->paramconnections[v2][v1] = real->pconnections[v1][v2];
        }
        des->updateBoundaries();


        for(size_t i = 0; i<des->paramsOrdered.size(); i=i+5){
                int v1 = atoi( des->paramsOrdered.substr(i,2).c_str() ) ;
                int v2 = atoi( des->paramsOrdered.substr(i+2,2).c_str() ) ;
                if(debug && AlgorithmJacobianRec::aDebug)  cout << i/5 << ": "
<< v1 << "-" << v2 << " " << des->paramconnections[v1][v2] << ",  ";
        }
        if(debug && AlgorithmJacobianRec::aDebug)  cout << "TB_EAB " <<
old_TB[0] << " " << old_TB[1] << " " << old_TB[2] << " " << old_eaB[0] << " " <<
old_eaB[1] << " " << old_eaB[2] << " " << endl;


        for(size_t i = 0; i<des->paramsOrdered.size(); i=i+5)
        {
                for(int k = 0; k<12; k++)	for(int j = 0; j<12; j++)
pconnections[k][j] = des->paramconnections[k][j];

                int v1 = atoi( des->paramsOrdered.substr(i,2).c_str() ) ;
                int v2 = atoi( des->paramsOrdered.substr(i+2,2).c_str() ) ;
                pconnections[v1][v2] += cayleySteps[i/5];
//Settings::Sampling::stepSize; pconnections[v2][v1] = pconnections[v1][v2];


                CartesianRealizer *relc = new CartesianRealizer(des, flip,
pconnections); if( !relc->getVolumeTest() ) //volume neg
                {
                        enteredNegVal = true;
                        delete relc; relc = NULL;
                        break;

//			pconnections[v1][v2] -= 2*Settings::Sampling::stepSize;
//try reverse direction
//			pconnections[v2][v1] = pconnections[v1][v2];
//			delete relc; relc = NULL;
//			relc = new CartesianRealizer(des, flip, pconnections);
//			if( relc->getVolumeTest() )
//				cayleySteps[i/5] =
-Settings::Sampling::stepSize; // SHOULD NOT IT BE   -
Settings::Sampling::stepSize ???
//			else //volume neg again
//			{
//
//				enteredNegVal = true;
//				delete relc; relc = NULL;
//				break;
//
////				if(debug && AlgorithmJacobianRec::aDebug) cout
<< "J volume negative when increased by stepsize" << endl;
////
////				//des->updateBoundaries();
////				double min = des->minconnections[v1][v2];
////				double max = des->maxconnections[v1][v2];
////				double prev = des->paramconnections[v1][v2];
////				double curr = pconnections[v1][v2];
////				if(debug && AlgorithmJacobianRec::aDebug) cout
<< "J " << i/5 << ": " << v1 << "-" << v2 << " max:" << max << " min:" << min <<
" prev:"  << prev << " curr:" << curr ;
////
////				for(int k = 0; k<12; k++)	for(int j = 0;
j<12; j++)	pconnections[k][j] = des->paramconnections[k][j];
//des->paramconnections[i][j];
////
////				if(abs(max-curr) < abs(curr-min) ){
////					pconnections[v1][v2]  = max ;
////					pconnections[v2][v1]  = max ;
////				}
////				else{
////					pconnections[v1][v2]  = min ;
////					pconnections[v2][v1]  = min ;
////				}
////
////				cayleySteps[i/5] = pconnections[v1][v2] - prev;
////				delete relc; relc = NULL;
////				relc = new CartesianRealizer(des, flip,
pconnections);
////
////				if(debug && AlgorithmJacobianRec::aDebug) cout
<< " now:" << pconnections[v1][v2] << endl;
//			}
                }

                if(debug && AlgorithmJacobianRec::aDebug)  cout << "TB_EAB after
walking on cayley " << i/5 << " is: " << relc->TB[0] << " " << relc->TB[1] << "
" << relc->TB[2] << " " << relc->eaB[0] << " " << relc->eaB[1] << " " <<
relc->eaB[2] << " " << endl;


                vector<double> diff_TB_eaB(6); /// 6 zero-initialized elements
                if( relc->getVolumeTest() ) //do not check angle now, maybe when
you change the direction of sampling using jacobian, the new realization will
have good angle
                {
                        Vector3d eaB = relc->eaB;
                        Vector3d TB = relc->TB;

                        Vector3d diff_TB = TB - old_TB;
                        Vector3d diff_eaB = eaB - old_eaB;
                        //if(debug && AlgorithmJacobianRec::aDebug)  cout << "J
" << eaB[0] << " " << old_eaB[0] << " " << diff_eaB[0] << " " ; for(int h=0;
h<3; h++) if( abs(diff_eaB[h]) > 180 )  diff_eaB[h] = -1 *
Utils::sign(diff_eaB[h]) * (360 - abs(diff_eaB[h]) );//to keep it between (-180,
180)
                        //for(int h=0; h<3; h++)   diff_eaB[h] =  ( (int)
diff_eaB[h]+180)%360 - 180 ;//to keep it between (-180, 180)
                        //if(debug && AlgorithmJacobianRec::aDebug)  cout <<
diff_eaB[0] << endl;

                        for(int h=0; h<3; h++) diff_TB_eaB[h] = diff_TB[h];
                        for(int h=0; h<3; h++) diff_TB_eaB[h+3] = diff_eaB[h];

                        delete relc; relc = NULL;
                }
                else{
                        //isCayleyBoundary[i/5] = true;
                        if(debug && AlgorithmJacobianRec::aDebug)  cout << i/5
<< " J OOOPS vol neg still.. then it is CayleyBoundary" << endl; delete relc;
relc = NULL; // delete in either case
                }

                diff_TB_eaB_list.push_back(diff_TB_eaB); //if volume neg, then
(0,,,0) is added to diff_TB_eaB for current params

//		pconnections[v1][v2] = real->pconnections[v1][v2];  //take it
back
//		pconnections[v2][v1] = real->pconnections[v1][v2];
        }


//double scale_back=1;
//if(!enteredNegVal && des->independent_directions.size() != 0)
//{
//	vector<double> abs_avg_diff_TB_eaB(6);
//	for(int h=0; h<6; h++)
//	{
//		for(int i=0; i < dim; i++)
//			abs_avg_diff_TB_eaB[h] += abs(diff_TB_eaB_list[i][h]);
//		abs_avg_diff_TB_eaB[h] /= dim;
//		abs_avg_diff_TB_eaB[h] /= Settings::Sampling::gridSteps[h];
//	}
//
//	double sum=0;
//	for(int i=0; i < dim; i++)
//	{
//		int d = des->independent_directions[i];
//		sum += abs_avg_diff_TB_eaB[d]; abs_avg_diff_TB_eaB[d];
//	}
//	scale_back =  (sum/dim)  /(1./dim); //pow(sum,1./dim);
//	if(debug && AlgorithmJacobianRec::aDebug) cout << "scale_back " <<
scale_back << endl;
//}

        MatrixXd JacInv;
        double cond;


bool toosmall = false;
for(int y=0;y<6;y++) if( abs(cayleySteps[y]) < 0.1 ) toosmall = true;

if(enteredNegVal && !toosmall  ) //i put toosmall condition to prevent loop in
case neg volume
{

        if(debug && AlgorithmJacobianRec::aDebug)  cout << "enteredNegVal during
compute Jacobian" << endl;
//	double oldStep = Settings::Sampling::stepSize;
//	if(oldStep > 0)
        if(cayleySteps[0] > 0)
                for(int y=0;y<6;y++) cayleySteps[y] = -cayleySteps[y];
//		Settings::Sampling::stepSize = -Settings::Sampling::stepSize;
//first try neg stepping, i.e. the other direction else for(int y=0;y<6;y++)
cayleySteps[y] = -cayleySteps[y]/2;
//		Settings::Sampling::stepSize = -Settings::Sampling::stepSize /
2; //if it doesnot help, try halfing stepsize

//	Settings::Sampling::stepSize = Settings::Sampling::stepSize / 2;


        JacInv =computeJacobian3(des, real,  flip, succeed, cayleySteps);
//	Settings::Sampling::stepSize = oldStep;

}
//else if(scale_back > 2 && cayleySteps[0] > 0.01  && cayleySteps[0] < 10)
//{
////	double oldStep = Settings::Sampling::stepSize;
////	Settings::Sampling::stepSize = Settings::Sampling::stepSize /
scale_back;
//	for(int y=0;y<6;y++) cayleySteps[y] = Settings::Sampling::stepSize /
scale_back;
//	JacInv =computeJacobian(des, real,  flip, succeed, cayleySteps);
////	Settings::Sampling::stepSize = oldStep;
//}
else if(!enteredNegVal)  //it may enter here when enteredNegVal (if recursive
computeJacobian keep entering NegVal), when this is the case
{

        if(debug && AlgorithmJacobianRec::aDebug)  cout<< "J cayley steps: ";
        for(int j=0; j<6; j++) if(debug && AlgorithmJacobianRec::aDebug)  cout
<< cayleySteps[j] << " "; if(debug && AlgorithmJacobianRec::aDebug)  cout<<endl;


//	//to fill out rest of the dimensions
//	//diff_TB_eaB.assign(6,0); //6 numbers with value 0
//	vector<double> diff_TB_eaB(6);
//	for(size_t j = dim; j<6; j++){
//		diff_TB_eaB_list.push_back(diff_TB_eaB);
//		//cayleySteps[j] = 0;
//	}



        if(des->independent_directions.size()==0)
        {
                MatrixXd A(6,dim);
                for(int i=0; i<6; i++)
                        for(int j=0; j<dim; j++)
                                A(i,j) = diff_TB_eaB_list[j][i]/cayleySteps[j];

                vector<int> indep = computeIndependentColumns(A, dim) ;
                des->independent_directions.clear();
                des->independent_directions = indep;

                if(dim == 5)
                {
                        des->independent_directions.clear();
                        des->independent_directions.push_back(0);
                        des->independent_directions.push_back(1);
                        des->independent_directions.push_back(2);
                        des->independent_directions.push_back(3);
                        des->independent_directions.push_back(4);
                }
        }



        MatrixXd Jac(dim,dim); //remove dependent columns by setting dim x dim

        for(int i=0; i < dim; i++)
        {
                for(int j=0; j < dim; j++){
                        int d = des->independent_directions[i];
                        Jac(i,j) = diff_TB_eaB_list[j][d];// / cayleySteps[j];
                }
        }

        if(debug && AlgorithmJacobianRec::aDebug)  cout <<  "Jac \n" << Jac <<
endl;

//	double scale[dim];
//	for(int j=0; j<dim; j++){
//		double row = 0;
//		for(int k=0; k<dim; k++)
//		{
//			row += Jac(j,k)*Jac(j,k);
//		}
//		row = sqrt(row);
//		int jI = des->independent_directions.at(j);
//		scale[j] = row/Settings::Sampling::gridSteps[jI];
//	}
//	for(int j=0; j<dim; j++) for(int k=0; k<dim; k++)  if(scale[j]!=0)
Jac(j,k) = Jac(j,k) / scale[j]; else Jac(j,k) = 0;


//		if(debug && AlgorithmJacobianRec::aDebug)  cout <<  "Jac \n" <<
Jac << endl;
//		double determ = Jac.determinant();
//		if(debug && AlgorithmJacobianRec::aDebug)  cout << "Determinant
" <<  determ ;


//		VectorXd sval =  Jac.base().svd().singularValues();
//		cond = sval(0) / sval( sval.size()-1 );
//		if(true)  cout << " cond: " << cond  << endl;
                JacInv = Jac.inverse();
//		if(debug && AlgorithmJacobianRec::aDebug)  cout <<  "JacIn \n"
<< JacInv << endl;
//		for(int j=0; j<dim; j++) for(int k=0; k<dim; k++)
if(scale[j]!=0)  JacInv(k,j) = JacInv(k,j) / scale[j]; else JacInv(k,j) =  0;
                if(debug && AlgorithmJacobianRec::aDebug)  cout <<  "JacInv \n"
<< JacInv << endl;

//		MatrixXd res(dim,dim);
//		res = Jac*JacInv;
//		if(debug && AlgorithmJacobianRec::aDebug)  cout << "Jac*JacInv
\n" << res << endl;

                bool entered=false;
                for(int i=0; i < dim; i++)
                {
                        //row_i of JacInv is multipliers for cayley_step_i in
different steps double max_i_above1 = 0, max_i_under1=0; for(int j=0; j < dim;
j++)
                        {
                                if( abs( JacInv(i,j)) > 1 &&   abs( JacInv(i,j))
> max_i_above1 ) max_i_above1 = abs( JacInv(i,j));

                                if( abs( JacInv(i,j)) < 1 &&   1./abs(
JacInv(i,j)) > max_i_under1 ) max_i_under1 = 1./abs( JacInv(i,j));
                        }
                        double scale = sqrt( max_i_under1/max_i_above1);
                        if( (scale < 0.5 || scale>2) && abs(cayleySteps[i]) >=
0.1  ){ entered = true; cayleySteps[i] = cayleySteps[i] / scale;
                        }
                }
                if(entered)
                        JacInv =computeJacobian3(des, real,  flip, succeed,
cayleySteps); else for(int j=0; j<dim; j++) for(int k=0; k<dim; k++)
if(cayleySteps[j]!=0)  JacInv(j,k) = JacInv(j,k) * cayleySteps[j];

                if(debug && AlgorithmJacobianRec::aDebug)  cout <<  "final
JacInv \n" << JacInv << endl;

//		JacobiSVD<MatrixXf> svd(Jac, ComputeThinU | ComputeThinV);
//		if(debug && AlgorithmJacobianRec::aDebug)  cout << "Its singular
values are:" << endl << svd.singularValues() << endl;
//		if(debug && AlgorithmJacobianRec::aDebug)  cout << "cond: " <<
svd.singularValues()(0) / svd.singularValues()(svd.singularValues().size()-1) <<
endl;
//Jac.MatrixBase()
                //Jac.jacobiSvd().singularValues();
                //MatrixBase a(Jac);
                //if(debug && AlgorithmJacobianRec::aDebug)  cout <<
"Jac.eigenvalues() " << Jac.eigenvalues() << endl;


}
else{
        succeed = false;

        MatrixXd JacInv(dim,dim);
        for(int i=0; i < dim; i++)
                for(int j=0; j < dim; j++)
                        JacInv(i,j) = 0;
        return JacInv;
//	MatrixXd JacInv = computeJacobianfromCartesian( des,  real,  flip);
//	return JacInv;
}


//		if(cond > 3 && Settings::Sampling::stepSize > 0.1  &&
Settings::Sampling::stepSize < 10)
//		{
//			double oldStep = Settings::Sampling::stepSize;
//			Settings::Sampling::stepSize =
Settings::Sampling::stepSize / 2; //determ;
//			JacInv =computeJacobian(des, real,  flip);
//			Settings::Sampling::stepSize = oldStep;
//		}

        return JacInv;

//------------------------------------------------

//	bool isBoundary[6]; for(int b=0; b<6; b++) isBoundary[b]= false;
//	if(isCayleyBoundary[dir])  //in case neg volume, no walking on that
direction (actually it is better to choose grid axis with minimum gradient that
corresponds the cayley parameter with negative volume) but this should work as
well
//	{
//		isBoundary[dir] = true;
//		continue;
//	}

//---------------------------------------------

//	for(int i = 0; i<12; i++)	for(int j = 0; j<12; j++)
pconnections[i][j] = des->paramconnections[i][j];
//
//	//walk along direction 'dir'
//	for(size_t p = 0; p<des->paramsOrdered.size(); p=p+5)
//	{
//		int v1 = atoi( des->paramsOrdered.substr(p,2).c_str() ) ;
//		int v2 = atoi( des->paramsOrdered.substr(p+2,2).c_str() ) ;
//
//		if(debug && AlgorithmJacobianRec::aDebug)  cout << v1 << "-" <<
v2 << " " << pconnections[v1][v2] ;
//		//pconnections[v1][v2] += des->direction[g] *
Settings::Sampling::stepSize * (diff_TB_eaB_list[p/5][g] /total_TB_eaB[g]) *
(gridSteps[g] /total_TB_eaB[g]);  //(gridSteps[g] /total_TB_eaB[g]) is for
having same step size with grid steps
//
//
//		int fac = 1;
//		if(dir > 2) fac =10;
//		pconnections[v1][v2] += des->direction[dir] * JacInv(p/5,dir)
* fac;
//
//		int isContact = des->getCurContacts().find(
des->paramsOrdered.substr(p,4) , 0);
//
//		//free cayley params
//		if( isContact < 0 && pconnections[v1][v2] < 0.8 *
des->mygetRadiusSum(v1, v2) ){ //tetrahedral boundary is handled by volume test,
but collison should be handled seperately
//			pconnections[v1][v2] = 0.8 * des->mygetRadiusSum(v1,
v2);
////			isBoundary[dir] = true;
//			if(debug && AlgorithmJacobianRec::aDebug)  cout << "f";
//		}
//
//		//contact for 6d
//		if( isContact>=0 && pconnections[v1][v2] < 0.8 *
des->mygetRadiusSum(v1, v2) ){ //tetrahedral boundary is handled by volume test,
but collison should be handled seperately
//			pconnections[v1][v2] = 0.8 * des->mygetRadiusSum(v1,
v2);
////			isBoundary[dir] = true;
//			if(debug && AlgorithmJacobianRec::aDebug)  cout << "c";
//		}
//		else if( isContact>=0 && pconnections[v1][v2] >
des->mygetRadiusSum(v1, v2) + Settings::Collision::contactUpperThreshold )
//		{
//			pconnections[v1][v2] =  des->mygetRadiusSum(v1, v2) +
Settings::Collision::contactUpperThreshold;
////			isBoundary[dir] = true;
//			if(debug && AlgorithmJacobianRec::aDebug)  cout << "b";
//		}
//
//		pconnections[v2][v1] = pconnections[v1][v2];
//
//		if(debug && AlgorithmJacobianRec::aDebug)  cout << " " <<
pconnections[v1][v2] << ", ";
//	}
//
//
//	CartesianRealizer *relg = new CartesianRealizer(des, flip,
pconnections);
//	relg->axis = dir;

}
*/

/*
MatrixXd AlgorithmJacobianRec::computeJacobian2(ConvexChart* des,
CartesianRealizer * real, int flip, bool & succeed)
{
        bool debug = false;
        debug = true;

        succeed = true;
        bool enteredNegVal = false;
        if(debug && AlgorithmJacobianRec::aDebug)  cout <<
"________computeJacobian " << endl; int dim = des->paramsOrdered.size()/5;
        double cayleySteps[6] =
{Settings::Sampling::stepSize,Settings::Sampling::stepSize,Settings::Sampling::stepSize,
Settings::Sampling::stepSize, Settings::Sampling::stepSize,
Settings::Sampling::stepSize};

        Vector3d old_eaB = real->eaB;
        Vector3d old_TB = real->TB;

        vector<vector<double> > diff_TB_eaB_list;

        for(size_t i = 0; i<des->paramsOrdered.size(); i=i+5){
                int v1 = atoi( des->paramsOrdered.substr(i,2).c_str() ) ;
                int v2 = atoi( des->paramsOrdered.substr(i+2,2).c_str() ) ;
                if(debug && AlgorithmJacobianRec::aDebug)  cout << i/5 << ": "
<< v1 << "-" << v2 << " " << des->paramconnections[v1][v2] << ",  "; } if(debug
&& AlgorithmJacobianRec::aDebug)  cout << endl;

        double pconnections[12][12];
        //for(int i = 0; i<12; i++)	for(int j = 0; j<12; j++)
pconnections[i][j] = real->pconnections[i][j]; //des->paramconnections[i][j];
        //for(int i = 0; i<12; i++)	for(int j = 0; j<12; j++)
des->paramconnections[i][j] = real->pconnections[i][j]; for(size_t i = 0;
i<des->paramsOrdered.size(); i=i+5){ int v1 = atoi(
des->paramsOrdered.substr(i,2).c_str() ) ; int v2 = atoi(
des->paramsOrdered.substr(i+2,2).c_str() ) ; des->paramconnections[v1][v2] =
real->pconnections[v1][v2]; //copy only parameters from real, not the rest of
the connections because they should be updated and changed later according to
updated params. if you copy the rest, they wont be updated in realization class
since there already exist a value for them des->paramconnections[v2][v1] =
des->paramconnections[v1][v2];
        }
        des->updateBoundaries();


//		for(int i = 0; i<12; i++)	for(int j = 0; j<12; j++)
if(debug && AlgorithmJacobianRec::aDebug)  cout << i << "-" << j << ": " <<
real->pconnections[i][j] << " ";
//		if(debug && AlgorithmJacobianRec::aDebug)  cout<< endl;
//		for(int i = 0; i<12; i++)	for(int j = 0; j<12; j++)
if(debug && AlgorithmJacobianRec::aDebug)  cout << i << "-" << j << ": " <<
des->paramconnections[i][j] << " ";
//		if(debug && AlgorithmJacobianRec::aDebug)  cout<< endl;

                                for(size_t i = 0; i<des->paramsOrdered.size();
i=i+5){ int v1 = atoi( des->paramsOrdered.substr(i,2).c_str() ) ; int v2 = atoi(
des->paramsOrdered.substr(i+2,2).c_str() ) ; if(debug &&
AlgorithmJacobianRec::aDebug)  cout << i/5 << ": " << v1 << "-" << v2 << " " <<
des->paramconnections[v1][v2] << ",  "; } //if(debug &&
AlgorithmJacobianRec::aDebug)  cout << endl; if(debug &&
AlgorithmJacobianRec::aDebug)  cout << "TB_EAB " << old_TB[0] << " " <<
old_TB[1] << " " << old_TB[2] << " " << old_eaB[0] << " " << old_eaB[1] << " "
<< old_eaB[2] << " " << endl;

//	bool isCayleyBoundary[dim]; for(int b=0; b<dim; b++)
isCayleyBoundary[b]= false; for(size_t i = 0; i<des->paramsOrdered.size();
i=i+5){

                for(int k = 0; k<12; k++)	for(int j = 0; j<12; j++)
pconnections[k][j] = des->paramconnections[k][j]; //des->paramconnections[i][j];

                int v1 = atoi( des->paramsOrdered.substr(i,2).c_str() ) ;
                int v2 = atoi( des->paramsOrdered.substr(i+2,2).c_str() ) ;
                pconnections[v1][v2] += Settings::Sampling::stepSize;
                pconnections[v2][v1] = pconnections[v1][v2];


                CartesianRealizer *relc = new CartesianRealizer(des, flip,
pconnections); if( !relc->getVolumeTest() ) //volume neg
                {
                        enteredNegVal = true;
                        delete relc; relc = NULL;
                        break;

//			pconnections[v1][v2] -= 2*Settings::Sampling::stepSize;
//try reverse direction
//			pconnections[v2][v1] = pconnections[v1][v2];
//			delete relc; relc = NULL;
//			relc = new CartesianRealizer(des, flip, pconnections);
//			if( relc->getVolumeTest() )
//				cayleySteps[i/5] =
-Settings::Sampling::stepSize; // SHOULD NOT IT BE   -
Settings::Sampling::stepSize ???
//			else //volume neg again
//			{
//				if(debug && AlgorithmJacobianRec::aDebug) cout
<< "J volume negative when increased by stepsize" << endl;
//
//				//des->updateBoundaries();
//				double min = des->minconnections[v1][v2];
//				double max = des->maxconnections[v1][v2];
//				double prev = des->paramconnections[v1][v2];
//				double curr = pconnections[v1][v2];
//				if(debug && AlgorithmJacobianRec::aDebug) cout
<< "J " << i/5 << ": " << v1 << "-" << v2 << " max:" << max << " min:" << min <<
" prev:"  << prev << " curr:" << curr ;
//
//				for(int k = 0; k<12; k++)	for(int j = 0;
j<12; j++)	pconnections[k][j] = des->paramconnections[k][j];
//des->paramconnections[i][j];
//
//				if(abs(max-curr) < abs(curr-min) ){
//					pconnections[v1][v2]  = max ;
//					pconnections[v2][v1]  = max ;
//				}
//				else{
//					pconnections[v1][v2]  = min ;
//					pconnections[v2][v1]  = min ;
//				}
//
//				cayleySteps[i/5] = pconnections[v1][v2] - prev;
//				delete relc; relc = NULL;
//				relc = new CartesianRealizer(des, flip,
pconnections);
//
//				if(debug && AlgorithmJacobianRec::aDebug) cout
<< " now:" << pconnections[v1][v2] << endl;
//			}
                }

                if(debug && AlgorithmJacobianRec::aDebug)  cout << "TB_EAB after
walking on cayley " << i/5 << " is: " << relc->TB[0] << " " << relc->TB[1] << "
" << relc->TB[2] << " " << relc->eaB[0] << " " << relc->eaB[1] << " " <<
relc->eaB[2] << " " << endl;


                vector<double> diff_TB_eaB(6); /// 6 zero-initialized elements
                if( relc->getVolumeTest() ) //do not check angle now, maybe when
you change the direction of sampling using jacobian, the new realization will
have good angle
                {
                        Vector3d eaB = relc->eaB;
                        Vector3d TB = relc->TB;

                        Vector3d diff_TB = TB - old_TB;
                        Vector3d diff_eaB = eaB - old_eaB;
                        //if(debug && AlgorithmJacobianRec::aDebug)  cout << "J
" << eaB[0] << " " << old_eaB[0] << " " << diff_eaB[0] << " " ; for(int h=0;
h<3; h++) if( abs(diff_eaB[h]) > 180 )  diff_eaB[h] = -1 *
Utils::sign(diff_eaB[h]) * (360 - abs(diff_eaB[h]) );//to keep it between (-180,
180)
                        //for(int h=0; h<3; h++)   diff_eaB[h] =  ( (int)
diff_eaB[h]+180)%360 - 180 ;//to keep it between (-180, 180)
                        //if(debug && AlgorithmJacobianRec::aDebug)  cout <<
diff_eaB[0] << endl;

                        for(int h=0; h<3; h++) diff_TB_eaB[h] = diff_TB[h];
                        for(int h=0; h<3; h++) diff_TB_eaB[h+3] = diff_eaB[h];

                        delete relc; relc = NULL;
                }
                else{
                        //isCayleyBoundary[i/5] = true;
                        if(debug && AlgorithmJacobianRec::aDebug)  cout << i/5
<< " J OOOPS vol neg still.. then it is CayleyBoundary" << endl; delete relc;
relc = NULL; // delete in either case
                }

                diff_TB_eaB_list.push_back(diff_TB_eaB); //if volume neg, then
(0,,,0) is added to diff_TB_eaB for current params

//		pconnections[v1][v2] = real->pconnections[v1][v2];  //take it
back
//		pconnections[v2][v1] = real->pconnections[v1][v2];
        }


double scale_back=1;
if(!enteredNegVal && des->independent_directions.size() != 0)
{
        vector<double> abs_avg_diff_TB_eaB(6);
        for(int h=0; h<6; h++)
        {
                for(int i=0; i < dim; i++)
                        abs_avg_diff_TB_eaB[h] += abs(diff_TB_eaB_list[i][h]);
                abs_avg_diff_TB_eaB[h] /= dim;
                abs_avg_diff_TB_eaB[h] /= Settings::Sampling::gridSteps[h];
        }

        double sum=0;
        for(int i=0; i < dim; i++)
        {
                int d = des->independent_directions[i];
                sum += abs_avg_diff_TB_eaB[d]; abs_avg_diff_TB_eaB[d];
        }
        scale_back =  (sum/dim)  /(1./dim); //pow(sum,1./dim);
        if(debug && AlgorithmJacobianRec::aDebug) cout << "scale_back " <<
scale_back << endl;
}

        MatrixXd JacInv;
        double cond;



if(enteredNegVal && abs(Settings::Sampling::stepSize) >= 0.1 ) //i put stepsize
condition to prevent loop in case neg volume
{
        if(debug && AlgorithmJacobianRec::aDebug)  cout << "enteredNegVal during
compute Jacobian" << endl; double oldStep = Settings::Sampling::stepSize;
        if(oldStep > 0)
                Settings::Sampling::stepSize = -Settings::Sampling::stepSize;
//first try neg stepping, i.e. the other direction else
                Settings::Sampling::stepSize = -Settings::Sampling::stepSize /
2; //if it doesnot help, try halfing stepsize JacInv =computeJacobian(des, real,
flip, succeed); Settings::Sampling::stepSize = oldStep;

}
else if(scale_back > 2 && Settings::Sampling::stepSize > 0.01  &&
Settings::Sampling::stepSize < 10)
{
        double oldStep = Settings::Sampling::stepSize;
        Settings::Sampling::stepSize = Settings::Sampling::stepSize /
scale_back; JacInv =computeJacobian(des, real,  flip, succeed);
        Settings::Sampling::stepSize = oldStep;
}
else if(!enteredNegVal)  //it may enter here when enteredNegVal (if recursive
computeJacobian keep entering NegVal), when this is the case
{

        if(debug && AlgorithmJacobianRec::aDebug)  cout<< "J cayley steps: ";
        for(int j=0; j<6; j++) if(debug && AlgorithmJacobianRec::aDebug)  cout
<< cayleySteps[j] << " "; if(debug && AlgorithmJacobianRec::aDebug)  cout<<endl;


//	//to fill out rest of the dimensions
//	//diff_TB_eaB.assign(6,0); //6 numbers with value 0
//	vector<double> diff_TB_eaB(6);
//	for(size_t j = dim; j<6; j++){
//		diff_TB_eaB_list.push_back(diff_TB_eaB);
//		//cayleySteps[j] = 0;
//	}

        MatrixXd A(6,dim);
        for(int i=0; i<6; i++)
                for(int j=0; j<dim; j++)
                        A(i,j) = diff_TB_eaB_list[j][i]/cayleySteps[j];


        if(des->independent_directions.size()==0)
        {
                vector<int> indep = computeIndependentColumns(A, dim) ;
                des->independent_directions.clear();
                des->independent_directions = indep;

                if(dim == 5)
                {
                        des->independent_directions.clear();
                        des->independent_directions.push_back(0);
                        des->independent_directions.push_back(1);
                        des->independent_directions.push_back(2);
                        des->independent_directions.push_back(3);
                        des->independent_directions.push_back(4);
                }
        }



        MatrixXd Jac(dim,dim); //remove dependent columns by setting dim x dim
        //for(int j=0; j<dim; j++) for(int k=0; k<dim; k++)  if( cayleySteps[k]
!= 0 ) Jac(j,k) = diff_TB_eaB_list[k][j] / cayleySteps[k]; else Jac(j,k) = 0;

//	Jac(0,0) = diff_TB_eaB_list[0][0] / cayleySteps[0];
//	Jac(0,1) = diff_TB_eaB_list[1][0] / cayleySteps[1];
//
//	Jac(1,0) = diff_TB_eaB_list[0][2] / cayleySteps[0];
//	Jac(1,1) = diff_TB_eaB_list[1][2] / cayleySteps[1];
//
//	des->independent_directions.clear();
//	des->independent_directions.push_back(0);
//	des->independent_directions.push_back(2);


        for(int i=0; i < dim; i++)
        {
                for(int j=0; j < dim; j++){
                        int d = des->independent_directions[i];
                        Jac(i,j) = diff_TB_eaB_list[j][d] / cayleySteps[j];
                }
        }



//	if(debug && AlgorithmJacobianRec::aDebug)  cout <<  "Jacc \n" << Jac <<
endl;

        double scale[dim];
        for(int j=0; j<dim; j++){
                double row = 0;
                for(int k=0; k<dim; k++)
                {
                        row += Jac(j,k)*Jac(j,k);
                }
                row = sqrt(row);
                int jI = des->independent_directions.at(j);
                scale[j] = row/Settings::Sampling::gridSteps[jI];
        }
        for(int j=0; j<dim; j++) for(int k=0; k<dim; k++)  if(scale[j]!=0)
Jac(j,k) = Jac(j,k) / scale[j]; else Jac(j,k) = 0;


//		if(debug && AlgorithmJacobianRec::aDebug)  cout <<  "Jac \n" <<
Jac << endl;
//		double determ = Jac.determinant();
//		if(debug && AlgorithmJacobianRec::aDebug)  cout << "Determinant
" <<  determ ;

//		if(isnan(determ))
//			return NULL;

//		VectorXd sval =  Jac.base().svd().singularValues();
//		cond = sval(0) / sval( sval.size()-1 );
//		if(true)  cout << " cond: " << cond  << endl;
                JacInv = Jac.inverse();
//		if(debug && AlgorithmJacobianRec::aDebug)  cout <<  "JacIn \n"
<< JacInv << endl; for(int j=0; j<dim; j++) for(int k=0; k<dim; k++)
if(scale[j]!=0)  JacInv(k,j) = JacInv(k,j) / scale[j]; else JacInv(k,j) =  0;
                if(debug && AlgorithmJacobianRec::aDebug)  cout <<  "JacInv \n"
<< JacInv << endl;



//		JacobiSVD<MatrixXf> svd(Jac, ComputeThinU | ComputeThinV);
//		if(debug && AlgorithmJacobianRec::aDebug)  cout << "Its singular
values are:" << endl << svd.singularValues() << endl;
//		if(debug && AlgorithmJacobianRec::aDebug)  cout << "cond: " <<
svd.singularValues()(0) / svd.singularValues()(svd.singularValues().size()-1) <<
endl;
//Jac.MatrixBase()
                //Jac.jacobiSvd().singularValues();
                //MatrixBase a(Jac);
                //if(debug && AlgorithmJacobianRec::aDebug)  cout <<
"Jac.eigenvalues() " << Jac.eigenvalues() << endl;


}
else{
        succeed = false;

        MatrixXd JacInv(dim,dim);
        for(int i=0; i < dim; i++)
                for(int j=0; j < dim; j++)
                        JacInv(i,j) = 0;
        return JacInv;
//	MatrixXd JacInv = computeJacobianfromCartesian( des,  real,  flip);
//	return JacInv;
}


//		if(cond > 3 && Settings::Sampling::stepSize > 0.1  &&
Settings::Sampling::stepSize < 10)
//		{
//			double oldStep = Settings::Sampling::stepSize;
//			Settings::Sampling::stepSize =
Settings::Sampling::stepSize / 2; //determ;
//			JacInv =computeJacobian(des, real,  flip);
//			Settings::Sampling::stepSize = oldStep;
//		}

        return JacInv;

//------------------------------------------------

//	bool isBoundary[6]; for(int b=0; b<6; b++) isBoundary[b]= false;
//	if(isCayleyBoundary[dir])  //in case neg volume, no walking on that
direction (actually it is better to choose grid axis with minimum gradient that
corresponds the cayley parameter with negative volume) but this should work as
well
//	{
//		isBoundary[dir] = true;
//		continue;
//	}

//---------------------------------------------

//	for(int i = 0; i<12; i++)	for(int j = 0; j<12; j++)
pconnections[i][j] = des->paramconnections[i][j];
//
//	//walk along direction 'dir'
//	for(size_t p = 0; p<des->paramsOrdered.size(); p=p+5)
//	{
//		int v1 = atoi( des->paramsOrdered.substr(p,2).c_str() ) ;
//		int v2 = atoi( des->paramsOrdered.substr(p+2,2).c_str() ) ;
//
//		if(debug && AlgorithmJacobianRec::aDebug)  cout << v1 << "-" <<
v2 << " " << pconnections[v1][v2] ;
//		//pconnections[v1][v2] += des->direction[g] *
Settings::Sampling::stepSize * (diff_TB_eaB_list[p/5][g] /total_TB_eaB[g]) *
(gridSteps[g] /total_TB_eaB[g]);  //(gridSteps[g] /total_TB_eaB[g]) is for
having same step size with grid steps
//
//
//		int fac = 1;
//		if(dir > 2) fac =10;
//		pconnections[v1][v2] += des->direction[dir] * JacInv(p/5,dir)
* fac;
//
//		int isContact = des->getCurContacts().find(
des->paramsOrdered.substr(p,4) , 0);
//
//		//free cayley params
//		if( isContact < 0 && pconnections[v1][v2] < 0.8 *
des->mygetRadiusSum(v1, v2) ){ //tetrahedral boundary is handled by volume test,
but collison should be handled seperately
//			pconnections[v1][v2] = 0.8 * des->mygetRadiusSum(v1,
v2);
////			isBoundary[dir] = true;
//			if(debug && AlgorithmJacobianRec::aDebug)  cout << "f";
//		}
//
//		//contact for 6d
//		if( isContact>=0 && pconnections[v1][v2] < 0.8 *
des->mygetRadiusSum(v1, v2) ){ //tetrahedral boundary is handled by volume test,
but collison should be handled seperately
//			pconnections[v1][v2] = 0.8 * des->mygetRadiusSum(v1,
v2);
////			isBoundary[dir] = true;
//			if(debug && AlgorithmJacobianRec::aDebug)  cout << "c";
//		}
//		else if( isContact>=0 && pconnections[v1][v2] >
des->mygetRadiusSum(v1, v2) + Settings::Collision::contactUpperThreshold )
//		{
//			pconnections[v1][v2] =  des->mygetRadiusSum(v1, v2) +
Settings::Collision::contactUpperThreshold;
////			isBoundary[dir] = true;
//			if(debug && AlgorithmJacobianRec::aDebug)  cout << "b";
//		}
//
//		pconnections[v2][v1] = pconnections[v1][v2];
//
//		if(debug && AlgorithmJacobianRec::aDebug)  cout << " " <<
pconnections[v1][v2] << ", ";
//	}
//
//
//	CartesianRealizer *relg = new CartesianRealizer(des, flip,
pconnections);
//	relg->axis = dir;

}
*/

// do not use adaptiveStepSize, it keeps walking from previous point. But if it
// encounters negVol point, it stops walking. but if we walk always from relmid
// with playing step size, we fix that issue. binary jump is necessary either we
// hit cayley boundary or jacobian is not that good.
//---
// Here the problem is: if relmid is cayley tetrahedral boundary point, whenever
// we try to walk, we keep getting vol neg samples then we are stuck! it doesnot
// stuck when it is sterics because it still walks on the steric boundary (non
// neg volume point). if we can find a way to walk on tetrahedral boundary, that
// it would be fine. CartesianRealizer*
// AlgorithmJacobianRec::binaryJump(ConvexChart* des, CartesianRealizer * relmid,
// int flip, int g, bool & tempBoundary, double grid_step)
//{
//	//CartesianRealizer* relj2 = jumpToOtherDirectionTwiceJacobian( des,
//relmid, flip, g, tempBoundary); //overrides tempBoundary to be false, if it
//can jump
//
//	cout << "binaryJump" << endl;
//	cout << "relmid TB_EAB " << relmid->TB[0] << " " << relmid->TB[1] << " "
//<< relmid->TB[2] << " " << relmid->eaB[0] << " " << relmid->eaB[1] << " " <<
//relmid->eaB[2] << " " << endl; 	cout << relmid->getVolumeTest() << endl;
//
//	//double gridSteps[6] = {.4,.4,.4, 4, 4, 4};
//	int gI = des->independent_directions.at(g);
//	//double grid_step = 2*des->direction[g]*gridSteps[gI];
//
//	MatrixXd JacInv = computeJacobian( des, relmid, flip);
//	bool hit;
//	CartesianRealizer* relj2 = walkByJacobian( des, relmid->pconnections, g,
//flip, JacInv, grid_step, hit);
//
//	int i=0;
//	//if(tempBoundary || des->isBoundary[g] )
//	//{
//
//		Vector3d old_eaB = relmid->eaB;
//		Vector3d old_TB = relmid->TB;
//
//		bool goodRatio = false;
//		double ratio = 1;
//		double walk=grid_step;
//		bool jump = true;
//		double oldStep = grid_step;
//		double newStep = grid_step;
//		double volPosStep = 0;
//		double volNegStep = 10*grid_step; //start from a high number and
//keep it decreasing as you have tight boundaries for tethadredral boundary
//		while(jump && i<30 ) //&& !hit (keep  walking even if hits, that
//will allow you to walk along the boundary)
//		{
//			i++;
//			if(relj2 != NULL){ //if( relj2->getVolumeTest()  ){ //&&
//!hit 				Vector3d new_eaB = relj2->eaB; 				Vector3d new_TB = relj2->TB; 				double change
//= 1; 				if(gI<3) 					change =  (new_TB[gI]-old_TB[gI])  ;  //abs 				else 					change =
//(new_eaB[gI-3]-old_eaB[gI-3]) ; //abs
//
//				ratio = change / grid_step;
//				if( ratio<1.1 && ratio>0.9 ) goodRatio = true;
//
//				walk = grid_step - change;
//
//				if(abs(volPosStep) < abs(newStep) )
//					volPosStep = newStep; //the highest step
//possible while staying in tetrahedral region.
//			}
//			else{
////				int sign = 1;
////				if(walk<0) sign = -1;
////				walk = - sign * abs(walk) / 2; //multiply with
///sign(walk) /				ratio = 2;
//
//				if( abs(volNegStep) > abs(newStep) )
//					volNegStep = newStep; //the smallest step that
//violets tetrahedral region.
//
//				double nstep = (volPosStep + volNegStep )/2;
//				ratio = newStep / nstep;
//				goodRatio = false;
//			}
//
//			if( ratio!= 0 && !goodRatio )  //(ratio < 0.9 || ratio
//> 1.1)
//			{
//				delete relj2;
//				oldStep = newStep;
////				newStep = oldStep + walk;
//				newStep = oldStep / ratio;
//				relj2 = walkByJacobian( des, relmid->pconnections, g,
//flip, JacInv, newStep, hit);
//			}
//			else
//				jump = false;
//		}
//
//	//}
////	else
//
//
//		tempBoundary = hit;
//
////	if(hit)
////		des->isBoundary[g] = true;
//
//		//des->isBoundary[0] = true; //keep the original value
//	//if(i>=30 || ratio == 0 )
//	if( !goodRatio )
//	{
//		delete relj2;
//		//tempBoundary = true;
//
////		if(hit)
////			des->isBoundary[g] = true;
//
//		return NULL;
//	}
//	else
//		return relj2; //base case, able to jump
//}

CayleyPoint* AlgorithmJacobianRec::jumpByPreviousDirection(
    ActiveConstraintGraph* cgK, vector<double> diff_point, ConvexChart* des,
    CayleyPoint* cayleyPointmid, int flip, int g, bool& tempBoundary,
    double grid_step, bool& negVol, AtlasNode* rnode) {
  bool debug = true;

  double out[6];
  cayleyPointmid->getPoint(out);
  Orientation* relmid = cayleyPointmid->getOrientations().front();

  // CartesianRealizer* relj2 = jumpToOtherDirectionTwiceJacobian( des,  relmid,
  // flip, g, tempBoundary); //overrides tempBoundary to be false, if it can jump
  negVol = false;
  if (debug && AlgorithmJacobianRec::aDebug)
    cout << "jumpByPreviousDirection" << endl;
  if (debug && AlgorithmJacobianRec::aDebug)
    cout << "relmid TB_EAB " << relmid->TB[0] << " " << relmid->TB[1] << " "
         << relmid->TB[2] << " " << relmid->eaB[0] << " " << relmid->eaB[1]
         << " " << relmid->eaB[2] << " " << endl;
  //	if(debug && AlgorithmJacobianRec::aDebug) cout <<
  //relmid->getVolumeTest() << endl;

  // double gridSteps[6] = {.4,.4,.4, 4, 4, 4};
  int gI = des->independent_directions.at(g);
  // double grid_step = 2*des->direction[g]*gridSteps[gI];

  double multiplier_of_diffpoint = 1;
  bool hit;
  //	CartesianRealizer* relj2 = walkByPreviousDirection( des,
  //relmid->pconnections, g,  flip, diff_point, multiplier_of_diffpoint, hit);
  CayleyPoint* cayleyPoint = walkByPreviousDirection(
      cgK, des, out, g, flip, diff_point, multiplier_of_diffpoint, hit, rnode);
  Orientation* relj2 = NULL;

  int i = 0;
  // if(tempBoundary || des->isBoundary[g] )
  //{

  Vector3d old_eaB = relmid->eaB;
  Vector3d old_TB = relmid->TB;

  //		double min_step_passes_target = 10*grid_step;
  double max_step_approaches_target = 0;

  //		double min_ratio_passes_target = 10;
  double max_ratio_approaches_target = 0;

  //		double step_closest_target = 0;
  //		double ratio_closest_target = 0;

  bool goodRatio = false;
  double ratio = 1;
  double walk = grid_step;
  bool jump = true;
  double oldStep = grid_step;
  double newStep = grid_step;
  double volPosStep = 0;
  double volNegStep =
      10 * grid_step;  // start from a high number and keep it decreasing as you
                       // have tight boundaries for tethadredral boundary
  while (jump && i < 20)  //&& !hit (keep  walking even if hits, that will allow
                          //you to walk along the boundary)
  {
    i++;
    if (cayleyPoint != NULL)  // if( relj2->getVolumeTest()  ){ //&& !hit
    {
      relj2 = cayleyPoint->getOrientations().front();
      Vector3d new_eaB = relj2->eaB;
      Vector3d new_TB = relj2->TB;
      double change = 1;
      if (gI < 3)
        change = (new_TB[gI] - old_TB[gI]);  // abs
      else
        change = (new_eaB[gI - 3] - old_eaB[gI - 3]);  // abs

      ratio = change / grid_step;
      if (ratio < 1.1 && ratio > 0.9) goodRatio = true;

      walk = grid_step - change;

      if (abs(volPosStep) < abs(newStep))
        volPosStep = newStep;  // the highest step possible while staying in
                               // tetrahedral region.

      //				if(ratio>1 &&
      //ratio<min_ratio_passes_target)
      //				{
      //					min_ratio_passes_target = ratio;
      //					min_step_passes_target =
      //newStep;
      //				}

      if (ratio < 1 && ratio > max_ratio_approaches_target) {
        max_ratio_approaches_target = ratio;
        max_step_approaches_target = newStep;
      }

    } else {
      //				int sign = 1;
      //				if(walk<0) sign = -1;
      //				walk = - sign * abs(walk) / 2;
      ////multiply with    sign(walk) 				ratio = 2;

      if (abs(volNegStep) > abs(newStep))
        volNegStep =
            newStep;  // the smallest step that violets tetrahedral region.

      double nstep = (volPosStep + volNegStep) / 2;
      ratio = newStep / nstep;
      goodRatio = false;
      negVol = true;
    }

    if (ratio != 0 && !goodRatio)  //(ratio < 0.9 || ratio > 1.1)
    {
      if (cayleyPoint != NULL) {
        cayleyPoint->trim_PointMultiD();
        delete cayleyPoint;
        cayleyPoint = NULL;
      }

      oldStep = newStep;
      //				newStep = oldStep + walk;
      newStep = oldStep / ratio;

      //				if(newStep/max_step_approaches_target <
      //1) //newStep < max_step_approaches_target
      //				{
      //					if(
      //min_step_passes_target/volNegStep < 1 ) //min_step_passes_target <
      //volNegStep 						newStep = (max_step_approaches_target +
      //min_step_passes_target)/2; 					else 						newStep = (max_step_approaches_target +
      //volNegStep)/2;
      //				}
      //
      //				if(newStep/min_step_passes_target > 1)
      ////newStep > min_step_passes_target //division handles sign differences
      //as well 					newStep = (min_step_passes_target +
      //max_step_approaches_target)/2;

      //				if(ratio < max_ratio_approaches_target
      //&& abs(oldStep)>abs(max_step_approaches_target)) //passed the peak
      //					newStep =
      //(oldStep+max_step_approaches_target)/2;

      //				relj2 = walkByPreviousDirection( des,
      //relmid->pconnections, g,  flip, diff_point, newStep, hit);
      cayleyPoint = walkByPreviousDirection(cgK, des, out, g, flip, diff_point,
                                            newStep, hit, rnode);
    } else
      jump = false;
  }

  //}
  //	else

  tempBoundary = hit;

  //	if(hit)
  //		des->isBoundary[g] = true;

  // des->isBoundary[0] = true; //keep the original value
  // if(i>=30 || ratio == 0 )
  if (!goodRatio) {
    if (cayleyPoint != NULL) {
      cayleyPoint->trim_PointMultiD();
      delete cayleyPoint;
      cayleyPoint = NULL;
    }

    // tempBoundary = true;

    //		if(hit)
    //			des->isBoundary[g] = true;

    return NULL;
  } else
    return cayleyPoint;  // base case, able to jump
}

// if it passes the peak then it will keep decreasing than the peak. At the end
// if you do not have binary search around the peak, you may loose 1 sample
// point.
CayleyPoint* AlgorithmJacobianRec::binaryJump(
    ActiveConstraintGraph* cgK, MatrixXd JacInv, ConvexChart* des,
    CayleyPoint* cayleyPointmid, int flip, int g, bool& tempBoundary,
    double grid_step, bool& negVol, AtlasNode* rnode) {
  bool debug = false;
  debug = true;

  double out[6];
  cayleyPointmid->getPoint(out);
  Orientation* relmid = cayleyPointmid->getOrientations().front();

  // CartesianRealizer* relj2 = jumpToOtherDirectionTwiceJacobian( des,  relmid,
  // flip, g, tempBoundary); //overrides tempBoundary to be false, if it can jump
  negVol = false;
  if (debug && AlgorithmJacobianRec::aDebug) cout << "binaryJump" << endl;
  if (debug && AlgorithmJacobianRec::aDebug)
    cout << "relmid TB_EAB " << relmid->TB[0] << " " << relmid->TB[1] << " "
         << relmid->TB[2] << " " << relmid->eaB[0] << " " << relmid->eaB[1]
         << " " << relmid->eaB[2] << " " << endl;
  // if(debug && AlgorithmJacobianRec::aDebug) cout << relmid->getVolumeTest()
  // << endl;

  // double gridSteps[6] = {.4,.4,.4, 4, 4, 4};
  int gI = des->independent_directions.at(g);
  // double grid_step = 2 * des -> direction[g] * gridSteps[gI];

  // MatrixXd JacInv = computeJacobian( des, relmid, flip);
  bool hit;
  double walkratio = 1.;
  // CartesianRealizer* relj2 = walkByJacobian( des, relmid->pconnections, g,
  // flip, JacInv, walkratio, hit);
  CayleyPoint* cayleyPoint =
      walkByJacobian(cgK, des, out, g, flip, JacInv, walkratio, hit, rnode);
  Orientation* relj2 = NULL;

  int i = 0;
  // if(tempBoundary || des -> isBoundary[g])
  // {

  Vector3d old_eaB = relmid->eaB;
  Vector3d old_TB = relmid->TB;

  // double min_step_passes_target = 10 * grid_step;
  double max_step_approaches_target = 0;

  // double min_ratio_passes_target = 10;
  double max_ratio_approaches_target = 0;

  // double step_closest_target = 0;
  // double ratio_closest_target = 0;

  bool goodRatio = false;
  bool goodEnoughRatio = false;
  double ratio = 1;
  double walk = grid_step;
  bool jump = true;
  double oldStep = grid_step;
  double newStep = grid_step;
  double volPosStep = 0;
  double volNegStep =
      10 * grid_step;  // start from a high number and keep it decreasing as you
                       // have tight boundaries for tethadredral boundary
  while (jump && i < 20)  // && !hit (keep walking even if hits, that will allow
                          // you to walk along the boundary)
  {
    i++;
    if (cayleyPoint != NULL)  // if(relj2 -> getVolumeTest()) { // && !hit
    {
      relj2 = cayleyPoint->getOrientations().front();
      Vector3d new_eaB = relj2->eaB;
      Vector3d new_TB = relj2->TB;
      double change = 1;
      if (gI < 3)
        change = (new_TB[gI] - old_TB[gI]);  // abs
      else
        change = (new_eaB[gI - 3] - old_eaB[gI - 3]);  // abs

      ratio = change / grid_step;
      if (ratio < 1.1 && ratio > 0.9) goodRatio = true;
      // if(ratio > 0.3) goodRatio = true;
      if (ratio < 1.3 && ratio > 0.7) goodEnoughRatio = true;

      // walkratio = grid_step / change;

      walk = grid_step - change;

      if (abs(volPosStep) < abs(newStep))
        volPosStep = newStep;  // the highest step possible while staying in
                               // tetrahedral region.

      // if(ratio > 1 && ratio < min_ratio_passes_target)
      // {
      // min_ratio_passes_target = ratio;
      // min_step_passes_target = newStep;
      // }

      if (ratio < 1 && ratio > max_ratio_approaches_target) {
        max_ratio_approaches_target = ratio;
        max_step_approaches_target = newStep;
      }

    } else {
      // int sign = 1;
      // if(walk < 0) sign = -1;
      // walk = -sign * abs(walk) / 2; // multiply with sign(walk)
      // ratio = 2;

      if (abs(volNegStep) > abs(newStep))
        volNegStep =
            newStep;  // the smallest step that violets tetrahedral region.

      double nstep = (volPosStep + volNegStep) / 2;
      ratio = newStep / nstep;
      goodRatio = false;
      goodEnoughRatio = false;
      negVol = true;
    }

    if (ratio != 0 && !goodRatio)  // (ratio < 0.9 || ratio > 1.1)
    {
      if (cayleyPoint != NULL) {
        // cayleyPoint -> trim_PointMultiD();
        delete cayleyPoint;
        cayleyPoint = NULL;
      }

      oldStep = newStep;
      // newStep = oldStep + walk;
      newStep = oldStep / ratio;

      walkratio = newStep / grid_step;

      // if(newStep / max_step_approaches_target < 1) // newStep <
      // max_step_approaches_target
      // {
      // if(min_step_passes_target / volNegStep < 1 ) // min_step_passes_target
      // < volNegStep newStep = (max_step_approaches_target +
      // min_step_passes_target) / 2;
      // else
      // newStep = (max_step_approaches_target + volNegStep) / 2;
      // }

      // if(newStep / min_step_passes_target > 1)  // newStep >
      // min_step_passes_target // division handles sign differences as well
      // newStep = (min_step_passes_target + max_step_approaches_target) / 2;

      // if(ratio < max_ratio_approaches_target && abs(oldStep) >
      // abs(max_step_approaches_target)) // passed the peak newStep = (oldStep
      // + max_step_approaches_target) / 2;

      // relj2 = walkByJacobian(des, relmid -> pconnections, g, flip, JacInv,
      // walkratio , hit); // newStep
      cayleyPoint =
          walkByJacobian(cgK, des, out, g, flip, JacInv, walkratio, hit, rnode);
    } else
      jump = false;
  }

  // }
  // else

  tempBoundary = hit;

  // if(hit)
  // des -> isBoundary[g] = true;

  // des -> isBoundary[0] = true; // keep the original value
  // if(i >= 30 || ratio == 0 )

  // if(!goodEnoughRatio)
  if (!goodRatio) {
    if (cayleyPoint != NULL) {
      cayleyPoint->trim_PointMultiD();
      delete cayleyPoint;
      cayleyPoint = NULL;
    }

    // tempBoundary = true;

    // if(hit)
    // des -> isBoundary[g] = true;

    return NULL;
  } else
    return cayleyPoint;  // base case, able to jump
}

/*
CayleyPoint* AlgorithmJacobianRec::binaryJump_Corner(
        ActiveConstraintGraph* cgK,
        MatrixXd JacInv,
        ConvexChart* des,
        CayleyPoint* cayleyPointmid,
        int flip,
        int prev_g,
        int g,
        bool & tempBoundary,
        double previous_step,
        double curr_step,
        bool & negVol)
{
        bool debug = false;
        negVol = false;

        double out[6];
        cayleyPointmid->getPoint(out);
        Orientation * relmid = cayleyPointmid->getOrientations().front();


        if(debug && AlgorithmJacobianRec::aDebug) cout << "binaryJump" << endl;
        if(debug && AlgorithmJacobianRec::aDebug) cout << "relmid TB_EAB " <<
relmid->TB[0] << " " << relmid->TB[1] << " " << relmid->TB[2] << " " <<
relmid->eaB[0] << " " << relmid->eaB[1] << " " << relmid->eaB[2] << " " << endl;
//	if(debug && AlgorithmJacobianRec::aDebug) cout <<
relmid->getVolumeTest() << endl;

        int gI = des->independent_directions.at(g);
        int gI_prev = des->independent_directions.at(prev_g);

        double walkratio1=1, walkratio2=1;

        bool hit;
//	CartesianRealizer* relj2 = walkByJacobian2( des, relmid->pconnections,
prev_g, g,  flip, JacInv, walkratio1, walkratio2, hit); CayleyPoint* cayleyPoint
= walkByJacobian2( cgK, des, out, prev_g, g,  flip, JacInv, walkratio1,
walkratio2, hit); Orientation* relj2 = NULL;

        int i=0;
                Vector3d old_eaB = relmid->eaB;
                Vector3d old_TB = relmid->TB;

                bool goodRatio = false;
                bool goodEnoughRatio = false;

                double ratio = 1;
                double ratio_prev=1, ratio_curr=1;
                bool jump = true;

                double new_previous_step=previous_step, new_curr_step=curr_step;
                double curr_prev_step = sqrt(curr_step*curr_step +
previous_step*previous_step); double newStep = curr_prev_step;

                double volPosStep = 0;
                double volNegStep = 10*curr_prev_step; //start from a high
number and keep it decreasing as you have tight boundaries for tethadredral
boundary while(jump && i<20 ) //&& !hit (keep  walking even if hits, that will
allow you to walk along the boundary)
                {
                        i++;
                        if(cayleyPoint != NULL) //if( relj2->getVolumeTest()  ){
//&& !hit
                        {
                                relj2 = cayleyPoint->getOrientations().front();
                                Vector3d new_eaB = relj2->eaB;
                                Vector3d new_TB = relj2->TB;
                                double change_curr = 1;
                                if(gI<3)
                                        change_curr =  (new_TB[gI]-old_TB[gI])
;  //abs else change_curr =  (new_eaB[gI-3]-old_eaB[gI-3]) ; //abs ratio_curr =
change_curr / curr_step;


                                double change_prev = 1;
                                if(gI_prev<3)
                                        change_prev =
(new_TB[gI_prev]-old_TB[gI_prev])  ;  //abs else change_prev =
(new_eaB[gI_prev-3]-old_eaB[gI_prev-3]) ; //abs ratio_prev = change_prev /
previous_step;

                                double change = sqrt(change_prev*change_prev+
change_curr*change_curr); ratio = change / curr_prev_step;

                                goodRatio=false;
                                if( ratio<1.2 && ratio>0.8 )
                                        if(ratio_prev>0 && ratio_curr>0) //to
make sure you have the right direction goodRatio = true;

                                goodEnoughRatio = false;
                                if( ratio<1.3 && ratio>0.7 )
                                        if(ratio_prev>0 && ratio_curr>0) //to
make sure you have the right direction goodEnoughRatio = true;


                                if(abs(volPosStep) < abs(newStep) )
                                        volPosStep = newStep; //the highest step
possible while staying in tetrahedral region.


                        }
                        else{

                                if( abs(volNegStep) > abs(newStep) )
                                        volNegStep = newStep; //the smallest
step that violets tetrahedral region.

                                double nstep = (volPosStep + volNegStep )/2;
                                ratio = newStep / nstep;
                                ratio_curr=ratio; ratio_prev=ratio;

                                goodRatio = false;
                                goodEnoughRatio = false;
                                negVol = true;
                        }

                        if( ratio!= 0 && !goodRatio )  //(ratio < 0.9 || ratio
> 1.1)
                        {
                                if(cayleyPoint != NULL)
                                {
                                        cayleyPoint->trim_PointMultiD();
                                        delete cayleyPoint; cayleyPoint = NULL;
                                }

                                new_previous_step = new_previous_step /
ratio_prev; new_curr_step = new_curr_step / ratio_curr;

                                newStep = sqrt(new_curr_step*new_curr_step +
new_previous_step*new_previous_step);



                                walkratio1 =  new_previous_step / previous_step;
                                walkratio2 =  new_curr_step / curr_step;

//				relj2 = walkByJacobian2( des,
relmid->pconnections, prev_g, g,  flip, JacInv, walkratio1, walkratio2, hit);
                                cayleyPoint = walkByJacobian2( cgK, des, out,
prev_g, g,  flip, JacInv, walkratio1, walkratio2 , hit);
                        }
                        else
                                jump = false;
                }


                tempBoundary = hit;

//	if( !goodEnoughRatio )
        if( !goodRatio )
        {
                if(cayleyPoint != NULL)
                {
                        cayleyPoint->trim_PointMultiD();
                        delete cayleyPoint; cayleyPoint = NULL;
                }
                return NULL;
        }
        else
                return cayleyPoint; //base case, able to jump
}

*/

/*
CartesianRealizer* AlgorithmJacobianRec::binaryJump_Corner_merged(MatrixXd
JacInv, ConvexChart* des, CartesianRealizer * relmid, int flip, int prev_g, int
g,  bool & tempBoundary, double previous_step, double curr_step, bool & negVol)
{
        bool debug = false;
        negVol = false;
        if(debug && AlgorithmJacobianRec::aDebug) cout << "binaryJump" << endl;
        if(debug && AlgorithmJacobianRec::aDebug) cout << "relmid TB_EAB " <<
relmid->TB[0] << " " << relmid->TB[1] << " " << relmid->TB[2] << " " <<
relmid->eaB[0] << " " << relmid->eaB[1] << " " << relmid->eaB[2] << " " << endl;
        if(debug && AlgorithmJacobianRec::aDebug) cout <<
relmid->getVolumeTest() << endl;

        int gI = des->independent_directions.at(g);
        int gI_prev = des->independent_directions.at(prev_g);

        double walkratio= 1;
        double walkratio1=1, walkratio2=1;

        bool hit;
        CartesianRealizer* relj2 = walkByJacobian( des, relmid->pconnections,
prev_g,  flip, JacInv, walkratio, hit);

        int i=0;
                Vector3d old_eaB = relmid->eaB;
                Vector3d old_TB = relmid->TB;

                bool goodRatio = false;
                double ratio = 1;
                double ratio_prev=1, ratio_curr=1;
                bool jump = true;

                double new_previous_step=previous_step, new_curr_step=curr_step;
                double curr_prev_step = sqrt(curr_step*curr_step +
previous_step*previous_step); double newStep = curr_prev_step;

                double volPosStep = 0;
                double volNegStep = 10*curr_prev_step; //start from a high
number and keep it decreasing as you have tight boundaries for tethadredral
boundary while(jump && i<20 ) //&& !hit (keep  walking even if hits, that will
allow you to walk along the boundary)
                {
                        i++;
                        if(relj2 != NULL){ //if( relj2->getVolumeTest()  ){ //&&
!hit

                                Vector3d new_eaB = relj2->eaB;
                                Vector3d new_TB = relj2->TB;
                                double change_curr = 1;
                                if(gI<3)
                                        change_curr =  (new_TB[gI]-old_TB[gI])
;  //abs else change_curr =  (new_eaB[gI-3]-old_eaB[gI-3]) ; //abs ratio_curr =
change_curr / curr_step;


                                double change_prev = 1;
                                if(gI_prev<3)
                                        change_prev =
(new_TB[gI_prev]-old_TB[gI_prev])  ;  //abs else change_prev =
(new_eaB[gI_prev-3]-old_eaB[gI_prev-3]) ; //abs ratio_prev = change_prev /
previous_step;

                                double change = sqrt(change_prev*change_prev+
change_curr*change_curr); ratio = change / curr_prev_step;

                                goodRatio=false;
                                if( ratio<1.2 && ratio>0.8 )
                                        if(ratio_prev>0 && ratio_curr>0) //to
make sure you have the right direction goodRatio = true;



                                if(abs(volPosStep) < abs(newStep) )
                                        volPosStep = newStep; //the highest step
possible while staying in tetrahedral region.


                        }
                        else{

                                if( abs(volNegStep) > abs(newStep) )
                                        volNegStep = newStep; //the smallest
step that violets tetrahedral region.

                                double nstep = (volPosStep + volNegStep )/2;
                                ratio = newStep / nstep;
                                ratio_curr=ratio; ratio_prev=ratio;

                                goodRatio = false;
                                negVol = true;
                        }

                        if( ratio!= 0 && !goodRatio )  //(ratio < 0.9 || ratio
> 1.1)
                        {
                                delete relj2;

                                new_previous_step = new_previous_step /
ratio_prev; new_curr_step = new_curr_step / ratio_curr;

                                newStep = sqrt(new_curr_step*new_curr_step +
new_previous_step*new_previous_step);



                                walkratio1 =  new_previous_step / previous_step;
                                walkratio2 =  new_curr_step / curr_step;
                                walkratio = (walkratio1+walkratio2)/2;

                                relj2 = walkByJacobian( des,
relmid->pconnections, prev_g,  flip, JacInv, walkratio, hit);

                        }
                        else
                                jump = false;
                }


                tempBoundary = hit;

        if( !goodRatio )
        {
                delete relj2;
                return NULL;
        }
        else
                return relj2; //base case, able to jump
}
*/

/*
CartesianRealizer* AlgorithmJacobianRec::binaryJumpp(ConvexChart* des,
CartesianRealizer * relmid, int flip, int g, bool & tempBoundary, double
grid_step)
{

        cout << "binaryJumpp" << endl;
        cout << "relmid TB_EAB " << relmid->TB[0] << " " << relmid->TB[1] << " "
<< relmid->TB[2] << " " << relmid->eaB[0] << " " << relmid->eaB[1] << " " <<
relmid->eaB[2] << " " << endl; cout << relmid->getVolumeTest() << endl;

        int gI = des->independent_directions.at(g);

        CartesianRealizer* relValid = relmid;
        bool succeed;
        MatrixXd JacInv = computeJacobian3( des, relmid, flip, succeed);
//computeJacobian bool hit; CartesianRealizer* relj2 = walkByJacobian( des,
relmid->pconnections, g,  flip, JacInv, grid_step, hit);

        int i=0;

        Vector3d old_eaB = relmid->eaB;
        Vector3d old_TB = relmid->TB;

        //old_TB[gI] as reference point. so it is accepted as location 0.
        double up_min_valid_loc=0;  //valid steps while trying to achieve best
ratio double down_min_valid_loc=0;

        double up_max_invalid_loc=10*grid_step;  //neg_vol boundaries.  you
cannot walk further than these boundaries double
down_max_invalid_loc=-10*grid_step; //neg_vol


        bool goodRatio = false;
        double ratio = 1;
        double walk=grid_step;
        bool jump = true;
        double oldStep = grid_step;
        double newStep = grid_step;
        double volPosStep = 0;
        double volNegStep = 10*grid_step; //start from a high number and keep it
decreasing as you have tight boundaries for tethadredral boundary double
lastPosStep = 0; while(jump && i<30)
        {
                bool enteredNegVal = false;
                i++;
                if(relj2 != NULL){ //if( relj2->getVolumeTest() ){
                        enteredNegVal = false;
                        relValid = relj2;
                        lastPosStep = newStep;

                        Vector3d new_eaB = relj2->eaB;
                        Vector3d new_TB = relj2->TB;
                        double change = 1;
                        if(gI<3)
                                change =  (new_TB[gI]-old_TB[gI])  ;  //abs
                        else
                                change =  (new_eaB[gI-3]-old_eaB[gI-3]) ; //abs

//			change = abs(change); //the change can be negative in
case you just passed the peak of the curve

                        ratio = change / grid_step;
                        if( ratio<1.1 && ratio>0.9 ) goodRatio = true;

                        walk = grid_step - change;

//			if(volPosStep < newStep)
//				volPosStep = newStep; //the highest step
possible while staying in tetrahedral region.


                        if(newStep>0 && newStep>up_min_valid_loc)
                                up_min_valid_loc=newStep;
                        else if(newStep<0 && newStep<down_min_valid_loc)
                                down_min_valid_loc=newStep;

//			if(change>0 && change>up_min_valid_loc)
//				up_min_valid_loc=change;

                }
                else{
//				int sign = 1;
//				if(walk<0) sign = -1;
//				walk = - sign * abs(walk) / 2; //multiply with
sign(walk)
//				ratio = 2;

//			if(volNegStep > newStep)
//				volNegStep = newStep; //the smallest step that
violets tetrahedral region.


                        if(newStep>0 && newStep<up_max_invalid_loc)
                                up_max_invalid_loc=newStep;
                        else if(newStep<0 && newStep>down_max_invalid_loc)
                                down_max_invalid_loc=newStep;



                        //double nstep = (volPosStep + volNegStep )/2;

                        double nstep;
//			if( newStep>0 )
//				nstep = (up_min_valid_loc + up_max_invalid_loc
)/2;
//			else
//				nstep = (down_min_valid_loc +
down_max_invalid_loc )/2; if( newStep>0 ) nstep = (lastPosStep +
up_max_invalid_loc )/2; else nstep = (lastPosStep + down_max_invalid_loc )/2;

                        ratio = newStep / nstep;
                        goodRatio = false;

                        delete relj2;
                        enteredNegVal = true;

                        walk = nstep - newStep;
                }

                if( ratio!= 0 && !goodRatio )  //(ratio < 0.9 || ratio > 1.1)
                {
//			delete relj2;
                        oldStep = newStep;
                        newStep = oldStep + walk;
//			newStep = oldStep / ratio;
                        bool succeed;
                        MatrixXd JacInv = computeJacobian3( des, relValid, flip,
succeed); //computeJacobian
//			double walkBy =  newStep-oldStep;
//			if(enteredNegVal)
//				walkBy = newStep;
//			double walkBy =  newStep-lastPosStep;

                        double walkBy = walk;
                        if(enteredNegVal)
                                walkBy =  newStep-lastPosStep;
                        relj2 = walkByJacobian( des, relValid->pconnections, g,
flip, JacInv, walkBy, hit);
                }
                else
                        jump = false;
        }


//	if(hit)
//		des->isBoundary[g] = true;

                //des->isBoundary[0] = true; //keep the original value
        //if(i>=30 || ratio == 0)
        if( !goodRatio )
        {
                delete relj2;
                tempBoundary = true;

//		if(hit)
//			des->isBoundary[g] = true;

                return NULL;
        }
        else
                return relj2; //base case, able to jump
}
*/

// CartesianRealizer* AlgorithmJacobianRec::binaryJumpp(ConvexChart* des,
// CartesianRealizer * relmid, int flip, int g, bool & tempBoundary, double
// grid_step)
//{
//
//	cout << "binaryJumpp" << endl;
//	cout << "relmid TB_EAB " << relmid->TB[0] << " " << relmid->TB[1] << " "
//<< relmid->TB[2] << " " << relmid->eaB[0] << " " << relmid->eaB[1] << " " <<
//relmid->eaB[2] << " " << endl; 	cout << relmid->getVolumeTest() << endl;
//
//	int gI = des->independent_directions.at(g);
//
//	CartesianRealizer* relValid = relmid;
//	MatrixXd JacInv = computeJacobian( des, relmid, flip);
//	bool hit;
//	CartesianRealizer* relj2 = walkByJacobian( des, relmid->pconnections, g,
//flip, JacInv, grid_step, hit);
//
//	int i=0;
//
//	Vector3d old_eaB = relmid->eaB;
//	Vector3d old_TB = relmid->TB;
//
//	bool goodRatio = false;
//	double ratio = 1;
//	double walk=grid_step;
//	bool jump = true;
//	double oldStep = grid_step;
//	double newStep = grid_step;
//	double volPosStep = 0;
//	double volNegStep = 10*grid_step; //start from a high number and keep it
//decreasing as you have tight boundaries for tethadredral boundary 	double
//lastPosStep = 0; 	while(jump && i<30)
//	{
//		bool enteredNegVal = false;
//		i++;
//		if(relj2 != NULL){ //if( relj2->getVolumeTest() ){
//			enteredNegVal = false;
//			relValid = relj2;
//			lastPosStep = newStep;
//
//			Vector3d new_eaB = relj2->eaB;
//			Vector3d new_TB = relj2->TB;
//			double change = 1;
//			if(gI<3)
//				change =  (new_TB[gI]-old_TB[gI])  ;  //abs
//			else
//				change =  (new_eaB[gI-3]-old_eaB[gI-3]) ; //abs
//
////			change = abs(change); //the change can be negative in case
///you just passed the peak of the curve
//
//			ratio = change / grid_step;
//			if( ratio<1.1 && ratio>0.9 ) goodRatio = true;
//
//			walk = grid_step - change;
//
//			if(volPosStep < newStep)
//				volPosStep = newStep; //the highest step possible
//while staying in tetrahedral region.
//		}
//		else{
////				int sign = 1;
////				if(walk<0) sign = -1;
////				walk = - sign * abs(walk) / 2; //multiply with
///sign(walk) /				ratio = 2;
//
//			if(volNegStep > newStep)
//				volNegStep = newStep; //the smallest step that violets
//tetrahedral region.
//
//			double nstep = (volPosStep + volNegStep )/2;
//			ratio = newStep / nstep;
//			goodRatio = false;
//
//			delete relj2;
//			enteredNegVal = true;
//
//			walk = nstep - newStep;
//		}
//
//		if( ratio!= 0 && !goodRatio )  //(ratio < 0.9 || ratio > 1.1)
//		{
////			delete relj2;
//			oldStep = newStep;
//			newStep = oldStep + walk;
////			newStep = oldStep / ratio;
//			MatrixXd JacInv = computeJacobian( des, relValid, flip);
////			double walkBy =  newStep-oldStep;
////			if(enteredNegVal)
////				walkBy = newStep;
////			double walkBy =  newStep-lastPosStep;
//
//			double walkBy = walk;
//			if(enteredNegVal)
//				walkBy =  newStep-lastPosStep;
//			relj2 = walkByJacobian( des, relValid->pconnections, g,
//flip, JacInv, walkBy, hit);
//		}
//		else
//			jump = false;
//	}
//
//
////	if(hit)
////		des->isBoundary[g] = true;
//
//		//des->isBoundary[0] = true; //keep the original value
//	//if(i>=30 || ratio == 0)
//	if( !goodRatio )
//	{
//		delete relj2;
//		tempBoundary = true;
//
////		if(hit)
////			des->isBoundary[g] = true;
//
//		return NULL;
//	}
//	else
//		return relj2; //base case, able to jump
//}

// CartesianRealizer* AlgorithmJacobianRec::binaryJumpp(ConvexChart* des,
// CartesianRealizer * relmid, int flip, int g, bool & tempBoundary, double
// grid_step)
//{
//
//	cout << "binaryJumpp" << endl;
//	cout << "relmid TB_EAB " << relmid->TB[0] << " " << relmid->TB[1] << " "
//<< relmid->TB[2] << " " << relmid->eaB[0] << " " << relmid->eaB[1] << " " <<
//relmid->eaB[2] << " " << endl; 	cout << relmid->getVolumeTest() << endl;
//
//	int gI = des->independent_directions.at(g);
//
//	CartesianRealizer* relValid = relmid;
//	MatrixXd JacInv = computeJacobian( des, relmid, flip);
//	CartesianRealizer* relj2 = walkByJacobian( des, relmid->pconnections, g,
//flip, JacInv, grid_step);
//
//	int i=0;
//
//	Vector3d old_eaB = relmid->eaB;
//	Vector3d old_TB = relmid->TB;
//
//	bool goodRatio = false;
//	double ratio = 1;
//	double walk=grid_step;
//	bool jump = true;
//	double oldStep = grid_step;
//	double newStep = grid_step;
//	double volPosStep = 0;
//	double volNegStep = 10*grid_step; //start from a high number and keep it
//decreasing as you have tight boundaries for tethadredral boundary 	double
//lastPosStep = 0; 	while(jump && i<30)
//	{
//		bool enteredNegVal = false;
//		i++;
//		if( relj2->getVolumeTest() ){
//			relValid = relj2;
//			lastPosStep = newStep;
//
//			Vector3d new_eaB = relj2->eaB;
//			Vector3d new_TB = relj2->TB;
//			double change = 1;
//			if(gI<3)
//				change =  (new_TB[gI]-old_TB[gI])  ;  //abs
//			else
//				change =  (new_eaB[gI-3]-old_eaB[gI-3]) ; //abs
//
////			change = abs(change); //the change can be negative in case
///you just passed the peak of the curve
//
//			ratio = change / grid_step;
//			if( ratio<1.1 && ratio>0.9 ) goodRatio = true;
//
//			walk = grid_step - change;
//
//			if(volPosStep < newStep)
//				volPosStep = newStep; //the highest step possible
//while staying in tetrahedral region.
//		}
//		else{
////				int sign = 1;
////				if(walk<0) sign = -1;
////				walk = - sign * abs(walk) / 2; //multiply with
///sign(walk) /				ratio = 2;
//
//			if(volNegStep > newStep)
//				volNegStep = newStep; //the smallest step that violets
//tetrahedral region.
//
//			double nstep = (volPosStep + volNegStep )/2;
//			ratio = newStep / nstep;
//			goodRatio = false;
//
//			delete relj2;
//			enteredNegVal = true;
//		}
//
//		if( ratio!= 0 && !goodRatio )  //(ratio < 0.9 || ratio > 1.1)
//		{
////			delete relj2;
//			oldStep = newStep;
////				newStep = oldStep + walk;
//			newStep = oldStep / ratio;
//			MatrixXd JacInv = computeJacobian( des, relValid, flip);
////			double walkBy =  newStep-oldStep;
////			if(enteredNegVal)
////				walkBy = newStep;
//			double walkBy =  newStep-lastPosStep;
//			relj2 = walkByJacobian( des, relValid->pconnections, g,
//flip, JacInv, walkBy);
//		}
//		else
//			jump = false;
//	}
//
//
//		//des->isBoundary[0] = true; //keep the original value
//	if(i>=30 || ratio == 0)
//	{
//		delete relj2;
//		tempBoundary = true;
//		return NULL;
//	}
//	else
//		return relj2; //base case, able to jump
//}

// list<CartesianRealizer*> AlgorithmJacobianRec::findRealizations(ConvexChart*
// des, map<int, double> stepratios){
//
//	list<CartesianRealizer*> output; //WHY DO WE KEEP WHOLE REALIZATION,
//DONT WE NEED JUST ORIENTATION PART? 	bool map = des->mapl; 	int msol=0; 	double
//a[6][4]; 	string ss; 	int dim = des->getParams().size()/5;
//
//	if( !map )
//		for(int i=0; i<8 ; i++ ){ //for each flip
//			int ratio = round( stepratios[i] );
//			for(int s=1; s<ratio; s++) //we need ratio times
//intermediate steps
//			{
//				CartesianRealizer *relz = new CartesianRealizer(des,
//i, s); 				if( relz->getVolumeTest() ){ 					if( !Settings::Collision::checkAngle ||
//!relz->angleViolated() ) //small anglee
////relz->angleBetweenTwoMolecularUnits() 						output.push_back( relz ); 					else{
//						//delete relz;
//						output.push_back( relz ); //child nodes
//can cause small angles so if this realization cause a new contact, then create
//the child not and keep it as witness at the child node. After that delete this
//realization from parent(this node), and make it orange.
//
//					}
//				}
//				else{ //volume negative
//					delete relz;
//					break; //in the next flips it will always get
//volume negative too
//				}
//			}
//		}
//	return output;
//}

/*
bool Algorithm::MySample(AtlasNode *rnode, bool dense, Orientation* orrw, bool
continu, bool breath, bool bret ){ if(this->verbose) cout << "Sample Started" <<
endl;

        ActiveConstraintGraph *cgK = rnode->getCG();

        if(this->verbose) cgK->printParticipant();

        int from = rnode->getID();   //atlas->getNodeNum(cgK); // from is the
node number of the incoming contact graph as stored in the Roadmap object.
        if(this->verbose)cout<< "from  " << from << endl;

        if(dense){ //it is now being refined which means it isn't a completed
graph anymore until refinement is finished rnode->setComplete(false); //later
this information will be saved to RoadMap.txt file
        }

        bool mynewstart = continu; //Flag of whether this is a first/new sample
or a continued sample.

        ConstraintCheck *detector = new ConstraintCheck();//Create a fresh
Collision detector if(this->verbose)cout << "Created detector" << endl;


        if(Settings::Sampling::dynamicStepSizeAmong)
        {
                double testStep = 0.4; //was 0.2 before
                bool userDefined_stericConstraint =
Settings::Collision::stericConstraint; Settings::Collision::stericConstraint =
false; //allow collision on parameters cgK->setStepSize(testStep); //increment
this stepsize, it is causing slowness ConvexChart *descv = new ConvexChart( cgK,
dense );//for volume computation, to make blue points as minimum as possible int
volume = 0; if( descv->mynewstartGrid( mynewstart) && ( !descv->mapl ||
Settings::General::runSolver)  ) for( ;   !descv->doneGrid()&&
!Settings::AtlasBuilding::stop;		descv->stepGrid1()  ) volume++;

                double sumradius = descv->mygetRadiusSum(0, 6);
                delete descv;



                //since we are not doing 6d sampling and forcing contact to be
exact distance
                //but grid having range for contact hence having sampling size
proportional with sumofradius, we should increase our sampling size per node
with that ratio double range = pow(sumradius+1, 3) - pow(sumradius*.8, 3);
                Settings::Collision::stericConstraint =
userDefined_stericConstraint;
                //tns:total_number_of_steps nspd:number_of_steps_per_dimension
s:stepsize
                //nspd = pow(tns, 0.2)
                //nspd1*s1=nspd2*s2
                //s2 = s1* nspd1/nspd2 = s1* pow(tns1, 0.2)/pow(tns2, 0.2) = s1*
pow(tns1/tns2, 0.2) = s1 / pow(tns2/tns1, 0.2)

                double voltimes = range*50000. / volume; //50000:testjose2
created 6.8 million //10000 gives 1.5 million samples in total   //1000
//100.000
                //1000 with linear sampling: gives 70.000 samples
                //80000 with linear sampling: gives 376.184 samples
                //1200000 take a month
                //600000 with linear sampling: takes 2 week: gives 139448496
samples, 101130295 inner samples
                //120000 running on willy now
                //10.000 reverse linear gives:  13699 samples, 12482 inner
samples
                //100.000 reverse linear gives:  there are files with more
than 4.5 GB! not enough RAM

                double pwr = pow( voltimes, 0.2);
                double stepsizeee = testStep / pwr;
        //	double stepsizeee = 0.355637; //NODE84~~~~~~~~
                cgK->setStepSize(stepsizeee);
        //	cout << " volume " << volume  << " sumradius " << sumradius <<
"range " << range  << " voltimes " << voltimes  << " pwr " << pwr << "
stepsizeee " << stepsizeee << endl;

        }




        ConvexChart *desc = new ConvexChart( cgK, dense );//Create a fresh
description if(this->verbose)cout << "Created description" << endl;

//	if( orrw!= NULL){
//		cout<< "entrance of mysample " << from << endl;
//		detector->test(cgK, orrw);
//	}


        snl->saveNode( rnode ); //save node after parameters is set (by creating
description object).


        if(Settings::OldOrTesting::just5d && cgK->getDim() != 5 ) //!=5   //JUST
TO TEST 5D NODES
        {
                delete detector;
                delete desc;
                rnode->setFoundGoodOrientation(true);
                rnode->setComplete(true);
                return true;
        }

        if(Settings::OldOrTesting::just345d && cgK->getDim() < 3 ) //!=5 //JUST
TO TEST 345D NODES
        {
                delete detector;
                delete desc;
                rnode->setFoundGoodOrientation(true);
                rnode->setComplete(true);
                return true;
        }


        if(!continu)
                this->snl->appendDimension(rnode );



        if(  (!continu && !dense) || cgK->getSpace().empty()){// a fresh start
                //this->snl->saveNode( this->atlas->getNode(from) );   //no need
to save anymore, already saved to roadmap by addnode function, and no need to
save contact info in the nodefile anymore

        }else if(!breath && !dense)//		 if a continue on a non-empty
graph
        {
                vector<int> con = rnode->getConnection();// con = the node
numbers connected to this one size_t dim = rnode->getDim();
                for(vector<int>::iterator it = con.begin(); it!=con.end(); it++)
                {
                        //ActiveConstraintGraph *cgKk =
atlas->getCgkofNode(*it); //cgKk is an actual pointer to the graph of a
connected node. AtlasNode *rnode_cgKK =  this->atlas->getNode(*it); if(
rnode_cgKK->getDim() < dim && !rnode_cgKK->isComplete()  ) // if it is a lesser
dimension and incomplete
                        {
                                //ActiveConstraintGraph *cgKk =
snl->loadNode(*it);  //roadmapdaki noda yuklese, sonra o nodun cgk pointerini
alsak
                                // cgKk is replaced with a copy generated from
the save file.
                                //rnode_cgKK->setCG( cgKk );
                                snl->loadNode( *it, rnode_cgKK->getCG() );
                                this->MySample(rnode_cgKK , false, NULL,
continu, false, false );
                        }
                }
        }

        bool noGoodOrientation=true;
        if( orrw!= NULL){
                PointMultiD* wp4 = desc->getWitnessPoint(cgK,orrw); //DO NOT
FORGET TO DELETE ORRW  //DO NOT DELETE ORRW, IT IS USED LATER AT THE PLACE IT IS
CREATED. DELETE A POINTER WHERE IT IS CREATED AFTER ITS JOB IS DONE.
//		double pt[6];
//		wp4->getPoint(pt);
//		cout << "[desc] ";
//		for(int i = 0; i < 6; i++){
//			cout << pt[i] << " ";
//		}
//		cout << endl;

                desc->setWitnessPoint(wp4);
                if( !Settings::Collision::checkAngle || !orrw->angleViolated()){
                        cgK->AddWitnessPoint(wp4); //  add it to the graph //
uncommented 2june12 it will be needed as a starting point in the description
class  //commented on 27may12, no need to keep it since it is saved to the file
in the next line
//			Start_t = clock();
                        this->snl->appendWitness(rnode, wp4); // save it too the
graph  //then it will be saved to file twice in next steps again?? no only space
points are saved at the next steps, there is no more witness and witness saving
job later.
//			End_t = clock();
//			time_appendSpacePoints += difftime(End_t, Start_t);
                        noGoodOrientation = false;
                }
        }

        if(Settings::General::reverseWitness){
                vector<CgMarker::Edge> edges = this->cgmarker.getEdge(cgK);
                cout << "[mark]" << *cgK << endl;
                cout << "[mark] getting witness " << edges.size() << endl;
                for(int i = 0; i < edges.size(); i++){
                        Orientation* ornt = edges[i].second; //DO NOT FORGET TO
DELETE ORNT PointMultiD* wp4 = desc->getWitnessPoint(cgK,ornt);
                        desc->setWitnessPoint(wp4);
                        if( !Settings::Collision::checkAngle ||
!orrw->angleViolated()){ cgK->AddWitnessPoint(wp4); //  add it to the graph, it
will be needed as a starting point in the description class
                                this->snl->appendWitness(rnode, wp4);// save it
too the graph noGoodOrientation = false;
                        }
                }



//		list<CartesianRealizer*> real = findRealizations(desc );




                if(cgK->getDim() == 0){
                        this->cgmarker.mark(cgK);
                        cout << "[mark] mark 0" << endl;

                }
        }


        //delete orrw;  //do not delete

        if(bret){
                if(this->verbose)cout << "Returning\n" << endl;
                //return true;
        }

//	if(from!=0)
//		return true;







//	if(Settings::OldOrTesting::just5d && cgK->getDim() != 5 ) //!=5 //JUST
TO TEST 5D NODES
//	{
//		delete detector;
//		delete desc;
//		rnode->setFoundGoodOrientation(true);
//		rnode->setComplete(true);
//		return true;
//	}
//
//	if(Settings::OldOrTesting::just345d && cgK->getDim() < 3 ) //!=5 //JUST
TO TEST 5D NODES
//	{
//		delete detector;
//		delete desc;
//		rnode->setFoundGoodOrientation(true);
//		rnode->setComplete(true);
//		return true;
//	}






        bret = breath;
        int noPoints = 0;
        int loop=0;


    if( desc->mynewstartGrid( mynewstart) && ( !desc->mapl ||
Settings::General::runSolver)  )
    {
//    	{// appending volume data to file
//    		this->snl->appendVolume(rnode );
//    	}


//    	list<CartesianRealizer*> oldreal;
//    	list<CartesianRealizer*> inbetween;
//    	list<CartesianRealizer*> newreal;
        list<CartesianRealizer*> real;
                for( ;   !desc->doneGrid()&& !Settings::AtlasBuilding::stop;
desc->stepGrid1()  )
                {


//--------------------------------------- numerical jacobian sampling
//			oldreal = newreal;
//			//if( !newreal.empty() ) oldreal.insert(oldreal.begin(),
newreal.begin(), newreal.end() );
//			newreal = findRealizations(desc );
//
//			if( !oldreal.empty() && !newreal.empty() )
//				inbetween = findRealizations( desc, oldreal,
newreal);
//
//			real.clear();
//			if( !newreal.empty() )    real.insert(real.begin(),
newreal.begin(), newreal.end() );
//			if( !inbetween.empty() )    real.insert(real.begin(),
inbetween.begin(), inbetween.end() );
//--------------------------------------------

                        real = findRealizations(desc );

//------------------------------------------------

                        //real = findRealizationsJacobian(desc, real, 1);

                        vector<double> out = desc->getPoint(); //should change
the point for each real PointMultiD* p4d = new PointMultiD(out); p4d->tetra(!
real.empty());

            int rr=0;
            if( !real.empty() )
            for(list<CartesianRealizer*>::iterator r = real.begin(); r !=
real.end()&& !Settings::AtlasBuilding::stop; r++)
                        {
                //cout << rr++ << endl;
//            	loop++;
//            	{
//            		Orientation* ort = (*r)->getOrienation();
//            		double fa[3][3],fb[3][3],ta[3][3],tb[3][3];
//            		ort->getFromTo(fa,fb,ta,tb);
// Stree::printMatrix(Stree::getSingleTransformFromTwo(fa,ta,fb,tb));
//            		cout <<
is_collide(cgK->helA->getStree(),cgK->helB->getStree(),fa,ta,fb,tb) << endl;
//
//            	}

//            	Start_t = clock();
                bool myvalid = valid(*r,detector);
//            	End_t = clock();
//				time_myvalid += difftime(End_t, Start_t);

                                if(myvalid)//checks against rest of the helix
                                {
                                        Orientation* ori_on_lattice =
(*r)->getOrienation();//19april2013

                                        //ActiveConstraintGraph *cgKprime = new
ActiveConstraintGraph( (*r)->getCG() );

                                        //list<pair<int,int> > contactList =
detector->getContacts(); int noEquality = (*r)->getVolumeNo(); noEquality = 0;

                                        if( (*r)->angleViolated() )  //1
realization is enough to set it bad angle.??? Because during saving it will
write each orientation with good angle, even the point shows yellow, bad angle.
                                                p4d->setBadAngle(true);
                                        else
                                                noGoodOrientation = false;

                                        bool addedOr = false;

                                        if(
Settings::RootNodeCreation::createChildren ) //JUST TO CREATE INTERIOR
                                        ////if( !contactList.empty() )
                                        {
                                                vector<double> pp=
desc->getPoint(); int flip = (*r)->getFlip(); CartesianRealizer* sr; bool
entered = false; while( !desc->stepGridContact()  &&
!Settings::AtlasBuilding::stop)//look up and down for one step in each direction
                                                {
                                                        sr = new
CartesianRealizer(desc, flip);//if needs maple, then it will always return first
root not the one with specified flip ! if( sr->getVolumeTest() &&
!valid(sr,detector)  ){ //collision  //ACTUALLY CHECK SHOULD BE DONE AROUND BLUE
REGION AS WELL, IF VOLUME NEGATIVE, SEARCH FOR BOUNDARY IN BETWEEN!!!
                                                                desc->stepGridBinary(false);
                                                                delete sr;
                                                                sr = new
CartesianRealizer(desc, flip );

//								cout << "before
while " << endl;
// detector->test(sr);

                                                                bool binvalid =
valid(sr,detector); //if it is not valid, i mean if there is a collision then it
should continue binary stepping

                                                                list<pair<int,int>
> contactList2 = detector->getFineContacts(flip); int binno=0;

                                                                //binary search
                                                                while(
(contactList2.empty() || !binvalid ) && binno <30){ // while(contactList2.size()
== 0 && binno <30)    //find new contact for new small threshold t binno++;
                                                                        //cout
<< "binno " << binno << endl; desc->stepGridBinary( binvalid ); delete sr; sr =
new CartesianRealizer(desc, flip  );
//
detector->test(sr);

                                                                        binvalid
= valid(sr,detector); contactList2 = detector->getFineContacts(flip); //get
contacts according to previous stepGrid
                                                                }//end of binary
search

                                                                //cout << "end
while "  << endl; if( valid(sr,detector) && contactList2.size() != 0  &&
(Settings::AtlasBuilding::ifBadAngleWitness_createChild || !sr->angleViolated())
) ///
                                                                {
//									cout <<
"inside if " << endl;
//
detector->test(sr);

                                                                        Orientation*
orie_on_boundary = sr->getOrienation();  //added on 22march 2013

                                                                        bool
firstPath = true; ActiveConstraintGraph *cgKprime[contactList2.size()]; for(int
c=0; c<contactList2.size(); c++)  //TRY ALL CONTACT LIST
                                                                        {
                                                                                cgKprime[c] = new ActiveConstraintGraph( (*r)->getCG() );
                                                                                //cgKprime[c]->setStepSize(stt->stepsize);
                                                                                //cgKprime[c]->setStepSize(); //commented on 30jan2014 ???
                                                                                entered = true;
                                                                                pair<int,int> contact = contactList2.front();
                                                                                cgKprime[c]->addContact(contact);
                                                                                contactList2.pop_front();

                                                                                if( cgKprime[c]->getK() > 6  || cgKprime[c]->isDependent() || ( cgKprime[c]->getK() == 6 &&  !cgKprime[c]->hasThirdAtom() ) )
                                                                                {
                                                                                        delete cgKprime[c]; cgKprime[c] = NULL;
                                                                                        contactList2.push_back( contact );
                                                                                        continue;
                                                                                }


// Start_t = clock();

                                                                                int nodeID;  // = this->atlas->getNodeNum(cgKprime[c]);
                                                                                int success = this->atlas->addNode(cgKprime[c], nodeID);
                                                                                AtlasNode* rnode_ID = atlas->getNode(nodeID);

// End_t = clock();
//
time_mapView += difftime(End_t, Start_t);


                                                                                bool contin = false;
                                                                                if(success == 0){   //(nodeID != -1){ //started maybe done
                                                                                        delete cgKprime[c]; cgKprime[c] = NULL;
                                                                                        cgKprime[c] = rnode_ID->getCG();  //atlas->getCgkofNode(nodeID);
                                                                                        if(!rnode_ID->isComplete() && !breath){ //not done
                                                                                        //	snl->loadNode(nodeID, cgKprime[c] );    // no need to load points
                                                                                                contin = true;
                                                                                        }
                                                                                }

                                                                                Orientation* wit_orr_toSend;
                                                                                wit_orr_toSend = sr->getOrienation() ;//added 21april2013
                                                                                wit_orr_toSend->addBoundary(from);

                                                                                if(  success == 1 || (contin && !breath ) )  //  ||   //if(cgKprime->getK() <= 6  ){
                                                                                {
//
Start_t = clock(); if(!contin) this->atlas->connect(from, nodeID);   //cgK,
cgKprime[c]);
//
End_t = clock();
//
time_mapView += difftime(End_t, Start_t);

                                                                                        if(c==0) //add once
                                                                                        {
//
Start_t = clock(); this->snl->appendSpacePoints( rnode );//from, cgK, noPoints);
                                                                                                noPoints=0;
//
End_t = clock();
//
time_appendSpacePoints += difftime(End_t, Start_t);
                                                                                        }
//
if(!contin)
//
this->snl->saveNode(  this->atlas->getNode(  nodeID )  );
                                                                                        //wit_orr_toSend = sr->getOrienation(); //commented 21april2013
                                                                                        //wit_orr_toSend->addBoundary(from);

//
cout<< "before call to sample" << endl;
//
detector->test(sr);

//
Start_t = clock(); if( this->MySample(rnode_ID, false, wit_orr_toSend, contin,
breath, bret) ) //cgKprime[c]  //if there is a realization with good angle at
the childs, noGoodOrientation = false;
//
End_t = clock();
//
time_recursive += difftime(End_t, Start_t);

                                                                                        ori_on_lattice->addBoundary( nodeID );  //19april2013
                                                                                        orie_on_boundary->addBoundary( nodeID );

//
if( Settings::Saving::saveBoundary )//false //|| Settings::OldOrTesting::just5d
//19april2013   //true:28nov12(since we do not save it to child)   //do not add
witness point to parent node!
// {
//
if( c==contactList2.size() ) //added on 22march2013 to add same orientation ONCE
at the last iteration
//
{
//
vector<double> outt = desc->getPoint(); //NOT SURE STILL NEED p4 INSTEAD OF p4d
//
PointMultiD* p4 = new PointMultiD( outt ); //maybe the reason is to save exact
position of point after binary search
//
//
if( !Settings::Collision::checkAngle || !orie_on_boundary->angleViolated()){
//
p4->addOrientation(orie_on_boundary); //DO NOT SAVE SAME ORIENTATION FOR EACH
BOUNDARY, ADD ALL BOUNDARIES TO ONE ORIENTATION THEN SAVE THAT ORIENTATION
//
addedOr = true;
//
}
//
else{
//
delete orie_on_boundary; //orange
//
p4->angleViolated() = true;
//
}
//
cgK->AddSamplePoint(p4); //should it be insertwitness
!!!!!!!!!!!!!!!!!!!!!!!!!!!!no, witness is the one saved at the child node
//
noPoints++;
//
}
// }


                                                                                        //add witness
                                                                                        //p4d->getBoundaries()
                                                                                }
                                                                                else   // if( !( nodeID == -1 || (contin && !breath) ) ){           //done
                                                                                {
                                                                                        this->atlas->connect(from, nodeID);   //cgK,cgKprime[c]);

                                                                                        ori_on_lattice->addBoundary( nodeID );  //19april2013
                                                                                        orie_on_boundary->addBoundary(nodeID);

//
if(Settings::Saving::saveBoundary )//false //|| Settings::OldOrTesting::just5d
//19april2013   //true:28nov12(since we do not save it to child)   //do not add
witness point to parent node!
// {
//
if( c==contactList2.size() ) //added on 22march2013 to add same orientation ONCE
at the last iteration
//
{
//
if( !Settings::Collision::checkAngle || !orie_on_boundary->angleViolated()){
//
p4d->addOrientation(orie_on_boundary);  //DO NOT SAVE SAME ORIENTATION FOR EACH
BOUNDARY, ADD ALL BOUNDARIES TO ONE ORIENTATION THEN SAVE THAT ORIENTATION
//
addedOr = true;
//
}
//
else{
//
delete orie_on_boundary; //orange
//
p4d->angleViolated() = true;
//
}
//
}
// }


                                                                                        //wit_orr_toSend = sr->getOrienation() ;  //commented 21april2013
                                                                                        if( false && (!Settings::Collision::checkAngle || !wit_orr_toSend->angleViolated()) ){  //do not add witness to child, it is alot of duplicates 28nov12
                                                                                                wit_orr_toSend->addBoundary(from);
                                                                                                PointMultiD* wp4 = desc->getWitnessPoint( cgKprime[c], wit_orr_toSend);//copy of wit_orr_toSend is added inside the finction SHOUL DELETE wit_orr_toSend????
                                                                                                //cgKprime[c]->AddWitnessPoint(wp4);  //EVEN DO NOT NEED TO ADD CGK, SAVING TO FILE SHOULD BE ENOUGH

//
Start_t = clock(); this->snl->appendWitness(rnode_ID, wp4);  //AFTER SAVING
SHOULD TRIM CGK  AGAIN?

                                                                                                wp4->trim_PointMultiD();  delete wp4;  //added 10july12

//
End_t = clock();
//
time_appendSpacePoints += difftime(End_t, Start_t);

                                                                                                rnode_ID->setFoundGoodOrientation(true);
//
if( !cgKprime[c]->hasAnyGoodOrientation() ){
//
cgKprime[c]->setFoundGoodOrientation();
//
this->snl->setFoundAGoodOrientation( this->atlas->getNode(nodeID) );
//
}
                                                                                        }

                                                                                        if( rnode_ID->hasAnyGoodOrientation() )//if the child has any good orientation, then parent should be displayed.
                                                                                                noGoodOrientation = false;


                                                                                }

// Start_t = clock(); createChildNodes(rnode_ID, contactList2, wit_orr_toSend,
this->atlas, this->snl, firstPath); firstPath = false;
// End_t = clock();
//
time_createchildnodes += difftime(End_t, Start_t);

                                                                                contactList2.push_back(contact);
                                                                                delete wit_orr_toSend;
                                                                        }//end
of FOR



                                                                        if(
Settings::Saving::saveBoundary )//false //|| Settings::OldOrTesting::just5d
//19april2013   //true:28nov12(since we do not save it to child)   //do not add
witness point to parent node!
                                                                        {
                                                                                vector<double> outt = desc->getPoint(); //NOT SURE STILL NEED p4 INSTEAD OF p4d
                                                                                PointMultiD* p4 = new PointMultiD( outt ); //maybe the reason is to save exact position of point after binary search

                                                                                if( !Settings::Collision::checkAngle || !orie_on_boundary->angleViolated()){
                                                                                        p4->addOrientation(orie_on_boundary); //DO NOT SAVE SAME ORIENTATION FOR EACH BOUNDARY, ADD ALL BOUNDARIES TO ONE ORIENTATION THEN SAVE THAT ORIENTATION
                                                                                        //addedOr = true; //commend it to let ori_on_lattice to be saved
                                                                                }
                                                                                else{
                                                                                        delete orie_on_boundary; //orange
                                                                                        p4->angleViolated() = true;
                                                                                }
                                                                                cgK->AddSamplePoint(p4); //should it be insertwitness !!!!!!!!!!!!!!!!!!!!!!!!!!!!no, witness is the one saved at the child node
                                                                                noPoints++;
                                                                        }


                                                                        delete
sr; break;//if found a collision in one direction then stop
                                                                }
                                                        }
//							if(sr != NULL)
                                                        delete sr;
                                                        desc->setInitialPoint(pp);
                                                }
//						if(!entered) //does not seem to
be able to be contact
//							cout <<
"UNVALIDDDDDDDDDDDDDDDD 0" << endl;

                                                desc->setInitialPoint(pp);
                                                desc->setDir();

                                        }
                                        if(false && ( noEquality ==
cgK->getVolumeNo()+1  && !addedOr )) //tetrahedra equality
                                        {
                                                cout << " TETRAHEDRA equality
"<< endl; ActiveConstraintGraph *cgKprimet = new ActiveConstraintGraph(
(*r)->getCG() ); for(int i =0; i<3; i++){ cgKprimet->setVolumeZero(i,
(*r)->isVolumeZero(i) );
                                                }
                                                cgKprimet->setVolumeNo(noEquality);
                                                //	if( !isnan(
(*r)->getIndicesVolumeZero(i)[0]) ) int inds[3][4];
                                                (*r)->getIndicesVolumeZero(inds);
                                                cgKprimet->setIndicesVolumeZero(
inds );


                                                int nodeID;  // =
this->atlas->getNodeNum(cgKprimet); int success =
this->atlas->addNode(cgKprimet, nodeID); AtlasNode* rnode_ID =
this->atlas->getNode(nodeID); bool contin = false; if(success == 0){ //started
maybe done delete cgKprimet; cgKprimet = NULL; cgKprimet = rnode_ID->getCG();
//atlas->getCgkofNode(nodeID); if(!rnode_ID->isComplete() ){ //not done
                                                                snl->loadNode(nodeID,
cgKprimet); contin = true;
                                                        }
                                                }

                                        if( ( success == 1 || (contin && !breath
) ) ) //not done or not completed  //|| cgKprime->isIndicesEqual( (*r)->getCG())
//cgKprime->getK() <= 6-noEquality &&
                                                {
                                                if(!contin)
                                                        this->atlas->connect(from,nodeID);
                                                this->snl->appendSpacePoints(
rnode );  //from, cgK, noPoints); noPoints=0;
//   							if(!contin)
//
this->snl->saveNode( this->atlas->getNode( nodeID ) );

                                                this->MySample(rnode_ID, false,
(*r)->getOrienation(), contin, breath, bret ); //cgKprimet

                                                        (*r)->setBoundary(
nodeID );//set boundary on CartesianRealizer to then be transferred to
orientation. cgKprimet->insertBoundary(from , (*r)->getOrienation() );

                                                }
                                        else  {  //done
                                                        this->atlas->connect(from,
nodeID);
                                                        //int newNodenum =
this->atlas->getNodeNum(cgKprimet);
                                                        (*r)->setBoundary(nodeID);
//newNodenum);//set boundary on CartesianRealizer to then be transferred to
orientation.

                                                                //ActiveConstraintGraph*
mapCgk = atlas->getCgkofNode(newNodenum); if(!cgKprimet->isWitness(from))
                                                        {
                                                                //double
positions[NO][3];
                                                                // the cgKprime
pointer is not the pointer within the roadmap
                                                                //(*r)->getPositions(positions);
                                                                Orientation* orr
= (*r)->getOrienation() ; orr->addBoundary(atlas->getNodeNum((*r)->getCG()));
                                                                PointMultiD* wp4
= desc->getWitnessPoint( cgKprimet ,orr); wp4->addOrientation(orr );
                                                                //cgKprimet->AddWitnessPoint(wp4);
//commented on 27may12 this->snl->appendWitness(rnode_ID, wp4);  //newNodenum
//								delete cgKprime;
//								cgKprime =
mapCgk;
                                                        }

                                        }

                                        //p4d->addOrientation((*r)->getOrienation());
                                }



                                        //COMMENTED ON 24NOV10 TO GET RID OF
COMPILE ERROR cgKprime, CHECK LATER TO MAKE SURE WHAT IS IT FOR !!!
//					int prime =
this->atlas->getNodeNum(cgKprime);
////					int now = this->atlas->getNodeNum(cgK);
//					if( prime != from && prime >= 0 &&
(*r)->getBoundary() < 0 ){
//
(*r)->setBoundary(this->atlas->getNodeNum(cgKprime));//set boundary on
CartesianRealizer to then be transferred to orientation.
//					}

                                        if(!addedOr){
                                                if(
Settings::Saving::saveOutGridOrientation || !Settings::Collision::checkAngle ||
!ori_on_lattice->angleViolated()){ p4d->addOrientation(ori_on_lattice);
                                                }
                                                else{
                                                        delete ori_on_lattice;
//orange p4d->angleViolated() = true;
                                                }
                                        }


                                } //endif(myvalid)

                        } //endif real iterator


//#ifdef WIN32
                    for(list<CartesianRealizer*>::iterator it = real.begin(); it
!= real.end(); it++){ delete (*it);
                    }
//#endif
                    real.clear();


//        	//---------------------- numerical jacobian sampling
//            for(list<CartesianRealizer*>::iterator it = inbetween.begin();it
!= inbetween.end();it++)
//           		delete (*it);
//            inbetween.clear();
//
//		    for(list<CartesianRealizer*>::iterator it =
oldreal.begin();it != oldreal.end();it++)
//					delete (*it);
//			oldreal.clear();
//			//------------------------------


                    cgK->AddSamplePoint(p4d);
                    noPoints++;
                    if( noPoints == 10000 ) //100
                    {
//		    	if(from==1)  cout <<"i am here 3"<<endl;
//		    	Start_t = clock();
                        this->snl->appendSpacePoints( rnode );
//		    	End_t = clock();
//		    	time_appendSpacePoints += difftime(End_t, Start_t);
//		    	if(from==1)  cout <<"i am here 4"<<endl;
                        noPoints=0;
                    }


                }
    }

                if(noPoints != 0){
//			Start_t = clock();
                        this->snl->appendSpacePoints( rnode );
//			End_t = clock();
//			time_appendSpacePoints += difftime(End_t, Start_t);
                noPoints=0;
                }

//		time_t end_mysample = clock();
//		int time_total = difftime(end_mysample, Start_mysample);
//		int time_computation = time_total - time_recursive -
time_appendSpacePoints - time_createchildnodes;
//		cout << "node " << from << " time_appendSpacePoints " <<
time_appendSpacePoints << " time_computation " << time_computation << "
time_createchildnodes " << time_createchildnodes << " ratio " <<
100.*time_appendSpacePoints/(time_computation+time_appendSpacePoints)<<  "
time_findReal " << time_findReal << " loop " << loop << " time_mapView "<<
time_mapView << " time_myvalid " << time_myvalid << " time_w " << time_w <<
endl;

                if(!noGoodOrientation ){
                        rnode->setFoundGoodOrientation(true);
//			cgK->setFoundGoodOrientation(); //this line doesnot
work, the below is doing the job I guess! cgK is not the real cgk in the
roadmap?
//			atlas->getCgkofNode(from)->setFoundGoodOrientation();
//
//			this->snl->setFoundAGoodOrientation( rnode );
                }


                if(!Settings::AtlasBuilding::stop && !breath){
//			cgK->setComplete(true);
////			int h = atlas->getNodeNum(cgK);
////			vector<AtlasNode*> nodes2 = atlas->getNodes();
////			nodes2[h]->getCG()->setComplete(true);
//			atlas->getCgkofNode(from)->setComplete(true);
//			this->snl->writeComplete( rnode);

                        rnode->setComplete(true);

                        cgK->trim();
//			atlas->getCgkofNode(from)->trim();  //commented on
27may12 since now i am sure that cgK is consistent with cgk of roadnode

                }

                delete detector;
                delete desc;

                //snl->saveNode( rnode );  //save it here because to have recent
info for the completeness and anygoodorientation

                return !noGoodOri
                * entation; //should delete this node and the path. // no do not
delete, keep it for reference not to sample it again


}
*/
