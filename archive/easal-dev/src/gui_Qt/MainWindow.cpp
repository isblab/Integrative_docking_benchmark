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
#include "MainWindow.h"

#include <thread>

MainWindow::MainWindow() {
  setWindowTitle("EASAL");
  setFocusPolicy(Qt::StrongFocus);
  this->resize(1200, 800);
  // INIT WIDGETS
  bondViewWidget = new BondViewWidget();
  realizationViewWidgetCornerVersion =
      new RealizationViewWidgetCornerVersion(&sharedData, bondViewWidget);
  sweepViewWidget = new SweepViewWidget(&sharedData);
  cayleySpaceViewWidget = new CayleySpaceViewWidget(&sharedData);
  inputWindowInstance = new InputWindow;
  statsWindowInstance = new StatsWindow();
  pseudoAtlasWindowInstance = new PseudoAtlasWindow();
  cayleyStatsWindowInstance = new CayleyStatisticsWindow(&sharedData);
  aboutWindowInstance = new AboutWindow();
  atlasViewWidget = new AtlasViewWidget(&sharedData);
  nodeParentChildernViewWidget = new NodeParentChildernViewWidget(&sharedData);
  stackedWidget = new QStackedWidget();
  tabWidget = new QTabWidget();
  realizationViewWidget = new RealisationViewWidget(&sharedData);
  pathfindingViewWidget = new PathfindingViewWidget3DPart(&sharedData);
  constrainSelectionWindow = new ConstraintSelectionWindow();
  setLayouts();
  //[ INIT WIDGETS]
  //[ HANDEL TIMER]
  updateTimeInMS = 100;
  timer = new QTimer(this);
  // start argument determins update rate
  timer->start(updateTimeInMS);
  connect(timer, SIGNAL(timeout()), this, SLOT(updateSlot()));
  //[~HANDEL TIMER]

  //[ INIT PATHFINDING DATA]
  for (int i = 0; i < 2; i++) {
    currentPathfindingSelectedAtlasNodeID[i] = -1;
    currentPathfindingSelectedCayleyPointType[i] = NoType;
    currentPathfindingSelectedCayleyPointID[i] = -1;
  }
  //[~INIT PATHFINDING DATA]
  //[ MISC]
  stopStartSampleButton.setText("Pause Sampling ");
  runModeSelectionComboBox.addItem("TreeSample");
  runModeSelectionComboBox.addItem("ForestSample");
  runModeSelectionComboBox.addItem("ForestSampleAS");
  runModeSelectionComboBox.addItem("BreadthFirst");
  runModeSelectionComboBox.addItem("RefineSampling");
  runModeSelectionComboBox.addItem("Constraint selection");
  runModeSelectionComboBox.setCurrentText("ForestSampleAS");

  findPathBetweenD0NodesButton.setText("Find path");
  pathFindingStatusLabel.setText("");

  pathfindingButtonText[START] = "Select start 0D node";
  pathfindingButtonText[END] = "Select end 0D node";

  pathfindingButton[START].setText(pathfindingButtonText[START]);
  pathfindingButton[END].setText(pathfindingButtonText[END]);

  pathfindingValueLabel[START].setText("Start D0 Node:");
  pathfindingValueLabel[END].setText("End D0 Node:");

  pathfindingFlipLabel[START].setText("Flip:");
  pathfindingFlipLabel[END].setText("Flip:");
  //    runModeSelectionComboBox.setVisible(false);
  runModeSelectionComboBox.setDisabled(true);

  connect(&stopStartSampleButton, SIGNAL(clicked(bool)), this,
          SLOT(stopStartSamplingButtonSlot()));
  connect(&findPathBetweenD0NodesButton, SIGNAL(clicked(bool)), this,
          SLOT(findPathBetweenD0NodesButtonSlot()));
  connect(&pathfindingButton[START], SIGNAL(clicked(bool)), this,
          SLOT(pathfindingStartNodeButtonSlot()));
  connect(&pathfindingButton[END], SIGNAL(clicked(bool)), this,
          SLOT(pathfindingEndNodeButtonSlot()));
  activeWidgetID = 1;

  //[~MISC]
}

void MainWindow::paStats() { pseudoAtlasWindowInstance->exec(); }

void MainWindow::atlasStats() { statsWindowInstance->exec(); }

void MainWindow::cayleyStats() { cayleyStatsWindowInstance->exec(); }

void MainWindow::openInputWindow() {
  Settings *sett = Settings::getInstance();
  inputWindowInstance->exec();
  if (sett->getRunMode() != Stopped) {
    createAction->setEnabled(false);
    processDataFromSettings();
    sharedData.isSamplingStarted = true;
  } else {
    createAction->setEnabled(false);
    processDataFromSettings();
  }
}

void MainWindow::openPaperURL() {
  QDesktopServices::openUrl(QUrl("https://arxiv.org/abs/1805.07450"));
}

void MainWindow::openUserGuide() {
  QDesktopServices::openUrl(
      QUrl("./doc/TOMS paper/User Guide/CompleteUserGuide.pdf"));
}

void MainWindow::openAboutMenu() { aboutWindowInstance->exec(); }

void MainWindow::stopStartSamplingButtonSlot() {
  Settings *sett = Settings::getInstance();
  if (sharedData.isSamplingStarted) {
    isSamplingStopped = !isSamplingStopped;

    if (isSamplingStopped) {
      stopStartSampleButton.setText("Resume Sampling");

      // Set current text with name of current sampling mode for dropdown that
      // allow to select sampling mode
      if (sett->getRunMode() == Modes::TreeSample) {
        runModeSelectionComboBox.currentText() = "TreeSample";
      } else if (sett->getRunMode() == Modes::ForestSample) {
        runModeSelectionComboBox.currentText() = "ForestSample";
      } else if (sett->getRunMode() == Modes::ForestSampleAS) {
        runModeSelectionComboBox.currentText() = "ForestSampleAS";
      } else if (sett->getRunMode() == Modes::BreadthFirst) {
        runModeSelectionComboBox.currentText() = "BreadthFirst";
      } else if (sett->getRunMode() == Modes::RefineSampling) {
        runModeSelectionComboBox.currentText() = "RefineSampling";
      }

      runModeSelectionComboBox.setDisabled(false);
      sett->AtlasBuilding.stop = true;
      sett->setRunMode(Modes::Stopped);
    } else {
      saveLoader->saveRoadMap(atlas);
      stopStartSampleButton.setText("Pause Sampling ");
      if (runModeSelectionComboBox.currentText() == "TreeSample") {
        sett->AtlasBuilding.breadthFirst = false;
        sett->setRunMode(Modes::TreeSample);
      } else if (runModeSelectionComboBox.currentText() == "ForestSample") {
        sett->AtlasBuilding.breadthFirst = false;
        sett->setRunMode(Modes::ForestSample);
      } else if (runModeSelectionComboBox.currentText() == "ForestSampleAS") {
        sett->AtlasBuilding.breadthFirst = false;
        sett->setRunMode(Modes::ForestSampleAS);
      } else if (runModeSelectionComboBox.currentText() == "BreadthFirst") {
        setSamplingNodeNumber(atlasViewWidget->getCurrentNode()->getID());
        sett->setRunMode(Modes::BreadthFirst);
      } else if (runModeSelectionComboBox.currentText() == "RefineSampling") {
        sett->setRunMode(Modes::RefineSampling);
      } else if (runModeSelectionComboBox.currentText() ==
                 "Constraint selection") {
        constrainSelectionWindow->exec();
        ActiveConstraintGraph *graph =
            constrainSelectionWindow->getConstraintGraph();
        if (graph != NULL) {
          graph->setStepSize(sett->Sampling.stepSize);
          int num = ThreadShare::atlas_view->findNodeNum(graph);
          if (num > -1) {
            delete graph;
            // atlasViewWidget->setCurrentNode(num);
          } else {
            ThreadShare::atlas_view->addNode(graph, num, ROOT_NODE_PARENT_ID);

            AtlasNode *rnode_ID = (*ThreadShare::atlas_view)[num];
            ThreadShare::save_loader->appendDimension(rnode_ID);
            setSamplingNodeNumber(rnode_ID->getID());
          }
        }
        sett->setRunMode(Modes::TreeSample);
      }

      sett->AtlasBuilding.stop = false;
      runModeSelectionComboBox.setDisabled(true);
    }
  }
}

void MainWindow::findPathBetweenD0NodesButtonSlot() {
  MessageWindow messageWindow("Error: Wrong Pathfinding Data");

  bool startNodeIDIsInt = true;
  bool endNodeIDIsInt = true;
  int startNodeID = sharedData.pathfindingSelectedAtlasNodeID[START];
  int startFlip = pathfindingFlipSelectionComboBox[0].currentText().toInt();
  int endNodeID = sharedData.pathfindingSelectedAtlasNodeID[END];
  int endFlip = pathfindingFlipSelectionComboBox[1].currentText().toInt();
  if (!startNodeIDIsInt) {
    messageWindow.addLineToMessage(QString("Start node is not integer"));
  } else if (startNodeID < 0 || startNodeID >= atlas->number_of_nodes()) {
    messageWindow.addLineToMessage(
        QString("Start node is outside allowed range(0 to ") +
        QString::number(atlas->number_of_nodes()));
  } else if ((*atlas)[startNodeID]->getDim() != 0) {
    messageWindow.addLineToMessage(
        QString("Start node is not dimension 0 node"));
  }
  if (!endNodeIDIsInt) {
    messageWindow.addLineToMessage(QString("End node is not integer"));
  } else if (endNodeID < 0 || endNodeID >= atlas->number_of_nodes()) {
    messageWindow.addLineToMessage(
        QString("End node is outside allowed range(0 to ") +
        QString::number(atlas->number_of_nodes()));
  } else if ((*atlas)[endNodeID]->getDim() != 0) {
    messageWindow.addLineToMessage(QString("End node is not dimension 0 node"));
  }

  if (!messageWindow.isEmpty()) {
    messageWindow.exec();
  } else {
    vector<PathfindingReturnTuple> path;

    PathFinder *path_finder = new PathFinder(atlas);
	
    // PathFinder *pf = new PathFinder(atlas,startNodeID, startFlip, endNodeID,
    // endFlip); std::vector<std::pair<CayleyPoint*, int> *> Atlaspath =
    // pf->find1DoFPath();

    // std::vector<std::pair<CayleyPoint*, int> *>
    // Atlaspath =  atlas->find1DPath(startNodeID, endNodeID);

    std::vector<std::tuple<int, int, int> > atlas_path = 
			path_finder->FindCayleyPath(startNodeID, 
            					startFlip, endNodeID, endFlip);
	/* =
    std::vector<std::tuple<int, CayleyPoint *, int> > Atlaspath; =
        atlas->find1DoFPath(
            startNodeID, endNodeID,
            pathfindingFlipSelectionComboBox[0].currentText().toInt(),
            pathfindingFlipSelectionComboBox[1].currentText().toInt());*/

    if (atlas_path.size() == 0) {
      cout << "Got an empty path" << endl;
      pathFindingStatusLabel.setText("Path not found.");
    } else {
      cout << "Got a path" << endl;
      pathFindingStatusLabel.setText("Path found.");
      for (auto& point: atlas_path) {
	    PathfindingReturnTuple tmpPathReturnTuple;
		auto& [atlas_node_id, cayley_point_id, flip] = point;
		AtlasNode* atlas_node = atlas->getNode(atlas_node_id);
        Orientation *tmpOrientation = atlas_node->getACR()->GetCayleyPointFromID(cayley_point_id)->getOrientation(flip);
        tmpPathReturnTuple.atlasNodeID = atlas_node_id;
        tmpPathReturnTuple.orientation = tmpOrientation;
        if (tmpOrientation != NULL) {
          path.push_back(tmpPathReturnTuple);
        }
        //(std::get<0>(*Atlaspath[i]))->printData();
        // cout<<(std::get<1>(*Atlaspath[i]));
        // cout<<endl;
      }
      pathfindingViewWidget->setPath(path);
    }

    // form path manually for testing purposes

    /*
    for(int i=startNodeID;i<=endNodeID;i++)
    {
        if((*atlas)[i]->getDim()==0 ||(*atlas)[i]->getDim()==1 )
        {
            //Get ACR to get orientations from hard drive
            ActiveConstraintRegion *tmpACR = new ActiveConstraintRegion();
            saveLoader->loadNode(i,tmpACR);
            for(int j=0;j<tmpACR->space.size();j++)
            {
                vector<Orientation*> tmpOrientations=
    tmpACR->space[j]->getOrientations(); for(int
    k=0;k<tmpOrientations.size();k++)
                {
                    PathfindingReturnTuple tmpPathReturnTouple;
                    Orientation *tmpOrientation = new
    Orientation(tmpOrientations[k]); tmpPathReturnTouple.atlasNodeID = i;
                    tmpPathReturnTouple.orientation = tmpOrientation;
                    path.push_back(tmpPathReturnTouple);
                }
            }
            delete tmpACR;
        }
    }*/
  }
}

void MainWindow::updateSlot() {
  Settings *sett = Settings::getInstance();
  if (sharedData.isSamplingStarted) {
    // Update info about sampling
    if (sharedData.atlasNodeWasSet) {
      currentSelectedNodeLabel.setText(
          "Selected Node = " +
          QString::number(sharedData.currentAtlasNode->getID()));
    } else {
      currentSelectedNodeLabel.setText("Selected Node = NONE");
    }

    currentSampledNodeLabel.setText("Sampled Node = " +
                                    QString::number(atlas->number_of_nodes()));
    currentRunModeLabel.setText("Run Mode = " +
                                samplingModesNames[sett->getRunMode()]);
    // Update pathfinding related data

  }
}

void MainWindow::pathfindingStartNodeButtonSlot() {
  if(sharedData.atlasNodeWasSet) {
    sharedData.pathfindingFlipVector[0].clear();
    cayleySpaceViewWidget->setAtlasNode(sharedData.currentAtlasNode);
    sharedData.pathfindingSelectedAtlasNodeID[0] = sharedData.currentAtlasNode->getID();
    sharedData.currentAtlasNode->getWitnessFlips(sharedData.pathfindingFlipVector[0]);
    sharedData.pathfindingSelectCayleyPointNow[START] = true;
    pathfindingButton[START].setText("Update start node");
	pathfindingFlipSelectionComboBox[0].clear();
    for (int j = 0; j < sharedData.pathfindingFlipVector[0].size(); j++) {
      pathfindingFlipSelectionComboBox[0].addItem(
      QString::number(sharedData.pathfindingFlipVector[0][j]));
    }
    pathfindingValueLabel[0].setText("Selected node: " + QString::number
										(sharedData.currentAtlasNode->getID())); 
  }
}

void MainWindow::pathfindingEndNodeButtonSlot() {
  if(sharedData.atlasNodeWasSet) {
    sharedData.pathfindingFlipVector[1].clear();
    cayleySpaceViewWidget->setAtlasNode(sharedData.currentAtlasNode);
    sharedData.pathfindingSelectedAtlasNodeID[1] = sharedData.currentAtlasNode->getID();
    sharedData.currentAtlasNode->getWitnessFlips(sharedData.pathfindingFlipVector[1]);
    sharedData.pathfindingSelectCayleyPointNow[END] = true;
    pathfindingButton[END].setText("Update end node");
	pathfindingFlipSelectionComboBox[1].clear();
    for (int j = 0; j < sharedData.pathfindingFlipVector[1].size(); j++) {
      pathfindingFlipSelectionComboBox[1].addItem(
      QString::number(sharedData.pathfindingFlipVector[1][j]));
    }
    pathfindingValueLabel[1].setText("Selected node: " + QString::number
										(sharedData.currentAtlasNode->getID())); 
  }
}

void MainWindow::processDataFromSettings() {
  Settings *sett = Settings::getInstance();
  // PredefinedInteractions *df = new PredefinedInteractions();
  // init_DistanceFinder_from_settings(df);

  // molecular unit A and B
  // muA->init_MolecularUnit_A_from_settings(df);
  // muB->init_MolecularUnit_B_from_settings(df);

  // SaveLoader
  saveLoader->setData(sett->Output.dataDirectory, sett->runTimeObjects.muA,
                      sett->runTimeObjects.muB);

  if (sett->Constraint.wholeCollision) {
    // reading the neighbour matrix
    ConstraintCheck::nei_matrix = Utils::getMatrixFromFileFromTo("nei.txt");
  } else {
    ConstraintCheck::nei_matrix = Utils::getIdentityMatrix();
  }

  atlas_builder->setData(sett->runTimeObjects.muA, sett->runTimeObjects.muB,
                         saveLoader, sett->runTimeObjects.df, atlas);
  sweepViewWidget->setMolecularUnits(sett->runTimeObjects.muA,
                                     sett->runTimeObjects.muB);
  sweepViewWidget->setSaveLoader(saveLoader);
  cayleySpaceViewWidget->setData(atlas, saveLoader);

  if (sett->Sampling.runSample) {
    sett->setRunMode(ForestSampleAS);
    atlas_builder->setup();
    cout << "AtlasBuilder Set up done." << endl;

  } else {
    saveLoader->loadMapView(atlas);
    cout << "Loads existing atlas." << endl;
  }

  ThreadShare::save_loader = saveLoader;
  ThreadShare::atlas_builder = atlas_builder;
  ThreadShare::atlas_view = atlas;

  *backEndThread = new std::thread(BackEnd_Thread, 0);
  cout << "BackEnd_Thread: Started." << endl;

  atlasViewWidget->setComponents(atlas, saveLoader);
  nodeParentChildernViewWidget->setComponents(atlas, saveLoader);

  realizationViewWidgetCornerVersion->setMolecularUnits(
      sett->runTimeObjects.muA, sett->runTimeObjects.muB);
  realizationViewWidgetCornerVersion->setSaveLoader(saveLoader);

  realizationViewWidget->setMolecularUnits(sett->runTimeObjects.muA,
                                           sett->runTimeObjects.muB);
  realizationViewWidget->setSaveLoader(saveLoader);

  pathfindingViewWidget->setData(atlas, saveLoader, sett->runTimeObjects.muA,
                                 sett->runTimeObjects.muB);

  constrainSelectionWindow->setData(atlas, sett->runTimeObjects.muA,
                                    sett->runTimeObjects.muB);
}

void MainWindow::setLayouts() {
  // Layout looks as vertical box with group box inside
  // Vertical box also contains menu bar
  mainLayout = new QVBoxLayout;
  // Group box contains atlas view and molecule widget
  mainGroupBox = new QWidget();
  // main layout belongs to main widget
  mainWidget = new QWidget();

  //[ HANDLE MENU]
  menuBar = new QMenuBar(mainWidget);

  //[FILE]
  fileMenu = new QMenu("File");
  createAction = new QAction(tr("Create"), this);
  connect(createAction, SIGNAL(triggered()), this, SLOT(openInputWindow()));
  closeAction = new QAction(tr("Close"), this);
  connect(closeAction, SIGNAL(triggered()), this, SLOT(close()));
  fileMenu->addAction(createAction);
  fileMenu->addAction(closeAction);

  //[HELP]
  helpMenu = new QMenu("Help");

  paperAction = new QAction(tr("EASAL Theory and ALgorithms"), this);
  helpMenu->addAction(paperAction);
  connect(paperAction, SIGNAL(triggered()), this, SLOT(openPaperURL()));

  userGuideAction = new QAction(tr("EASAL User/Developer Guide"), this);
  helpMenu->addAction(userGuideAction);
  connect(paperAction, SIGNAL(triggered()), this, SLOT(openUserGuide()));

  aboutAction = new QAction(tr("About"), this);
  helpMenu->addAction(aboutAction);
  connect(aboutAction, SIGNAL(triggered()), this, SLOT(openAboutMenu()));

  //[STATISTICS]
  statsMenu = new QMenu("Statistics");
  basicStatisticsMenu = new QMenu("Basic Statistics");
  atlasStatsAction = new QAction(tr("Atlas Statistics"), this);
  connect(atlasStatsAction, SIGNAL(triggered()), this, SLOT(atlasStats()));
  nodeStatsAction = new QAction(tr("Node Statistics"), this);
  dimensionStatsAction = new QAction(tr("Sampling Time"), this);
  cayleyStatsAction = new QAction(tr("Cayley Statistics"), this);
  connect(cayleyStatsAction, SIGNAL(triggered()), this, SLOT(cayleyStats()));
  basicStatisticsMenu->addAction(atlasStatsAction);
  basicStatisticsMenu->addAction(nodeStatsAction);
  basicStatisticsMenu->addAction(dimensionStatsAction);
  basicStatisticsMenu->addAction(cayleyStatsAction);

  gridStatisticsMenu = new QMenu("Grid Statistics");
  epsilonCoverageAction = new QAction(tr("Epsilon Coverage"), this);
  gridCoverageAction = new QAction(tr("Grid Coverage"), this);
  gridStatisticsMenu->addAction(epsilonCoverageAction);
  gridStatisticsMenu->addAction(gridCoverageAction);

  statsMenu->addMenu(basicStatisticsMenu);
  statsMenu->addMenu(gridStatisticsMenu);

  //[Atlas]
  atlasMenu = new QMenu("Atlas");
  basinMenu = new QMenu("Basin");
  basinVolumeAction = new QAction("Basin Volume");
  basinRegionsAction = new QAction("Basin Regions");
  basinMenu->addAction(basinVolumeAction);
  basinMenu->addAction(basinRegionsAction);
  atlasMenu->addMenu(basinMenu);
  paAction = new QAction(tr("Pseudo Atlas"), this);
  atlasMenu->addAction(paAction);
  connect(paAction, SIGNAL(triggered()), this, SLOT(paStats()));

  // Add the menus to the menu bar
  menuBar->setFixedHeight(22);
  menuBar->addMenu(fileMenu);
  menuBar->addMenu(statsMenu);
  menuBar->addMenu(atlasMenu);
  menuBar->addMenu(helpMenu);

  //[~HANDLE MENU]
  //[ HANDLE STACKEDWIDGET]
  tabWidget->addTab(atlasViewWidget, "Atlas");
  tabWidget->addTab(nodeParentChildernViewWidget, "Node neighbors");
  tabWidget->addTab(sweepViewWidget, "Sweep");
  tabWidget->addTab(cayleySpaceViewWidget, "Cayley space");
  tabWidget->addTab(realizationViewWidget, "Realization");
  tabWidget->addTab(pathfindingViewWidget, "Pathfinding");
  tabWidget->setCurrentWidget(atlasViewWidget);
  //[~HANDLE STACKEDWIDGET]

  mainGroupBoxLayout = new QGridLayout;

  mainGroupBoxLayout->setContentsMargins(1, 1, 1, 1);
  atlasViewWidget->setContentsMargins(1, 1, 1, 1);
  nodeParentChildernViewWidget->setContentsMargins(1, 1, 1, 1);
  sweepViewWidget->setContentsMargins(1, 1, 1, 1);
  cayleySpaceViewWidget->setContentsMargins(1, 1, 1, 1);
  realizationViewWidgetCornerVersion->setContentsMargins(1, 1, 1, 1);
  tabWidget->setContentsMargins(0, 0, 0, 0);
  mainWidget->setContentsMargins(1, 1, 1, 1);
  menuBar->setContentsMargins(1, 1, 1, 1);
  mainLayout->setContentsMargins(1, 1, 1, 1);
  mainGroupBox->setContentsMargins(1, 1, 1, 1);
  this->setContentsMargins(1, 1, 1, 1);

  // Set interface layout (the one that contains all buttons and active
  // constrain graph)
  int interfaceLayoutRow = 0;
  interfaceLayout.addWidget(bondViewWidget, interfaceLayoutRow, 0, 3, 1);
  interfaceLayout.addWidget(&currentSelectedNodeLabel, interfaceLayoutRow, 1, 1,
                            1);
  interfaceLayoutRow++;
  interfaceLayout.addWidget(&currentSampledNodeLabel, interfaceLayoutRow, 1, 1,
                            1);
  interfaceLayoutRow++;
  interfaceLayout.addWidget(&currentRunModeLabel, interfaceLayoutRow, 1, 1, 1);
  interfaceLayoutRow++;
  interfaceLayout.addWidget(&findPathBetweenD0NodesButton, interfaceLayoutRow,
                            0, 1, 2);
  interfaceLayoutRow++;
  interfaceLayout.addWidget(&pathfindingButton[START], interfaceLayoutRow, 0, 1,
                            1);
  interfaceLayout.addWidget(&pathfindingFlipLabel[START], interfaceLayoutRow, 1,
                            1, 1);
  interfaceLayoutRow++;
  interfaceLayout.addWidget(&pathfindingValueLabel[START], interfaceLayoutRow,
                            0, 1, 1);
  interfaceLayout.addWidget(&pathfindingFlipSelectionComboBox[START],
                            interfaceLayoutRow, 1, 1, 1);
  interfaceLayoutRow++;
  interfaceLayout.addWidget(&pathfindingButton[END], interfaceLayoutRow, 0, 1,
                            1);
  interfaceLayout.addWidget(&pathfindingFlipLabel[END], interfaceLayoutRow, 1,
                            1, 1);
  interfaceLayoutRow++;
  interfaceLayout.addWidget(&pathfindingValueLabel[END], interfaceLayoutRow, 0,
                            1, 1);
  interfaceLayout.addWidget(&pathfindingFlipSelectionComboBox[END],
                            interfaceLayoutRow, 1, 1, 1);
  interfaceLayoutRow++;
  interfaceLayout.addWidget(&pathFindingStatusLabel,
                            interfaceLayoutRow, 0, 1, 1);

  interfaceBox.setFixedWidth(370);
  interfaceBox.setLayout(&interfaceLayout);

  mainGroupBoxLayout->addWidget(tabWidget, 0, 0, 10, 2);
  mainGroupBoxLayout->addWidget(realizationViewWidgetCornerVersion, 0, 3, 4, 1);
  mainGroupBoxLayout->addWidget(&interfaceBox, 6, 3, 3, 1);
  mainGroupBoxLayout->addWidget(&stopStartSampleButton, 4, 3, 1, 1);
  mainGroupBoxLayout->addWidget(&runModeSelectionComboBox, 5, 3, 1, 1);
  mainGroupBoxLayout->setSpacing(1);
  mainLayout->setSpacing(1);
  mainGroupBox->setLayout(mainGroupBoxLayout);

  mainLayout->addWidget(menuBar);
  mainLayout->addWidget(mainGroupBox);

  mainWidget->setLayout(mainLayout);
  setCentralWidget(mainWidget);
}

void MainWindow::init_DistanceFinder_from_settings(PredefinedInteractions *df) {
  Settings *sett = Settings::getInstance();
  vector<vector<string> > data = readData(sett->DistanceData.file);

  if (data.size() == 0) return;
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

  df->assign(output);
  delete output;
}

void MainWindow::keyPressEvent(QKeyEvent *event) {
  if (event->key() == Qt::Key_F1) {
    std::cout << "F1" << std::endl;
    stackedWidget->setCurrentWidget(atlasViewWidget);
  }
  if (event->key() == Qt::Key_F2) {
    std::cout << "F2" << std::endl;
    stackedWidget->setCurrentWidget(nodeParentChildernViewWidget);

    nodeParentChildernViewWidget->setAtlasNode(
        atlasViewWidget->getCurrentNode());
  }
  if (event->key() == Qt::Key_F3) {
    std::cout << "F3" << std::endl;
    stackedWidget->setCurrentWidget(sweepViewWidget);

    sweepViewWidget->setAtlasNode(atlasViewWidget->getCurrentNode());
  }
  if (event->key() == Qt::Key_F4) {
    std::cout << "F4" << std::endl;
    stackedWidget->setCurrentWidget(cayleySpaceViewWidget);
    cayleySpaceViewWidget->setAtlasNode(atlasViewWidget->getCurrentNode());
  }
  if (event->key() == Qt::Key_F5) {
    std::cout << "F5" << std::endl;
    stackedWidget->setCurrentWidget(realizationViewWidget);
  }
}

void MainWindow::setDependenices(int argc, char **argv, MolecularUnit *muA,
                                 MolecularUnit *muB, SaveLoader *save_loader,
                                 AtlasBuilder *atlas_builder, Atlas *atlas_view,
                                 std::thread *backEndThread) {
  this->argc = argc;
  this->argv = argv;
  this->muA = muA;
  this->muB = muB;
  this->saveLoader = save_loader;
  this->atlas_builder = atlas_builder;
  this->backEndThread = &backEndThread;
  this->atlas = atlas_view;
}
