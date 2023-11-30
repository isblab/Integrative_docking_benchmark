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
#ifndef ATLASWINDOW_H
#define ATLASWINDOW_H

#include <easalcore/Atlas.h>
#include <easalcore/AtlasBuilder.h>
#include <easalcore/SaveLoader.h>
#include <easalcore/Settings.h>
#include <easalcore/readIn.h>

#include <QtWidgets>
#include <iostream>
#include <thread>
//#include <easalcore/Thread.h>
#include <easalcore/ThreadShare.h>
#include <easalcore/Thread_BackEnd.h>
#include <easalcore/PathFinder.h>

#include "AboutWindow.h"
#include "AtlasViewWidget.h"
#include "BondViewWidget.h"
#include "CayleySpaceViewWidget.h"
#include "CayleyStatisticsWindow.h"
#include "ConstraintSelectionWindow.h"
#include "InputWindow.h"
#include "MessageWindow.h"
#include "NodeParentChildernViewWidget.h"
#include "PathfindingViewWidget3DPart.h"
#include "PseudoAtlasWindow.h"
#include "RealisationViewWidget.h"
#include "RealizationViewWidgetCornerVersion.h"
#include "SharedDataGUI.h"
#include "StatsWindow.h"
#include "SweepViewWidget.h"

class InputWindow;

class MainWindow : public QMainWindow {
  Q_OBJECT
 public:
  MainWindow();

  //[ COMPONENTS]
  InputWindow *inputWindowInstance;
  StatsWindow *statsWindowInstance;
  PseudoAtlasWindow *pseudoAtlasWindowInstance;
  CayleyStatisticsWindow *cayleyStatsWindowInstance;
  AboutWindow *aboutWindowInstance;

  RealizationViewWidgetCornerVersion *realizationViewWidgetCornerVersion;
  AtlasViewWidget *atlasViewWidget;
  SweepViewWidget *sweepViewWidget;
  QWidget *mainWidget;
  BondViewWidget *bondViewWidget;
  NodeParentChildernViewWidget *nodeParentChildernViewWidget;
  CayleySpaceViewWidget *cayleySpaceViewWidget;
  QStackedWidget *stackedWidget;
  QTabWidget *tabWidget;
  RealisationViewWidget *realizationViewWidget;
  ConstraintSelectionWindow *constrainSelectionWindow;
  PathfindingViewWidget3DPart *pathfindingViewWidget;

  QLayout *mainLayout;
  QGridLayout *mainGroupBoxLayout;
  QGridLayout interfaceLayout;
  QWidget *mainGroupBox;
  QWidget interfaceBox;

  QLineEdit pathfindingLineEdit[2];
  QPushButton pathfindingButton[2];
  QLabel pathfindingValueLabel[2];
  QLabel pathfindingFlipLabel[2];
  QComboBox pathfindingFlipSelectionComboBox[2];
  QPushButton findPathBetweenD0NodesButton;
  QLabel pathFindingStatusLabel;

  SharedDataGUI sharedData;
  //[~COMPONENTS]
  //[ DEPENDENCIES]
  int argc;
  char **argv;
  MolecularUnit *muA;
  MolecularUnit *muB;
  SaveLoader *saveLoader;
  AtlasBuilder *atlas_builder;
  Atlas *atlas;
  std::thread **backEndThread;
  //[~DEPENDENCIES]
  //[ MENU]
  QMenuBar *menuBar;
  //[MENU ITEMS]
  //[FILE]
  QMenu *fileMenu;
  QAction *createAction;
  QAction *closeAction;

  //[EDIT]
  QMenu *editMenu;

  //[Statistics]
  QMenu *statsMenu;
  //[Statistics SUB MENU ITEMS]
  QMenu *basicStatisticsMenu;
  QAction *atlasStatsAction;
  QAction *nodeStatsAction;
  QAction *dimensionStatsAction;
  QAction *cayleyStatsAction;
  QMenu *gridStatisticsMenu;
  QAction *epsilonCoverageAction;
  QAction *gridCoverageAction;

  //[HELP]
  QMenu *helpMenu;
  QAction *paperAction;
  QAction *userGuideAction;
  QAction *aboutAction;

  //[ATLAS]
  QMenu *atlasMenu;
  QAction *paAction;
  QMenu *basinMenu;
  QAction *basinVolumeAction;
  QAction *basinRegionsAction;

  QPushButton stopStartSampleButton;
  QComboBox runModeSelectionComboBox;
  QLabel currentSampledNodeLabel;
  QLabel currentSelectedNodeLabel;
  QLabel currentRunModeLabel;
  //[~MENU]
  //[ TIMER]
  QTimer *timer;
  double updateTimeInMS;
  //[~TIMER]
  QString samplingModesNames[6] = {"Stopped",      "TreeSample",
                                   "ForestSample", "ForestSampleAS",
                                   "BreadthFirst", "RefineSampling"};
  void setDependenices(int argc, char **argv, MolecularUnit *muA,
                       MolecularUnit *muB, SaveLoader *saveLoader,
                       AtlasBuilder *atlas_builder, Atlas *atlas,
                       std::thread *backEndThread);

 private Q_SLOTS:
  void openInputWindow();
  void stopStartSamplingButtonSlot();
  void findPathBetweenD0NodesButtonSlot();
  void atlasStats();
  void cayleyStats();
  void paStats();
  void openPaperURL();
  void openUserGuide();
  void openAboutMenu();
  /**
   * @brief updateSamplingStatus updates currentSampledNodeLabel,
   * currentSelectedNodeLabel, currentRunModeLabel each updateTimeInMS
   */
  void updateSlot();
  void pathfindingStartNodeButtonSlot();
  void pathfindingEndNodeButtonSlot();

 private:
  bool isSamplingStopped = false;
  /**
   * @brief activeWidgetID: int that hold ID of current widget in stackedWidget
   */
  int activeWidgetID;
  /**
   * @brief processDataFromSettings: is called after dependencies are set to
   * start back end with correct data and settings
   */
  void processDataFromSettings();
  /**
   * @brief setLayouts: wrapper for code that determines how all widgets are
   * arranged
   */
  void setLayouts();
  /**
   * @brief init_DistanceFinder_from_settings: this function used to be in
   * main.cpp and when I was making QT GUI I was forced to move it here
   * @param df
   */
  void init_DistanceFinder_from_settings(PredefinedInteractions *df);
  /**
   * @brief keyPressEvent: Qt function that we override to handle keyborad input
   * @param event
   */
  void keyPressEvent(QKeyEvent *event) Q_DECL_OVERRIDE;

  //[ PATHFINDING VARS]
  int currentPathfindingSelectedAtlasNodeID[2];
  PointType currentPathfindingSelectedCayleyPointType[2];
  int currentPathfindingSelectedCayleyPointID[2];
  QString pathfindingButtonText[2];
  //[~PATHFINDING VARS]
};

#endif  // ATLASWINDOW_H
