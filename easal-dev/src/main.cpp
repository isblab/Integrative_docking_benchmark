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
#ifdef CAF
#include <chrono>

#include "caf/all.hpp"
#include "easalcore/Actors.h"
#include "easalcore/AtlasBuilderActor.h"
using namespace caf;
#endif

#include <easalcore/Atlas.h>
#include <easalcore/AtlasBuilder.h>
#include <easalcore/SaveLoader.h>
#include <easalcore/Settings.h>
#include <easalcore/Thread_BackEnd.h>
#include <easalcore/readIn.h>
#include <glog/logging.h>
#include <gflags/gflags.h>

//#include <absl/flags/flag.h>
//#include <absl/flags/parse.h>
//#include <absl/strings/string_view.h>

#ifdef QT
// gui related includes
#include <gui_Qt/MainWindow.h>

#include <QtWidgets>
#endif

// standard
#include <iostream>
#include <map>

#ifdef SERVER
#include <mongocxx/client.hpp>
#include <mongocxx/instance.hpp>
#include <mongocxx/uri.hpp>
#endif

#ifdef GPERF
#include <gperftools/profiler.h>
#endif

using namespace std;

//ABSL_FLAG(std::string, settings, "settings.ini", "Settings file");
DEFINE_string(settings_file, "settings.ini", "settings file");
DEFINE_int32(max_threads, -1, "No. of threads");
DEFINE_string(log_file, "", "Log File");
DEFINE_string(output_dir, "./", "Path where all the output is generated");


/*
 * Need to write a meaningful comment descibing what main does.
 */
int main(int argc, char **argv) {
  
  // Parse command line flags
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  if (FLAGS_log_file != "") {
    google::SetLogSymlink(google::GLOG_INFO, FLAGS_log_file.c_str()); 
  }
  if (FLAGS_log_dir == "") {
      FLAGS_log_dir = ".";
  }
  google::InitGoogleLogging(argv[0]);

#ifdef GPERF
  ProfilerStart("testmain.prof");
#endif

  SaveLoader *save_loader;
  AtlasBuilder *atlas_builder;
  Atlas *atlas = new Atlas();

#ifdef SERVER
  mongocxx::instance inst{};
  mongocxx::client *mongoDBClient = new mongocxx::client({mongocxx::uri{}});
#endif

  // load the settings
  LOG(INFO) << "Loading settings from \"" << FLAGS_settings_file << "\".";

  cout << "Loading settings from \"" << FLAGS_settings_file << "\"." << endl;
  Settings *sett = Settings::getInstance();
  sett->load(FLAGS_settings_file.c_str());
  if (FLAGS_output_dir.back() != '/') {
    FLAGS_output_dir += "/";
  }
  sett->Output.dataDirectory = FLAGS_output_dir + sett->Output.dataDirectory;
  sett->save(sett->Output.dataDirectory.c_str());
  LOG(INFO) << "Output.dataDirectory: " << sett->Output.dataDirectory;

#ifdef CAF
  init_global_meta_objects<caf::id_block::EASAL>();
  core::init_global_meta_objects();
  actor_system_config cfg;
  if (FLAGS_max_threads > 0) {
    cfg.set("scheduler.max-threads", FLAGS_max_threads);
  }
  actor_system system{cfg};
  LOG(INFO) << "Number of threads used: " << system.scheduler().num_workers();
  sett->ActorSystem.sys = &system;
#endif

  // SaveLoader
#ifdef SERVER
  save_loader =
      new SaveLoader(sett->Output.dataDirectory, sett->runTimeObjects.muA,
                     sett->runTimeObjects.muB, mongoDBClient);
  save_loader->setSessionID(sett->Output.sessionID);
#else
  LOG(INFO) << sett->Output.dataDirectory;
  save_loader =
      new SaveLoader(sett->Output.dataDirectory, sett->runTimeObjects.muA,
                     sett->runTimeObjects.muB);
#endif

  // TODO: This needs to be looked into.
  if (sett->Constraint.wholeCollision) {
    // reading the neighbour matrix
    ConstraintCheck::nei_matrix = Utils::getMatrixFromFileFromTo("nei.txt");
  } else {
    ConstraintCheck::nei_matrix = Utils::getIdentityMatrix();
  }
  sett->setSaveLoader(save_loader);

#ifndef QT  // Just run the EASAL backend
#ifdef CAF
  cout << "CAF Code is being called" << endl;
  sett->createAtlasBuilderActor(atlas);
  cout << "Finished Creating AB" << endl;

  // sett->ActorSystem.AB = sett->ActorSystem.sys->spawn(typedAtlasBuilder,
  // save_loader, atlas);
  anon_send(sett->ActorSystem.AB, Start_v);

  // std::this_thread::sleep_for(std::chrono::seconds(180));
#else
  atlas_builder =
      new AtlasBuilder(sett->runTimeObjects.muA, sett->runTimeObjects.muB,
                       save_loader, sett->runTimeObjects.df, atlas);
  if (sett->Sampling.runSample) {
    atlas_builder->setup();
    cout << "Thread_Main: AtlasBuilder Set up done." << endl;
  } else {
    save_loader->loadMapView(atlas);
    cout << "Thread_Main: Loads existing atlas." << endl;
  }
  cout << "Thread_Main: Calling atlas_builder->startAtlasBuilding()." << endl;
  atlas_builder->startAtlasBuilding();

  cout << "Thread_Main: Calling this->save_loader->saveRoadMap(this->atlas)."
       << endl;
  save_loader->saveRoadMap(atlas);
  cout << "Thread_Main: Finishes and Exits.." << endl;
#endif
#else  // Run EASAL with the QT GUI

  // QApplication is necessary to render all windows. They will be rendered only
  // after app.exec() call;
  QApplication app(argc, argv);
  //[ CONTEXT SELECTION] Mystical code below searches for supported OpenGl
  //context
  // We need to do this to render 3D stuff
  QSurfaceFormat fmt;
  fmt.setDepthBufferSize(24);
  fmt.setStencilBufferSize(8);
  if (QCoreApplication::arguments().contains(QStringLiteral("--multisample")))
    fmt.setSamples(4);

  // Request OpenGL 3.3 compatibility or OpenGL ES 3.0.
  if (QOpenGLContext::openGLModuleType() == QOpenGLContext::LibGL) {
    // Requesting 3.3 compatibility context
    fmt.setVersion(3, 3);
    fmt.setProfile(QSurfaceFormat::CompatibilityProfile);
#if defined(__APPLE__) && defined(__MACH__)
    fmt.setProfile(QSurfaceFormat::CoreProfile);
#endif
  } else {
    // Requesting 3.0 context
    fmt.setVersion(3, 0);
  }
  QSurfaceFormat::setDefaultFormat(fmt);
  //[ SETTING GUI ARGUMENTS]
  sett->Arg.argc = argc;
  sett->Arg.argv = argv;
  MainWindow *mainWindow;
  // muA = new MolecularUnit();
  // muB = new MolecularUnit();
  save_loader = new SaveLoader();
  mainWindow = new MainWindow();
#ifndef CAF
  atlas_builder = new AtlasBuilder();
  std::thread *backEndThread;
  //[~SETTING GUI ARGUMENTS]

  mainWindow->setDependenices(argc, argv, sett->runTimeObjects.muA,
                              sett->runTimeObjects.muB, save_loader,
                              atlas_builder, atlas, backEndThread);
#endif
  mainWindow->show();
  // exec() runs till GUI is closed
  app.exec();
#endif

#ifdef GPERF
  ProfilerStop();
#endif

  return 0;
}
