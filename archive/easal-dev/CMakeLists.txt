cmake_minimum_required(VERSION 3.19.7)
#project(easal VERSION 1.0.0 LANGUAGES CXX)
project(easal LANGUAGES CXX)

# Set CXX standard to C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find includes in corresponding build directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)

# ON - Write sampling information to mongodb
# OFF - Write sampling information to files
option(WITH_MONGO "Build using MongoDB as backend" OFF)

# ON - Build the GUI code with the backend
# OFF - Build just the backend
option(WITH_QT "Build with QT frontend" OFF)

# ON - Build the parallel version which uses CAF
# OFF - Build the single threaded version
option(WITH_CAF "Build with CAF backend" OFF)

# ON - Build the easal test suite
# OFF - Do not build the easal test suite
option(WITH_TEST "Build the test code" OFF)

# ON - Print Debug messages
# OFF - Do not print Debug messages
option(WITH_VERBOSE "Build with VERBOSE" OFF)

# ON - Print Debug messages
# OFF - Do not print Debug messages
option(WITH_GPERF "Build with CPU profiling using Gperf" OFF)

# ON - Build the Jacobian sampling code
# OFF -Build the grid sampling code
option(WITH_JACOBIAN_SAMPLING "Build with Jacobian Sampling" OFF)

if(NOT CMAKE_BUILD_TYPE)
 set(CMAKE_BUILD_TYPE Release)
endif()
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

if(APPLE)
    set(LIB_EXTENSION dylib)
else()
    set(LIB_EXTENSION so) 
endif() 

if(WITH_JACOBIAN_SAMPLING)
	add_definition(-DJACOBIAN_SAMPLING)
endif()

if(WITH_VERBOSE)
    add_definitions(-DVERBOSE)
endif()

if(WITH_CAF)

    list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR})
    add_definitions(-DCAF)
    message("Building with CAF backend")

    find_package(CAF COMPONENTS core io REQUIRED)
    include_directories(${CAF_INCLUDE_DIRS})
endif()

if(WITH_MONGO)
    add_definitions(-DSERVER)
    message("Building with MongoDB as backend")
    
    find_package(PkgConfig REQUIRED)
    pkg_check_modules(JSONCPP jsoncpp)
    link_libraries(${JSONCPP_LIBRARIES})
    
    find_package(libbsoncxx REQUIRED)
    include_directories(${LIBBSONCXX_INCLUDE_DIRS})
    
    find_package(libmongocxx REQUIRED)
    include_directories(${LIBMONGOCXX_INCLUDE_DIRS})
endif()

# Currently builds both the sqlite version and the flat file version.
# Will be updated at a later date to write only sqlite, at which
# point, WITH_FLAT_FILE will be it's own target.
if(WITH_SQLITE)
  find_package (SQLite3)
  add_definitions(-DSQLITE)
  include_directories(${SQLite3_INCLUDE_DIRS})
endif()


#Make sure we always use pthreads
set(THREADS_PREFER_PTHREAD_FLAG ON)

# Find the QtWidgets library
if(WITH_QT)
    # Instruct CMake to run moc automatically when needed
    set(CMAKE_AUTOMOC ON)
    # Create code from a list of Qt designer ui files
    set(CMAKE_AUTOUIC ON)
    
    add_definitions(-DQT)
    find_package(Qt5 COMPONENTS Core Widgets REQUIRED )
endif()


#########Download and install external dependencies#################
include(ExternalProject)
# Eigen
set(EIGEN ${CMAKE_SOURCE_DIR}/include/Eigen)

ExternalProject_Add(
    Eigen
	PREFIX ${EIGEN}
	GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
	GIT_TAG 53a7864c48a2b2fc226c69e9184d7cc8e6b36b6d
	CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
)


# GFLAGS
set(GFLAGS ${CMAKE_SOURCE_DIR}/include/gflags)

ExternalProject_Add(
    gflags
    PREFIX ${GFLAGS}
    GIT_REPOSITORY https://github.com/gflags/gflags.git
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR> -DBUILD_SHARED_LIBS=ON
    BUILD_COMMAND make
    INSTALL_COMMAND make install
)

# GLOG
set(GLOG ${CMAKE_SOURCE_DIR}/include/glog)

ExternalProject_Add(
    glog
    PREFIX ${GLOG}
	GIT_REPOSITORY https://github.com/google/glog.git
	CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR> -DWITH_GTEST=OFF
    	BUILD_COMMAND make
    	INSTALL_COMMAND make install
)

# Find pthreads
find_package(Threads REQUIRED)

# Find Boost
find_package( Boost 1.40 COMPONENTS program_options REQUIRED )
include_directories(${Boost_INCLUDE_DIR})

# Add SimpleIni
include_directories(include)

# Find JsonCpp
#find_package(JsonCpp REQUIRED)

include_directories(${GLOG}/include)

include_directories(${GFLAGS}/include)

include_directories(${EIGEN}/include/eigen3)

include_directories(src)

add_executable(easal "")

add_subdirectory(src)

add_dependencies(easal glog gflags Eigen)

# Use the Widgets module from Qt 5
target_link_libraries(easal 
    LINK_PRIVATE
            ${Boost_LIBRARIES} 
            ${CMAKE_THREAD_LIBS_INIT} 
            Threads::Threads
	    ${GFLAGS}/lib/libgflags${CMAKE_SHARED_LIBRARY_SUFFIX}
	    ${GLOG}/lib64/libglog${CMAKE_SHARED_LIBRARY_SUFFIX})

if(WITH_QT)
    target_link_libraries(easal
        LINK_PRIVATE
		easalqtgui
            )
else()
    target_link_libraries(easal
        LINK_PRIVATE
            easalcore)
endif()

if(WITH_MONGO)
    target_link_libraries(easal
        LINK_PRIVATE
            ${LIBMONGOCXX_LIBRARIES}
            ${LIBBSONCXX_LIBRARIES} 
            ${JSONCPP_LIBRARIES})
endif()

if(WITH_CAF)
    target_link_libraries(easal
       LINK_PRIVATE
       ${CAF_LIBRARIES})
endif()

if(WITH_TEST)
	find_package(GTest REQUIRED)
	add_subdirectory(test)
endif()

if(WITH_SQLITE)
    target_link_libraries(easal
       LINK_PRIVATE
       ${SQLite3_LIBRARIES})
endif()


if(WITH_GPERF)
    
    list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR})
    add_definitions(-DGPERF)
    message("Building with Gperf profiling")

    find_package(Gperftools REQUIRED)
    include_directories(${GPERFTOOLS_INCLUDE_DIR})
    
    target_link_libraries(easal
       LINK_PRIVATE
       ${GPERFTOOLS_LIBRARIES})

endif()
