# Install script for directory: /home/muskaan/easal-dev/include/Eigen/src/Eigen/Eigen

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/home/muskaan/easal-dev/include/Eigen")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "0")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xDevelx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/eigen3/Eigen" TYPE FILE FILES
    "/home/muskaan/easal-dev/include/Eigen/src/Eigen/Eigen/Cholesky"
    "/home/muskaan/easal-dev/include/Eigen/src/Eigen/Eigen/CholmodSupport"
    "/home/muskaan/easal-dev/include/Eigen/src/Eigen/Eigen/Core"
    "/home/muskaan/easal-dev/include/Eigen/src/Eigen/Eigen/Dense"
    "/home/muskaan/easal-dev/include/Eigen/src/Eigen/Eigen/Eigen"
    "/home/muskaan/easal-dev/include/Eigen/src/Eigen/Eigen/Eigenvalues"
    "/home/muskaan/easal-dev/include/Eigen/src/Eigen/Eigen/Geometry"
    "/home/muskaan/easal-dev/include/Eigen/src/Eigen/Eigen/Householder"
    "/home/muskaan/easal-dev/include/Eigen/src/Eigen/Eigen/IterativeLinearSolvers"
    "/home/muskaan/easal-dev/include/Eigen/src/Eigen/Eigen/Jacobi"
    "/home/muskaan/easal-dev/include/Eigen/src/Eigen/Eigen/LU"
    "/home/muskaan/easal-dev/include/Eigen/src/Eigen/Eigen/MetisSupport"
    "/home/muskaan/easal-dev/include/Eigen/src/Eigen/Eigen/OrderingMethods"
    "/home/muskaan/easal-dev/include/Eigen/src/Eigen/Eigen/PaStiXSupport"
    "/home/muskaan/easal-dev/include/Eigen/src/Eigen/Eigen/PardisoSupport"
    "/home/muskaan/easal-dev/include/Eigen/src/Eigen/Eigen/QR"
    "/home/muskaan/easal-dev/include/Eigen/src/Eigen/Eigen/QtAlignedMalloc"
    "/home/muskaan/easal-dev/include/Eigen/src/Eigen/Eigen/SPQRSupport"
    "/home/muskaan/easal-dev/include/Eigen/src/Eigen/Eigen/SVD"
    "/home/muskaan/easal-dev/include/Eigen/src/Eigen/Eigen/Sparse"
    "/home/muskaan/easal-dev/include/Eigen/src/Eigen/Eigen/SparseCholesky"
    "/home/muskaan/easal-dev/include/Eigen/src/Eigen/Eigen/SparseCore"
    "/home/muskaan/easal-dev/include/Eigen/src/Eigen/Eigen/SparseLU"
    "/home/muskaan/easal-dev/include/Eigen/src/Eigen/Eigen/SparseQR"
    "/home/muskaan/easal-dev/include/Eigen/src/Eigen/Eigen/StdDeque"
    "/home/muskaan/easal-dev/include/Eigen/src/Eigen/Eigen/StdList"
    "/home/muskaan/easal-dev/include/Eigen/src/Eigen/Eigen/StdVector"
    "/home/muskaan/easal-dev/include/Eigen/src/Eigen/Eigen/SuperLUSupport"
    "/home/muskaan/easal-dev/include/Eigen/src/Eigen/Eigen/UmfPackSupport"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xDevelx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/eigen3/Eigen" TYPE DIRECTORY FILES "/home/muskaan/easal-dev/include/Eigen/src/Eigen/Eigen/src" FILES_MATCHING REGEX "/[^/]*\\.h$")
endif()

