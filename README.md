# EASAL
Steps:

cd easal-dev

sudo dnf install boost-devel 

sudo dnf install freeglut 

sudo dnf install freeglut-devel 

sudo dnf install binutils-gold 

sudo dnf groupinstall "Development Tools" 

sudo dnf install qt5-qtbase-devel qt5-qdbusviewer qtchooser 

sudo dnf install eigen3-devel

sudo dnf install glog-devel

sudo dnf install qt-devel

sudo dnf upgrade cmake 

cd include
mkdir Eigen
mkdir include
cp -r /usr/include/eigen3 ./

cd ../../
mkdir glog
mkdir include
cp -r /usr/include/glog/* ./

cd ../../../

mkdir build

cd build

cmake -DWITH_QT=ON -DWITH_CAF=OFF ..

make

