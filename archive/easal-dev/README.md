```
    _________   _____ ___    __
   / ____/   | / ___//   |  / /
  / __/ / /| | \__ \/ /| | / /
 / /___/ ___ |___/ / ___ |/ /___
/_____/_/  |_/____/_/  |_/_____/


```

Description
===========
This software is the open-source project EASAL: Efficient Atlasing, Analysis and Search of Molecular Assembly Landscapes. For a detailed explanation of the theory behind the concepts in this program, refer to the papers referenced in Publications.

This project leverages the convexity of the Cayley parameter space for atlasing and sampling molecular assembly configuration spaces. With respect to many performance measures the method is provably superior to the typical Monte Carlo approach. In particular it achieves comparable coverage of the assembly configuration space far more efficiently and permits flexible tradeoff between accuracy and efficiency. Hybridizing monte carlo with EASAL promises the best of both worlds.

This project should be viewed primarily as mathematical software. For routine use by molecular scientists it would require significant improvement in its user friendliness.

A web version of the software is hosted at http://ufo-host.cise.ufl.edu/

Requirements
============
- Operating system: 
	- Ubuntu 16.04 or higher, 
	- Fedora 23 or higher, 
	- OSX Darwin or higher.
- C++ compiler
	- g++ version 7.3.x or higher
	- clang++ version 3.3 or higher
- CMAKE
	- cmake version 3.10.x or higher
- The following requirements are for the optional QT GUI included
	- OpenGL libraries
	- QT library (Version 5.9 or higher)
	- Nvidia Graphics card (with driver version-340 or above)

This software uses the QT framework for the GUI (https://www.qt.io/).

Building
========
- Prerequisites
    - Boost Libraries
        - $ sudo apt-get install libboost-all-dev
    - Install GLUT (Optional, required only for GUI)
	- $ sudo apt-get install freeglut3 freeglut3-dev
	- $ sudo apt-get install binutils-gold
    - Install QT (Optional, required only for GUI)
	- Follow the install instructions at the following link https://doc.qt.io/qt-5/gettingstarted.html
        - $ export CMAKE_PREFIX_PATH=<path_to_qt_installation> (example /usr/include/x86_64-linux-gnu/qt5).
    - Install GCC and CMAKE
        - $ sudo apt install build-essential
	- $ sudo apt install cmake
    - Install QT
        - $ sudo apt install qtbase5-dev qt5-qmake qtbase5-dev-tools qtchooser

- EASAL uses cmake to build its sources
	- $ mkdir build
	- $ cd build
	- $ cmake -DWITH_QT=ON . .
	- $ make
	
- Complete list of cmake flags include
	- WITH_QT - For the GUI
	- WITH_CAF (under development) - For parallel sampling
	- WITH_TEST (under development) - Using google test framework to test
	- WITH_MONGO (under development) - Using mongo db instead of files to store sampled data

- EASAL BackEnd (without GUI)
    - Follow the build process described above with the cmake flag -DWITH_QT set to OFF.
	- To run EASAL run ‘build/easal’ from the main directory.

- EASAL with QT GUI
	- Follow the build process described above with the cmake flag -DWITH_QT set to ON.
	- To run EASAL run ‘build/easal’ from the main directory.

Resources
=========
1. Input Files
	The 'files' directory in the root folder contains all the example input molecular data.
2. Include libraries
    the 'include' directory in the root folder contains all the required includes for the project. 

#Usage

- Command-line arguments
    - `-settings <settings_file>` Gives the settings file from which to load the input.
- Main View
    - Tabs: Allow users to switch between the different views listed below.
    - Right Panel
        - Sweep view of the currently selected node.
        - Sampling mode conrols
            - This has controls for stopping/re-starting sampling, choosing a different root node, sampling in a BFS fashion, and refining the sampling. 
        - Contact graph of the currently selected node, and sampling information.

##Tabs
- Atlas
    - This tab shows the main atlas.
    - Clicking on a node here displays the corresponding sweep view and the active constraint graphs in the right pane.
- Node Neighbors
    -  For a selected node, this tab shows its ancestors and descedants in the atlas.
- Sweep 
    - Shows the Cartesian Sweep of the currently selected node.
    - Controls at the bottom allow for selection of particular flips to see.
    - The user can step through all possible flips by pressing the up and down arrow keys.
    - Video controls at the bottom allow for animation of all the realizations and their flips.
- Cayley Space 
    - Shows the Cayley Space of the currently selected node.
    - For nodes of dimensions higher than 3, sliders at the bottom right corner allow users to step through each dimension.
    - Allows the user to inspect the boundaries of each node.
    - Controls at the bottom allow users to filter and inspect Cayley points with different properties.
- Realization
    -  Shows realizations in the currently selected nodes.
    -  Controls at the bottom allow for filtering realizations by flip numbers.
      
##Right Panel

- Constraint Selection Dialogue Box
    - Lets users choose a node based on the participating contacts.
    - If it doesn't already exist, it is created and sampled. Otherwise, the node is selected and 

##Menu
- File
    -Create
        - Opens the Input Window.
        - The fields are pre-populated with entries from the settings file. The user can edit these before starting sampling.
        - You will need to point the program to:
            - The two input point sets.
            - Distance data for the participating atoms.
            - A directory in which to store the atlas. If this directory already contains an atlas you will be prompted to either overwrite the old data or to continue where the other was left off. 

- Statistics
    - Basic Statistics
        - Atlas Statistics - Shows a dimensional breakdown of the number of nodes.
        - Node Statistics - Shows Cayley statistics for the selected node.
    - Grid Statistics
        - Grid Coverage
        - Epsilon Coverage

- Atlas
    - Basin
        -Basin Volume
        -Basin Regions
    - Pseudo Atlas - Generates the pseudo atlas and compares the number of nodes in the pseudo atlas to the atlas.
 
Example
=======
- Run EASAL from the command line by running the following command.
    `$build/easal' 
- Click on File->Create from the menu.
- In the input window, select the following either using the browse option or by entering the text in the text box provided
    - Data for Molecule A - files/A.pdb
    - Data for Molecule B - files/B.pdb
    - Distance Data - files/source_files/union computed desired distances.txt
    - Data Directory - data/

- The user can edit the data for Molecular Unit A, Molecular Unit B and Distance Data loaded from files by clicking on 'Set Data' and editing the values in the new pop-up that appears.

- Enter the values for Bonding Thresholds and step size . The user can either change the values for all these or just continue with the default values.

- The user can then click on the 'Advanced Options' button to set advanced user inputs. This opens a new pop-up where the user can enter either enter the data or choose to accept the default values.

- Click on 'Accept'. This opens the Atlas View.

- In the Atlas View, we initially see a root node at the center of a 3D grid. As and when more nodes are discovered, they are populated on the Atlas.

- In this view, the user can control how the sampling proceeds by using any of the controls on the left side of the atlas. The different controls available are
    - Stop Sampling.
    - Constraint Selection Dialogue Box.
    - Sample the Current Tree.
	- Sample all Incomplete Trees.
	- Auto-Solve.
    - BFS Sampling.
    - Refine Sampling.

- Clicking on a particular node does the following
    - Loads the Active Constraint Graph at the bottom left corner.
    - Loads the Cartesian realization for that node in the top right corner if the sampling for the node is complete.

- Pressing the space bar after selecting a node takes the user to the 'Cayley Space View' of that node.

- In the Cayley Space View, initially green points are shown which correspond to all the realizable points in the Cayley space.

- Clicking on the Red square at the bottom, shows the points that have collision. Clicking the red also shows two more options viz, cyan and pink. The Cyan points represent the points which have angle collision and the pink points represent points which have distance collision. 

- Clicking on the Blue square shows all the points sampled including, the good, the collision and the unrealizable.

- Clicking on the Boundaries shows the boundary points and the user can step through them along each dimension.

- Pressing the space bar here takes the user to the 'Realization View'

- Pressing 'v' on the keyboard generates the sweep view of the molecule.

- Once the sweep view is generated, the user can use the up and down arrow keys to view all the flips of the molecule.

- Clicking on boundaries shows the different realization along different boundaries.

- Clicking on the video controls at the bottom and clicking the play button on it animates and shows all the possible realizations of the molecule.

- Once sampling is completed, we find the number of paths between every pair of 0D and 1D nodes and write it to the path_matrix.txt file.

License
=========
EASAL is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. The GNU GPL license can be found at  <http://www.gnu.org/licenses/>.

Authors
=======
- James Pence
- Aysegul Ozkan
- Rahul Prabhu
- Troy Baker
- Amit Verma

Contact
=======
Meera Sitharam, CISE @ UF

Publications
============
 - Aysegul Ozkan and Meera Sitharam. 2011. ``EASAL: Efficient Atlasing, Analysis and Search of Molecular Assembly Landscapes''. In Proceedings of the ISCA 3rd International Conference on Bioinformatics and Computational Biology (BICoB-2011).
 - Aysegul Ozkan and Meera Sitharam. 2014. ``Best of Both Worlds: Uniform sampling in Cartesian and Cayley Molecular Assembly Configuration Space''. (2014). (on arxiv).
 - Aysegul Ozkan, Ruijin Wu, Jorg Peters, and Meera Sitharam. 2014b. ``Efficient Atlasing and Sampling of Assembly Free Energy Landscapes using EASAL: Stratification and Convexification via Customized Cayley Parametrization.'' (2014). (on arxiv).
 - Aysegul Ozkan, Rahul Prabhu, Troy Baker, James Pence, Jorg Peters, and Meera Sitharam ``Algorithm 990: Efficient Atlasing and Search of Configuration Spaces of Point-Sets Constrained by Distance Intervals'', ACM Transactions on Mathematical Software (TOMS), Volume 44 Issue 4, August 2018. 
Article No. 48 
