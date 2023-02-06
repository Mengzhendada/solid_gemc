# solid_gemc

overview
--------
This repository is for [SoLID](https://solid.jlab.org) simulation based on [GEMC](https://gemc.jlab.org).

It includes all codes of defining and running simulation, digitizing and analyzing results.

It should be run with container from [solid_release](https://github.com/JeffersonLab/solid_release) which has prerequisites included. 

Read the gemc simulation general note https://solid.jlab.org/wiki/index.php/Gemc_simulation_general_note and solid software wiki https://solid.jlab.org/wiki/index.php/Software

code structure
--------
* "geometry", detector definition
* "field" field definition
* "script" example script to run simulation
* "digitization" code to digitize simulation output 
* "analysis" code to analyze results and do some studies
* "source" source code of solid_gemc which add solid specific hit process routine and load gemc as library
* "mod" some modification of original gemc code to fit solid need
* "prototype" code for solid prototype simulation

development rule and how to contribute
--------
* Use "master" branch as devel version and github release as production version
* most small development will be done by "fork->modify->submit pull request" mode.
* some major collaborative development can be done with "branch->modify->submit pull request" mode as needed 

To keep the fork in sync, do this on your computer
* git clone https://github.com/[your_github_name]/solid_gemc.git
* git remote add upstream https://github.com/JeffersonLab/solid_gemc.git
* git fetch upstream
* git checkout master
* git merge upstream/master
* git push origin masters
