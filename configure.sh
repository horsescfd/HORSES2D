#!/bin/bash
#
#///////////////////////////////////////////////////////////////////////////////////////////////////////
#
#    HORSES2D - A high-order discontinuous Galerkin spectral element solver.
#    Copyright (C) 2017  Juan Manzanero Torrico (juan.manzanero@upm.es)
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#////////////////////////////////////////////////////////////////////////////////////////////////////////
#
echo '       //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// '
echo '      |                                  _    _    ____    _____     _____   ______    _____                                   |'
echo '      |                                 | |  | |  / __ \  |  __ \   / ____| |  ____|  / ____|                                  |'
echo '      |                                 | |__| | | |  | | | |__) | | (___   | |__    | (___                                    |'
echo '      |                                 |  __  | | |  | | |  _  /   \___ \  |  __|    \___ \                                   |'
echo '      |                                 | |  | | | |__| | | | \ \   ____) | | |____   ____) |  				       |'
echo '      |                                 |_|  |_|  \____/  |_|  \_\ |_____/  |______| |_____/                                   |'
echo '      |                                                                                                                        |'
echo '      |                                  +m:y:-`                                       `-:y:m+                                 |'
echo '      |                          .:::odmMMMNMms:-                                   -:smMNMMMmdo:::.                           |'
echo '      |                          `+mMMNMMMMNNNd+:                                   :+dNNNMMMMNMMm+`                           |'
echo '      |                        /shdNMNymMMMMshMmdNmd:                           :dmNdmMhsMMMMmyNMNdhs/                         |'
echo '      |                       ooosMMMNMMMNmNhNMMMNmhd-                         -dhmNMMMNhNmNMMMNMMMsooo                        |'
echo '      |                      `-/+yhMMmMMMhhMhss+/+o+.                           .+o+/+sshMhhMMMmMMhy+/-`                       |'
echo '      |                          yMNMhMMMMNMs                                           sMNMMMMhMNMy                           |'
echo '      |                        :ydmMMhMMMMy+M/                                         /M+yMMMMhMMmdy:                         |'
echo '      |                       .+s+omMhMMMMMy+Ms`                                     `sM+yMMMMMhMmo+s+.                        |'
echo '      |                         +ymNMNmMMMMMmdMm+`                                 `+mMdmMMMMMmNMNmy+                          |'
echo '      |           :-           //oymMMmMMMMMMMmNMd-       /yys:       :syy/       -dMNmMMMMMMMmMMmyo//           -:            |'
echo '      |          - -y.          ./oymMMMMMMMMhddhNMNs..+yNMMMMm-     -mMMMMNy+..sNMNhddhMMMMMMMMmyo/.          .y- -           |'
echo '      |          +o.-N.       `-:+yhdNMMMMMsMMMMMMMMMMmMMMmd/msd     dsm/dmMMMmMMMMMMMMMMsMMMMMNdhy+:-`       .N-.o+           |'
echo '      |         : :NsMo           oyo+-hMMMMMMMMMMMMMM+NMd+`.+mM/   /Mm+.`+dMN+MMMMMMMMMMMMMMh-+oyo           oMsN: :          |'
echo '      |         /y.yMM:          `.   :hoMMMMMMMMMMMMMMmMNMMMMmsN` `NsmMMMMNMmMMMMMMMMMMMMMMoh:   .`          :MMy.y/          |'
echo '      |          oNMMh               oh+MMMMMMMMMMMMMMMMNdNNyNMdms smdMNyNNdNMMMMMMMMMMMMMMMM+ho               hMMNo           |'
echo '      |          oMMM-             /dyhMMMMMMMMMMMMMMMMNs/. `mmsMm mMsmm` ./sNMMMMMMMMMMMMMMMMhyd/             -MMMo           |'
echo '      |          dNym           -oNMMMMMMMMMMMMMMMMNoo/.    odh-mN Nm-hdo    ./ooNMMMMMMMMMMMMMMMMNo-           myNd           |'
echo '      |         `Mdsh        `+mMMMMMMMMMMMMMMMh+Mh-       -NM/ `. .` /MN-       -hM+hMMMMMMMMMMMMMMMm+`        hsdM`          |'
echo '      |         `NN:N.      oNMMMMMMNoMMMMMMMhoyh:       .oNMd`       `dMNo.       :hyohMMMMMMMoNMMMMMMNo      .N:NN`          |'
echo '      |          /Mm/do`  .mMMMMMMMMMMMMMMMmhyo.         NNM+           +MNN         .oyhmMMMMMMMMMMMMMMMm.  `od/mM/           |'
echo '      |           .smNMMdymMMMMMMMMNmMMMMMhy:            sy-             -ys            :yhMMMMMmNMMMMMMMMmydMMNms.            |'
echo '      |             -NhMMdMhmMMMMMMMyMMMhhMM+                                           +MMhhMMMyMMMMMMMmhMdMMhN-              |'
echo '      |             /yhNd-hNMMMMMMMMmymhmMMd                                             dMMmhmymMMMMMMMMNh-dNhy/              |'
echo '      |             :-sm/-`dMMMMMMMMhh/MMMm.                                             .mMMM/hhMMMMMMMMd`-/ms-:              |'
echo '      |               `yy  `oNMMMMMMMmsMMm.                                               .mMMsmMMMMMMMNo`  yy`                |'
echo '      |                 y`   `omMMMMMMsMy`    -/:/+++o.                       .o+++/:/-    `yMsMMMMMMmo`   `y                  |'
echo '      |                         sMMhmoNNhysyhmMMmNMMho`                       `ohMMNmMMmhysyhNNomhMMs                          |'
echo '      |                         oMMN:yMmyo/-.``                                       ``.-/oymMy:NMMo                          |'
echo '      |                        oNMN.                                                             .NMNo                         |'
echo '      |                        -sNmo.             High Order Spectral Element Solver            .omNs-                         |'
echo '      |                           /dyo.                                                       .oyd/                            |'
echo '      |                            `oMNd-      2D Compressible Navier-Stokes equations      -dNMo`                             |'
echo '      |                              :dMNy`                                               `yNMd:                               |'
echo '      |                               -ymms                                               smmy-                                |'
echo '      |                                                                                                                        |'
echo '       \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ '
echo ' '
echo ' '

# Look for folders named "SETUP"
FILES=$(find . -name "SETUP")

# Get the current script absolute path
script_path=$(pwd)

# Include the absolute path of the solver in the Makefile.template
sed "1s/.*/HORSES_DIR=$(echo $PWD/Solver | sed 's_/_\\/_g')/" ./Utils/ProblemFile/Makefile.problemFile > ./TestCases/Makefile

echo " 	This script performs the following:"
echo "	 -> Replaces old binary links with new ones."
echo "	 -> Copies the newest version of the Problem file f90 file"
echo "	 -> Copies the newest version of the Problem file Makefile"
echo "	 -> Generates the cases folder structure."
echo " 		"

# Loop in all folders to configure the environment
for i in $FILES; do

#	Get the path (removing the SETUP folder) in both LINUX and MacOS
	FILE_PATH=${i%\/SETUP}	
	FILE_PATH=${FILE_PATH%\/SETUP/}	

#	Print the status
	echo "      ** Configuring directory " $FILE_PATH

# 	Remove old binary links
	rm -f $FILE_PATH/HORSES2D.NS
	rm -f $FILE_PATH/HORSES2D.Euler
	rm -f $FILE_PATH/horses2d.ns
	rm -f $FILE_PATH/horses2d.euler
	rm -f "$FILE_PATH/horses2d (Case Conflict"*
	rm -f "$FILE_PATH/horses2d (case conflict"*
	rm -f "$FILE_PATH/horses2d ("*"conflicted copy"*
	rm -f "$FILE_PATH/HORSES2D ("*"conflicted copy"*

#	Generate new binary links
	ln -s $script_path/Solver/bin/HORSES2D.NS $FILE_PATH
	ln -s $script_path/Solver/bin/HORSES2D.Euler $FILE_PATH

#	Copy the Makefile.template
	cp $script_path/TestCases/Makefile $FILE_PATH/SETUP/Makefile.template

#	Copy the ProblemFileTemplate.f90
	cp $script_path/Solver/src/IO/ProblemFile.f90 $FILE_PATH/SETUP/ProblemFileTemplate.f90

# 	Generate the remaining folders
	mkdir -p $FILE_PATH/RESULTS
	mkdir -p $FILE_PATH/MESH
done

rm ./TestCases/Makefile
