#!/bin/bash
FILES=$(find . -name "SETUP")
script_path=$(/bin/pwd)
cp ./TestCases/Makefile.problemFile ./TestCases/Makefile
for i in $FILES; do
	FILE_PATH=${i//\/SETUP/}	
	/bin/rm $FILE_PATH/HORSES2D.NS
	/bin/rm $FILE_PATH/HORSES2D.Euler
	/bin/ln -v -s $script_path/Solver/bin/HORSES2D.NS $FILE_PATH
	/bin/ln -v -s $script_path/Solver/bin/HORSES2D.Euler $FILE_PATH
	/usr/bin/sed -i .case "1s/.*/HORSES_DIR=$(echo $PWD/Solver | sed 's_/_\\/_g')/" ./TestCases/Makefile
	cp -v $script_path/TestCases/Makefile.case $FILE_PATH/SETUP/Makefile
done

/bin/rm ./TestCases/Makefile.case
/bin/rm ./TestCases/Makefile
