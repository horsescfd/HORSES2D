#!/bin/bash
FILES=$(find . -name "SETUP")
script_path=$(pwd)
sed "1s/.*/HORSES_DIR=$(echo $PWD/Solver | sed 's_/_\\/_g')/" ./TestCases/Makefile.problemFile > ./TestCases/Makefile
for i in $FILES; do
	FILE_PATH=${i%\/SETUP}	
	FILE_PATH=${FILE_PATH%\/SETUP/}	
	rm $FILE_PATH/HORSES2D.NS
	rm $FILE_PATH/HORSES2D.Euler
	rm $FILE_PATH/horses2d.ns
	rm $FILE_PATH/horses2d.euler
	rm "$FILE_PATH/horses2d (Case Conflict"*
	ln -v -s $script_path/Solver/bin/HORSES2D.NS $FILE_PATH
	ln -v -s $script_path/Solver/bin/HORSES2D.Euler $FILE_PATH
	cp -v $script_path/TestCases/Makefile $FILE_PATH/SETUP/Makefile.template
done

rm ./TestCases/Makefile
