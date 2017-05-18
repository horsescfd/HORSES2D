#!/bin/bash
FILES=$(find . -name "SETUP")
script_path=$(/bin/pwd)
for i in $FILES; do
	PATH=${i//\/SETUP/}	
	/bin/rm $PATH/HORSES2D.NS
	/bin/rm $PATH/HORSES2D.Euler
	/bin/ln -v -s $script_path/Solver/bin/HORSES2D.NS $PATH
	/bin/ln -v -s $script_path/Solver/bin/HORSES2D.Euler $PATH
done
