#!/bin/bash

shopt -s globstar
for file in **/*.ipynb ; do
	echo "clearing output cells from $file"
	jupyter nbconvert --clear-output --inplace $file
done
