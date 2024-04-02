#!/bin/bash

####################################################
#
#  Script updates all existing packages in given 
#  working copy of capd library
#  Invoke from root directory of capd library
#
####################################################

echo "[CAPD core package]"
svn update
for dir in capdDynSys capdDynSys4 capdExtHom capdRedHom; do    
    if [ -e $dir ]; then  
        echo "[package $dir]"
        svn update $dir;
#    else
# echo "skipping $dir"
    fi
done

