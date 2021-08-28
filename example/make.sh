#!/bin/bash
#changedir.sh
cd ..
make clean
make -j12
cd -
make clean 
make -j12