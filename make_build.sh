#!/bin/bash

sudo rm -rf build
mkdir build
cp compile.sh build
cd build/
chmod +x compile.sh
./compile.sh
