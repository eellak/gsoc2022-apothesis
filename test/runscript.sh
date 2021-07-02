#!/bin/bash

PATH=$(pwd)
REL_PATH=$PATH/$1
pushd $REL_PATH

~/Apothesis/src/build/Apothesis .
popd
