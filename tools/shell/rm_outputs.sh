#!/usr/bin/env bash

shopt -s extglob
cd ./test
rm -rf !(inputs)
cd ../
