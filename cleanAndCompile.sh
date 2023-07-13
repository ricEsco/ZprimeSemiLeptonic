#!/bin/bash

scram b clean
make clean
scram b -j1
make -j1
