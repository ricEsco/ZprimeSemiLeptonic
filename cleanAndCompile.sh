#!/bin/bash

scram b clean
make clean
scram b -j4
make -j4
