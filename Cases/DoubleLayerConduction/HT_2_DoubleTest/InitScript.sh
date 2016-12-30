#!/bin/bash
unset FOAM_SIGFPE

blockMesh

mkdir -p 0
cp -r A/* 0/

setFields

interThermalFoam2

