#!/bin/bash
echo "This is a test script to be used as tutorial"
echo "Convert a single smile to 3D molecular formats ready for docking"
python3 ../../smile2dock.py -i "c1ccccc1" -o benzene > test.log
