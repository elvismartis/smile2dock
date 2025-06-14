#!/bin/bash
echo "This is a test script to be used as tutorial"
echo "Convert multiple smiles to 3D molecular formats ready for docking"
python3 ../../smile2dock.py -i test.smi -o ./ > test.log
