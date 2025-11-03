#!/bin/bash
echo "This is a test script to be used as tutorial"
echo "Convert multiple smiles to 3D molecular formats ready for docking"
echo "Calculate Tanimoto similarity with a reference structure represented as smiles"
python3 ../../smile2dock_v3.py -i test.smi -o ./ --reference "c1ccccc1" > test.log
