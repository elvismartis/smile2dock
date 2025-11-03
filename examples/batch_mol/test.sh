#!/bin/bash
echo "This is a test script to be used as tutorial"
echo "Convert a single smile to 3D molecular formats ready for docking"
python3 ../../smile2dock_v3.py -i test.smi -o ./ --reference "CCC(=O)CC(O)C" --protonate --ph_min 1.0 --ph_max 7.4 > test.log

