#!/bin/bash
echo "This is a test script to be used as tutorial"
echo "Convert a single smile to 3D molecular formats ready for docking"
python3 ../../smile2dock_v2.py -i "COC1=C(C=C2C(=C1)CC(C2=O)CC3CCN(CC3)CC4=CC=CC=C4)OC" -o something --reference "CCC(=O)CC(O)C" --protonate --ph_min 7.4 --ph_max=7.4 > test.log

