#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Test suite for SMILES-to-3D converters
"""

import unittest
import os
import sys
import logging
from rdkit import Chem
import tempfile
from typing import Dict, List, Optional

# Suppress non-critical logging during tests
logging.getLogger().setLevel(logging.ERROR)

class TestSMILES2Dock(unittest.TestCase):
    def setUp(self):
        """Set up test cases"""
        self.test_smiles = {
            'simple': 'CC(=O)O',  # Acetic acid
            'complex': 'CC1=C(C(=O)NC(=N1)C2=CC=CC=C2)C(=O)O',  # More complex molecule
            'invalid': 'XXX'  # Invalid SMILES
        }
        self.temp_dir = tempfile.mkdtemp()
        
    def test_v3_basic_conversion(self):
        """Test basic SMILES conversion in v3"""
        sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
        import smile2dock_v3 as sd3
        
        # Test valid SMILES
        mol, props, variants = sd3.smiles_to_3d(
            self.test_smiles['simple'],
            output_base=os.path.join(self.temp_dir, "test_simple")
        )
        self.assertIsNotNone(mol)
        self.assertIsNotNone(props)
        self.assertTrue(os.path.exists(os.path.join(self.temp_dir, "test_simple.pdb")))
        
        # Test invalid SMILES
        mol, props, variants = sd3.smiles_to_3d(
            self.test_smiles['invalid'],
            output_base=os.path.join(self.temp_dir, "test_invalid")
        )
        self.assertIsNone(mol)
        
    def test_v3_protonation(self):
        """Test protonation functionality in v3"""
        sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
        import smile2dock_v3 as sd3
        
        mol, props, variants = sd3.smiles_to_3d(
            self.test_smiles['simple'],
            output_base=os.path.join(self.temp_dir, "test_protonation"),
            protonate=True,
            ph_min=7.0,
            ph_max=7.4
        )
        self.assertIsNotNone(variants)
        if variants:  # Type check for variants
            self.assertGreater(len(variants), 0)
        
    def test_v3_gnn_basic(self):
        """Test basic functionality of GNN version"""
        sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
        import smile2dock_v3_GNNimplicitsolvent as sd3gnn
        
        # Test without GNN (should work even if GNN not available)
        mol, props, variants = sd3gnn.smiles_to_3d(
            self.test_smiles['simple'],
            output_base=os.path.join(self.temp_dir, "test_gnn_simple")
        )
        self.assertIsNotNone(mol)
        self.assertIsNotNone(props)
        
    @unittest.skipIf(not hasattr(sys.modules.get('smile2dock_v3_GNNimplicitsolvent', None), 'HAS_GNN'),
                     "GNNImplicitSolvent not available")
    def test_v3_gnn_solvation(self):
        """Test GNN-based solvation (only if GNNImplicitSolvent is available)"""
        sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
        import smile2dock_v3_GNNimplicitsolvent as sd3gnn
        
        mol, props, variants = sd3gnn.smiles_to_3d(
            self.test_smiles['simple'],
            output_base=os.path.join(self.temp_dir, "test_gnn_solvation"),
            use_gnn=True,
            solvent="water"
        )
        self.assertIsNotNone(mol)
        self.assertIsNotNone(props)
        
    def test_file_generation(self):
        """Test generation of all output file formats"""
        sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
        import smile2dock_v3 as sd3
        
        output_base = os.path.join(self.temp_dir, "test_formats")
        mol, props, variants = sd3.smiles_to_3d(
            self.test_smiles['simple'],
            output_base=output_base
        )
        
        for fmt in ['pdb', 'mol2', 'sdf', 'pdbqt']:
            self.assertTrue(
                os.path.exists(f"{output_base}.{fmt}"),
                f"Failed to generate {fmt} file"
            )
            
    def test_property_calculation(self):
        """Test molecular property calculations"""
        sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
        import smile2dock_v3 as sd3
        
        mol, props, variants = sd3.smiles_to_3d(self.test_smiles['simple'])
        
        required_properties = [
            "Molecular Weight",
            "Crippen_LogP",
            "H-Bond Donors",
            "H-Bond Acceptors",
            "TPSA",
            "Rotatable Bonds"
        ]
        
        # Check if props is not None before accessing
        if props:
            for prop in required_properties:
                self.assertIn(prop, props)
                self.assertIsInstance(props[prop], (int, float))
            
    def tearDown(self):
        """Clean up temporary files"""
        import shutil
        shutil.rmtree(self.temp_dir)

if __name__ == '__main__':
    unittest.main(verbosity=2)