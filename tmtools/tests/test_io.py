import unittest

import numpy.testing as nptest

from ..io import get_residue_data, get_structure
from ..testing import get_pdb_path, get_mmcif_path

from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue
from Bio.PDB.Atom import Atom


class TestIO(unittest.TestCase):
    def test_get_structure(self):
        # Given - an example structure
        protein_id = "2gtl"

        # When - load the example PDB file
        pdb_path = get_pdb_path(protein_id)
        pdb_structure = get_structure(pdb_path)
        self.assertEqual(pdb_structure.id, "2gtl")

        # When - load the example MMCIF file
        mmcif_path = get_mmcif_path(protein_id)
        mmcif_structure = get_structure(mmcif_path, format="mmcif")
        self.assertEqual(mmcif_structure.id, "2gtl")

        # Then - both files contain the same structure
        self.assertEqual(pdb_structure, mmcif_structure)

    def test_get_residue_coordinates(self):
        # Given
        pdb = get_pdb_path("2gtl")
        structure = get_structure(pdb)
        chain0 = next(structure.get_chains())

        # When
        coords, seq = get_residue_data(chain0)

        # Then
        self.assertEqual(coords.shape, (147, 3))
        nptest.assert_array_almost_equal(coords[0], [14.58, 114.133, 44.707], decimal=5)
        nptest.assert_array_almost_equal(
            coords[-1], [18.708, 134.427, 44.33], decimal=5
        )

        self.assertEqual(len(seq), coords.shape[0])
        self.assertEqual(seq[:5], "DCCSY")
        self.assertEqual(seq[-5:], "AKDLP")
    
    def test_get_residue_data(self):
        # Given
        chain = Chain('A')
        ca = Atom('CA', coord=np.array([1.0, 0.0, 0.0]), bfactor=1.0, occupancy=1.0, altloc='', fullname='CA', serial_number=1, element='C')
        
        r = Residue(id=(' ', 1, ' '), resname='VAL', segid=' ')
        r.add(ca)
        chain.add(r)

        r = Residue(id=(' ', 2, ' '), resname='PYL', segid=' ')
        r.add(ca)
        chain.add(r)

        r = Residue(id=(' ', 3, ' '), resname='SEC', segid=' ')
        r.add(ca)
        chain.add(r)
        
        r = Residue(id=(' ', 4, ' '), resname='UNK', segid=' ')
        r.add(ca)
        chain.add(r)

        # When 
        coords, seq = get_residue_data(chain)

        # Then
        self.assertEqual(seq, "VOUX")

        
