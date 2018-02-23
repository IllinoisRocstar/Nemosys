import unittest
import json
import os
import sys

class TestPyNemosys(unittest.TestCase):

    def testCreateMeshBase(self):
        sys.path.append('/Nemosys/python/')
        print('LD_LIBRARY_PATH=', os.environ['LD_LIBRARY_PATH'])
        from pyNemosys import meshBase, diffMesh

        datadir_path = '/Nemosys/testing/test_data/test_pyNemosys/'

        with open(os.path.join(datadir_path, 'transfer2.json'), 'r') as jsonfile:
            transfer2 = json.load(jsonfile)

        source = meshBase.Create(os.path.join(datadir_path, transfer2[0]['Mesh File Options']['Input Mesh Files']['Source Mesh']))
        target = meshBase.Create(os.path.join(datadir_path, transfer2[0]['Mesh File Options']['Input Mesh Files']['Target Mesh']))
        source.transfer(target, transfer2[0]['Transfer Options']['Method'], transfer2[0]['Transfer Options']['Array Names'])
        target.write('/Nemosys/testing/output/result.vtu')

        self.assertEqual(diffMesh(meshBase.Create(os.path.join(datadir_path, 'driver_test.vtu')), target), 0)

if __name__ == '__main__':
    unittest.main()
