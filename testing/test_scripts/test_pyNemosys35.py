import unittest
import json
import os
import sys
import tempfile
from shutil import copy2

class TestPyNemosys(unittest.TestCase):

    # def setUp(self):
    #     self.CMAKE_CURRENT_SOURCE_DIR = '/Nemosys/testing/'
    #     self.datadir_path = os.path.join(self.CMAKE_CURRENT_SOURCE_DIR, 'test_data/test_pyNemosys/')

    def testMeshBase(self):
        from pyNemosys import meshBase, diffMesh

        with tempfile.TemporaryDirectory() as tmpdirname:
            os.chdir(tmpdirname)
            path = '/Nemosys/testing/test_data/test_pyNemosys/transfer/'
            for f in os.listdir(path):
                copy2(os.path.join(path, f), tmpdirname)

            with open('transfer.json', 'r') as jsonfile:
                inputjson = json.load(jsonfile)

            source_file   = inputjson['Mesh File Options']['Input Mesh Files']['Source Mesh']
            target_file   = inputjson['Mesh File Options']['Input Mesh Files']['Target Mesh']
            method_name   = inputjson['Transfer Options']['Method']
            array_names   = inputjson['Transfer Options']['Array Names']
            output_file   = inputjson['Mesh File Options']['Output Mesh File']
            source = meshBase.Create(source_file)
            target = meshBase.Create(target_file)
            source.transfer(target, method_name, array_names)
            target.write('meshbase_test.vtu')

            gold_output_file = 'gold_transfer_test.vtu'

            self.assertEqual(diffMesh(meshBase.Create(gold_output_file), meshBase.Create(output_file)), 0)
            self.assertEqual(diffMesh(meshBase.Create(gold_output_file), target), 0)


    def testTransferDriver(self):
        from pyNemosys import meshBase, TransferDriver, diffMesh
        with tempfile.TemporaryDirectory() as tmpdirname:
            os.chdir(tmpdirname)
            path = '/Nemosys/testing/test_data/test_pyNemosys/transfer/'
            for f in os.listdir(path):
                copy2(os.path.join(path, f), tmpdirname)

            with open('transfer.json', 'r') as jsonfile:
                inputjson = json.load(jsonfile)

            source_file   = inputjson['Mesh File Options']['Input Mesh Files']['Source Mesh']
            target_file   = inputjson['Mesh File Options']['Input Mesh Files']['Target Mesh']
            method_name   = inputjson['Transfer Options']['Method']
            array_names   = inputjson['Transfer Options']['Array Names']
            output_file   = inputjson['Mesh File Options']['Output Mesh File']
            check_quality = inputjson['Transfer Options']['Check Transfer Quality'] in ['true', 'True']
            gold_output_file = 'gold_transfer_test.vtu'

            TransferDriver(source_file, target_file, method_name, array_names, output_file, check_quality)
            self.assertEqual(diffMesh(meshBase.Create(gold_output_file), meshBase.Create(output_file)), 0)


    def testTransferDriver_load_json_obj(self):
        from pyNemosys import meshBase, TransferDriver, diffMesh
        with tempfile.TemporaryDirectory() as tmpdirname:
            os.chdir(tmpdirname)
            path = '/Nemosys/testing/test_data/test_pyNemosys/transfer/'
            for f in os.listdir(path):
                copy2(os.path.join(path, f), tmpdirname)

            with open('transfer.json', 'r') as jsonfile:
                inputjson = json.load(jsonfile)

            output_file   = inputjson['Mesh File Options']['Output Mesh File']
            gold_output_file = 'gold_transfer_test.vtu'

            td = TransferDriver.readJSON(inputjson)
            self.assertEqual(diffMesh(meshBase.Create(gold_output_file), meshBase.Create(output_file)), 0)


    def testTransferDriver_load_json_filename(self):
        from pyNemosys import meshBase, TransferDriver, diffMesh
        with tempfile.TemporaryDirectory() as tmpdirname:
            os.chdir(tmpdirname)
            path = '/Nemosys/testing/test_data/test_pyNemosys/transfer/'
            for f in os.listdir(path):
                copy2(os.path.join(path, f), tmpdirname)

            output_file   = 'transfer_test.vtu'
            gold_output_file = 'gold_transfer_test.vtu'

            td = TransferDriver.readJSON('transfer.json')
            self.assertEqual(diffMesh(meshBase.Create(gold_output_file), meshBase.Create(output_file)), 0)


    def testNemDriver_readJSON(self):
        from pyNemosys import NemDriver
        #TODO: this


    def testRefineDriver_value(self):
        from pyNemosys import RefineDriver, diffMesh, meshBase
        with tempfile.TemporaryDirectory() as tmpdirname:
            os.chdir(tmpdirname)
            path = '/Nemosys/testing/test_data/test_pyNemosys/refine/'
            for f in os.listdir(path):
                copy2(os.path.join(path, f), tmpdirname)

            with open('refine_value.json', 'r') as jsonfile:
                inputjson = json.load(jsonfile)

            _mesh = inputjson["Mesh File Options"]["Input Mesh File"]
            ofname = inputjson["Mesh File Options"]["Output Mesh File"]
            method = inputjson["Refinement Options"]["Refinement Method"]
            arrayName = inputjson["Refinement Options"]["Array Name"]
            dev_mult = inputjson["Refinement Options"]["StdDev Multiplier"]
            maxIsmin = inputjson["Refinement Options"]["Max Is Min for Scaling"] in [1, 'true', 'True']
            edgescale = 0.0

            gold_output_file = 'gold_refined_beam_value.vtu'

            refdrvobj = RefineDriver(_mesh, method, arrayName, dev_mult, maxIsmin,
                                     edgescale, ofname)
            self.assertEqual(diffMesh(meshBase.Create(gold_output_file), meshBase.Create(ofname)), 0)


    def testRefineDriver_uniform(self):
        from pyNemosys import RefineDriver, diffMesh, meshBase
        with tempfile.TemporaryDirectory() as tmpdirname:
            os.chdir(tmpdirname)
            path = '/Nemosys/testing/test_data/test_pyNemosys/refine/'
            for f in os.listdir(path):
              copy2(os.path.join(path, f), tmpdirname)

            with open('refine_uniform.json', 'r') as jsonfile:
                inputjson = json.load(jsonfile)

            _mesh = inputjson["Mesh File Options"]["Input Mesh File"]
            ofname = inputjson["Mesh File Options"]["Output Mesh File"]
            method = inputjson["Refinement Options"]["Refinement Method"]
            edgescale = inputjson["Refinement Options"]["Edge Scaling"]

            gold_output_file = 'gold_refined_beam_uniform.vtu'

            refdrvobj = RefineDriver(_mesh, method, edgescale, ofname)
            self.assertEqual(diffMesh(meshBase.Create(gold_output_file), meshBase.Create(ofname)), 0)


    def testRefineDriver_readJSON_obj(self):
        from pyNemosys import RefineDriver, diffMesh, meshBase

        with tempfile.TemporaryDirectory() as tmpdirname:
            os.chdir(tmpdirname)
            path = '/Nemosys/testing/test_data/test_pyNemosys/refine/'
            for f in os.listdir(path):
              copy2(os.path.join(path, f), tmpdirname)

            with open('refine_value.json', 'r') as jsonfile:
                inputjson = json.load(jsonfile)

            output_file   = inputjson['Mesh File Options']['Output Mesh File']
            gold_output_file = 'gold_refined_beam_value.vtu'

            refdrvobj = RefineDriver.readJSON(inputjson)

            self.assertEqual(diffMesh(meshBase.Create(gold_output_file), meshBase.Create(output_file)), 0)


    def testRefineDriver_readJSON_filename(self):
        from pyNemosys import RefineDriver, diffMesh, meshBase

        with tempfile.TemporaryDirectory() as tmpdirname:
            os.chdir(tmpdirname)
            path = '/Nemosys/testing/test_data/test_pyNemosys/refine/'
            for f in os.listdir(path):
              copy2(os.path.join(path, f), tmpdirname)

            output_file   = 'refined_beam.vtu'
            gold_output_file = 'gold_refined_beam_value.vtu'

            refdrvobj = RefineDriver.readJSON('refine_value.json')

            self.assertEqual(diffMesh(meshBase.Create(gold_output_file), meshBase.Create(output_file)), 0)


    def testMeshGenDriver_readJSON_obj(self):
        from pyNemosys import MeshGenDriver, diffMesh, meshBase

        with tempfile.TemporaryDirectory() as tmpdirname:
            os.chdir(tmpdirname)
            path = '/Nemosys/testing/test_data/test_pyNemosys/meshGen/'
            for f in os.listdir(path):
              copy2(os.path.join(path, f), tmpdirname)

            with open('meshGen.json', 'r') as jsonfile:
                inputjson = json.load(jsonfile)

            output_file   = inputjson['Mesh File Options']['Output Mesh File']
            gold_output_file = 'gold_hinge.vtu'

            mshgndrvobj = MeshGenDriver.readJSON(inputjson)

            self.assertEqual(diffMesh(meshBase.Create(gold_output_file), meshBase.Create(output_file)), 0)


    def testMeshGenDriver_readJSON_filename(self):
        from pyNemosys import MeshGenDriver, diffMesh, meshBase

        with tempfile.TemporaryDirectory() as tmpdirname:
            os.chdir(tmpdirname)
            path = '/Nemosys/testing/test_data/test_pyNemosys/meshGen/'
            for f in os.listdir(path):
              copy2(os.path.join(path, f), tmpdirname)

            output_file   = 'hinge.vtu'
            gold_output_file = 'gold_hinge.vtu'

            mshgndrvobj = MeshGenDriver.readJSON('meshGen.json')

            self.assertEqual(diffMesh(meshBase.Create(gold_output_file), meshBase.Create(output_file)), 0)


    def testMeshQualityDriver(self):
        from pyNemosys import MeshQualityDriver, meshBase

        with tempfile.TemporaryDirectory() as tmpdirname:
            os.chdir(tmpdirname)
            path = '/Nemosys/testing/test_data/test_pyNemosys/meshQuality/'
            for f in os.listdir(path):
              copy2(os.path.join(path, f), tmpdirname)

            with open('MeshQuality.json', 'r') as jsonfile:
                inputjson = json.load(jsonfile)

            _mesh = inputjson["Input Mesh File"]
            ofname = inputjson["Output File"]
            gold_output_file = 'gold_meshQual.txt'

            qualdrvobj = MeshQualityDriver(_mesh, ofname)

            self.assertTrue(qualityfile_equals(ofname, gold_output_file, 1.0e-8))


def qualityfile_equals(filename1, filename2, epsilon):
    with open(filename1, 'r') as file1:
        with open(filename2, 'r') as file2:
            for lineNum in range(7):
                line1 = file1.readline().split()
                line2 = file2.readline().split()
                if lineNum > 2:
                    for valNum in range(6):
                        if valNum > 0:
                            if not abs(float(line1[valNum]) - float(line2[valNum])) <= epsilon:
                                return False

    return True





if __name__ == '__main__':
    unittest.main()
