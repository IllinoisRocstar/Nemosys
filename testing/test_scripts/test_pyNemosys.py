import unittest
import json
import os
import sys
import tempfile
from shutil import copy2
import argparse
from inspect import currentframe, getframeinfo

topsrcdir = "@CMAKE_CURRENT_BINARY_DIR@"

# TODO: Replace the tests so they don't rely on JSON

class TestPyNemosys(unittest.TestCase):

    # def setUp(self):
    #     self.CMAKE_CURRENT_SOURCE_DIR = '/Nemosys/testing/'
    #     self.datadir_path = os.path.join(self.CMAKE_CURRENT_SOURCE_DIR, 'test_data/test_pyNemosys/')

    def testMeshBase(self):
        frameinfo = getframeinfo(currentframe())
        print(frameinfo.filename + '-' + str(frameinfo.lineno))
        from pyNemosys import meshBase, diffMesh, TransferDriver, TransferBase

        # tmpdirname = tempfile.mkdtemp()
        # os.chdir(tmpdirname)
        # path = topsrcdir + '/test_data/test_pyNemosys/transfer/'
        # for f in os.listdir(path):
        #     copy2(os.path.join(path, f), tmpdirname)
        path = topsrcdir + '/test_data/test_pyNemosys/transfer/'
        os.chdir(path)

        with open('transfer.json', 'r') as jsonfile:
            inputjson = json.load(jsonfile)

        source_file = inputjson['Mesh File Options']['Source Mesh File']
        target_file = inputjson['Mesh File Options']['Target Mesh File']
        method_name = inputjson['Transfer Options']['Method']
        array_names = [str(i) for i in inputjson['Transfer Options']['Array Names']]
        output_file = inputjson['Mesh File Options']['Output Mesh File']
        source = meshBase.Create(source_file)
        target = meshBase.Create(target_file)
        print(source, target, method_name, array_names)
        for name in array_names:
            print(name)
        transfer = TransferDriver.CreateTransferObject(source, target, method_name)
        transfer.transferPointData(source.getArrayIDs(array_names))
        target.write('transfer_test.vtu')

        gold_output_file = 'gold_transfer_test.vtu'

        self.assertEqual(diffMesh(meshBase.Create(gold_output_file), meshBase.Create(output_file)), 0)
        self.assertEqual(diffMesh(meshBase.Create(gold_output_file), target), 0)

    def testTransferDriver(self):
        frameinfo = getframeinfo(currentframe())
        print(frameinfo.filename + '-' + str(frameinfo.lineno))
        from pyNemosys import meshBase, TransferDriver, diffMesh

        # tmpdirname = tempfile.mkdtemp()
        # os.chdir(tmpdirname)
        # path = topsrcdir + '/test_data/test_pyNemosys/transfer/'
        # for f in os.listdir(path):
        #     copy2(os.path.join(path, f), tmpdirname)
        path = topsrcdir + '/test_data/test_pyNemosys/transfer/'
        os.chdir(path)

        with open('transfer.json', 'r') as jsonfile:
            inputjson = json.load(jsonfile)

        source_file = inputjson['Mesh File Options']['Source Mesh File']
        target_file = inputjson['Mesh File Options']['Target Mesh File']
        method_name = inputjson['Transfer Options']['Method']
        array_names = [str(i) for i in inputjson['Transfer Options']['Array Names']]
        output_file = inputjson['Mesh File Options']['Output Mesh File']
        check_quality = inputjson['Transfer Options']['Check Transfer Quality'] in ['true', 'True']
        gold_output_file = 'gold_transfer_test.vtu'

        opts = TransferDriver.Opts(method_name, check_quality)
        opts.arrayNames = array_names

        driver = TransferDriver(TransferDriver.Files(source_file, target_file, output_file), opts)
        driver.execute()
        self.assertEqual(diffMesh(meshBase.Create(gold_output_file), meshBase.Create(output_file)), 0)

    def testRefineDriver_value(self):
        frameinfo = getframeinfo(currentframe())
        print(frameinfo.filename + '-' + str(frameinfo.lineno))
        from pyNemosys import SizeFieldRefineDriver, diffMesh, meshBase

        # tmpdirname = tempfile.mkdtemp()
        # os.chdir(tmpdirname)
        # path = topsrcdir + '/test_data/RefinementTest/'
        # for f in os.listdir(path):
        #     copy2(os.path.join(path, f), tmpdirname)
        path = topsrcdir + '/test_data/RefinementTest/'
        os.chdir(path)

        with open('refine_value.json', 'r') as jsonfile:
            inputjson = json.load(jsonfile)

        _mesh = inputjson["Mesh File Options"]["Input Mesh File"]
        ofname = inputjson["Mesh File Options"]["Output Mesh File"]
        method = inputjson["Refinement Options"]["Refinement Method"]
        arrayName = inputjson["Refinement Options"]["Array Name"]
        dev_mult = inputjson["Refinement Options"]["StdDev Multiplier"]
        maxIsmin = inputjson["Refinement Options"]["Max Is Min for Scaling"] in [1, 'true', 'True']
        transferData = inputjson["Refinement Options"]["Transfer Data"] in [1, 'true', 'True']
        print(maxIsmin)

        gold_output_file = 'gold_refined_beam_value.vtu'

        opts = SizeFieldRefineDriver.Opts(SizeFieldRefineDriver.Opts.Method.GRADIENT
                                          if method == "gradient"
                                          else SizeFieldRefineDriver.Opts.Method.VALUE,
            arrayName, dev_mult, maxIsmin, transferData)

        refdrvobj = SizeFieldRefineDriver(SizeFieldRefineDriver.Files(_mesh, ofname), opts)
        refdrvobj.execute()
        # self.assertEqual(diffMesh(meshBase.Create(gold_output_file), meshBase.Create(ofname)), 0)

    def testRefineDriver_uniform(self):
        frameinfo = getframeinfo(currentframe())
        print(frameinfo.filename + '-' + str(frameinfo.lineno))
        from pyNemosys import UniformRefineDriver, diffMesh, meshBase

        # tmpdirname = tempfile.mkdtemp()
        # os.chdir(tmpdirname)
        # path = topsrcdir + '/test_data/RefinementTest/'
        # for f in os.listdir(path):
        #     copy2(os.path.join(path, f), tmpdirname)
        path = topsrcdir + '/test_data/RefinementTest/'
        os.chdir(path)

        with open('refine_uniform.json', 'r') as jsonfile:
            inputjson = json.load(jsonfile)

        _mesh = inputjson["Mesh File Options"]["Input Mesh File"]
        ofname = inputjson["Mesh File Options"]["Output Mesh File"]
        method = inputjson["Refinement Options"]["Refinement Method"]
        edgescale = inputjson["Refinement Options"]["Edge Scaling"]
        transferData = bool(inputjson["Refinement Options"]["Transfer Data"])

        gold_output_file = 'gold_refined_beam_uniform.vtu'

        refdrvobj = UniformRefineDriver(UniformRefineDriver.Files(_mesh, ofname),
                                        UniformRefineDriver.Opts(transferData, edgescale))
        refdrvobj.execute()
        self.assertEqual(diffMesh(meshBase.Create(gold_output_file), meshBase.Create(ofname)), 0)

    def testMeshGenDriver(self):
        frameinfo = getframeinfo(currentframe())
        print(frameinfo.filename + '-' + str(frameinfo.lineno))
        import pyNemosys
        if not hasattr(pyNemosys, "NetgenMeshGenDriver"):
            print("Nothing to test. Netgen not enabled")
            return
        from pyNemosys import NetgenMeshGenDriver, netgenParams, diffMesh, meshBase

        # tmpdirname = tempfile.mkdtemp()
        # os.chdir(tmpdirname)
        # path = topsrcdir + '/test_data/test_pyNemosys/meshGen/'
        # for f in os.listdir(path):
        #     copy2(os.path.join(path, f), tmpdirname)
        path = topsrcdir + '/test_data/test_pyNemosys/meshGen/'
        os.chdir(path)

        with open('meshGen.json', 'r') as jsonfile:
            inputjson = json.load(jsonfile)

        in_file = inputjson["Mesh File Options"]["Input Geometry File"]
        output_file = inputjson['Mesh File Options']['Output Mesh File']
        gold_output_file = 'gold_hinge.vtu'

        params_json = inputjson["Mesh Generation Options"]["Netgen Parameters"]
        params = netgenParams()
        # We can do this because they happen to be named the same
        for param in ["uselocalh", "maxh", "fineness", "grading", "elementsperedge", "elementspercurve",
                      "closeedgeenable", "closeedgefact", "second_order", "meshsize_filename", "quad_dominated",
                      "optvolmeshenable", "optsteps_2d", "optsteps_3d", "invert_tets", "invert_trigs", "check_overlap",
                      "check_overlapping_boundary", "refine_with_geom", "refine_without_geom"]:
            if param in params_json:
                setattr(params, param, params_json[param])

        mshgndrvobj = NetgenMeshGenDriver(NetgenMeshGenDriver.Files(in_file, output_file), params)
        mshgndrvobj.execute()

        refMesh = meshBase.Create(gold_output_file)
        newMesh = meshBase.Create(output_file)

        # self.assertEqual(diffMesh(meshBase.Create(gold_output_file), meshBase.Create(output_file)), 0)

        # Due to Netgen algorithm changing slightly between different versions,
        # we are not going to compare the meshes point-by-point, cell-by-cell.
        # The test will check the number of cells and points and ensure they
        # are within a 0.5% tolerance.

        if os.name == 'nt':
            divisor = 39  # 2.56% = 1 / 39
        else:
            divisor = 200  # 0.5% = 1 / 200

        refpoints = refMesh.getNumberOfPoints()
        refcells = refMesh.getNumberOfCells()
        newpoints = newMesh.getNumberOfPoints()
        newcells = newMesh.getNumberOfCells()

        self.assertEqual(
            refpoints - refpoints / divisor <= newpoints <= refpoints + refpoints / divisor
            and refcells - refcells / divisor <= newcells <= refcells + refcells / divisor,
            True)

    def testMeshQualityDriver(self):
        frameinfo = getframeinfo(currentframe())
        print(frameinfo.filename + '-' + str(frameinfo.lineno))
        from pyNemosys import CheckMeshQualDriver, meshBase

        # tmpdirname = tempfile.mkdtemp()
        # os.chdir(tmpdirname)
        # path = topsrcdir + '/test_data/test_pyNemosys/meshQuality/'
        # for f in os.listdir(path):
        #     copy2(os.path.join(path, f), tmpdirname)
        path = topsrcdir + '/test_data/test_pyNemosys/meshQuality/'
        os.chdir(path)

        with open('MeshQuality.json', 'r') as jsonfile:
            inputjson = json.load(jsonfile)

        _mesh = inputjson["Mesh File Options"]["Input Mesh File"]
        ofname = inputjson["Mesh File Options"]["Output File"]
        gold_output_file = 'gold_meshQual.txt'

        qualdrvobj = CheckMeshQualDriver(CheckMeshQualDriver.Files(_mesh, ofname))
        qualdrvobj.execute()

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
    parser = argparse.ArgumentParser()
    parser.add_argument('unittest_args', nargs='*')

    args = parser.parse_args()

    # Now set the sys.argv to the unittest_args (leaving sys.argv[0] alone)
    sys.argv[1:] = args.unittest_args
    unittest.main()
