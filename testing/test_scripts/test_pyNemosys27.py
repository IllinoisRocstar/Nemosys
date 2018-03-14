import unittest
import json
import os
import sys
import tempfile
from shutil import copy2
import argparse
from inspect import currentframe, getframeinfo

global topsrcdir 

class TestPyNemosys(unittest.TestCase):

    # def setUp(self):
    #     self.CMAKE_CURRENT_SOURCE_DIR = '/Nemosys/testing/'
    #     self.datadir_path = os.path.join(self.CMAKE_CURRENT_SOURCE_DIR, 'test_data/test_pyNemosys/')

    def testMeshBase(self):
        frameinfo = getframeinfo(currentframe())
        print (str(frameinfo.filename)+'-'+str(frameinfo.lineno))
        from pyNemosys import meshBase, diffMesh
       # tmpdirname=tempfile.mkdtemp()
       # os.chdir(tmpdirname)
       # path = topsrcdir + '/test_data/test_pyNemosys/transfer/'
       # for f in os.listdir(path):
       #     copy2(os.path.join(path, f), tmpdirname)
        path =  topsrcdir + '/test_data/test_pyNemosys/transfer/'
        os.chdir(path)

        with open('transfer.json', 'r') as jsonfile:
            inputjson = json.load(jsonfile)

        source_file   = str(inputjson['Mesh File Options']['Input Mesh Files']['Source Mesh'])
        target_file   = str(inputjson['Mesh File Options']['Input Mesh Files']['Target Mesh'])
        method_name   = str(inputjson['Transfer Options']['Method'])
        array_names   = inputjson['Transfer Options']['Array Names']
        array_names = [str(i) for i in array_names]
        output_file   = str(inputjson['Mesh File Options']['Output Mesh File'])
        source = meshBase.Create(source_file)
        target = meshBase.Create(target_file)
        source.transfer(target, method_name, array_names)
        target.write('transfer_test.vtu')

        gold_output_file = 'gold_transfer_test.vtu'

        self.assertEqual(diffMesh(meshBase.Create(gold_output_file), meshBase.Create(output_file)), 0)
        self.assertEqual(diffMesh(meshBase.Create(gold_output_file), target), 0)


    def testTransferDriver(self):
        frameinfo = getframeinfo(currentframe())
        print (str(frameinfo.filename)+'-'+str(frameinfo.lineno))
        from pyNemosys import meshBase, TransferDriver, diffMesh
        #tmpdirname=tempfile.mkdtemp()
        #os.chdir(tmpdirname)
        #path = topsrcdir + '/test_data/test_pyNemosys/transfer/'
        #for f in os.listdir(path):
        #    copy2(os.path.join(path, f), tmpdirname)
        path =  topsrcdir + '/test_data/test_pyNemosys/transfer/'
        os.chdir(path)

        with open('transfer.json', 'r') as jsonfile:
            inputjson = json.load(jsonfile)

        source_file   = str(inputjson['Mesh File Options']['Input Mesh Files']['Source Mesh'])
        target_file   = str(inputjson['Mesh File Options']['Input Mesh Files']['Target Mesh'])
        method_name   = str(inputjson['Transfer Options']['Method'])
        array_names   = inputjson['Transfer Options']['Array Names']
        array_names = [str(i) for i in array_names]
        output_file   = str(inputjson['Mesh File Options']['Output Mesh File'])
        check_quality = inputjson['Transfer Options']['Check Transfer Quality'] in ['true', 'True']
        gold_output_file = 'gold_transfer_test.vtu'

        TransferDriver(source_file, target_file, method_name, array_names, output_file, check_quality)
        self.assertEqual(diffMesh(meshBase.Create(gold_output_file), meshBase.Create(output_file)), 0)


    def testTransferDriver_load_json_obj(self):
        frameinfo = getframeinfo(currentframe())
        print (str(frameinfo.filename)+'-'+str(frameinfo.lineno))
        from pyNemosys import meshBase, TransferDriver, diffMesh
        #tmpdirname = tempfile.mkdtemp()
        #os.chdir(tmpdirname)
        #path =  topsrcdir + '/test_data/test_pyNemosys/transfer/'
        #for f in os.listdir(path):
        #    copy2(os.path.join(path, f), tmpdirname)
        path =  topsrcdir + '/test_data/test_pyNemosys/transfer/'
        os.chdir(path)

        with open('transfer.json', 'r') as jsonfile:
            inputjson = json.load(jsonfile)

        output_file   = str(inputjson['Mesh File Options']['Output Mesh File'])
        gold_output_file = 'gold_transfer_test.vtu'

        td = TransferDriver.readJSON(inputjson)
        self.assertEqual(diffMesh(meshBase.Create(gold_output_file), meshBase.Create(output_file)), 0)


    def testTransferDriver_load_json_filename(self):
        frameinfo = getframeinfo(currentframe())
        print (str(frameinfo.filename)+'-'+str(frameinfo.lineno))
        from pyNemosys import meshBase, TransferDriver, diffMesh
        #tmpdirname = tempfile.mkdtemp()
        #os.chdir(tmpdirname)
        #path =  topsrcdir + '/test_data/test_pyNemosys/transfer/'
        #for f in os.listdir(path):
        #    copy2(os.path.join(path, f), tmpdirname)
        path =  topsrcdir + '/test_data/test_pyNemosys/transfer/'
        os.chdir(path)

        output_file   = 'transfer_test.vtu'
        gold_output_file = 'gold_transfer_test.vtu'

        td = TransferDriver.readJSON('transfer.json')
        self.assertEqual(diffMesh(meshBase.Create(gold_output_file), meshBase.Create(output_file)), 0)


    def testNemDriver_readJSON(self):
        frameinfo = getframeinfo(currentframe())
        print (str(frameinfo.filename)+'-'+str(frameinfo.lineno))
        from pyNemosys import NemDriver
        #TODO: this


    def testRefineDriver_value(self):
        frameinfo = getframeinfo(currentframe())
        print (str(frameinfo.filename)+'-'+str(frameinfo.lineno))
        from pyNemosys import RefineDriver, diffMesh, meshBase
        #tmpdirname = tempfile.mkdtemp()
        #os.chdir(tmpdirname)
        #path =  topsrcdir + '/test_data/test_pyNemosys/refine/'
        #for f in os.listdir(path):
        #    copy2(os.path.join(path, f), tmpdirname)
        path =  topsrcdir + '/test_data/test_pyNemosys/refine/'
        os.chdir(path)


        with open('refine_value.json', 'r') as jsonfile:
            inputjson = json.load(jsonfile)

        _mesh = str(inputjson["Mesh File Options"]["Input Mesh File"])
        ofname = str(inputjson["Mesh File Options"]["Output Mesh File"])
        method = str(inputjson["Refinement Options"]["Refinement Method"])
        arrayName = str(inputjson["Refinement Options"]["Array Name"])
        dev_mult = inputjson["Refinement Options"]["StdDev Multiplier"]
        maxIsmin = inputjson["Refinement Options"]["Max Is Min for Scaling"] in [1, 'true', 'True']
        transferData = inputjson["Refinement Options"]["Transfer Data"] in [1, 'true', 'True']
        print(maxIsmin)
        edgescale = 0.0

        gold_output_file = 'gold_refined_beam_value.vtu'

        refdrvobj = RefineDriver(_mesh, method, arrayName, dev_mult, maxIsmin,
                                 edgescale, ofname, transferData)
        #self.assertEqual(diffMesh(meshBase.Create(gold_output_file), meshBase.Create(ofname)), 0)


    def testRefineDriver_uniform(self):
        frameinfo = getframeinfo(currentframe())
        print (str(frameinfo.filename)+'-'+str(frameinfo.lineno))
        from pyNemosys import RefineDriver, diffMesh, meshBase
        #tmpdirname = tempfile.mkdtemp()
        #os.chdir(tmpdirname)
        #path =  topsrcdir + '/test_data/test_pyNemosys/refine/'
        #for f in os.listdir(path):
        #  copy2(os.path.join(path, f), tmpdirname)
        path =  topsrcdir + '/test_data/test_pyNemosys/refine/'
        os.chdir(path)

        with open('refine_uniform.json', 'r') as jsonfile:
            inputjson = json.load(jsonfile)

        _mesh = str(inputjson["Mesh File Options"]["Input Mesh File"])
        ofname = str(inputjson["Mesh File Options"]["Output Mesh File"])
        method = str(inputjson["Refinement Options"]["Refinement Method"])
        edgescale = inputjson["Refinement Options"]["Edge Scaling"]
        transferData = inputjson["Refinement Options"]["Transfer Data"] in [1, 'true', 'True']

        gold_output_file = 'gold_refined_beam_uniform.vtu'

        refdrvobj = RefineDriver(_mesh, method, edgescale, ofname, transferData)
        self.assertEqual(diffMesh(meshBase.Create(gold_output_file), meshBase.Create(ofname)), 0)


    def testRefineDriver_readJSON_obj(self):
        frameinfo = getframeinfo(currentframe())
        print (str(frameinfo.filename)+'-'+str(frameinfo.lineno))
        from pyNemosys import RefineDriver, diffMesh, meshBase

        #tmpdirname = tempfile.mkdtemp()
        #os.chdir(tmpdirname)
        #path =  topsrcdir + '/test_data/test_pyNemosys/refine/'
        #for f in os.listdir(path):
        #  copy2(os.path.join(path, f), tmpdirname)
        path =  topsrcdir + '/test_data/test_pyNemosys/refine/'
        os.chdir(path)

        with open('refine_value.json', 'r') as jsonfile:
            inputjson = json.load(jsonfile)

        output_file   = str(inputjson['Mesh File Options']['Output Mesh File'])
        gold_output_file = 'gold_refined_beam_value.vtu'

        refdrvobj = RefineDriver.readJSON(inputjson)
        #self.assertEqual(diffMesh(meshBase.Create(gold_output_file), meshBase.Create(output_file)), 0) 

    def testRefineDriver_readJSON_filename(self):
        frameinfo = getframeinfo(currentframe())
        print (str(frameinfo.filename)+'-'+str(frameinfo.lineno))
        from pyNemosys import RefineDriver, diffMesh, meshBase

        #tmpdirname = tempfile.mkdtemp()
        #os.chdir(tmpdirname)
        #path =  topsrcdir + '/test_data/test_pyNemosys/refine/'
        #for f in os.listdir(path):
        #  copy2(os.path.join(path, f), tmpdirname)
        path =  topsrcdir + '/test_data/test_pyNemosys/refine/'
        os.chdir(path)

        output_file   = 'refined_beam1.vtu'
        gold_output_file = 'gold_refined_beam_value.vtu'

        refdrvobj = RefineDriver.readJSON('refine_value.json')
        #self.assertEqual(diffMesh(meshBase.Create(gold_output_file), meshBase.Create(output_file)), 0)


    def testMeshGenDriver_readJSON_obj(self):
        frameinfo = getframeinfo(currentframe())
        print (str(frameinfo.filename)+'-'+str(frameinfo.lineno))
        from pyNemosys import MeshGenDriver, diffMesh, meshBase

        #tmpdirname = tempfile.mkdtemp()
        #os.chdir(tmpdirname)
        #path =  topsrcdir + '/test_data/test_pyNemosys/meshGen/'
        #for f in os.listdir(path):
        #  copy2(os.path.join(path, f), tmpdirname)
        path =  topsrcdir + '/test_data/test_pyNemosys/meshGen/'
        os.chdir(path)

        with open('meshGen.json', 'r') as jsonfile:
            inputjson = json.load(jsonfile)

        output_file   = str(inputjson['Mesh File Options']['Output Mesh File'])
        gold_output_file = 'gold_hinge.vtu'

        mshgndrvobj = MeshGenDriver.readJSON(inputjson)

        self.assertEqual(diffMesh(meshBase.Create(gold_output_file), meshBase.Create(output_file)), 0)


    def testMeshGenDriver_readJSON_filename(self):
        frameinfo = getframeinfo(currentframe())
        print (str(frameinfo.filename)+'-'+str(frameinfo.lineno))
        from pyNemosys import MeshGenDriver, diffMesh, meshBase

        #tmpdirname = tempfile.mkdtemp()
        #os.chdir(tmpdirname)
        #path =  topsrcdir + '/test_data/test_pyNemosys/meshGen/'
        #for f in os.listdir(path):
        #  copy2(os.path.join(path, f), tmpdirname)
        path =  topsrcdir + '/test_data/test_pyNemosys/meshGen/'
        os.chdir(path)

        output_file   = 'hinge.vtu'
        gold_output_file = 'gold_hinge.vtu'

        mshgndrvobj = MeshGenDriver.readJSON('meshGen.json')

        self.assertEqual(diffMesh(meshBase.Create(gold_output_file), meshBase.Create(output_file)), 0)


    def testMeshQualityDriver(self):
        frameinfo = getframeinfo(currentframe())
        print (str(frameinfo.filename)+'-'+str(frameinfo.lineno))
        from pyNemosys import MeshQualityDriver, meshBase

        #tmpdirname = tempfile.mkdtemp()
        #os.chdir(tmpdirname)
        #path =  topsrcdir + '/test_data/test_pyNemosys/meshQuality/'
        #for f in os.listdir(path):
        #  copy2(os.path.join(path, f), tmpdirname)
        path =  topsrcdir + '/test_data/test_pyNemosys/meshQuality/'
        os.chdir(path)

        with open('MeshQuality.json', 'r') as jsonfile:
            inputjson = json.load(jsonfile)

        _mesh = str(inputjson["Input Mesh File"])
        ofname = str(inputjson["Output File"])
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
    parser = argparse.ArgumentParser()
    parser.add_argument('--topsrcdir', default='.')
    parser.add_argument('unittest_args', nargs='*')

    args = parser.parse_args()
    global topsrcdir
    topsrcdir = args.topsrcdir

    # Now set the sys.argv to the unittest_args (leaving sys.argv[0] alone)
    sys.argv[1:] = args.unittest_args
    unittest.main()
