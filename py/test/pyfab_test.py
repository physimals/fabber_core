import unittest
import os
import sys
import numpy as np

d = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
sys.path.insert(0, os.path.join(d))
from fabber import FabberLib, FabberExec, FabberRunData

# Misc numbers for generating test data
c0=0.4
c1=0.8
c2=-0.3

class TestFabberLib(unittest.TestCase):

    def setUp(self):
        self.fab = FabberLib()
        self.progress = []

    def quad_data(self, x, y, z, t):
        t=t+1
        return c0 + c1*t + c2*t*t

    def test_get_no_lib(self):
        raised=False
        try:
            FabberLib(fabber_lib="junk")
        except:
            raised=True
        self.assertTrue(raised)

    def test_get_no_modellib(self):
        raised=False
        try:
            FabberLib(models_lib="junk")
        except:
            raised=True
        self.assertTrue(raised)

    def test_get_models(self):
        models = self.fab.get_models()
        self.assertTrue("poly" in models)

    def test_get_methods(self):
        methods = self.fab.get_methods()
        self.assertTrue("vb" in methods)
        self.assertTrue("nlls" in methods)
        self.assertTrue("spatialvb" in methods)

    def test_get_options(self):
        opts, desc = self.fab.get_options()
        opts = [opt["name"] for opt in opts]
        self.assertTrue("save-model-fit" in opts)

    def test_get_options_vb(self):
        opts, desc = self.fab.get_options(method="vb")
        opts = [opt["name"] for opt in opts]
        self.assertTrue("PSP_byname<n>" in opts)

    def test_get_options_poly(self):
        opts, desc = self.fab.get_options(model="poly")
        opts = [opt["name"] for opt in opts]
        self.assertTrue("degree" in opts)

    def test_get_model_params(self):
        rundata = FabberRunData()
        rundata["model"] = "poly"
        rundata["degree"] = "2"
        params = self.fab.get_model_params(rundata)
        self.assertTrue("c0" in params)
        self.assertTrue("c1" in params)
        self.assertTrue("c2" in params)

    def test_model_evaluate(self):
        rundata = FabberRunData()
        rundata["model"] = "poly"
        rundata["degree"] = "2"
        ret = self.fab.model_evaluate(rundata, {"c0":1, "c1":2, "c2":3}, 5)
        self.assertEquals(5, len(ret))
        for t in range(5):
            self.assertEquals(1 + 2*(t+1) + 3*(t+1)*(t+1), ret[t])

    def test_run_no_mask(self):
        data = np.fromfunction(self.quad_data, (3,3,3,3))
        rundata = FabberRunData()
        rundata["model"] = "poly"
        rundata["degree"] = "2"
        rundata["save-mean"] = ""
        run = self.fab.run_with_data(rundata, {"data" : data})
        self.assertAlmostEqual(run.data["mean_c0"][0,0,0], c0, delta=0.1)
        self.assertAlmostEqual(run.data["mean_c1"][0,0,0], c1, delta=0.1)
        self.assertAlmostEqual(run.data["mean_c2"][0,0,0], c2, delta=0.1)

    def progress_cb(self, voxel, nvoxels):
        self.progress.append((voxel, nvoxels))

    def test_run_cb(self):
        data = np.fromfunction(self.quad_data, (3,3,3,3))
        rundata = FabberRunData()
        rundata["model"] = "poly"
        rundata["degree"] = "2"
        rundata["save-mean"] = ""
        run = self.fab.run_with_data(rundata, {"data" : data}, progress_cb=self.progress_cb)

        self.assertEquals(0, self.progress[0][0])
        last = self.progress[-1]
        self.assertTrue(last[0] > 0)
        self.assertEquals(last[0], last[1])

    def test_run_mask(self):
        data = np.fromfunction(self.quad_data, (3,3,3,3))
        mask = np.zeros((3, 3, 3))
        mask[0,1,2] = 1
        rundata = FabberRunData()
        rundata["model"] = "poly"
        rundata["degree"] = "2"
        rundata["save-mean"] = ""
        run = self.fab.run_with_data(rundata, {"data" : data}, mask)
        self.assertEqual(run.data["mean_c0"][0,0,0], 0)
        self.assertEqual(run.data["mean_c1"][0,0,0], 0)
        self.assertEqual(run.data["mean_c2"][0,0,0], 0)
        self.assertAlmostEqual(run.data["mean_c0"][0,1,2], c0, delta=0.1)
        self.assertAlmostEqual(run.data["mean_c1"][0,1,2], c1, delta=0.1)
        self.assertAlmostEqual(run.data["mean_c2"][0,1,2], c2, delta=0.1)

    def test_run_empty_mask(self):
        data = np.fromfunction(self.quad_data, (3,3,3,3))
        mask = np.zeros((3, 3, 3))
        rundata = FabberRunData()
        rundata["model"] = "poly"
        rundata["degree"] = "2"
        rundata["save-mean"] = ""
        run = self.fab.run_with_data(rundata, {"data" : data}, mask)
        self.assertEqual(run.data["mean_c0"].max(), 0)
        self.assertEqual(run.data["mean_c1"].max(), 0)
        self.assertEqual(run.data["mean_c2"].max(), 0)

class TestFabberExec(unittest.TestCase):

    def setUp(self):
        # FIXME should do this properly using tempfiles containing generated
        self.datadir = os.path.abspath("%s/../../test/" % os.path.dirname(__file__))
        #print(self.datadir)
        self.fab = FabberExec()

    def quad_data(self, x, y, z, t):
        t=t+1
        return c0 + c1*t + c2*t*t

    def test_get_no_lib(self):
        raised=False
        try:
            FabberExec(ex="junk")
        except:
            raised=True
        self.assertTrue(raised)

    def test_get_no_modellib(self):
        raised=False
        try:
            FabberExec(models_lib="junk")
        except:
            raised=True
        self.assertTrue(raised)

    def test_get_models(self):
        models = self.fab.get_models()
        self.assertTrue("poly" in models)

    def test_get_methods(self):
        methods = self.fab.get_methods()
        self.assertTrue("vb" in methods)
        self.assertTrue("nlls" in methods)
        self.assertTrue("spatialvb" in methods)

    def test_get_options(self):
        opts, desc = self.fab.get_options()
        opts = [opt["name"] for opt in opts]
        self.assertTrue("save-model-fit" in opts)

    def test_get_options_vb(self):
        opts, desc = self.fab.get_options(method="vb")
        opts = [opt["name"] for opt in opts]
        self.assertTrue("PSP_byname<n>" in opts)

    def test_get_options_poly(self):
        opts, desc = self.fab.get_options(model="poly")
        opts = [opt["name"] for opt in opts]
        self.assertTrue("degree" in opts)

    # def test_get_model_params(self):
    #     rundata = FabberRunData()
    #     rundata["model"] = "poly"
    #     rundata["degree"] = "2"
    #     params = self.fab.get_model_params(rundata)
    #     self.assertTrue("c0" in params)
    #     self.assertTrue("c1" in params)
    #     self.assertTrue("c2" in params)

    def test_run_no_mask(self):
        rundata = FabberRunData()
        rundata["model"] = "poly"
        rundata["degree"] = "2"
        # rundata["save-mean"] = "" FIXME causes error due to compatibility options
        rundata["data"] = os.path.join(self.datadir, "test_data_small")
        run = self.fab.run(rundata)
        #print(run.log)
        #self.assertAlmostEqual(run.data["mean_c0"][0,0,0], c0, delta=0.1)
        #self.assertAlmostEqual(run.data["mean_c1"][0,0,0], c1, delta=0.1)
        #self.assertAlmostEqual(run.data["mean_c2"][0,0,0], c2, delta=0.1)

    def test_run_mask(self):
        rundata = FabberRunData()
        rundata["model"] = "poly"
        rundata["degree"] = "2"
        # rundata["save-mean"] = "" FIXME causes error due to compatibility options
        rundata["data"] = os.path.join(self.datadir, "test_data_small")
        rundata["mask"] = os.path.join(self.datadir, "test_mask_small")
        run = self.fab.run(rundata)
        #self.assertAlmostEqual(run.data["mean_c0"][0,0,0], c0, delta=0.1)
        #self.assertAlmostEqual(run.data["mean_c1"][0,0,0], c1, delta=0.1)
        #self.assertAlmostEqual(run.data["mean_c2"][0,0,0], c2, delta=0.1)

if __name__ == '__main__':
    unittest.main()
