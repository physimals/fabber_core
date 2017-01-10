import unittest
import os
import sys
import numpy as np

sys.path.append(os.path.join(os.environ["FSLDIR"], "lib/python"))
from pyfab.fabber import FabberLib
from pyfab.model import FabberRunData

# Misc numbers for generating test data
c0=0.4
c1=0.8
c2=-0.3

class TestFabberLib(unittest.TestCase):

    def setUp(self):
        self.fab = FabberLib()

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

if __name__ == '__main__':
    unittest.main()
