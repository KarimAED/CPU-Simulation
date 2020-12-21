# -*- coding: utf-8 -*-
"""
Created on Sun Dec 20 09:37:43 2020

@author: Karim

This uses unittest to perform proper typing assertions
"""
import unittest


import system as sys

class MaterialTest(unittest.TestCase):
    """
    Test that all of the attributes must have the right types.
    """
    
    def test_name(self):
        with self.assertRaises(TypeError):
            sys.Material(1, 0, 0)
    
    def test_power(self):
        with self.assertRaises(TypeError):
            sys.Material("Air", "not_a_number", 0)
    
    def test_cond(self):
        with self.assertRaises(TypeError):
            sys.Material("Air", 0, "not_a_number")



class HeatSinkTest(unittest.TestCase):
    """
    Validate the width and height calculations in the heat sink initialiser.
    """
    
    def valid_width_height(self):
        l = 10
        dist = 2
        ct = 10
        bs_ht = 4
        
        hs = sys.HeatSink(l, dist, ct, bs_ht)
        
        self.assertEqual(hs.height, bs_ht+l)
        self.assertEqual(hs.width, (dist+1)*ct)


class SystemTest(unittest.TestCase):
    
    def min_res_check(self):
        """
        Confirm restriction to minimum resolution.
        """
        with self.assertRaises(ValueError):
            m1 = sys.Material("Air")
            m2 = sys.Material("Air")
            m3 = sys.Material("Air")
            mats = [m1, m2, m3] 
            sys.System("Test", mats, resolution=1)
    
    def valid_mat_check(self):
        """
        Confirm proper material array properties. (len, e.g.)
        """
        with self.assertRaises(ValueError):
            sys.System("Needs Materials", [])
        
        with self.assertRaises(ValueError):
            m1 = sys.Material("Air")
            m2 = sys.Material("Air")
            m3 = sys.Material("Air")
            m4 = sys.Material("Air")
            mats = [m1, m2, m3, m4] 
            sys.System("Too many Materials", mats)
        
        with self.assertRaises(ValueError):
            m1 = sys.Material("Air")
            m2 = sys.Material("Air")
            m3 = sys.Material("Air")
            mats = [m1, m2, m3]
            hs = sys.HeatSink()
            sys.System("Wrong no. Materials for HeatSink", mats, heat_sink=hs)
            

def suite():
    """
    Combines all of the tests into a test suite.

    Returns
    -------
    suite : unittest.TestSuite
        The test suite for the system tests.

    """
    suite = unittest.TestSuite()
    suite.addTest(MaterialTest("test_name"))
    suite.addTest(MaterialTest("test_power"))
    suite.addTest(MaterialTest("test_cond"))
    suite.addTest(HeatSinkTest("valid_width_height"))
    suite.addTest(SystemTest("min_res_check"))
    suite.addTest(SystemTest("valid_mat_check"))
    return suite

def check_mat_hs_sys():
    """
    The function for the system tests.

    """
    print("\n\nTest: Material, HeatSink and System dataclasses work properly")
    r = unittest.TextTestRunner()
    return not r.run(suite()).errors # return error free