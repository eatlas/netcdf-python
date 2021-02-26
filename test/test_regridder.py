import regrid.builder as builder
import regrid.regridder as Regridder
import math
import numpy as np
import unittest


class TestRegridder(unittest.TestCase):

    def test_calculate_weighted_mean(self):
        values = [1.0, 2.0, 3.0, 4.0]
        weight = [3.0, 2.0, 1.0, math.nan]
        expected_result = (values[0] * weight[0] + values[1] * weight[1] + values[2] * weight[2]) / \
                          (weight[0] + weight[1] + weight[2])
        self.assertEqual(Regridder._calculate_weighted_mean(values, weight), expected_result)

    def test_regridded_1_depth(self):
        """
        Test the Regridder with an input grid that has no depth dimension.
        """
        # Define test data.
        input_latitude_array = np.array([1.0, 2.0, 3.0])
        input_longitude_array = np.array([4.0, 5.0, 6.0])
        salt = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]])

        # Retrieve and verify the properties of the Regridder.
        regridder = builder.build_from_input_grid(input_latitude_array, input_longitude_array, resolution=0.5)
        self.assertEqual(regridder.latitude_count, 4)
        self.assertIn(1.0, regridder.latitude_array)
        self.assertIn(1.5, regridder.latitude_array)
        self.assertIn(2.0, regridder.latitude_array)
        self.assertIn(2.5, regridder.latitude_array)
        self.assertEqual(regridder.longitude_count, 4)
        self.assertIn(4.0, regridder.longitude_array)
        self.assertIn(4.5, regridder.longitude_array)
        self.assertIn(5.0, regridder.longitude_array)
        self.assertIn(5.5, regridder.longitude_array)

        output = regridder.regrid(salt)
