import regrid.builder as builder
import math
import numpy as np
import unittest


class TestRegridBuilder(unittest.TestCase):

    def test_InputGridProxy_1d(self):
        """
        Test the InputGridProxy class with a 1-dimensional input array.
        """

        # Define test data.
        input_grid_proxy = builder.InputGridProxy()
        input_latitude_array = np.array([1.0, 2.0, 3.0])
        input_longitude_array = np.array([4.0, 5.0, 6.0])
        input_grid_proxy.populate(input_latitude_array, input_longitude_array)

        # Perform tests.
        self.assertEqual(input_grid_proxy.latitude_min, 1.0)
        self.assertEqual(input_grid_proxy.latitude_max, 3.0)
        self.assertEqual(input_grid_proxy.longitude_min, 4.0)
        self.assertEqual(input_grid_proxy.longitude_max, 6.0)

        self.assertEqual(len(input_grid_proxy.index_to_location_map),
                         len(input_latitude_array) * len(input_longitude_array))
        self.assertEqual(input_grid_proxy.index_to_location_map[0], (1.0, 4.0))
        self.assertEqual(input_grid_proxy.index_to_location_map[2], (1.0, 6.0))
        self.assertEqual(input_grid_proxy.index_to_location_map[4], (2.0, 5.0))
        self.assertEqual(input_grid_proxy.index_to_location_map[6], (3.0, 4.0))
        self.assertEqual(input_grid_proxy.index_to_location_map[8], (3.0, 6.0))

    def test_InputGridProxy_2d(self):
        """
        Test the InputGridProxy class with a 2-dimensional input array.
        """

        # Definte test data.
        input_grid_proxy = builder.InputGridProxy()
        input_latitude_array = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]])
        input_longitude_array = np.array([[8.0, 7.0, 6.0], [5.0, 4.0, 3.0], [2.0, 1.0, 0.0]])
        input_grid_proxy.populate(input_latitude_array, input_longitude_array)

        # Perform tests.
        self.assertEqual(input_grid_proxy.latitude_min, 1.0)
        self.assertEqual(input_grid_proxy.latitude_max, 9.0)
        self.assertEqual(input_grid_proxy.longitude_min, 0.0)
        self.assertEqual(input_grid_proxy.longitude_max, 8.0)

        self.assertEqual(len(input_grid_proxy.index_to_location_map),
                         len(input_latitude_array) * len(input_longitude_array))
        self.assertEqual(input_grid_proxy.index_to_location_map[0], (1.0, 8.0))
        self.assertEqual(input_grid_proxy.index_to_location_map[2], (3.0, 6.0))
        self.assertEqual(input_grid_proxy.index_to_location_map[4], (5.0, 4.0))
        self.assertEqual(input_grid_proxy.index_to_location_map[6], (7.0, 2.0))
        self.assertEqual(input_grid_proxy.index_to_location_map[8], (9.0, 0.0))

    def test_calculate_distance(self):
        """
        Test the _calculate_distance function.
        """

        # Define test data
        input_grid_proxy = builder.InputGridProxy()
        input_latitude_array = np.array([1.0, 2.0, 3.0])
        input_longitude_array = np.array([4.0, 5.0, 6.0])
        input_grid_proxy.populate(input_latitude_array, input_longitude_array)

        # Perform tests.
        self.assertEqual(input_grid_proxy.calculate_distance(1, 4, 0, input_grid_proxy), 0.0)
        self.assertEqual(input_grid_proxy.calculate_distance(2, 4, 0, input_grid_proxy), 1.0)
        self.assertEqual(input_grid_proxy.calculate_distance(3, 6, 8, input_grid_proxy), 0.0)
        self.assertEqual(input_grid_proxy.calculate_distance(3, 5, 8, input_grid_proxy), 1.0)

    def test_find_closest_four(self):
        """
        Test the _find_closest_four function.
        """

        # Define test data
        input_grid_proxy = builder.InputGridProxy()
        input_latitude_array = np.array([1.0, 2.0, 3.00])
        input_longitude_array = np.array([4.0, 5.0, 6.0])
        input_grid_proxy.populate(input_latitude_array, input_longitude_array)

        # Perform tests.

        # Using the middle cell grid as the point of interest, find the closest four of all of the other cells. This
        # will be the middle cell, and 3 of the 4 closest cells that are equal distance.
        closest_four = builder._find_closest_four(2, 5, input_grid_proxy, [0, 1, 2, 3, 4, 5, 6, 7, 8])
        self.assertEqual(len(closest_four), 4)
        self.assertEqual(closest_four[0].index, 4)
        self.assertEqual(closest_four[0].distance, 0.0)
        self.assertIn(closest_four[1].index, [1, 3, 5, 7])
        self.assertEqual(closest_four[1].distance, 1.0)
        self.assertIn(closest_four[2].index, [1, 3, 5, 7])
        self.assertEqual(closest_four[2].distance, 1.0)
        self.assertIn(closest_four[3].index, [1, 3, 5, 7])
        self.assertEqual(closest_four[3].distance, 1.0)

        # Using the left middle cell grid as the point of interest, find the closest four of all of the other cells.
        closest_four = builder._find_closest_four(2, 4, input_grid_proxy, [0, 1, 2, 3, 4, 5, 6, 7, 8])
        self.assertEqual(len(closest_four), 4)
        self.assertEqual(closest_four[0].index, 3)
        self.assertEqual(closest_four[0].distance, 0.0)
        self.assertIn(closest_four[1].index, [0, 4, 6])
        self.assertEqual(closest_four[1].distance, 1.0)
        self.assertIn(closest_four[2].index, [0, 4, 6])
        self.assertEqual(closest_four[2].distance, 1.0)
        self.assertIn(closest_four[3].index, [0, 4, 6])
        self.assertEqual(closest_four[3].distance, 1.0)

    def test_InputGridProxy_geo_search(self):
        """
        Test the _find_within_box function.
        """

        # Define test data. Note that one latitude is an outlier.
        input_grid_proxy = builder.InputGridProxy()
        input_latitude_array = np.array([1.0, 2.0, 4.0])
        input_longitude_array = np.array([4.0, 5.0, 6.0])
        input_grid_proxy.populate(input_latitude_array, input_longitude_array)

        # Perform tests.

        # Too restricted box, expecting none to be found.
        found_indexes = input_grid_proxy.geo_search(50, 50, 0.5)
        self.assertEqual(len(found_indexes), 0)

        # Small box, expecting matching index.
        found_indexes = input_grid_proxy.geo_search(2, 5, 0.5)
        self.assertEqual(len(found_indexes), 1)
        self.assertEqual(found_indexes[0], 4)

        # Larger box, but small enough to ignore outlier latitude.
        found_indexes = input_grid_proxy.geo_search(2, 5, 1.1)
        self.assertEqual(len(found_indexes), 6)
        self.assertIn(found_indexes[0], [0, 1, 2, 3, 4, 5])
        self.assertIn(found_indexes[1], [0, 1, 2, 3, 4, 5])
        self.assertIn(found_indexes[2], [0, 1, 2, 3, 4, 5])
        self.assertIn(found_indexes[3], [0, 1, 2, 3, 4, 5])
        self.assertIn(found_indexes[4], [0, 1, 2, 3, 4, 5])
        self.assertIn(found_indexes[5], [0, 1, 2, 3, 4, 5])

        # Large box that covers all coordinates.
        found_indexes = input_grid_proxy.geo_search(2, 5, 10)
        self.assertEqual(len(found_indexes), 9)
        self.assertIn(found_indexes[0], [0, 1, 2, 3, 4, 5, 6, 7, 8])
        self.assertIn(found_indexes[1], [0, 1, 2, 3, 4, 5, 6, 7, 8])
        self.assertIn(found_indexes[2], [0, 1, 2, 3, 4, 5, 6, 7, 8])
        self.assertIn(found_indexes[3], [0, 1, 2, 3, 4, 5, 6, 7, 8])
        self.assertIn(found_indexes[4], [0, 1, 2, 3, 4, 5, 6, 7, 8])
        self.assertIn(found_indexes[5], [0, 1, 2, 3, 4, 5, 6, 7, 8])
        self.assertIn(found_indexes[6], [0, 1, 2, 3, 4, 5, 6, 7, 8])
        self.assertIn(found_indexes[7], [0, 1, 2, 3, 4, 5, 6, 7, 8])
        self.assertIn(found_indexes[8], [0, 1, 2, 3, 4, 5, 6, 7, 8])

    def test_build_regular_output_grid(self):
        """
        Test the _build_regular_grid function.
        """

        # Define test data
        input_grid_proxy = builder.InputGridProxy()
        input_latitude_array = np.array([1.0, 4.0, 7.0])
        input_longitude_array = np.array([2.0, 5.0, 8.0])
        input_grid_proxy.populate(input_latitude_array, input_longitude_array)

        # Build the output grid.
        output_latitude_array, output_longitude_array = builder._build_regular_output_grid(
            input_grid_proxy=input_grid_proxy, resolution=2.0, number_of_decimal_places=2)
        self.assertEqual(len(output_latitude_array), 3)
        self.assertIn(1.0, output_latitude_array)
        self.assertIn(3.0, output_latitude_array)
        self.assertIn(5.0, output_latitude_array)
        self.assertNotIn(7.0, output_latitude_array)
        self.assertEqual(len(output_longitude_array), 3)
        self.assertIn(2.0, output_longitude_array)
        self.assertIn(4.0, output_longitude_array)
        self.assertIn(6.0, output_longitude_array)
        self.assertNotIn(8.0, output_longitude_array)

    def test_build_output_grid_index_to_closest_input_cells(self):
        """
        Test the _build_output_grid_index_to_closest_input_cells function.
        """

        # Define test data. Note that one latitude is an outlier.
        input_grid_proxy = builder.InputGridProxy()
        input_latitude_array = np.array([1.0, 2.0, 4.0])
        input_longitude_array = np.array([2.0, 3.0, 5.0])
        input_grid_proxy.populate(input_latitude_array, input_longitude_array)

        # Build the output grid.
        output_latitude_array, output_longitude_array = builder._build_regular_output_grid(
            input_grid_proxy=input_grid_proxy, resolution=1.0, number_of_decimal_places=2)

        # Build the mappings between output grid cells and closest input grid cells.
        output_grid_index_to_closest_input_cells = builder._build_output_grid_index_to_closest_input_cells(
            input_grid_proxy, output_latitude_array, output_longitude_array, search_box_increment_size=1,
            max_search_box_increments=3, min_found_cell_count=4)

        self.assertEqual(len(output_grid_index_to_closest_input_cells), 9)

        closest_input_cells = output_grid_index_to_closest_input_cells[0]
        self.assertEqual(len(closest_input_cells), 4)
        self.assertEqual(closest_input_cells[0].index, 0)
        self.assertEqual(closest_input_cells[0].distance, 0.0)
        self.assertIn(closest_input_cells[1].index, [1, 3])
        self.assertEqual(closest_input_cells[1].distance, 1.0)
        self.assertIn(closest_input_cells[2].index, [1, 3])
        self.assertEqual(closest_input_cells[2].distance, 1.0)
        self.assertEqual(closest_input_cells[3].index, 4)
        self.assertEqual(closest_input_cells[3].distance, math.sqrt(2))

        closest_input_cells = output_grid_index_to_closest_input_cells[8]
        self.assertEqual(len(closest_input_cells), 4)
        self.assertIn(closest_input_cells[0].index, [4, 5, 7, 8])
        self.assertEqual(closest_input_cells[0].distance, math.sqrt(2))
        self.assertIn(closest_input_cells[1].index, [4, 5, 7, 8])
        self.assertEqual(closest_input_cells[1].distance, math.sqrt(2))
        self.assertIn(closest_input_cells[2].index, [4, 5, 7, 8])
        self.assertEqual(closest_input_cells[2].distance, math.sqrt(2))
        self.assertIn(closest_input_cells[3].index, [4, 5, 7, 8])
        self.assertEqual(closest_input_cells[3].distance, math.sqrt(2))
