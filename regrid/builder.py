from typing import List

import logging
import math
from collections import namedtuple
import numpy as np
from . import regridder
from multiprocessing import Process, Manager

#
# NOTES:
#
#   - Depending on the dataset, the latitude and longitude can be either a 1-dimensional array or a 2-dimensional
#     array. Generally a 1-dimensional array means the underlying grid is rectilinear in nature, and each value
#     (either latitude or longitude) is unique in the array. However, 2-dimensional arrays generally mean the underlying
#     grid is non-rectilinear, allowing the grid locations to be anything/anywhere, such as curvilinear.
#

# Define the logger for this file.
logger = logging.getLogger('regrid')

# Constant determining the number of decimal places to use for precision when converting from the original dimensions
# to longitude/latitude.
DEFAULT_NUMBER_OF_DECIMAL_PLACES = 6

DEFAULT_RESOLUTION = 0.03

# A NamedTuple that represents a coordinate (latitude/longitude).
Coordinate = namedtuple('Coordinate', ['latitude', 'longitude'])

# A NamedTuple that contains an index (location in an array) and a distance between this index point and some other
# known point.
IndexWithDistance = namedtuple('IndexWithDistance', [
    'index', 'distance'
])


# ----------------------------------------------------------------------------------------------------------------------

class InputGridProxy:
    """
    A proxy object for the input grid. The primary role of this object is to provide a mechanism to perform a geo search
    for input cells that can be used to calculate the value of an output cell during the regridding process.
    """

    index_to_location_map: List[Coordinate] = []
    """
    A map between the two possible representations of each input cell within the input grid. Each input cell has a geo
    location identified by the lat/lon, and an index within a 1-dimensional representation (flattened) of the input 
    grid. This allows the identification of the lat/lon from the index. This map is populated by the "populate()"
    function.
    """

    latitude_value_to_index_maps = []
    """
    Finding the input cells nearest an output cell requires a geo search (ie: use of lat/lon). This property maps every
    latitude value in the input grid to it's corresponding index(es) in the 1-dimensional representation (flattened) of
    the input grid. For example, latitude -21.3 can be found in cells 37, 93, 1345, 9593, 10583. This map is populated
    by the "populate()" method. In order to speed up searching, the mapping is split across smaller groups to reduce the
    number of compares.
    """

    longitude_value_to_index_maps = []
    """
    Finding the input cells nearest an output cell requires a geo search (ie: use of lat/lon). This property maps every
    longitude value in the input grid to it's corresponding index(es) in the 1-dimensional representation (flattened) of
    the input grid. For example, longitude -21.3 can be found in cells 37, 93, 1345, 9593, 10583. This map is populated
    by the "populate()" method. In order to speed up searching, the mapping is split across smaller groups to reduce the
    number of compares.
    """

    latitude_min = None
    latitude_max = None
    longitude_min = None
    longitude_max = None

    def __init__(self):

        # Initialise all class/instance variables.
        self.index_to_location_map = []
        self.latitude_value_to_index_maps = None
        self.longitude_value_to_index_maps = None
        self.latitude_min = None
        self.latitude_max = None
        self.longitude_min = None
        self.longitude_max = None

    def populate(self, input_latitude_array, input_longitude_array):
        """
        Populate the various maps used for geo search based on the input grid as defined by the input arrays. Note that
        the input arrays can be either 1- (regular grid) or 2-dimensional (eg: curvilinear grid).

        :param input_latitude_array: the 1- or 2-dimensional numpy containing latitude values from the input grid.
        :param input_longitude_array: the 1- or 2-dimensional numpy array containing longitude values from the input
        grid.
        """

        logger.debug('Building map of input grid')

        self.latitude_min = None
        self.latitude_max = None
        self.longitude_min = None
        self.longitude_max = None

        self.index_to_location_map = []
        _latitude_value_to_index_map = []
        _longitude_value_to_index_map = []

        # Flatten the input grid, depending on the number of dimensions of the input coordinates.
        if len(input_latitude_array.shape) == 1:

            # 1-dimensional input array, so regular grid.
            logger.debug('Inputs are 1-dimensional.')

            num_of_cells = len(input_latitude_array) * len(input_longitude_array)
            index = 0
            for latitude in input_latitude_array:
                for longitude in input_longitude_array:
                    if index % 10000 == 0:
                        logger.debug(f'{index} of {num_of_cells}')

                    # Position-dependent map, requiring that every input cell is stored.
                    self.index_to_location_map.append(Coordinate(latitude, longitude))

                    # The rest of the calculations are not position-dependent, so ignore NaN values.
                    if not math.isnan(latitude):

                        # Record the link between the latitude and the index as a tuple. This can be hidden behind the
                        # .isnan() call as position is not important, so we can ignore any NaNs.
                        _latitude_value_to_index_map.append((latitude, index))

                        if self.latitude_min is None or latitude < self.latitude_min:
                            self.latitude_min = latitude
                        if self.latitude_max is None or latitude > self.latitude_max:
                            self.latitude_max = latitude

                    if not math.isnan(longitude):

                        # Record the link between the longitude and the index as a tuple. This can be hidden behind the
                        # .isnan() call as position is not important, so we can ignore any NaNs.
                        _longitude_value_to_index_map.append((longitude, index))

                        if self.longitude_min is None or longitude < self.longitude_min:
                            self.longitude_min = longitude
                        if self.longitude_max is None or longitude > self.longitude_max:
                            self.longitude_max = longitude

                    index += 1

        else:

            # 2-dimensional input array, so non-regular grid, such as curvilinear.
            logger.debug('Inputs are 2-dimensional.')

            num_of_cells = input_latitude_array.shape[0] * input_longitude_array.shape[1]
            index = 0
            for i in range(input_latitude_array.shape[0]):
                for j in range(input_latitude_array.shape[1]):
                    if index % 10000 == 0:
                        logger.debug(f'{index} of {num_of_cells}')
                    latitude = input_latitude_array[i, j]
                    longitude = input_longitude_array[i, j]

                    # Position-dependent map, requiring that every input cell is stored.
                    self.index_to_location_map.append(Coordinate(latitude, longitude))

                    # The rest of the calculations are not position-dependent, so ignore NaN values.
                    if not math.isnan(latitude):

                        # Record the link between the latitude and the index as a tuple. This can be hidden behind the
                        # .isnan() call as position is not important, so we can ignore any NaNs.
                        _latitude_value_to_index_map.append((latitude, index))

                        if self.latitude_min is None or latitude < self.latitude_min:
                            self.latitude_min = latitude
                        if self.latitude_max is None or latitude > self.latitude_max:
                            self.latitude_max = latitude

                    if not math.isnan(longitude):

                        # Record the link between the longitude and the index as a tuple. This can be hidden behind the
                        # .isnan() call as position is not important, so we can ignore any NaNs.
                        _longitude_value_to_index_map.append((longitude, index))

                        if self.longitude_min is None or longitude < self.longitude_min:
                            self.longitude_min = longitude
                        if self.longitude_max is None or longitude > self.longitude_max:
                            self.longitude_max = longitude

                    index += 1

        # Private helper method to slice the geo search maps (for lat and lon) into smaller maps to minimise compares
        # when searching.
        def _slice_value_to_index_map(arr, field_name):
            _output = []

            # Sort the input map by value before slicing into sub-maps. The map contains tuples, where the first value
            # is the latitude or longitude, and the second value is the index.
            _sorted_map = np.sort(np.array(arr, dtype=[(field_name, float), ('index', int)]), order=field_name)

            # Slice the input map into smaller maps. These numbers are chosen arbitrarily.
            map_size = len(_sorted_map)
            sub_map_size = int(map_size / 20) if map_size > 200 else map_size
            start_index = 0
            while start_index < map_size - 1:
                end_index = start_index + sub_map_size - 1
                if end_index < map_size - 1:
                    # Although this is meant to be the strict size of the sub map, include any further cells that have
                    # the same value. This makes searching easier later.
                    end_value = _sorted_map[end_index][0]
                    while end_index < map_size and _sorted_map[end_index][0] == end_value:
                        end_index += 1
                if end_index > map_size:
                    end_index = map_size - 1

                _output.append(
                    (_sorted_map[start_index], _sorted_map[end_index], _sorted_map[start_index:end_index + 1])
                )
                start_index = end_index + 1

            return _output

        # Separate the value to index maps into smaller maps for quicker searching.
        self.latitude_value_to_index_maps = _slice_value_to_index_map(_latitude_value_to_index_map, 'latitude')
        self.longitude_value_to_index_maps = _slice_value_to_index_map(_longitude_value_to_index_map, 'longitude')

    def geo_search(self, latitude, longitude, search_box_width):
        """
        Perform a search on the input grid to find which cells fall within the search box, and then return the indexes
        of those cells. Note that the return values have meaning only within the 1-dimensional array equivalent of the
        input grid.

        :param latitude: the search latitude.
        :param longitude: the search longitude.
        :param search_box_width: the width of the box (in decimal degrees) to search within.
        :return: an array of indexes to cells in the input grid that fall within the search area.
        """

        # Calculate the search parameters for convenience.
        search_latitude_start = latitude - search_box_width
        search_latitude_end = latitude + search_box_width
        search_longitude_start = longitude - search_box_width
        search_longitude_end = longitude + search_box_width

        # Find the indexes of the latitude array that fall within the search box.
        input_latitude_indexes = set()

        # Loop through the sub-maps of the latitude to index mappings.
        num_index_maps = len(self.latitude_value_to_index_maps)
        for latitude_value_to_index_map in self.latitude_value_to_index_maps:

            # Find the start and end latitudes for the sub-map. The submap contains a tuple, where the 1st value is
            # the first latitude in the sub-map, the 2nd value is the last latitude in the sub-map, and the 3rd value is
            # the actual sub-map. Each entry in the sub-map is also a tuple, where the 1st value is the latitude and the
            # 2nd value is the corresponding index.
            start_value = latitude_value_to_index_map[0][0]
            end_value = latitude_value_to_index_map[1][0]

            # Only search this sub-map if it contains latitudes of interest, or if this is the only sub-map.
            if num_index_maps == 1 or start_value <= search_latitude_start <= end_value or \
                    start_value <= search_latitude_end <= end_value:
                for value_and_index in latitude_value_to_index_map[2]:
                    input_latitude = value_and_index[0]
                    input_index = value_and_index[1]
                    if search_latitude_start <= input_latitude <= search_latitude_end:
                        input_latitude_indexes.add(input_index)

        # Find the indexes of the input longitude array that fall within the search box.
        input_longitude_indexes = set()

        # Loop through the sub-maps of the longitude to index mappings.
        num_index_maps = len(self.longitude_value_to_index_maps)
        for longitude_value_to_index_map in self.longitude_value_to_index_maps:

            # Find the start and end longitudes for the sub-map. The submap contains a tuple, where the 1st value is
            # the first longitude in the sub-map, the 2nd value is the last longitude in the sub-map, and the 3rd value
            # is the actual sub-map. Each entry in the sub-map is also a tuple, where the 1st value is the longitude and
            # the 2nd value is the corresponding index.
            start_value = longitude_value_to_index_map[0][0]
            end_value = longitude_value_to_index_map[1][0]

            # Only search this sub-map if it contains latitudes of interest, or if this is the only sub-map.
            if num_index_maps == 1 or start_value < search_longitude_start <= end_value or \
                    start_value <= search_longitude_end <= end_value:
                for value_and_index in longitude_value_to_index_map[2]:
                    input_longitude = value_and_index[0]
                    input_index = value_and_index[1]
                    if search_longitude_start <= input_longitude <= search_longitude_end:
                        input_longitude_indexes.add(input_index)

        # Retain only those indexes that are in both lists.
        found_indexes = []
        for index in input_latitude_indexes:
            if index in input_longitude_indexes:
                found_indexes.append(index)
        return found_indexes

    def calculate_distance(self, output_latitude, output_longitude, input_cell_index, input_grid_proxy):
        """
        Calculate the distance between the output grid cell (identified by 'output_latitude/output_longitude') and the
        input grid cell (identified by 'input_cell_index').

        :param output_latitude: the latitude of the output grid cell.
        :param output_longitude: the longitude of the output grid cell.
        :param input_cell_index: the index of the input grid cell.
        :return: the distance in degree decimals.
        """
        input_latitude = self.index_to_location_map[input_cell_index].latitude
        input_longitude = self.index_to_location_map[input_cell_index].longitude
        return math.sqrt(pow((output_latitude - input_latitude), 2) + pow((output_longitude - input_longitude), 2))


# ----------------------------------------------------------------------------------------------------------------------

# Helper utility for taking an input array (either latitude or longitude), and returning the min, max, and
# indexToValues and valueToIndexes lists.

def build_from_input_grid(input_latitude_array, input_longitude_array, resolution=DEFAULT_RESOLUTION):
    """
    Build the regridder based on the input latitudes, longitudes and resolution. This function assumes that each
    input dimension (latitude/longitude) does NOT contain unique values, but a map of the index (the position within
    the array) to the corresponding dimensional value. For example, to identify the geo-location of the 5th cell of
    data, it is necessary to get the latitude from the 5th cell of the latitude array, and the longitude from the
    5th cell of the longitude array. This type of input dimension format is used for curvilinear grids, or grids
    that are random (not regular).

    :param latitudes: this is an array of latitudes. There
    :param longitudes:
    :param resolution:
    :return:
    """

    # Build a proxy for easy access to details of the input grid.
    input_grid_proxy = InputGridProxy()
    input_grid_proxy.populate(input_latitude_array, input_longitude_array)

    # Build the output grid based on the bounds of the input grid.
    # Populate the latitudes and longitudes for the output grid.
    output_latitude_array, output_longitude_array = _build_regular_output_grid(input_grid_proxy,
                                                                               resolution,
                                                                               DEFAULT_NUMBER_OF_DECIMAL_PLACES)

    # Build a map linking each output grid cell to the closest input grid cells for calculating weighted means when
    # regridding.
    output_grid_index_to_closest_input_cells = \
        _build_output_grid_index_to_closest_input_cells(input_grid_proxy,
                                                        output_latitude_array, output_longitude_array,
                                                        search_box_increment_size=resolution / 2,
                                                        max_search_box_increments=3,
                                                        min_found_cell_count=4)

    return regridder.Regridder(output_latitude_array, output_longitude_array, output_grid_index_to_closest_input_cells)


# ----------------------------------------------------------------------------------------------------------------------

def _build_regular_output_grid(input_grid_proxy,
                               resolution=DEFAULT_RESOLUTION,
                               number_of_decimal_places=DEFAULT_NUMBER_OF_DECIMAL_PLACES):
    """
    Build a regular output grid based on the bounds of the input grid, the desired resolution, and the number of
    decimal places.

    :param input_grid_proxy: the proxy object representing the input grid.
    :param resolution: the resolution of the output grid, used to determine distance between cells.
    :param number_of_decimal_places: the number of decimal places to maintain in the geo location values.
    :return: A tuple containing the latitude and longitude arrays for the output grid.
    """
    latitude_count = int(math.ceil((input_grid_proxy.latitude_max - input_grid_proxy.latitude_min) / resolution))
    output_latitude_array = []
    for index in range(latitude_count):
        output_latitude_array.append(
            round(input_grid_proxy.latitude_min + (resolution * index), number_of_decimal_places))

    longitude_count = int(math.ceil((input_grid_proxy.longitude_max - input_grid_proxy.longitude_min) / resolution))
    output_longitude_array = []
    for index in range(longitude_count):
        output_longitude_array.append(
            round(input_grid_proxy.longitude_min + (resolution * index), number_of_decimal_places))

    return output_latitude_array, output_longitude_array


# ----------------------------------------------------------------------------------------------------------------------

# Define a helper method to do the work, allowing for multi-processing.
def _populate_output_cell(input_grid_proxy, output_latitude, output_longitude, search_box_increment_size,
                          max_search_box_increments, min_found_cell_count, output_grid_index_to_closest_input_cells_map,
                          output_cell_index, size):
    """
    Helper function for execution in a sub-process to calculate the closest input cells to a given output cell.
    """

    # Provide feedback on operation.
    if output_cell_index % 100 == 0:
        logger.debug(f'Count {output_cell_index} of {size}')

    # Find the closest input pixels/cells for the current output latitude/longitude. The search area (box)
    # starts small and expands until at least the specified minimum number of input grid cells are found, or
    # until the number of increments exceeds the specified maximum.
    search_box_width = 0.0
    search_count = 1
    found_input_cells_indexes = []
    while (len(found_input_cells_indexes) < min_found_cell_count) and (
            search_count <= max_search_box_increments):
        search_box_width += search_box_increment_size
        found_input_cells_indexes = input_grid_proxy.geo_search(output_latitude,
                                                                output_longitude,
                                                                search_box_width)
        search_count += 1

    # If no cells/pixels were found, this likely means the output pixel/grid is outside the bounds of the input
    # grid (eg: on land or something), so ignore. Otherwise, identify and cache the four (4) closest input
    # pixels/cells along with their distances to the output pixel/cell.
    if len(found_input_cells_indexes) > 0:
        output_grid_index_to_closest_input_cells_map[output_cell_index] = _find_closest_four(output_latitude,
                                                                                             output_longitude,
                                                                                             input_grid_proxy,
                                                                                             found_input_cells_indexes)
    else:
        output_grid_index_to_closest_input_cells_map[output_cell_index] = []


def _build_output_grid_index_to_closest_input_cells(input_grid_proxy,
                                                    output_latitude_array, output_longitude_array,
                                                    search_box_increment_size,
                                                    max_search_box_increments=3,
                                                    min_found_cell_count=4):
    """
    Build a map linking each output grid cell (latitude x longitude) to the the closest input grid cells based on the
    maximum search box size and the minimum number of cells to find.
    :param input_grid_proxy: a proxy object representing the input grid.
    :param output_latitude_array: the 1-dimensional array containing latitude values of the output grid.
    :param output_longitude_array: the 1-dimensional array containing longitude values of the output grid.
    :param search_box_increment_size: the size by which the search box is increased when searching for matching input
     grid cells.
    :param max_search_box_increments: the maximum number of increments in search box size before the search is
    abandoned.
    :param min_found_cell_count: the minimum number of input cells to find for each output cell. Once this number of
    cells is found, searching is abandoned and the found cells are used.
    :return:
    """

    logger.debug('Building map of input grid -> output grid')

    # This function uses a small helper function to perform the calculation so the load can be shared across processes.
    _processes = []

    # Initialise the output map so each sub-process can write directly to it.
    size = len(output_latitude_array) * len(output_longitude_array)
    output_grid_index_to_closest_input_cells_map = Manager().list()
    for i in range(size):
        output_grid_index_to_closest_input_cells_map.append(None)

    # Loop through every cell in the output grid.
    output_cell_index = 0
    for output_latitude in output_latitude_array:
        for output_longitude in output_longitude_array:

            # Create a sub-process for calculating the closest input cells for each output cell.
            process = Process(
                target=_populate_output_cell,
                args=(input_grid_proxy, output_latitude, output_longitude, search_box_increment_size,
                      max_search_box_increments, min_found_cell_count, output_grid_index_to_closest_input_cells_map,
                      output_cell_index, size,)
            )
            _processes.append(process)
            process.start()

            output_cell_index += 1

    # Once all sub-processes have started, wait until they finish.
    for process in _processes:
        process.join()

    return output_grid_index_to_closest_input_cells_map


# ----------------------------------------------------------------------------------------------------------------------

def _find_closest_four(output_latitude, output_longitude, input_grid_proxy, found_input_cells_indexes):
    """
    Find the four (4) closest cells in 'found_input_cells_indexes' to the specified 'output_latitude/output_longitude',
    using the properties of 'input_grid_proxy' which represents the input grid.

    :return: the closest cells, up to four (4), sorted by distance.
    """

    # Build a list of distances.
    unsorted_list = []
    for input_cell_index in found_input_cells_indexes:
        distance = input_grid_proxy.calculate_distance(output_latitude, output_longitude, input_cell_index, input_grid_proxy)
        unsorted_list.append(IndexWithDistance(input_cell_index, distance))

    # Sort the list by distance.
    sorted_list = sorted(unsorted_list, key=lambda entry: entry[1])

    # Return only the top (closest) four pixels to the coordinate.
    if len(sorted_list) > 4:
        return sorted_list[0:4]
    else:
        return sorted_list

