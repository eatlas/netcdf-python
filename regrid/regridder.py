import math
import numpy as np

# Constant identifying the weight value to use to represent infinite weight when distance is zero.
INFINITE_WEIGHT = 999.0

def _calculate_weighted_mean(values, weights):
    """
    Utility method to perform a weighted mean calculation based on the input arrays.
    :param values: an array of values.
    :param weights: an array of weights for the corresponding values.
    :return: the weighted mean, or NaN if no real values included.
    """
    sum_values = 0.0
    sum_weights = 0.0
    is_non_zero = False
    for index in range(len(values)):

        # Only include if both value and weight are numbers.
        value = values[index]
        weight = weights[index]
        if not math.isnan(value) and not math.isnan(weight):
            sum_values += value * weight
            sum_weights += weight
            is_non_zero = True

    if is_non_zero:
        return sum_values / sum_weights
    else:
        return math.nan


class Regridder:
    """
    Utility class to convert a location on an input grid to the corresponding location on an output regular grid. This
    class is normally instantiated by 'builder.py'.
    """

    latitude_count = None
    """
    Cache the number of cells on the latitude dimension of the output grid.
    """

    latitude_array = None
    """
    Cached reference to the latitude array of the output grid.
    """

    longitude_count = None
    """
    Cache the number of cells on the longitude dimension of the output grid.
    """

    longitude_array = None
    """
    Cached reference to the longitude array of the output grid.
    """

    output_grid_index_to_closest_input_cells = None
    """
    Cached reference to the mapping of output cells to the closest corresponding input cells. This allows for the 
    calculation of weighted means for doing the regridding.
    """

    def __init__(self, latitude_array, longitude_array, output_grid_index_to_closest_input_cells):
        self.latitude_array = latitude_array
        self.latitude_count = len(latitude_array)
        self.longitude_array = longitude_array
        self.longitude_count = len(longitude_array)
        self.output_grid_index_to_closest_input_cells = output_grid_index_to_closest_input_cells

    def regrid(self, input_data):
        """
        Regrid the input data to the output grid configured for this instance of the regridder.

        :param input_data: a 2- (lat/lon) or 3-dimensional (lat/lon/depth) array of data in the same format used by the
        latitude/longitude data used to instantiate this regridder instance.
        :return:
        """

        nan_counter = 0
        value_counter = 0
        i = 0
        while i < np.shape(input_data)[0]:
            j = 0
            while j < np.shape(input_data)[1]:
                if math.isnan(input_data[i,j]):
                    nan_counter += 1
                else:
                    value_counter += 1
                j += 1
            i += 1

        # Ensure the input grid is either 2-dimensional (lat/lon) or 3-dimensional (lat/lon/depth).
        input_shape = np.shape(input_data)
        num_dimensions = len(input_shape)
        if num_dimensions != 2 and num_dimensions != 3:
            raise IndexError(f'Invalid number of dimensions. Expected 2 or 3; found {num_dimensions}.')

        # Repackage the data into a 1-dimensional array for latitude/longitude, and ensure there is a 3rd dimension
        # (height/depth) for consistent processing.
        repackaged_input_data = []
        if num_dimensions == 3:
            return np.reshape(
                self._regrid(
                    np.reshape(input_data,(input_shape[0] * input_shape[1], input_shape[2]))
                ),
                (self.latitude_count, self.longitude_count, input_shape[2])
            )
        else:
            return np.reshape(
                np.squeeze(
                    self._regrid(
                        np.reshape(input_data,(input_shape[0] * input_shape[1], 1))
                    ),
                    axis=1
                ),
                (self.latitude_count, self.longitude_count)
            )

    def _regrid(self, repackaged_input_data):

        nan_counter = 0
        value_counter = 0
        i = 0
        while i < np.shape(repackaged_input_data)[0]:
            j = 0
            while j < np.shape(repackaged_input_data)[1]:
                if math.isnan(repackaged_input_data[i,j]):
                    nan_counter += 1
                else:
                    value_counter += 1
                j += 1
            i += 1

        # Build the output result.
        num_layers = np.shape(repackaged_input_data)[1]
        output_data = np.empty((self.latitude_count * self.longitude_count, num_layers))

        # Loop through every cell in the output grid to populate it.
        for latitude_index in range(self.latitude_count):
            for longitude_index in range(self.longitude_count):

                # Calculate the position (index) being processed within the output grid.
                output_grid_index = latitude_index * self.longitude_count + longitude_index

                # Identify the input grid indexes that are used to calculate the output grid value.
                indexes_with_distances = self.output_grid_index_to_closest_input_cells[output_grid_index]

                # Process each packaged layer in turn.
                for layer_index in range(num_layers):

                    if len(indexes_with_distances) == 0:
                        # Handle scenario for no data.
                        output_data[output_grid_index, layer_index] = math.nan
                    else:
                        # Calculate weighted mean.
                        values = []
                        weights = []
                        for index_with_distance in indexes_with_distances:
                            values.append(repackaged_input_data[index_with_distance.index, layer_index])
                            weights.append(
                                INFINITE_WEIGHT if index_with_distance.distance == 0.0 else 1/index_with_distance.distance
                            )
                        weighted_mean = _calculate_weighted_mean(values, weights)
                        output_data[output_grid_index, layer_index] = weighted_mean

        return output_data
