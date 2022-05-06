import numpy as np

def import_bin_data(filename):
	# Imports the data from filename and organizes
	# it for analysis.
	
	np.set_printoptions(threshold=np.inf)
	file_data = np.genfromtxt(filename,dtype=float,delimiter=',',missing_values='none',filling_values=np.nan)
	sorted_file_data = file_data[np.lexsort((file_data[:,1],file_data[:,0]))]
	sorted_file_data[:,0] += 1.0 ;  # Shifts the 0 radial value to 1
	radial_bin_width,theta_bin_width = get_bin_size_information(filename)
	#print(radial_bin_width,theta_bin_width)
	sorted_file_data[:,0] *= radial_bin_width
	sorted_file_data[:,1] *= theta_bin_width
	print(sorted_file_data)
	return sorted_file_data,radial_bin_width,theta_bin_width

def get_num_bins(array_name):
	# Gets the number of bins in a 
	# column of the array

	bins_in_first_column = int(array_name[-1,0])
	j = 0
	for i in array_name:
		if i[0] == array_name[0,0]:
			j += 1

	
	print("Number of bins is: {}, {}".format(bins_in_first_column,j))
	return bins_in_first_column,j

def get_bin_size_information(filename):
	# Takes the file as input and
	# reads the comments to determine
	# the bin sizes for use in data
	# analysis.
	
	with open(filename,"r") as f:
		for line in f:
			if line.startswith('#'):
				bin_sizes = line.split()

	radial_bin_size = float(bin_sizes[1])
	theta_bin_size = float(bin_sizes[2])

	return radial_bin_size,theta_bin_size

def nearest_neighbor_test(array_name):
	# Takes an array as input and tests
	# each row to see if it has at least 
	# one neighbor in each dimension. Returns
	# indices of rows that have both nearest
	# neighbors.

	# Get number of bins in each direction
	bins_in_first_column, bins_in_second_column = get_num_bins(array_name)

	
	indices_with_neighbors = []
	index = 0
	for i in array_name:
		print("Index is {}".format(index))
		radial_ahead = index + bins_in_second_column
		radial_behind = index - bins_in_second_column
		theta_ahead = index + 1
		theta_behind = index - 1
		# Handles the radial exceptions of the first and last radial bins which can only have nearest neighbors
		# in one direction.
		if i[0] == 1.0:
			# If the radial dimension nearest neighbor is empty, don't add to the list
			if np.isnan(array_name[radial_ahead,2]):
				print("Exception 1")
				index += 1
				continue
		if i[0] == float(bins_in_first_column):
			if np.isnan(array_name[radial_behind,2]):
				print("Exception 2")
				index += 1
				continue

		# If both ahead and behind radial bins are empty, don't add to the list
		if i[0] != 1.0 and i[0] != float(bins_in_first_column):
			if np.isnan(array_name[radial_ahead,2]) and np.isnan(array_name[radial_behind,2]):
				print("Exception 3")
				index += 1
				continue

		# These exceptions handle the periodic nature of the polar plot
		# Resets theta_ahead and theta_behind appropriately
		if i[1] == 0.0:
			theta_behind = index + bins_in_second_column
		if i[1] == float(bins_in_second_column - 1 ):
			theta_ahead = index - bins_in_second_column


		if np.isnan(array_name[theta_ahead,2]) and np.isnan(array_name[theta_behind,2]):
			print("Exception 4")
			index +=1
			continue

		print(i)
		indices_with_neighbors.append(index)
		index += 1

	print(indices_with_neighbors)
	return indices_with_neighbors

def vector_elements_for_normal_to_surface(r, radial_difference, radial_step, theta, theta_difference, theta_step):
	# Identifies the nearest neighbor cells for 
	# computation of the derivatives to calculate
	# the membrane surface normal.
	# Preferentially uses the central difference 
	# method, followed by the forward/backward
	# difference method. These are unit vectors that
	# have been normalized

	c1 = 1/r
	c2 = 1/r**2
	radial_derivative = radial_difference / radial_step
	theta_derivative = theta_difference / theta_step
	normalization_factor = np.sqrt( 1 + (c2 * theta_derivative**2) + radial_derivative**2  )
	neg_hx = ( c1 * np.sin(theta) * theta_derivative - np.cos(theta) * radial_derivative ) / normalization_factor
	neg_hy = ( -1 * ( np.sin(theta) * radial_derivative + c1 * np.cos(theta) * theta_derivative) ) / normalization_factor
	z = 1 / normalization_factor

	return neg_hx,neg_hy,z

def get_surface_normals(filename):
	# Takes an array of 2 dimensional data
	# representing the Monge surface and 
	# calculates the surface normals for 
	# each bin.

	array_name,radial_bin_width,theta_bin_width = import_bin_data(filename)

	# Gets the indices of the array that have nearest neighbors
	# to calculate the derivative
	indices_with_neighbors = nearest_neighbor_test(array_name)

	bins_in_first_column,bins_in_second_column = get_num_bins(array_name)
	number_of_entries = bins_in_first_column * bins_in_second_column
	normal_vectors = np.zeros((number_of_entries,5))

	j = 0
	for i in array_name:
		#print("The last row of ith line is {}".format(i[2]))
		normal_vectors[j,0] = i[0]
		normal_vectors[j,1] = i[1]
		if np.isnan(i[2]):
			normal_vectors[j,2] = np.nan
			normal_vectors[j,3] = np.nan
			normal_vectors[j,4] = np.nan
			print(normal_vectors[j])
			j += 1
			continue
		elif j not in indices_with_neighbors:
			normal_vectors[j,2] = np.nan
			normal_vectors[j,3] = np.nan
			normal_vectors[j,4] = np.nan
			print(normal_vectors[j])
			j += 1
			continue
		# Tests whether there is forward and backward neighbors
		# for the cell for the numerical derivative later.
		radial_ahead = j + bins_in_second_column
		radial_behind = j - bins_in_second_column
		radial_step = 1.0
		theta_ahead = j + 1
		theta_behind = j - 1
		theta_step = 1.0

		# The if-else loop calculates the radial difference for all
		# contingencies (I think)
		if i[0] == 1.0:
			r = i[0]
			radial_difference =  array_name[radial_ahead,2] - i[2]
		elif i[0] == float(bins_in_first_column):
			r = i[0]
			radial_difference = i[2] - array_name[radial_behind,2]
		elif np.isfinite(array_name[radial_ahead,2]) and np.isfinite(array_name[radial_behind,2]):
			r = i[0]
			radial_difference = array_name[radial_ahead,2] - array_name[radial_behind,2]
			radial_step = 2.0
		elif np.isfinite(array_name[radial_ahead,2]) and np.isnan(array_name[radial_behind,2]):
			r = i[0]
			radial_difference = array_name[radial_ahead,2] - i[2]
		elif np.isnan(array_name[radial_ahead,2]) and np.isfinite(array_name[radial_behind,2]):
			r = i[0]
			radial_difference = i[2] - array_name[radial_behind,2]
		
		# Corrects theta_ahead and theta_behind for periodicity
		if i[1] == 0.0:
			theta_behind = j + bins_in_second_column
		elif i[1] == float(bins_in_second_column - 1):
			theta_ahead = j - bins_in_second_column
		
		if np.isfinite(array_name[theta_ahead,2]) and np.isfinite(array_name[theta_behind,2]):
			theta = i[1]
			theta_difference = array_name[theta_ahead,2] - array_name[theta_behind,2]
			theta_step = 2.0
		elif np.isfinite(array_name[theta_ahead,2]) and np.isnan(array_name[theta_behind,2]):
			theta = i[1]
			theta_difference = array_name[theta_ahead,2] - i[1]
		elif np.isnan(array_name[theta_ahead,2]) and np.isfinite(array_name[theta_behind,2]):
			theta = i[1]
			theta_difference = i[1] - array_name[theta_behind,2]

		normal_x, normal_y, normal_z = vector_elements_for_normal_to_surface(r,radial_difference,radial_step,theta,theta_difference,theta_step)
		normal_vectors[j,2] = normal_x
		normal_vectors[j,3] = normal_y
		normal_vectors[j,4] = normal_z
		print(normal_vectors[j])
		j += 1


#surface = import_bin_data('/data/marcario/elicNanodisc/elic110A/mongeSurface.upper.txt')
#nearest_neighbor_test(surface)
get_surface_normals('/data/marcario/elicNanodisc/elic110A/mongeSurface.upper.txt')
