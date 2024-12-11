extends Resource

class_name Tensor

const G_CONST : float = 6.6743e-11
const S_O_L : float = 299792458
const SMALL_VALUE : float = 1e-10

const CHR_KEY_ARR : Array = [["t", "t", "t"], ["r", "r", "r"], ["phi", "phi", "phi"], ["theta", "theta", "theta"]]
const MET_KEY_ARR : Array = [["t", "t"], ["r", "r"], ["phi", "phi"], ["theta", "theta"]]
const FOUR_VEC_KEY_ARR : Array = [["t"], ["r"], ["phi"], ["theta"]]



var components : Dictionary = {}

func ZERO(key_array : Array):
	set_multiple_components(key_array, [0, 0, 0, 0])

func set_component(indices : Array, value : float):
	components[indices] = value

func set_multiple_components(indices : Array, value_array : Array):
	var value_index : float = 0.0
	for mu in indices:
		components[mu] = value_array[value_index]
		value_index += 1

func assign_all_tensor_components_to_other(other_tensor : Tensor):
	components = other_tensor.components

func get_components(indices : Array) -> float:
	return components[indices]

func get_spatial_components() -> Vector3:
	var spatial_comps = Vector3(components[FOUR_VEC_KEY_ARR[1]], components[FOUR_VEC_KEY_ARR[2]], components[FOUR_VEC_KEY_ARR[3]])
	
	return spatial_comps

func metric_to_minkowski():
	var metric_comp_value_array : Array = [-1, 1, 1, 1]
	
	set_multiple_components(MET_KEY_ARR, metric_comp_value_array)

func additive_on_component(indicies : Array, value : float):
	set_component(indicies, get_components(indicies) + value)

func add(tensor_a: Tensor, tensor_b: Tensor) -> Tensor:
	var result = Tensor.new()
	for key in tensor_a.components.keys():
		if tensor_b.components.has(key):
			result.set_component([key], tensor_a.get_components([key]) + tensor_b.get_components([key]))
	return result

func subtract(tensor_a: Tensor, tensor_b: Tensor) -> Tensor:
	var result = Tensor.new()
	for key in tensor_a.components.keys():
		if tensor_b.components.has(key):
			result.set_component(key, tensor_a.get_components(key) - tensor_b.get_components(key))
	return result

func scale(tensor: Tensor, scalar: float) -> Tensor:
	var result = Tensor.new()
	result.ZERO(tensor.components.keys())
	for key in tensor.components.keys():
		result.set_component(key, tensor.get_components(key) * scalar)
	return result

func dot_product(tensor_a: Tensor, tensor_b: Tensor, metric: Tensor = null) -> float:
	var result = 0.0
	for key in tensor_a.components.keys():
		var term = tensor_a.get_components([key]) * tensor_b.get_components([key])
		if metric != null:
			term *= metric.get_components([key])
		result += term
	return result

func multiply(tensor_a: Tensor, tensor_b: Tensor) -> Tensor:
	var result = Tensor.new()
	for key in tensor_a.components.keys():
		if tensor_b.components.has(key):
			result.set_components([key], tensor_a.get_components([key]) * tensor_b.get_components([key]))
	return result

func magnitude_with_metric(metric: Tensor) -> float:
	var magnitude = 0.0
	
	# Loop through the tensor's components
	for index in components.keys():
		# Use the metric tensor to calculate the inner product (dot product)
		# e.g., for each component, multiply by the corresponding metric component
		var metric_component = metric.get_components([index[0], index[0]])
		magnitude += components[index] * metric_component * components[index]
	
	# Take the square root of the sum to get the magnitude
	return sqrt(magnitude)

func normalize_four_velocity(velocity: Tensor, metric: Tensor) -> Tensor:
	# Compute the Minkowski norm: v^μ v_μ using the metric tensor
	var norm = 0.0
	for key in velocity.components.keys():
		var metric_component = metric.get_components([key[0], key[0]])
		norm += metric_component * velocity.get_components(key) ** 2
	
	# Ensure norm is -1 (timelike); correct the time component to normalize
	if norm != 0:
		var correction_factor = 1.0 / sqrt(abs(norm))  # Adjust for normalization
		velocity.set_component(["t"], velocity.get_components(["t"]) * correction_factor)
	
	return velocity

func norm_squared_spatial() -> float:
	# Assuming spatial components are stored under the keys "x", "y", "z"
	var r = self.get_components(["r"])
	var phi = self.get_components(["phi"])
	var theta = self.get_components(["theta"])
	
	# Compute the squared norm of the spatial components
	return r * r + phi * phi + theta * theta



func gaussian_elimination_invert() -> Tensor:
	"""
	Performs Gaussian elimination to invert the metric tensor.
	Assumes the metric tensor only has diagonal components.
	"""
	var diag_keys = ["t", "r", "phi", "theta"]  # List of diagonal keys
	var n = diag_keys.size()  # Number of spacetime dimensions
	
	# Initialize augmented matrix with metric components and identity
	var augmented_matrix = []
	for i in range(n):
		var row = []
		for j in range(n):
			if i == j:
				row.append(get_components([diag_keys[i], diag_keys[j]]))  # Diagonal element
			else:
				row.append(0.0)  # Off-diagonal elements are zero
		augmented_matrix.append(row)
	
	# Perform Gaussian elimination to invert the matrix
	for i in range(n):
		# Find the pivot row with the largest absolute value in column i
		var pivot_row = i
		for j in range(i + 1, n):
			if abs(augmented_matrix[j][i]) > abs(augmented_matrix[pivot_row][i]):
				pivot_row = j
		
		# Swap the pivot row with the current row
		if pivot_row != i:
			var temp = augmented_matrix[i]
			augmented_matrix[i] = augmented_matrix[pivot_row]
			augmented_matrix[pivot_row] = temp
		
		# Normalize the pivot row
		var pivot = augmented_matrix[i][i]
		for j in range(n):
			augmented_matrix[i][j] /= pivot
		
		# Eliminate the column below and above the pivot
		for j in range(n):
			if j != i:
				var factor = augmented_matrix[j][i]
				for k in range(n):
					augmented_matrix[j][k] -= factor * augmented_matrix[i][k]
	
	# Create a new Tensor for the inverse
	var inverse_tensor = Tensor.new()
	var inverse_components = {}
	
	# Fill in the inverse tensor components
	for i in range(n):
		for j in range(n):
			if i == j:  # Only diagonal components exist
				inverse_components[[diag_keys[i], diag_keys[j]]] = augmented_matrix[i][j]
	
	# Set the components in the inverse tensor
	inverse_tensor.set_multiple_components(inverse_components.keys(), inverse_components.values())
	
	return inverse_tensor


func inverse() -> Tensor:
	var inverse_tensor = Tensor.new()
	var inverse_components = {}
	
	# Loop through the keys defined in MET_KEY_ARR
	for row in MET_KEY_ARR:
		for col in MET_KEY_ARR:
			# Access components using the current row and column indices
			var g_ij = get_components(row)
			
			# Diagonal components are inverted directly
			var inverse_value = 1.0 / g_ij if row == col else 0.0
			
			# Store the inverse value with 2-string keys
			inverse_components[row] = inverse_value
			
	# Set the calculated inverse components in the new tensor
	inverse_tensor.set_multiple_components(inverse_components.keys(), inverse_components.values())
	
	return inverse_tensor

func clone() -> Tensor:
	var cloned_tensor : Tensor = Tensor.new()
	cloned_tensor.set_multiple_components(components.keys(), components.values())
	return cloned_tensor
