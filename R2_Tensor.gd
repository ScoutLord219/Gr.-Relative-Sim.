extends Tensor

class_name R2_Tensor

var within_singularity : bool

var max_error : float = 0.0

#region Constant functions : accept no parameters and always return the same result

func R2_ZERO() -> void:
	for mu in ARR_IND:
		for nu in ARR_IND:
			set_value(mu, 0.0, nu)

func MINKOWSKI() -> void:
	R2_ZERO()
	for mu in ARR_IND:
		for nu in ARR_IND:
			set_value(mu, 0.0, nu)
	set_value("t", -1.0, "t")

func duplicate_tensor() -> R2_Tensor:
	var tensor_duplicate : R2_Tensor = R2_Tensor.new()
	
	for mu in ARR_IND:
		for nu in ARR_IND:
			tensor_duplicate.set_value(mu, get_value(mu, nu), nu)
	
	return tensor_duplicate

#endregion

#region Parameter Functions

func schwarz(mass : float, radial_distance : float, theta_angle : float) -> R2_Tensor:
	var schwarz_radius = (1 - (2 * G_CONST * mass) / radial_distance)
	
	var new_metric : R2_Tensor = R2_Tensor.new() 
	
	for mu in ARR_IND:
		for nu in ARR_IND:
			new_metric.set_value(mu, 0.0, nu)
	
	new_metric.set_value("t", -schwarz_radius, "t")
	new_metric.set_value("r", schwarz_radius ** -1, "r")
	new_metric.set_value("theta", radial_distance ** 2, "theta")
	new_metric.set_value("phi", (radial_distance ** 2) * sin(theta_angle) ** 2, "phi")
	
	return new_metric

func time_reversed_schwarz(mass : float, radial_distance : float, theta_angle : float) -> R2_Tensor:
	var schwarz_radius = (1 - (2 * G_CONST * mass) / radial_distance)
	
	var new_metric : R2_Tensor = R2_Tensor.new() 
	
	for mu in ARR_IND:
		for nu in ARR_IND:
			new_metric.set_value(mu, 0.0, nu)
	
	new_metric.set_value("t", schwarz_radius, "t")
	new_metric.set_value("r", -schwarz_radius ** -1, "r")
	new_metric.set_value("theta", radial_distance ** 2, "theta")
	new_metric.set_value("phi", (radial_distance ** 2) * sin(theta_angle) ** 2, "phi")
	
	return new_metric

func tensor_addition_rank2(metric_1 : R2_Tensor, metric_2 : R2_Tensor) -> R2_Tensor:
	var new_metric : R2_Tensor = R2_Tensor.new()
	
	for mu in ARR_IND:
		for nu in ARR_IND:
			var component_sum : float = metric_1.get_value(mu, nu) + metric_2.get_value(mu, nu)
			new_metric.set_value(mu, component_sum, nu)
	return new_metric

func subtract_R2_tensors(tensor_a: R2_Tensor, tensor_b: R2_Tensor) -> R2_Tensor:
	# Create a new tensor for the result
	var result_tensor = R2_Tensor.new()
	for mu in ARR_IND:
		for nu in ARR_IND:
			var difference = tensor_a.get_value(mu, nu) - tensor_b.get_value(mu, nu)
			result_tensor.set_value(mu, difference, nu)
	return result_tensor


func initialize_matrix(size: int) -> Array:
	var matrix = []
	for mu in range(size):
		matrix.append([])
		for nu in range(size):
			matrix[mu].append(0.0)
	return matrix

func matrix_R2_multiplication(matrix: Array, vector: R2_Tensor) -> R2_Tensor:
	#Perform matrix-vector multiplication.
	#
	#Parameters:
	#- matrix: A nested Array representing the matrix (e.g., 4x4 or larger for Jacobian).
	#- vector: A Tensor object representing the vector.
	#
	#Returns:
	#- A new Tensor object as the result of the multiplication.
	
	var result = R2_Tensor.new()  # Initialize an empty Tensor for the result
	
	for mu in ARR_IND:  # Loop over the output indices
		var sum_value = 0.0
		for nu in ARR_IND:  # Loop over the input indices
			sum_value += matrix[ARR_IND.find(mu)][ARR_IND.find(nu)] * vector.get_value(mu, nu)
			result.set_value(mu, sum_value, nu)
	
	return result

func tensor_dot_rank2(tensor_a: R2_Tensor, tensor_b: R2_Tensor, contraction_index: String) -> R2_Tensor:
	#Compute the contraction (dot product) of two rank-2 tensors along a specified index.
	#
	#Parameters:
	#- tensor_a: A rank-2 Tensor object.
	#- tensor_b: A rank-2 Tensor object.
	#- contraction_index: The index over which to contract (shared dimension, e.g., 't', 'r', etc.).
	#
	#Returns:
	#- A new rank-2 Tensor representing the result of the contraction.
	#
	var result = R2_Tensor.new()  # Resultant rank-2 tensor

	for mu in ARR_IND:  # Iterate over first index of result
		for nu in ARR_IND:  # Iterate over second index of result
			var sum_value = 0.0
			
			for sigma in ARR_IND:  # Contraction index (e.g., summing over 't', 'r', ...)
				sum_value += tensor_a.get_value(mu, sigma) * tensor_b.get_value(sigma, nu)
			
			result.set_value(mu, sum_value, nu)
	
	return result

func tensor_norm_rank2(tensor: R2_Tensor) -> float:
	#Compute the Frobenius norm of a rank-2 tensor.
	#Parameters:
	#- tensor: R2_Tensor, the rank-2 tensor to compute the norm for.
#
	#Returns:
	#- The Frobenius norm (a float).
	
	var norm_squared = 0.0
	for mu in ARR_IND:
		for nu in ARR_IND:
			var value = tensor.get_value(mu, nu)
			norm_squared += value * value

	return sqrt(norm_squared)

func tensor_scale_rank2(tensor: R2_Tensor, scale: float) -> R2_Tensor:
#
	#Scale all components of a rank-2 tensor by a scalar.
	#Parameters:
	#- tensor: R2_Tensor, the rank-2 tensor to scale.
	#- scale: float, the scalar multiplier.
#
	#Returns:
	#- A new scaled R2_Tensor.
	#
	var scaled_tensor = R2_Tensor.new()
	for mu in ARR_IND:
		for nu in ARR_IND:
			scaled_tensor.set_value(mu, tensor.get_value(mu, nu) * scale, nu)

	return scaled_tensor

#endregion

#region Christoffel Symbols

func christoffel_symbols(metric : R2_Tensor, inv_metric : R2_Tensor, delta : float) -> Dictionary:
	var gamma : R3_Tensor = R3_Tensor.new()
	
	var derivative_array : Array = []
	var term1_derivatives : R3_Tensor = R3_Tensor.new()
	var term2_derivatives : R3_Tensor = R3_Tensor.new() 
	var term3_derivatives : R3_Tensor = R3_Tensor.new()
	
	for mu in ARR_IND:
		for nu in ARR_IND:
			for lambda in ARR_IND:
				# Compute Christoffel symbols using the metric
				var term1 = compute_metric_derivative(metric, mu, nu)
				var term2 = compute_metric_derivative(metric, mu, lambda)
				var term3 = compute_metric_derivative(metric, nu, lambda)
				
				term1_derivatives.set_R3_value(mu, nu, lambda, term1)
				term2_derivatives.set_R3_value(mu, lambda, nu, term2)
				term3_derivatives.set_R3_value(nu, lambda, mu, term3)
				
				var sum_terms = term1 + term2 - term3
				
				# Multiply by 0.5 * inverse_metric
				gamma.set_R3_value(mu, nu, lambda, 0.5 * inv_metric.get_value(lambda, nu) * sum_terms)
	
	derivative_array.append(term1_derivatives)
	derivative_array.append(term2_derivatives)
	derivative_array.append(term3_derivatives)
	
	return {"Christoffel" : gamma, "Derivatives" : derivative_array}

func compute_metric_derivative(metric: R2_Tensor, mu: String, nu: String) -> float:
	var h = EPI  # Step size for finite difference
	var forward_shift = metric.duplicate_tensor()
	var backward_shift = metric.duplicate_tensor()
	
	# Perturb the specified coordinate direction
	forward_shift.self_additive(mu, nu, h)
	backward_shift.self_additive(mu, nu, -h)
	
	# Compute finite difference
	var value_forward = forward_shift.get_value(mu, nu)
	var value_backward = backward_shift.get_value(mu, nu)
	return (value_forward - value_backward) / (2 * h)

#endregion

#region Inverse Calculation

# All this for the inverse :\

func invert_tensor_with_fallback(tensor: R2_Tensor) -> R2_Tensor:
	var tensor_matrix : Array = initialize_matrix(4)
	for i in range(4):
		for j in range(4):
			tensor_matrix[i][j] = tensor.get_value(ARR_IND[i], ARR_IND[j])
	
	var determinant = compute_determinant(tensor_matrix)
	if abs(determinant) < EPI:
		print("Determinant is zero or too small; returning null")
		return null
	else:
		return invert_tensor(tensor)

func invert_tensor(metric: R2_Tensor) -> R2_Tensor:
	# Create a new R2_Tensor to store the inverse
	var inv_metric = R2_Tensor.new()

	# Convert the R2_Tensor into a 4x4 nested array for inversion
	var matrix : Array = initialize_matrix(4)
	var indices = metric.ARR_IND

	# Populate the 4x4 array from the R2_Tensor components
	for i in range(4):
		for j in range(4):
			matrix[i][j] = metric.get_value(indices[i], indices[j])
	
	# Compute the determinant (needed for inversion)
	var determinant = compute_determinant(matrix)
	if abs(determinant) < metric.EPI:
		push_error("Metric R2_Tensor is singular and cannot be inverted.")
		return null
	
	# Calculate the adjugate matrix (transpose of cofactor matrix)
	var adjugate = compute_adjugate(matrix)

	# Compute the inverse by dividing adjugate by the determinant
	for i in range(4):
		for j in range(4):
			matrix[i][j] = adjugate[i][j] / determinant
	
	# Convert the 4x4 array back to the R2_Tensor structure
	for i in range(4):
		for j in range(4):
			inv_metric.set_value(indices[i], matrix[i][j], indices[j])
	
	return inv_metric

func compute_determinant(matrix: Array) -> float:
	return (
		matrix[0][0] * (
			matrix[1][1] * (matrix[2][2] * matrix[3][3] - matrix[2][3] * matrix[3][2]) -
			matrix[1][2] * (matrix[2][1] * matrix[3][3] - matrix[2][3] * matrix[3][1]) +
			matrix[1][3] * (matrix[2][1] * matrix[3][2] - matrix[2][2] * matrix[3][1])
		) -
		matrix[0][1] * (
			matrix[1][0] * (matrix[2][2] * matrix[3][3] - matrix[2][3] * matrix[3][2]) -
			matrix[1][2] * (matrix[2][0] * matrix[3][3] - matrix[2][3] * matrix[3][0]) +
			matrix[1][3] * (matrix[2][0] * matrix[3][2] - matrix[2][2] * matrix[3][0])
		) +
		matrix[0][2] * (
			matrix[1][0] * (matrix[2][1] * matrix[3][3] - matrix[2][3] * matrix[3][1]) -
			matrix[1][1] * (matrix[2][0] * matrix[3][3] - matrix[2][3] * matrix[3][0]) +
			matrix[1][3] * (matrix[2][0] * matrix[3][1] - matrix[2][1] * matrix[3][0])
		) -
		matrix[0][3] * (
			matrix[1][0] * (matrix[2][1] * matrix[3][2] - matrix[2][2] * matrix[3][1]) -
			matrix[1][1] * (matrix[2][0] * matrix[3][2] - matrix[2][2] * matrix[3][0]) +
			matrix[1][2] * (matrix[2][0] * matrix[3][1] - matrix[2][1] * matrix[3][0])
		)
	)

func compute_adjugate(matrix: Array) -> Array:
	var adjugate = initialize_matrix(4)
	for i in range(4):
		for j in range(4):
			# Compute the minor matrix (3x3 submatrix) by removing row i and column j
			var minor = compute_minor(matrix, i, j)
			# Compute the determinant of the minor
			var minor_det = compute_determinant_3x3(minor)
			# Apply the cofactor sign alternation
			adjugate[j][i] = (-1 if (i + j) % 2 != 0 else 1) * minor_det
	return adjugate

func compute_determinant_3x3(matrix: Array) -> float:
	return (
		matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1]) -
		matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0]) +
		matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0])
	)

func compute_minor(matrix: Array, row: int, col: int) -> Array:
	var minor = []
	for i in range(4):
		if i == row:
			continue
		var minor_row = []
		for j in range(4):
			if j == col:
				continue
			minor_row.append(matrix[i][j])
		minor.append(minor_row)
	return minor

# I am astounded at how complex that turns out to be
#endregion



#region Compute Metric Adjustment

func compute_metric_adjustment(metric: R2_Tensor, einstein_function : Callable, residual_error: R2_Tensor, epsilon: float = 1e-6) -> R2_Tensor:
	# Step 1: Compute the Jacobian Matrix
	var jacobian = compute_jacobian_matrix(metric, residual_error, einstein_function, epsilon)
	
	# Step 2: Invert the Jacobian Matrix
	var jacobian_inverse = invert_matrix(jacobian)
	if jacobian_inverse == null:
		push_error("Jacobian matrix inversion failed. Metric adjustment cannot be computed.")
		return null
	
	# Step 3: Compute the metric adjustment
	var delta_metric = matrix_R2_multiplication(jacobian_inverse, residual_error)
	
	return delta_metric


func compute_jacobian_matrix(metric: R2_Tensor, initial_einstein : R2_Tensor, einstein_function : Callable, epsilon: float) -> Array:
	var jacobian : Array = initialize_matrix(4)
	
	for mu in range(4):
		for nu in range(4):
			# Create a perturbed metric tensor
			var perturbed_metric = metric.duplicate_tensor()

			# Apply a small positive perturbation
			perturbed_metric.self_additive(ARR_IND[mu], ARR_IND[nu], epsilon)

			# Compute the Einstein tensor for the perturbed metric
			var perturbed_einstein = einstein_function.call(perturbed_metric)
			
			if perturbed_einstein == null:
				break

			# Compute the finite difference approximation for the Jacobian matrix entry
			jacobian[mu][nu] = (perturbed_einstein.get_value(ARR_IND[mu], ARR_IND[nu]) - 
								initial_einstein.get_value(ARR_IND[mu], ARR_IND[nu])) / epsilon

	return jacobian

func invert_matrix(matrix: Array) -> Array:
	var determinant = compute_determinant_jacobian(matrix)
	if abs(determinant) < EPI:
		push_error("Matrix determinant is zero or near-zero. Cannot invert.")
		return []
	
	var adjugate = compute_adjugate_jacobian(matrix)
	var inverse = []
	
	for i in range(4):
		inverse.append([])
		for j in range(4):
			inverse[i].append(adjugate[i][j] / determinant)
	
	return inverse

func compute_determinant_jacobian(matrix: Array) -> float:
	if len(matrix) != 4 or len(matrix[0]) != 4:
		push_error("Determinant calculation only supports 4x4 matrices.")
		return 0.0
	
	var determinant = 0.0
	for j in range(4):
		determinant += matrix[0][j] * compute_minor_jacobian(matrix, 0, j)
	return determinant

func compute_adjugate_jacobian(matrix: Array) -> Array:
	var adjugate = []
	for i in range(4):
		adjugate.append([])
		for j in range(4):
			var minor = compute_minor_jacobian(matrix, i, j)
			adjugate[i].append(((-1) ** (i + j)) * minor)
	return adjugate

func compute_minor_jacobian(matrix: Array, row: int, col: int) -> float:
	var minor_matrix = []
	for i in range(4):
		if i == row:
			continue
		var minor_row = []
		for j in range(4):
			if j == col:
				continue
			minor_row.append(matrix[i][j])
		minor_matrix.append(minor_row)
	return compute_determinant_3x3_jacobian(minor_matrix)

func compute_determinant_3x3_jacobian(matrix: Array) -> float:
	if len(matrix) != 3 or len(matrix[0]) != 3:
		push_error("Determinant calculation only supports 3x3 matrices.")
		return 0.0

	# Using the determinant formula for 3x3 matrices:
	var determinant = (
		matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1]) -
		matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0]) +
		matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0])
	)

	return determinant

func compute_residual_error(einstein_tensor: R2_Tensor) -> R2_Tensor:
	var residual : R2_Tensor = einstein_tensor.duplicate_tensor()
	for mu in ARR_IND:
		for nu in ARR_IND:
			var value = einstein_tensor.get_value(mu, nu)
			residual.set_value(mu, value, nu)  # Use raw value, without abs
	return residual


#endregion
