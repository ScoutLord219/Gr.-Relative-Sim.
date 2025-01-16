extends Area3D

class_name Cell

const EPI : float = 1e-6


@export var time_step : float = 0.1

@onready var pivot = $Pivot

signal progress
signal metric_computed

var s_coordinates : Vector3
var t_coordinate : float

var affecting_masses : Array[Massive_Object]

var metric_history : Dictionary = {}

var metric : R2_Tensor
var inverse_metric : R2_Tensor = R2_Tensor.new() 
var metric_derivatives : Array

var christoffel_symbols : R3_Tensor = R3_Tensor.new() 
var riemann : R4_Tensor = R4_Tensor.new()
var ricci_tensor : R2_Tensor = R2_Tensor.new() 
var ricci_scalar : float = 0.0
var einstein_tensor : R2_Tensor = R2_Tensor.new() 

var lapse_function : float = 0.0


var flat_metric : bool = false


var max_error : float = 0.0


var spatial_ind : Array = ["r", "theta", "phi"]

#region Initialize

func _ready():
	Spacetime.progress_simulation.connect(progress_simulation.bind())
	Spacetime.evolve_simulation.connect(evolve_system.bind())

func init(grid_spacing : float):
	scale.x = grid_spacing
	scale.z = grid_spacing
	
	define_spherical_coordinates()

func define_spherical_coordinates():
	var x = global_transform.origin.x
	var y = global_transform.origin.y
	var z = global_transform.origin.z
	
	var r : float = sqrt((x ** 2) + (y ** 2) + (z ** 2))
	var theta : float
	var phi : float 
	
	if r == 0:
		theta = 0
		phi = 0
	else:
		theta = acos(z/r)  
		phi = atan2(y, x)
	
	s_coordinates = Vector3(r, theta, phi)
	#print(s_coordinates)

#endregion



#region Main Loop



func progress_simulation():
	accurately_approximate_metric_tensor()

func accurately_approximate_metric_tensor():
	if metric == null:
		assume_initial_metric()
	
	for i in range(25):
		#loop goes here
		compute_geometric_terms()
		
		var residual_error : R2_Tensor = metric.compute_residual_error(einstein_tensor)
		
		#print(einstein_tensor.components)
		
		
		if metric.max_error <= EPI:
			#print("METRIC : Accurately Approximated")
			break
		
		elif  metric.within_singularity || inverse_metric == null || flat_metric:
			#print("Flat Metric")
			break
		
		else:
			var metric_adjustment = metric.compute_metric_adjustment(metric, edify_einstein_tensor, residual_error, EPI)
			
			if metric_adjustment == null:
				#print("Null Adjustment")
				break
			
			metric = adjust_metric_tensor(metric, metric_adjustment)
			metric_history[t_coordinate] = metric
		#print("Max Error : ", metric.max_error, " : Iterations : ", i)
	
	scale_cell_mesh(ricci_scalar)
	
	print("Ricci Scalar : ", ricci_scalar, " : Distance from Center : ", abs(s_coordinates.x))
	
	metric_computed.emit()

func adjust_metric_tensor(input_metric: R2_Tensor, delta_metric : R2_Tensor, alpha: float = 0.1) -> R2_Tensor:
	var updated_metric : R2_Tensor = input_metric.duplicate_tensor()

	# Apply the adjustment with damping
	for mu in Tensor.ARR_IND:
		for nu in Tensor.ARR_IND:
			var current_value = input_metric.get_value(mu, nu)
			var adjustment = alpha * delta_metric.get_value(mu, nu)
			updated_metric.set_value(mu, current_value + adjustment, nu)

	return updated_metric

#endregion

#region Define Geometry
func invalid_geometry():
	christoffel_symbols.R3_ZERO()
	riemann.R4_ZERO()
	ricci_tensor.R2_ZERO()
	ricci_scalar = 0.0
	einstein_tensor.R2_ZERO() 

func assume_initial_metric():
	metric = R2_Tensor.new()
	
	if affecting_masses.size() == 0:
		metric.MINKOWSKI()
		flat_metric = true
		print("No Masses")
		return
	
	var black_hole
	var white_hole
	
	for mass in affecting_masses:
		if mass.is_in_group("Black Hole"):
			black_hole = mass
		elif mass.is_in_group("White Hole"):
			white_hole = mass
	
	
	var dist_sqr_to_bh : float = sqrt(s_coordinates.distance_squared_to(black_hole.s_coordinates) + black_hole.radius)
	var dist_sqr_to_wh : float = sqrt(s_coordinates.distance_squared_to(white_hole.s_coordinates) + white_hole.radius) 
	
	#print(dist_sqr_to_wh, "Dist to White")
	#print(dist_sqr_to_bh, "Dist to Black")
	
	var black_hole_metric : R2_Tensor = R2_Tensor.new()
	var white_hole_metric : R2_Tensor = R2_Tensor.new() 
	
	black_hole_metric.components = black_hole_metric.schwarz(black_hole.mass, s_coordinates.x, s_coordinates.y).components
	white_hole_metric.components = white_hole_metric.time_reversed_schwarz(white_hole.mass, sqrt(dist_sqr_to_wh), white_hole.s_coordinates.y - s_coordinates.y).components
	
	#print(dist_sqr_to_bh - dist_sqr_to_wh)
	
	if dist_sqr_to_bh < (dist_sqr_to_wh + 5):
		metric.components = black_hole_metric.components
		if  dist_sqr_to_bh < black_hole.radius:
			metric.within_singularity = true
			metric.MINKOWSKI()
			flat_metric = true
		print("Black Metric")
		
	elif dist_sqr_to_wh < (dist_sqr_to_bh + 5):
		metric.components = white_hole_metric.components
		if  dist_sqr_to_wh < white_hole.radius:
			metric.within_singularity = true
			metric.MINKOWSKI()
			flat_metric = true
		print("White Metric")
	
	elif abs(dist_sqr_to_bh - dist_sqr_to_wh) <= 15:
		metric.components = metric.tensor_addition_rank2(black_hole_metric, white_hole_metric).components
		print("Combined Metric")
	
	else:
		metric.MINKOWSKI()
		flat_metric = true
	



func calculate_christoffel(input_metric : R2_Tensor) -> R3_Tensor:
	var christoffel_and_derivatives : Dictionary = {}
	inverse_metric = input_metric.invert_tensor_with_fallback(input_metric)
	
	#print("Metric Determinant : ", metric.compute_determinant(matrix), " : Location : ", s_coordinates.x)
	if inverse_metric == null || input_metric.within_singularity:
		#print("Invalid Metric - Likely Within Singularity")
		christoffel_symbols.R3_ZERO()
		return null
	else:
		christoffel_and_derivatives = input_metric.christoffel_symbols(input_metric, inverse_metric, 0.1)
		christoffel_symbols = christoffel_and_derivatives["Christoffel"]
		metric_derivatives = christoffel_and_derivatives["Derivatives"]
		#print(christoffel_symbols)
	return christoffel_symbols

func rouse_riemann_tensor(input_christoffel_symbols : R3_Tensor) -> R4_Tensor:
	riemann = input_christoffel_symbols.compute_riemann_tensor(input_christoffel_symbols)
	return riemann

func rect_ricci_tensor(input_riemann : R4_Tensor) -> R2_Tensor:
	# Contract the Riemann tensor over the first and third indices
	ricci_tensor = input_riemann.contract_R4_tensor(input_riemann)
	return ricci_tensor

func reckon_ricci_scalar(input_inverse_metric : R2_Tensor, input_ricci_tensor : R2_Tensor) -> float:
	if input_inverse_metric == null:
		print("Inverse Metric NULL - How'd you get here?")
		return 1.0
	
	var r_s = 0.0
	
	#Contraction : multiply values of i_m and r_t together, and then add up all those products
	for mu in Tensor.ARR_IND:
		for nu in Tensor.ARR_IND:
			r_s += input_inverse_metric.get_value(mu, nu) * input_ricci_tensor.get_value(mu, nu)
	
	ricci_scalar = r_s
	
	return r_s

func edify_einstein_tensor(input_metric : R2_Tensor) -> R2_Tensor:
	var calculated_et : R2_Tensor = R2_Tensor.new()
	
	
	var input_christoffel : R3_Tensor = calculate_christoffel(input_metric)
	if input_christoffel == null:
		print("Invalid Christoffel")
		return null
	
	var input_riemann : R4_Tensor = rouse_riemann_tensor(input_christoffel)
	if input_christoffel == null:
		print("Invalid Riemann")
		return null
		
	var input_ricci_tensor : R2_Tensor = rect_ricci_tensor(input_riemann)
	if input_christoffel == null:
		print("Invalid Ricci Tensor")
		return null
		
	var input_ricci_scalar : float = reckon_ricci_scalar(inverse_metric, input_ricci_tensor)
	if input_christoffel == null:
		print("Invalid Ricci Scalar")
		return null
	
	for mu in Tensor.ARR_IND:
		for nu in Tensor.ARR_IND:
			var ricci_component = input_ricci_tensor.get_value(mu, nu)
			var metric_component = input_metric.get_value(mu, nu)
			
			var einstein_component = ricci_component - 0.5 * metric_component * input_ricci_scalar
			calculated_et.set_value(mu, einstein_component, nu)
	
	einstein_tensor = calculated_et
	return einstein_tensor

func compute_geometric_terms():
	if metric.within_singularity || inverse_metric == null || flat_metric:
		#print("Point in Singularity - Halting Calculations : ")
		invalid_geometry()
		return
	
	var einstein_check = edify_einstein_tensor(metric)
	
	if metric.within_singularity || inverse_metric == null || flat_metric || einstein_check == null:
		#print("Point in Singularity - Halting Calculations : ")
		invalid_geometry()
		return
	#print("GEOMETRIC TERMS : IDENTIFIED")
	#YIPEPEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE

#endregion

#region Compute Evolved Values

func evolve_system():
	var spatial_matrix : Array = []
	var extrinsic_curvature : R3_Tensor = R3_Tensor.new()
	var partial_t : R3_Tensor = R3_Tensor.new()
	var lapse : float = 0.0
	var shift_vector : Dictionary
	
	#spatial_matrix = extract_spatial_metric(metric)
	
	
	t_coordinate += time_step
	print("cells evolving")

func calculate_extrinsic_curvature(
		spatial_metric: R3_Tensor,
		partial_t_spatial_metric: R3_Tensor,
		lapse: float,
		shift_vector: Vector3,
		christoffel_symbols: R3_Tensor
	) -> R3_Tensor:
	# Initialize the Extrinsic Curvature Tensor
	var extrinsic_curvature = R3_Tensor.new()
	extrinsic_curvature.R3_ZERO()

	# Calculate the covariant derivatives of the shift vector
	var covariant_shift = R3_Tensor.new()
	covariant_shift.R3_ZERO()
	for i in Tensor.ARR_IND:
		for j in Tensor.ARR_IND:
			for k in Tensor.ARR_IND:
				var term = shift_vector[j] * christoffel_symbols.get_R3_value(k, i, j)
				covariant_shift.set_R3_value(k, i, j, term)
	
	# Compute the extrinsic curvature tensor
	for i in Tensor.ARR_IND:
		for j in Tensor.ARR_IND:
			var d_ij = partial_t_spatial_metric.get_R3_value("0", i, j)  # Time derivative of g_ij
			var cov_shift = covariant_shift.get_R3_value("0", i, j) + covariant_shift.get_R3_value("0", j, i)
			
			# Combine terms according to the formula
			var k_ij = (d_ij - cov_shift) / (2.0 * lapse)
			extrinsic_curvature.set_R3_value("0", i, j, k_ij)

	return extrinsic_curvature

func extract_spatial_metric(metric: R2_Tensor) -> R2_Tensor:
	var spatial_metric = R2_Tensor.new()
	spatial_metric.R3_ZERO()
	
	for i in ["x", "y", "z"]:
		for j in ["x", "y", "z"]:
			spatial_metric.set_R3_value(i, j, "0", metric.get_R3_value(i, j, "0"))
	
	return spatial_metric

func calculate_partial_t_spatial_metric(
		metric_history: Array,  # Array of past metrics for finite difference
		dt: float               # Time step size
	) -> R3_Tensor:
	if metric_history.size() < 2:
		push_error("Insufficient metric history for finite difference.")
		return null
	
	var recent_metric = metric_history[-1]
	var previous_metric = metric_history[-2]
	
	var partial_t_metric = R3_Tensor.new()
	partial_t_metric.R3_ZERO()
	
	for i in ["x", "y", "z"]:
		for j in ["x", "y", "z"]:
			var delta = (recent_metric.get_R3_value(i, j, "0") - previous_metric.get_R3_value(i, j, "0")) / dt
			partial_t_metric.set_R3_value(i, j, "0", delta)
	
	return partial_t_metric

func compute_lapse(
		initial_alpha: float,
		extrinsic_curvature: R3_Tensor,
		metric: R3_Tensor,
		dt: float,
		max_iterations: int = 10,
		tolerance: float = 1e-5
	) -> float:
	# Start with the initial guess for lapse
	var alpha = initial_alpha
	
	for iteration in range(max_iterations):
		# Compute the trace of the extrinsic curvature K = g^{ij} K_{ij}
		var trace_k = 0.0
		for i in ["x", "y", "z"]:
			for j in ["x", "y", "z"]:
				var gij = metric.get_R3_value(i, j, "0")
				var kij = extrinsic_curvature.get_R3_value(i, j, "0")
				trace_k += gij * kij
		
		# Update lapse using 1+log slicing condition
		var new_alpha = alpha - 2.0 * alpha * trace_k * dt
		
		# Check convergence
		if abs(new_alpha - alpha) < tolerance:
			break
		
		alpha = new_alpha
	
	return max(alpha, 0.01)  # Ensure lapse stays positive

func compute_shift_vector(
		initial_beta: Dictionary,
		metric: R3_Tensor,
		christoffel_symbols: R3_Tensor,
		extrinsic_curvature: R3_Tensor,
		eta: float,
		dt: float,
		max_iterations: int = 10,
		tolerance: float = 1e-5
	) -> Dictionary:
	# Initialize the shift vector and auxiliary variable
	var beta = initial_beta.duplicate()
	var B = { "x": 0.0, "y": 0.0, "z": 0.0 }
	
	for iteration in range(max_iterations):
		# Compute the Christoffel symbols contracted with the spatial metric
		var gamma_i = {}
		for i in ["x", "y", "z"]:
			gamma_i[i] = 0.0
			for j in ["x", "y", "z"]:
				for k in ["x", "y", "z"]:
					var gij = metric.get_R3_value(i, j, "0")
					var gamma_jk = christoffel_symbols.get_R3_value(i, j, k)
					gamma_i[i] += gij * gamma_jk
		
		# Update auxiliary variable B^i
		for i in ["x", "y", "z"]:
			B[i] += dt * (eta * gamma_i[i] - B[i])
		
		# Update shift vector beta^i
		var new_beta = {}
		for i in ["x", "y", "z"]:
			new_beta[i] = beta[i] + dt * B[i]
		
		# Check convergence
		var max_change = 0.0
		for i in ["x", "y", "z"]:
			max_change = max(max_change, abs(new_beta[i] - beta[i]))
		
		if max_change < tolerance:
			break
		
		beta = new_beta
	
	return beta

func vec_to_dictionary(input_vector : Vector3) -> Dictionary:
	var new_dictionary : Dictionary = {"x" : input_vector.x, "y" : input_vector.y, "z" : input_vector.z}
	return new_dictionary

#endregion



#region Godot Specific

func scale_cell_mesh(scalar : float):
	var new_pos_y = 0 - scalar
	
	pivot.scale.y = scalar

func append_mass(mass : Massive_Object):
	affecting_masses.append(mass)

func remove_mass(mass : Massive_Object):
	affecting_masses.erase(mass)

#endregion
