extends Node3D

class_name Cell


const SMALL_VALUE = 1e-6 #so we dont encounter any division by 0 errors

@onready var pivot = $Pivot
@onready var mesh = $Pivot/Mesh


signal values_computed(four_position : Tensor, metric : Tensor, christoffel : Tensor)
signal initializing_finished

var space_time
var update_components : bool = false

var affecting_bodies : Array[Body]

var four_position_components = Tensor.new()
var metric_tensor = Tensor.new()
var schwar_christ_symbs = Tensor.new()

var effective_r : float
var sum_of_masses : float
#only t, r, and phi are really important; theta is not needed - unless we need it for calculations >:( - because we'll be limiting spacial movement to a 2D plane

func init_started(cell_dimensions : Vector2, space_time_parent : SPACETIME) -> int:
	init(cell_dimensions, space_time_parent)
	await initializing_finished
	
	return 1

func init(cell_dimensions : Vector2, space_time_parent : SPACETIME):
	mesh.scale.x = cell_dimensions.x
	mesh.scale.z = cell_dimensions.y
	
	space_time = space_time_parent
	
	var letter_array = ["t", "r", "phi", "theta"]
	
	for mu in letter_array:
		for nu in letter_array:
			if nu == mu:
				metric_tensor.set_components([mu, nu], 0.0)
			for lambda in letter_array:
				schwar_christ_symbs.set_components([mu, nu, lambda], 0.0)
	
	define_flat_space_time_metric()
	
	initializing_finished.emit()
	update_components = true

func warp(object):
	affecting_bodies.append(object)
	sum_of_masses += object.mass
	calculate_schwar_christ_symbs()
	update_components = true

func unwarp(object):
	affecting_bodies.remove_at(affecting_bodies.find(object))
	sum_of_masses -= object.mass
	calculate_schwar_christ_symbs()
	update_components = false

func _physics_process(_delta):
	if update_components:
		define_flat_space_time_metric()
		calculate_schwar_christ_symbs()

func calculate_effective_r(masses : Array, metric : Tensor) -> float:
	#finds harmonic mean of distances to all bodies - good for modeling tidal forces and localized phenomena apparently
	
	if masses.size() == 0:
		return 1
	
	var distance_sum = 0.0
	for body in masses:
		var distance = four_position_components.metric_based_distance(body.new_four_position, metric) + body.radius
		if distance > 0:
			distance_sum += 1.0 / distance
	
	if distance_sum > 0:
		effective_r = masses.size() / distance_sum
	else:
		effective_r = 0.0
		print("approaching singularity")
	return effective_r

func calculate_proper_time() -> float:
	if metric_tensor.components.size() == 0:
		return 0.0
	var g_tt = metric_tensor.get_components(["t", "t"])
	for body in affecting_bodies:
		var distance = four_position_components.metric_based_distance(body.new_four_position, metric_tensor) + body.radius
		g_tt -= 2.0 * space_time.G_CONST * body.mass / max(distance, 1e-10)
	return g_tt

func define_metric_tensor_components(defined_sum_of_masses : float):
	# initialize flat spacetime to begin estimates
	#lowkey add some kinda check so that if theres already a viable value for metric tensor then just skip this step
	define_flat_space_time_metric()
	cartesian_to_four_position_multi(global_transform.origin, affecting_bodies, defined_sum_of_masses, metric_tensor, calculate_proper_time())
	
	var number_of_iterations : int = 3
	var calculated_effective_r : float
	
	for iteration in range(number_of_iterations):
		calculated_effective_r = calculate_effective_r(affecting_bodies, metric_tensor)
	
		metric_tensor.set_components(["t", "t"], metric_tensor.get_components(["t", "t"]) + (2 * space_time.G_CONST * defined_sum_of_masses) / calculated_effective_r)
		metric_tensor.set_components(["r", "r"], metric_tensor.get_components(["r", "r"]) + (2 * space_time.G_CONST * defined_sum_of_masses) / (calculated_effective_r * (1 - (2 * space_time.G_CONST * defined_sum_of_masses) / calculated_effective_r)))
		metric_tensor.set_components(["phi", "phi"], metric_tensor.get_components(["phi", "phi"]) + calculated_effective_r ** 2 * sin(four_position_components.get_components(["theta"])) ** 2)
		metric_tensor.set_components(["theta", "theta"], metric_tensor.get_components(["theta", "theta"]) + calculated_effective_r ** 2)

func define_flat_space_time_metric():
	metric_tensor.set_components(["t", "t"], -1)  # Start with flat spacetime
	metric_tensor.set_components(["r", "r"], 1)  # Start with flat spacetime
	metric_tensor.set_components(["phi", "phi"], 1)  # Initialize angular components
	metric_tensor.set_components(["theta", "theta"], 1)

func cartesian_to_four_position(reference_point : Vector3, cartesian_coords: Vector3, metric: Tensor, time: float) -> Tensor:
	# Convert Cartesian to spherical coordinates
	var r = reference_point.distance_to(cartesian_coords)
	var phi = atan2(cartesian_coords.y, cartesian_coords.x)  # azimuthal angle
	var theta = acos(cartesian_coords.z / r)  # polar angle

	# Create the spatial position in spherical coordinates
	var spatial_position = Tensor.new()
	spatial_position.set_components(["r"], r)
	spatial_position.set_components(["phi"], phi)
	spatial_position.set_components(["theta"], theta)

	# Apply metric adjustments for curvature
	var corrected_position = Tensor.new()
	for key in spatial_position.components.keys():
		var g_ij = metric.get_components([key, key])  # Assuming diagonal metric for simplicity
		corrected_position.set_components([key], sqrt(abs(g_ij)) * spatial_position.get_components([key]))

	# Combine with time component to form the four-position
	var four_position = Tensor.new()
	four_position.set_components(["t"], time)  # Set time component
	for key in corrected_position.components.keys():
		four_position.set_components([key], corrected_position.get_components([key]))

	return four_position

func calculate_schwar_christ_symbs():
	#reset christoffel symbols for re-evaluation
	if affecting_bodies.size() == 0:
		return 
	
	define_metric_tensor_components(sum_of_masses)
	
	var r_averaged_twixt_masses : float = calculate_effective_r(affecting_bodies, metric_tensor)
	
	var metric_tensor_rr = metric_tensor.get_components(["r","r"])
	
	var array_of_body_positions = []
	
	
	if metric_tensor_rr == 0:
		#print("Error: metric_tensor[\"rr\"] is zero. Adjust calculations.")
		print("Metric RR of 0")
		return # Or handle gracefully
	
		
	for current_body in affecting_bodies:
		array_of_body_positions.append(current_body.global_transform.origin)
		
		var gm_over_r = (space_time.G_CONST * current_body.mass) / max(r_averaged_twixt_masses, 1e-10)
		var one_minus_two_gm_over_r = 1 - (2 * gm_over_r)
		
		var shared_t_symbol_value = (gm_over_r) * one_minus_two_gm_over_r ** -1
		schwar_christ_symbs.set_components(["t", "r", "t"], schwar_christ_symbs.get_components(["t", "r", "t"]) + shared_t_symbol_value)
		schwar_christ_symbs.set_components(["t", "t", "r"], schwar_christ_symbs.get_components(["t", "t", "r"]) + shared_t_symbol_value)
		
		schwar_christ_symbs.set_components(["r", "t", "t"], schwar_christ_symbs.get_components(["r", "t", "t"]) + (gm_over_r ** -2) * one_minus_two_gm_over_r)
		schwar_christ_symbs.set_components(["r", "r", "r"], schwar_christ_symbs.get_components(["r", "r", "r"]) + -(gm_over_r ** -2) * (one_minus_two_gm_over_r) ** -1)
		schwar_christ_symbs.set_components(["r", "phi", "phi"], schwar_christ_symbs.get_components(["r", "phi", "phi"]) + -metric_tensor_rr * one_minus_two_gm_over_r)
		schwar_christ_symbs.set_components(["r", "theta", "theta"], schwar_christ_symbs.get_components(["r", "theta", "theta"]) + -metric_tensor_rr * sin(metric_tensor.get_components(["theta", "theta"])) ** 2 * one_minus_two_gm_over_r)
		
		var shared_theta_and_phi_symbol_value = 1/metric_tensor_rr
		schwar_christ_symbs.set_components(["theta", "phi", "r"], schwar_christ_symbs.get_components(["theta", "phi", "r"]) + shared_theta_and_phi_symbol_value)
		schwar_christ_symbs.set_components(["theta", "r", "phi"], schwar_christ_symbs.get_components(["theta", "r", "phi"]) + shared_theta_and_phi_symbol_value)
		schwar_christ_symbs.set_components(["theta", "theta", "theta"], schwar_christ_symbs.get_components(["theta", "theta", "theta"]) + -sin(metric_tensor.get_components(["theta", "theta"])) * cos(metric_tensor.get_components(["theta", "theta"])))
		
		var shared_phi_symbol_value = cos(metric_tensor.get_components(["theta", "theta"]))/sin(metric_tensor.get_components(["theta", "theta"]))
		schwar_christ_symbs.set_components(["phi", "theta", "r"], schwar_christ_symbs.get_components(["phi", "theta", "r"]) + shared_theta_and_phi_symbol_value)
		schwar_christ_symbs.set_components(["phi", "theta", "phi"], schwar_christ_symbs.get_components(["phi", "theta", "phi"]) + shared_theta_and_phi_symbol_value)
		schwar_christ_symbs.set_components(["phi", "r", "theta"], schwar_christ_symbs.get_components(["phi", "r", "theta"]) + shared_phi_symbol_value)
		schwar_christ_symbs.set_components(["phi", "phi", "theta"], schwar_christ_symbs.get_components(["phi", "phi", "theta"]) + shared_phi_symbol_value)
	
	four_position_components.components = cartesian_to_four_position_multi(global_transform.origin, array_of_body_positions, affecting_bodies, metric_tensor, calculate_proper_time())

	values_computed.emit(four_position_components, metric_tensor, schwar_christ_symbs)
	print("values computed")

func scale_cell(value : float):
	pivot.scale.y = value

func cartesian_to_four_position_multi(cartesian_coords: Vector3, reference_points: Array, masses: Array, metric: Tensor, time: float) -> Tensor:
	# Compute effective radial distance
	var r_effective = calculate_effective_r(masses, metric)

	# Continue with the rest of the original function using r_effective
	var relative_position = cartesian_coords - reference_points[0]  # Example for orientation
	var phi = atan2(relative_position.y, relative_position.x)
	var theta = acos(relative_position.z / r_effective)

	var spatial_position = Tensor.new()
	spatial_position.set_components(["r"], r_effective)
	spatial_position.set_components(["phi"], phi)
	spatial_position.set_components(["theta"], theta)

	# Apply metric adjustments
	var corrected_position = Tensor.new()
	for key in spatial_position.components.keys():
		var g_ij = metric.get_components([key, key])  # Diagonal metric for simplicity
		corrected_position.set_components([key], sqrt(abs(g_ij)) * spatial_position.get_components([key]))

	var four_position = Tensor.new()
	four_position.set_components(["t"], time)
	for key in corrected_position.components.keys():
		four_position.set_components([key], corrected_position.get_components([key]))

	return four_position
