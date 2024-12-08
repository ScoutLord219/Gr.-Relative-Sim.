extends CharacterBody3D

class_name Body

const SMALL_VALUE : float = 1e-6

@export var mass : float
@export var radius : float
@export var space_time : SPACETIME

@onready var cell_updater_shape = $"Cell Updater_Body Checker/Cell Updater Shape"

var cells_initialized : bool

var just_computed_four_position = Tensor.new()

var new_four_position = Tensor.new()
var previous_four_position = Tensor.new()

var four_velocity = Tensor.new()

var current_spacetime_metric = Tensor.new()
var current_christoffel = Tensor.new()

var affected_cells : Array
var current_closest_cell : Cell 

var proper_time : float = 0.0

var previous_velocity_position_distances : float
var difference_req_for_vel_recalc : float

func _init():
	var letter_array = ["t", "r", "phi", "theta"]
	
	for mu in letter_array:
		previous_four_position.set_components([mu], 0.0)
		new_four_position.set_components([mu], 0.0)
		just_computed_four_position.set_components([mu], 0.0)
		four_velocity.set_components([mu], 0.0)
		for nu in letter_array:
			if nu == mu:
				current_spacetime_metric.set_components([mu, nu], 0.0)
			for lambda in letter_array:
				current_christoffel.set_components([mu, nu, lambda], 0.0)

func on_cells_initialized():
	cells_initialized = true

func _ready():
	space_time.all_cells_loaded.connect(on_cells_initialized)
	
	await space_time.all_cells_loaded
	
	cell_updater_shape.scale *= radius * 2
	var test = Tensor.new()
	
	test.set_components(["r"], global_transform.origin.x)
	test.set_components(["phi"], global_transform.origin.y)
	test.set_components(["theta"], global_transform.origin.z)
	print(test.components)
	
	wait_for_cell_data(test)
	
	previous_four_position.components = just_computed_four_position.components
	previous_four_position.set_components(["t"], proper_time)

func find_closest_cells(four_position : Tensor) -> Array:
	var array_by_cell_distance : Array = []
	
	if current_closest_cell != null:
		current_closest_cell.disconnect("values_computed", initialize_cell_data)
	else:
		print(affected_cells)
	
	if affected_cells.size() > 0:
		var closest_cell = affected_cells[0]
		
		if closest_cell.four_position_components.components.size() == 0:
			print("No four_position on closest cell")
			affected_cells
		
		for cells in affected_cells:
			if four_position.metric_based_distance(cells.four_position_components, current_spacetime_metric) < four_position.metric_based_distance(closest_cell.four_position_components, current_spacetime_metric):
				closest_cell = cells
				array_by_cell_distance.push_front(closest_cell)
		current_closest_cell = closest_cell
	else:
		current_closest_cell = null
		print("No Cell Array @: ", self)
	
	if current_closest_cell != null:
		current_closest_cell.connect("values_computed", initialize_cell_data)
	
	return array_by_cell_distance

func initialize_cell_data(four_position_components : Tensor, metric_tensor : Tensor, schwar_christ_symbs : Tensor):
	just_computed_four_position.components = four_position_components.components
	current_spacetime_metric.components = schwar_christ_symbs.components
	current_christoffel.components = metric_tensor.components

func wait_for_cell_data(four_pos : Tensor) -> bool:
	if current_closest_cell == null:
		return false
	print("stopped")
	await current_closest_cell.values_computed
	print("moving")
	return true

func _physics_process(delta):
	if not cells_initialized:
		return
	
	if await wait_for_cell_data(previous_four_position):
		proper_time += current_closest_cell.calculate_proper_time() * delta
		
		new_four_position.components = just_computed_four_position.components
		new_four_position.set_components(["t"], proper_time)
		
		define_four_velocity_of_massive_body(previous_four_position, new_four_position)
		progress_along_geodesic(new_four_position, four_velocity, proper_time)
		
		move_and_slide()
		
		wait_for_cell_data(new_four_position)
		
		previous_four_position.components = new_four_position.components

# RK4 Integration for position and velocity
func runge_kutta_step(_position: Tensor, new_velocity: Tensor, delta: float) -> Dictionary:
	
	if new_velocity or _position == null:
		return {"position" : Vector3.ZERO, "velocity" : Vector3.ZERO}
	
	var acceleration = define_acceleration_of_massive_body(_position, new_velocity)
	
	print(acceleration)
	
	new_velocity.set_components(["t"], new_velocity.get_components(["t"]) + acceleration["t"] * delta)
	new_velocity.set_components(["r"], new_velocity.get_components(["r"]) + acceleration["r"] * delta)
	new_velocity.set_components(["phi"], new_velocity.get_components(["phi"]) + acceleration["phi"] * delta)
	new_velocity.set_components(["theta"], new_velocity.get_components(["theta"]) + acceleration["theta"] * delta)

	# k1 values
	var k1_velocity = acceleration  # a(t, y)
	var k1_position = new_velocity
	
	
	var altered_k1_pos : Tensor = _position
	for k1_keys in _position.keys():
		altered_k1_pos[k1_keys] = _position[k1_keys] + k1_position[k1_keys] * delta / 2
	
	var altered_k1_vel : Tensor = new_velocity
	for k1_keys in new_velocity.keys():
		altered_k1_vel[k1_keys] = new_velocity[k1_keys] + k1_velocity[k1_keys] * delta / 2
	
	# k2 values
	var k2_velocity = define_acceleration_of_massive_body(altered_k1_pos, altered_k1_vel)
	var k2_position = altered_k1_pos
	
	
	var altered_k2_pos : Tensor = _position
	for k2_keys in _position.keys():
		altered_k2_pos[k2_keys] = _position[k2_keys] + k2_position[k2_keys] * delta / 2
	
	var altered_k2_vel : Tensor = new_velocity
	for k2_keys in new_velocity.keys():
		altered_k2_vel[k2_keys] = new_velocity[k2_keys] + k2_velocity[k2_keys] * delta / 2
	
	# k3 values
	var k3_velocity = define_acceleration_of_massive_body(altered_k2_pos, altered_k2_vel)
	var k3_position = altered_k2_pos
	
	
	var altered_k3_pos : Tensor = _position
	for k3_keys in _position.keys():
		altered_k3_pos[k3_keys] = _position[k3_keys] + k3_position[k3_keys] * delta
	
	var altered_k3_vel : Tensor = new_velocity
	for k3_keys in new_velocity.keys():
		altered_k3_vel[k3_keys] = new_velocity[k3_keys] + k2_velocity[k3_keys] * delta
	
	# k4 values
	var k4_velocity = define_acceleration_of_massive_body(altered_k3_pos, altered_k3_vel)
	var k4_position = altered_k3_pos
	
	# Update position and velocity
	var new_position_tens : Tensor = _position
	var new_velocity_tens : Tensor = new_velocity
	
	for letter_keys in _position.keys():
		new_position_tens.set_components([letter_keys], _position[letter_keys] + delta / 6 * (k1_position[letter_keys] + 2 * k2_position[letter_keys] + 2 * k3_position[letter_keys] + k4_position[letter_keys]))
		new_velocity_tens.set_components([letter_keys], new_velocity[letter_keys] + delta / 6 * (k1_velocity[letter_keys] + 2 * k2_velocity[letter_keys] + 2 * k3_velocity[letter_keys] + k4_velocity[letter_keys]))
	
	var new_position_cart : Vector3 = spherical_to_cartesian_coords(new_position_tens)
	var new_velocity_cart : Vector3 = spherical_to_cartesian_coords(new_velocity_tens)
	return {"position": new_position_cart, "velocity": new_velocity_cart}

func progress_along_geodesic(current_four_position : Tensor, current_four_velocity : Tensor, affine_parameter : float):
	var result = runge_kutta_step(current_four_position, current_four_velocity, affine_parameter)
	
	for values in result["position"]:
		if values == NAN or INF:
			result["position"] = Vector3.ZERO
			print("invalid position detected")
	
	for values in result["velocity"]:
		if values == NAN or INF:
			result["velocity"] = Vector3.ZERO
			print("invalid velocity detected")
	
	result["position"].y = 0
	result["velocity"].y = 0
	
	
	velocity = result["velocity"]

func define_four_velocity_of_massive_body(old_pos : Tensor, new_pos : Tensor):
	if previous_velocity_position_distances == new_pos.metric_based_distance(old_pos, current_spacetime_metric): 
		four_velocity = four_velocity
		return
	
	previous_velocity_position_distances = new_pos.metric_based_distance(old_pos, current_spacetime_metric)
	
	for indices in new_pos.components.keys():
		four_velocity.set_components(indices, new_pos.subtract(new_pos, old_pos).get_components(indices) / proper_time)

func define_acceleration_of_massive_body(_position : Tensor, new_velocity : Tensor) -> Tensor:
	var four_acceleration = Tensor.new()
	var spatial_coord_vector = Vector3(_position.get_components(["r"]), _position.get_components(["phi"]), _position.get_components(["theta"]))
	
	var interpolated_christoffel_symbols = Tensor.new()
	var interpolated_metric = Tensor.new()  #Interpolate this in a similar fashion to how you interpolated christoffel symbols
	var array_of_cells : Array = []
	var limited_array_of_nearest_cells : Array = []
	var limit_on_number_of_neighbor_cells : int = 8
	
	if array_of_cells.size() == 0:
		return four_acceleration
	
	for cell_index in range(array_of_cells.size()):
		if cell_index > limit_on_number_of_neighbor_cells - 1:
			break
		limited_array_of_nearest_cells.append(array_of_cells[cell_index])
	
	interpolated_christoffel_symbols = interpolate_christoffel_symbols(spatial_coord_vector, limited_array_of_nearest_cells)
	interpolated_metric = interpolate_metric_tensor_components(spatial_coord_vector, limited_array_of_nearest_cells)
	
	var equation_value_prior_to_negative_symbol : float = 0
	
	for mu_key in four_acceleration.keys():
		equation_value_prior_to_negative_symbol = 0
		for nu_key in interpolated_christoffel_symbols[mu_key].keys():
			for lambda_key in interpolated_christoffel_symbols[mu_key][nu_key].keys():
				var term_value = 0 
				if proper_time > 0: #may have to substitute with simply a realyl small value if errors keep happening
					term_value = interpolated_christoffel_symbols[mu_key][nu_key][lambda_key] * new_velocity[nu_key] * new_velocity[lambda_key] / proper_time
				else:
					term_value = 0
				equation_value_prior_to_negative_symbol += term_value
		four_acceleration[mu_key] = -(equation_value_prior_to_negative_symbol)
	
	four_acceleration = ensure_orthogonality(four_acceleration, four_velocity, interpolated_metric)
	
	four_acceleration = normalize_acceleration(four_acceleration, interpolated_metric)
	
	return four_acceleration

func normalize_acceleration(accel_dict : Tensor, interpolated_met : Tensor) -> Tensor:
		# Compute norm of four-acceleration
	var norm = 0.0
	for mu_key in accel_dict.keys():
		var term_value : float
		term_value = interpolated_met[mu_key + mu_key] * accel_dict[mu_key] ** 2
		norm += term_value
	norm = sqrt(abs(norm))
		
	# Limit to speed of light
	var c = 1.0  # Assuming natural units
	if norm > c:
		for mu_key in accel_dict.keys():
			accel_dict[mu_key] *= c / norm
	return accel_dict

func ensure_orthogonality(four_acceleration: Tensor, four_vel: Tensor, metric: Tensor) -> Tensor:
	# Calculate inner product g_{\mu\nu} u^\mu a^\nu
	var inner_product = 0.0
	for mu_key in four_vel.components.keys():
		inner_product += metric[mu_key + mu_key] * four_vel[mu_key] * four_acceleration[mu_key]
		
	# Calculate norm of the four-velocity g_{\mu\nu} u^\mu u^\nu
	var velocity_norm = 0.0
	for mu_key in four_vel.components.keys():
		velocity_norm += metric[mu_key + mu_key] * four_vel[mu_key] * four_vel[mu_key]
	
	# Ensure orthogonality
	var orthogonal_acceleration = Tensor.new()
	for mu_key in four_acceleration.components.keys():
		orthogonal_acceleration.set_components([mu_key], four_acceleration.get_components([mu_key]) - (inner_product / velocity_norm) * four_vel.get_components([mu_key]))
		
	return orthogonal_acceleration

func interpolate_christoffel_symbols(_position: Vector3, neighboring_cells: Array) -> Tensor:
	var interpolated_symbols = Tensor.new()  # Store interpolated Christoffel symbols
	var total_weight = 0.0
	var weights = []
	var epsilon = 1e-6  # Small value to avoid division by zero
			
	# Step 1: Calculate weights based on distance
	for cell in neighboring_cells:
		var distance = new_four_position.metric_based_distance(cell.four_position.components, current_spacetime_metric)
		var weight = 1.0 / max(distance, epsilon)
		weights.append(weight)
		total_weight += weight
			
	# Step 2: Normalize weights
	for i in range(weights.size()):
		weights[i] /= total_weight
			
	# Step 3: Interpolate Christoffel symbols
	for i in range(neighboring_cells.size()):
		var cell = neighboring_cells[i]
		var weight = weights[i]
		var cell_christoffel = cell.schwar_christ_symbs
			
		for mu in cell_christoffel.keys():
			for nu in cell_christoffel.components.keys():
				for lambda in cell_christoffel.components.keys():
					interpolated_symbols.set_components([mu][nu][lambda], (interpolated_symbols.get_components([mu, nu, lambda])) + weight * cell_christoffel[mu][nu][lambda])
			
	return interpolated_symbols

func interpolate_metric_tensor_components(four_position : Tensor, neighboring_cells : Array) -> Tensor:
	var interpolated_metric_components = Tensor.new()  # Store interpolated Christoffel symbols
	var total_weight = 0.0
	var weights = []
	var epsilon = 1e-6  # Small value to avoid division by zero
			
	# Step 1: Calculate weights based on distance
	for cell in neighboring_cells:
		var distance = four_position.metric_based_distance(cell.four_position_components, current_spacetime_metric)
		var weight = 1.0 / max(distance, epsilon)
		weights.append(weight)
		total_weight += weight
			
	# Step 2: Normalize weights
	for i in range(weights.size()):
		weights[i] /= total_weight
			
	# Step 3: Interpolate Christoffel symbols
	for i in range(neighboring_cells.size()):
		var cell = neighboring_cells[i]
		var weight = weights[i]
		var cell_metric = cell.metric_tensor_components
			
		for mu in cell_metric.components.keys():
			interpolated_metric_components.set_components([mu], interpolated_metric_components.get_components([mu]) + weight * cell_metric[mu])
			
	return interpolated_metric_components

func spherical_to_cartesian_coords(components : Tensor) -> Vector3:
	var spatial_coordinates_in_cart : Vector3 = Vector3.ZERO
	
	spatial_coordinates_in_cart.x = components.get_components(["r"]) * sin(components.get_components(["theta"])) * cos(components.get_components(["phi"]))
	spatial_coordinates_in_cart.y = 0
	spatial_coordinates_in_cart.z = components.get_components(["r"]) * cos(components.get_components(["theta"]))
	
	if spatial_coordinates_in_cart.is_finite():
		return spatial_coordinates_in_cart
	else:
		print("Invalid cartesian coordinates: ")
		return Vector3.ZERO  # Return a default value if invalid


func _on_cell_updater_area_entered(area):
	if area.is_in_group("Cell"):
		area.warp(self)
		affected_cells.append(area)

func _on_cell_updater_area_exited(area):
	if area.is_in_group("Cell"):
		area.unwarp(self)
		affected_cells.remove_at(affected_cells.find(area))
