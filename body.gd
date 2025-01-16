extends CharacterBody3D

class_name Massive_Object

@onready var mesh_pivot = $"Mesh Pivot"
@onready var cell_updater_shape = $"Cell Updater_Body Checker/Cell Updater Shape"

@export var mass : float
@export var is_static : bool
var radius : float

var s_coordinates : Vector3
var t_coordinate : float = 0.0
var s_velocity : Vector3
var christoffels_at_position : R3_Tensor = R3_Tensor.new()

var lapse : float = 0.0

var closest_cells : Array = []


func _ready():
	define_spherical_coordinates()
	radius = schwarzschild_radius()
	mesh_pivot.scale *= radius
	
	cell_updater_shape.scale *= 100.0
	
	#Spacetime.evolve_simulation.connect(evolve_self.bind())

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

func schwarzschild_radius() -> float:
	var _metric : R2_Tensor = R2_Tensor.new() 
	
	return (2 * _metric.G_CONST * mass) / _metric.S_O_L ** 2

# Converts spherical coordinates to Cartesian coordinates.
# Inputs:
# - r: float, radial distance
# - theta: float, polar angle (0 to π)
# - phi: float, azimuthal angle (0 to 2π)
# Output:
# - Vector3 containing (x, y, z) Cartesian coordinates
func spherical_to_cartesian(r: float, theta: float, phi: float) -> Vector3:
	var x = r * sin(theta) * cos(phi)
	var y = r * sin(theta) * sin(phi)
	var z = r * cos(theta)
	return Vector3(x, y, z)

func _input(event):
	if event.is_action_pressed("evolve") && is_static == false:
		evolve_self()


func evolve_self():
	print("EVOLVING")
	define_closest_cells()
	
	var spatial_coords_of_closest_cells_simple : Array = []
	var spatial_coords_of_closest_cells_as_dictionary : Array = []
	var christoffel_of_closest_cells : Dictionary = {}
	
	var previous_position : Vector3 = s_coordinates
	
	#chose 8 because the closest points in like basically all directions are gonna amount to 8 points, but ill make it 16 :\
	for i in range(closest_cells.size()):
		var current_cell = closest_cells[i]
		var current_spatial_coord : Vector3 = current_cell.s_coordinates
		var current_christoffel : R3_Tensor = current_cell.christoffel_symbols
		
		spatial_coords_of_closest_cells_simple.append(current_spatial_coord)
		spatial_coords_of_closest_cells_as_dictionary.append({"x" : current_spatial_coord.x, "y" : current_spatial_coord.y, "z" : current_spatial_coord.z})
		christoffel_of_closest_cells[current_spatial_coord] = current_cell.christoffel_symbols
	
	lapse = lapse_function()
	
	s_coordinates = interpolate_spatial_coordinates(spatial_coords_of_closest_cells_simple, s_coordinates)
	
	s_velocity = define_velocity(previous_position, s_coordinates, lapse)
	
	
	var new_definitions : Dictionary = rk4_integration(s_coordinates, s_velocity, lapse, christoffel_of_closest_cells, spatial_coords_of_closest_cells_as_dictionary)
	
	s_coordinates = new_definitions["_position"]
	s_velocity = new_definitions["_velocity"]
	
	global_transform.origin = spherical_to_cartesian(s_coordinates.x, s_coordinates.y, s_coordinates.z)
	t_coordinate += lapse

func vec_to_dictionary(input_vector : Vector3) -> Dictionary:
	var new_dictionary : Dictionary = {"x" : input_vector.x, "y" : input_vector.y, "z" : input_vector.z}
	return new_dictionary

#region Geodesic Calculations

func interpolate_christoffel_symbols(
		christoffel: Dictionary,
		target_coords: Dictionary,
		grid: Array
	) -> R3_Tensor:
	# Initialize the interpolated tensor
	var interpolated : R3_Tensor = R3_Tensor.new()
	interpolated.R3_ZERO()
	
	# Find bounds for interpolation
	var bounds = {}
	for coord in target_coords.keys():
		var lower_bound = null
		var upper_bound = null
		for point in grid:
			if point[coord] <= target_coords[coord]:
				if lower_bound == null or point[coord] > lower_bound:
					lower_bound = point[coord]
			if point[coord] >= target_coords[coord]:
				if upper_bound == null or point[coord] < upper_bound:
					upper_bound = point[coord]
		if lower_bound == null or upper_bound == null:
			push_error("Failed to find valid bounds for interpolation.")
			return null
		bounds[coord] = [lower_bound, upper_bound]
	
	# Trilinear interpolation for each component of the Christoffel symbols
	for mu in Tensor.ARR_IND:
		for nu in Tensor.ARR_IND:
			for lambda_ in Tensor.ARR_IND:
				var corners = []
				var weights = []
				
				# Collect corner values and compute weights
				for x in bounds["x"]:
					for y in bounds["y"]:
						for z in bounds["z"]:
							var corner_coords = { "x": snapped(x, 6), "y": snapped(y, 6), "z": snapped(z, 6) }
							var found_key = false
							for key in christoffel.keys():
								if is_near_equal_dict(vec_to_dictionary(key), corner_coords):
									found_key = true
									break
							if not found_key:
								push_error("Christoffel symbol not defined at grid point: %s" % str(corner_coords))
								return null
							
							var corner_tensor = christoffel[corner_coords]
							var corner_value = corner_tensor.get_R3_value(lambda_, mu, nu)
							corners.append(corner_value)
							
							# Weight is inversely proportional to distance
							var weight = (
								(abs(x - target_coords["x"]) if x != target_coords["x"] else 1.0) *
								(abs(y - target_coords["y"]) if y != target_coords["y"] else 1.0) *
								(abs(z - target_coords["z"]) if z != target_coords["z"] else 1.0)
							)
							weights.append(weight)
				
				# Normalize weights and interpolate
				var total_weight = weights.sum()
				if total_weight == 0.0:
					push_error("Total weight is zero during interpolation.")
					return null
				var interpolated_value = 0.0
				for i in range(len(corners)):
					interpolated_value += corners[i] * (weights[i] / total_weight)
				
				# Set the interpolated value
				interpolated.set_R3_value(lambda_, mu, nu, interpolated_value)
	
	return interpolated

func is_near_equal_dict(dict1: Dictionary, dict2: Dictionary, tolerance: float = 1e-6) -> bool:
	for key in dict1.keys():
		if abs(dict1[key] - dict2[key]) > tolerance:
			return false
	return true


# Interpolates spatial coordinates using linear interpolation
# grid_points: Array of discrete spacetime points
# obj_position: Vector3 representing the massive object's current _position
func interpolate_spatial_coordinates(grid_points: Array, obj_position: Vector3) -> Vector3:
	# Identify the bounding grid points
	var lower_bound = Vector3()
	var upper_bound = Vector3()
	for point in grid_points:
		if point.x <= obj_position.x and point.y <= obj_position.y and point.z <= obj_position.z:
			lower_bound = point
		elif point.x > obj_position.x or point.y > obj_position.y or point.z > obj_position.z:
			upper_bound = point
			break
	# Perform linear interpolation
	var t_x = (obj_position.x - lower_bound.x) / (upper_bound.x - lower_bound.x)
	var t_y = (obj_position.y - lower_bound.y) / (upper_bound.y - lower_bound.y)
	var t_z = (obj_position.z - lower_bound.z) / (upper_bound.z - lower_bound.z)
	
	return Vector3(
		lower_bound.x + t_x * (upper_bound.x - lower_bound.x),
		lower_bound.y + t_y * (upper_bound.y - lower_bound.y),
		lower_bound.z + t_z * (upper_bound.z - lower_bound.z)
	)

func define_velocity(initial_position : Vector3, new_position, lapse : float) -> Vector3:
	var new_velocity : Vector3
	
	new_velocity = Vector3(
		(new_position.x - initial_position.x) / lapse, 
		(new_position.y - initial_position.y) / lapse, 
		(new_position.z - initial_position.z) / lapse)
	
	return new_velocity

# Adds a sum function to arrays for convenience
func sum(array: Array) -> float:
	var total = 0.0
	for element in array:
		total += element
	return total


# Returns a constant lapse function value
func lapse_function() -> float:
	return 0.1  # Replace with the real lapse function when ready

# Computes the acceleration based on Christoffel symbols
# _position: Vector3, current _position
# _velocity: Vector3, current _velocity
# christoffel_symbols: Dictionary, interpolated Christoffel symbols
# lapse: float, lapse function value
func compute_acceleration(_position: Vector3, _velocity: Vector3, christoffel_symbols: Dictionary, lapse: float) -> Vector3:
	var acceleration = Vector3()
	for mu in Tensor.ARR_IND:
		var sum_term = 0.0
		for nu in Tensor.ARR_IND:
			for lambda in Tensor.ARR_IND:
				var gamma = christoffel_symbols[mu + nu + lambda]
				sum_term += gamma * _velocity[nu.to_int()] * _velocity[lambda.to_int()]
		acceleration[mu.to_int()] = -sum_term / lapse
	return acceleration

# Performs RK4 integration to update the _position and _velocity of a massive object
# Inputs:
# - _position: Vector3, the object's current _position
# - _velocity: Vector3, the object's current _velocity
# - timestep: float, the timestep for integration
# - christoffel_grid: Dictionary containing the Christoffel symbols at grid points
# - grid_points: Array of discretized spacetime points
# - mass: float, the mass of the object (for later use, if needed)
func rk4_integration(_position: Vector3, _velocity: Vector3, input_lapse: float, christoffel_grid: Dictionary, grid_points: Array, mass: float = 1.0) -> Dictionary:
	# Interpolation helpers
	var christoffel_symbols = interpolate_christoffel_symbols(christoffel_grid, vec_to_dictionary(_position), grid_points)
	
	if christoffel_symbols == null:
		print("Null Christoffel")
		return {"_position": _position, "_velocity": _velocity}
	
	# Initial state
	var k1_position = _velocity
	var k1_velocity = compute_acceleration(_position, _velocity, christoffel_symbols, lapse)
	
	# Step 2: Compute k2
	var mid_position_k1 = _position + k1_position * (input_lapse / 2.0)
	var mid_velocity_k1 = _velocity + k1_velocity * (input_lapse / 2.0)
	var christoffel_k2 = interpolate_christoffel_symbols(christoffel_grid, vec_to_dictionary(mid_position_k1), grid_points)
	var k2_position = mid_velocity_k1
	var k2_velocity = compute_acceleration(mid_position_k1, mid_velocity_k1, christoffel_k2, lapse)
	
	# Step 3: Compute k3
	var mid_position_k2 = _position + k2_position * (input_lapse / 2.0)
	var mid_velocity_k2 = _velocity + k2_velocity * (input_lapse / 2.0)
	var christoffel_k3 = interpolate_christoffel_symbols(christoffel_grid, vec_to_dictionary(mid_position_k2), grid_points)
	var k3_position = mid_velocity_k2
	var k3_velocity = compute_acceleration(mid_position_k2, mid_velocity_k2, christoffel_k3, lapse)
	
	# Step 4: Compute k4
	var end_position_k3 = _position + k3_position * input_lapse
	var end_velocity_k3 = _velocity + k3_velocity * input_lapse
	var christoffel_k4 = interpolate_christoffel_symbols(christoffel_grid, vec_to_dictionary(end_position_k3), grid_points)
	var k4_position = end_velocity_k3
	var k4_velocity = compute_acceleration(end_position_k3, end_velocity_k3, christoffel_k4, lapse)
	
	# Final update
	var new_position = _position + input_lapse / 6.0 * (k1_position + 2.0 * k2_position + 2.0 * k3_position + k4_position)
	var new_velocity = _velocity + input_lapse / 6.0 * (k1_velocity + 2.0 * k2_velocity + 2.0 * k3_velocity + k4_velocity)	
	return {"_position": new_position, "_velocity": new_velocity}

#endregion 

#region Functions Reliant on Godot Architechture 

func define_closest_cells():
	#the culling
	for cell in closest_cells:
		if global_transform.origin.distance_squared_to(cell.global_transform.origin) > 4:
			closest_cells.erase(cell)
	#the sorting
	closest_cells.sort_custom(the_distance_sorter)

func the_distance_sorter(a, b) -> bool:
	if a.global_transform.origin.distance_squared_to(global_transform.origin) < b.global_transform.origin.distance_squared_to(global_transform.origin):
		return true
	return false

func _on_cell_updater_body_checker_area_entered(area):
	if area.is_in_group("Cell"):
		area.append_mass(self)
		closest_cells.append(area)
func _on_cell_updater_body_checker_area_exited(area):
	if area.is_in_group("Cell"):
		area.remove_mass(self)
		closest_cells.erase(area)

#endregion
