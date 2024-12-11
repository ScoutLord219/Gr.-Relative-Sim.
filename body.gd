extends CharacterBody3D

class_name Body

@export var mass : float
@export var radius : float

@onready var four_position : Tensor = Tensor.new()
@onready var four_velocity : Tensor = Tensor.new()

@onready var current_relevant_metric : Tensor = Tensor.new()

@onready var christoffel : Tensor = Tensor.new()

@onready var cell_updater_shape = $"Cell Updater_Body Checker/Cell Updater Shape"

var affecting_bodies_and_christoffels = {}




var proper_time : float = 0.0

# Integration step size
var step_size : float = 0.01


func _ready():
	var coord_to_four_pos : Tensor = translate_to_four_position(Vector3.ZERO, global_transform.origin)
	
	four_position.set_multiple_components(coord_to_four_pos.components.keys(), coord_to_four_pos.components.values())
	four_velocity.set_multiple_components(Tensor.FOUR_VEC_KEY_ARR, [1, 0, 0, 0])
	current_relevant_metric.metric_to_minkowski()
	christoffel.ZERO(Tensor.CHR_KEY_ARR)
	
	cell_updater_shape.scale *= radius
	# Initialize the Minkowski metric (flat spacetime)
	var minkowski_metric = current_relevant_metric

	# Normalize the velocity using the Minkowski metric
	four_velocity = four_velocity.normalize_four_velocity(four_velocity, minkowski_metric)
	
	# Combine initial position and velocity for debugging or initialization
	#print("Initial Four-Position:", four_position)
	#print("Normalized Four-Velocity:", four_velocity)

func _physics_process(delta):
	#account for multiple massive objects & define effective r
	var coord_to_four_pos : Tensor
	if affecting_bodies_and_christoffels.keys().size() == 0:
		coord_to_four_pos = translate_to_four_position(Vector3.ZERO, global_transform.origin)
	else:
		#put a function here that finds the effecive_four_position using the calculate_effetive_r method, and just use that as the ref point
		coord_to_four_pos = translate_to_four_position(translate_spherical_to_cartesian(calculate_effective_reference_point(four_position, affecting_bodies_and_christoffels.keys(), current_relevant_metric, affecting_bodies_and_christoffels.values())), translate_spherical_to_cartesian(four_position))
		#coord_to_four_pos = translate_to_four_position()
	
	# Flat spacetime projection (initial guess)
	var projected_position = project_flat_spacetime_trajectory(four_position, four_velocity, step_size, step_size, 50, coord_to_four_pos.get_components(["r"]), affecting_bodies_and_christoffels.keys())
	
	#for four_pos in projected_position:
		#print("Time Position: ", four_pos.get_components(["t"]))
		#pass
	# Metric tensor at current position
	var metric_tensor = calculate_metric_tensor(four_position, affecting_bodies_and_christoffels.keys())
	var metric_derivatives = calculate_metric_derivatives(metric_tensor, four_position, affecting_bodies_and_christoffels.keys(), coord_to_four_pos.get_components(["r"]))
	#print(metric_derivatives)
	
	# Solve geodesic equation using Runge-Kutta
	integrate_geodesic(four_position, four_velocity, metric_tensor, affecting_bodies_and_christoffels.keys())
	
	# Update proper time
	proper_time += step_size



func calculate_effective_reference_point(
		body_position: Tensor,
		input_affecting_bodies: Array,
		metric_tensor: Tensor,  # 4x4 Metric tensor for curved spacetime
		Christoffel_symbols: Array  # Array of Christoffel symbols for each body
	) -> Tensor:
	"""
	Calculates an effective reference point considering the curvature effects of multiple bodies
	under general relativity, incorporating the metric tensor and Christoffel symbols.

	Parameters:
	- body_position: Tensor representing the current four-position of the body (reference point).
	- affecting_bodies_and_christoffels.keys(): Array of bodies affecting the reference point, each with mass, position, and velocity.
	- metric_tensor: The 4x4 tensor representing the spacetime metric.
	- Christoffel_symbols: Array containing the Christoffel symbols for each body.

	Returns:
	- Tensor representing the effective reference four-position.
	"""
	var effective_reference = Tensor.new()
	effective_reference.ZERO(Tensor.FOUR_VEC_KEY_ARR)

	var total_mass = 0.0
	var weighted_position = Tensor.new()
	weighted_position.ZERO(Tensor.FOUR_VEC_KEY_ARR)

	# 1. Loop over each affecting body and calculate their contribution
	for other_body in affecting_bodies_and_christoffels.keys():
		var other_position: Tensor = other_body.get("four_position")
		var other_mass: float = other_body.get("mass")

		# Calculate the spatial separation vector
		var spatial_separation = Tensor.new()
		for mu in ["r", "phi", "theta"]:
			var separation_component = other_position.get_components([mu]) - body_position.get_components([mu])
			spatial_separation.set_component([mu], separation_component)
		var distance_squared = spatial_separation.norm_squared_spatial()

		# Avoid division by zero or very small distances
		if distance_squared < Tensor.SMALL_VALUE:
			continue  # Skip singularities (overlapping bodies)

		var distance = sqrt(distance_squared)

		# 2. Gravitational potential approximation (1/r^2 law)
		var gravitational_influence = other_mass / distance_squared

		# 3. Weight the position by gravitational influence and sum
		for mu in ["r", "phi", "theta"]:
			var weighted_component = spatial_separation.get_components([mu]) * gravitational_influence
			weighted_position.additive_on_component([mu], weighted_component)

		# Accumulate total mass for normalization
		total_mass += gravitational_influence

	# 4. Normalize the weighted position by the total gravitational influence
	if total_mass > 0.0:
		for mu in ["r", "phi", "theta"]:
			effective_reference.set_component([mu], weighted_position.get_components([mu]) / total_mass)

	# 5. Compute the time component based on the influence of the masses
	effective_reference.set_component(["t"], body_position.get_components(["t"]))

	# 6. Calculate the four-velocity and four-acceleration (via geodesic equation)
	var four_velocity = Tensor.new()
	# Here, compute the four-velocity using the geodesic equation, i.e., 
	# proper time derivatives, and incorporating Christoffel symbols for curvature.
	# The geodesic equation can be numerically integrated for this purpose.

	# 7. Return the effective reference four-position, accounting for relativistic corrections
	return effective_reference


func project_flat_spacetime_trajectory( 
		initial_four_position: Tensor,
		initial_four_velocity: Tensor, 
		input_proper_time: float, 
		delta_tau: float, 
		num_steps: int,
		effective_r : float,
		input_affecting_bodies: Array
	) -> Array:
	"""
	Projects the trajectory of a body in flat or slightly curved spacetime.

	Parameters:
	- initial_four_position: Tensor with 4 components [t, x, y, z].
	- initial_four_velocity: Tensor with 4 components [v_t, v_x, v_y, v_z].
	- input_proper_time: Initial proper time.
	- delta_tau: Step size for proper time.
	- num_steps: Number of steps to project the trajectory.
	- affecting_bodies_and_christoffels.keys(): Array of nearby bodies influencing the trajectory (optional).

	Returns:
	- trajectory: Array of Tensor spacetime points [(t, x, y, z), ...].
	"""
	var _position: Tensor = initial_four_position.clone()
	var _velocity: Tensor = initial_four_velocity.clone()
	var current_tau: float = input_proper_time
	
	# Normalize the initial velocity for Minkowski metric
	var norm = -_velocity.get_components(["t"]) ** 2 + _velocity.get_components(["r"]) ** 2 + \
			   _velocity.get_components(["phi"]) ** 2 + _velocity.get_components(["theta"]) ** 2
	if abs(norm + 1.0) > Tensor.SMALL_VALUE:  # Ensure normalization for timelike trajectories
		_velocity = normalize_velocity_minkowski(_velocity)
		
	var trajectory = []  # Store projected points
	
	#var trajectory_line =  Line2D.new()
	#trajectory_line.add_point(Vector2(_position["r"]))
	
	for step in range(num_steps):
		# Append a deep copy of the current position to the trajectory
		var position_copy: Tensor = Tensor.new()
		position_copy.assign_all_tensor_components_to_other(_position)
		trajectory.append(position_copy)
		
		# Calculate the effective acceleration if affecting bodies are present
		var effective_acceleration: Tensor = Tensor.new()
		
		
		if affecting_bodies_and_christoffels.keys().size() > 0:
			effective_acceleration = calculate_effective_acceleration(_position, _velocity, current_relevant_metric, input_affecting_bodies)
		else:
			effective_acceleration.ZERO(Tensor.FOUR_VEC_KEY_ARR)
		
		# Update velocity based on acceleration (if any)
		for mu in Tensor.FOUR_VEC_KEY_ARR:
			var acc_component = effective_acceleration.get_components(mu)
			_velocity.additive_on_component(mu, acc_component * delta_tau)
		
		# Update position using the updated velocity
		for mu in Tensor.FOUR_VEC_KEY_ARR:
			var velocity_component = _velocity.get_components(mu)
			_position.additive_on_component(mu, velocity_component * delta_tau)
		
		# Re-normalize velocity to prevent numerical drift
		_velocity = normalize_velocity_minkowski(_velocity)
		
		# Update proper time
		current_tau += delta_tau
	
	return trajectory

func project_flat_spacetime(_position, _velocity, delta):
	# Simple Minkowski metric projection
	return _position + _velocity * delta

func normalize_velocity_minkowski(_velocity):
	#"""
	#Normalize a _velocity vector in Minkowski spacetime (timelike normalization).
	#"""
	var minkowski_norm = -_velocity.get_components(["t"]) ** 2 + _velocity.get_components(["r"]) ** 2 + _velocity.get_components(["phi"]) ** 2 + _velocity.get_components(["theta"]) ** 2
	if minkowski_norm < 0:
		var _scale = 1.0 / sqrt(-minkowski_norm)
		_velocity.scale(_velocity, _scale)
	return _velocity


#func integrate_geodesic(position, velocity, christoffel, step):
	## Use Runge-Kutta or another numerical solver for integration
	#return RungeKuttaSolver.solve(position, velocity, christoffel, step)



func calculate_metric_tensor(input_four_position: Tensor, input_affecting_bodies: Array) -> Tensor:
	# Start with the Minkowski metric (flat spacetime)
	var metric_tensor : Tensor = Tensor.new()
	
	metric_tensor.metric_to_minkowski()
	
	var effective_r = input_four_position.get_components(["r"])
	
	
	var gm_over_r = (Tensor.G_CONST * calculate_effective_mass(input_affecting_bodies)) / effective_r
	var one_minus_two_gm_over_r = 1 - 2 * gm_over_r
	# Loop through each body
	for body in input_affecting_bodies:
		# Calculate the effective radial distance (r) to the body
		# Ensure a minimum radial distance to avoid singularities
		# Compute Schwarzschild metric components for this body
		# Adjust the metric components
		metric_tensor.additive_on_component(["t", "t"], -(1 - one_minus_two_gm_over_r))  # Time-time component
		metric_tensor.additive_on_component(["r", "r"], 1 / one_minus_two_gm_over_r)      # Radial component
		metric_tensor.additive_on_component(["phi", "phi"], effective_r ** 2 * sin(input_four_position.get_components(["phi"])) ** 2)  # Angular phi component
		metric_tensor.additive_on_component(["theta", "theta"], effective_r ** 2)                # Angular theta component
	
	return metric_tensor

func calculate_effective_r(input_four_position: Tensor, metric : Tensor, input_affecting_bodies : Array) -> float:
	# Compute the spatial distance (ignoring time component) between the two positions
	var total_effective_r : float = 0.0
	
	
	for body in input_affecting_bodies:
		var spatial_distance = input_four_position.subtract(input_four_position, body.four_position).magnitude_with_metric(metric)
		
		# Compute Schwarzschild-like metric (this is a simplified version for clarity)
		# Adjust distance based on the body's mass and the radial distance to avoid singularities
		var r = max(spatial_distance, Tensor.SMALL_VALUE)  # Avoid division by zero or very small r
		var gm_over_r = (Tensor.G_CONST * body.mass) / r
		var schwarzschild_factor = 1 - 2 * gm_over_r
		
		# Ensure the Schwarzschild factor is valid (avoid division by zero or negative roots)
		schwarzschild_factor = max(schwarzschild_factor, Tensor.SMALL_VALUE)
		# Effective radial distance using the full metric
		total_effective_r += r * (schwarzschild_factor) ** -0.5
		
	return total_effective_r

func calculate_effective_mass(input_affecting_bodies : Array) -> float : 
	var sum_of_masses : float = 0.0
	for body in input_affecting_bodies:
		sum_of_masses += body.mass
	return sum_of_masses

func calculate_metric_change(
		input_four_position: Tensor, 
		mu: String, 
		nu: String, 
		coordinate: String, 
		input_affecting_bodies: Array
	) -> float:
	#"""
	#Approximates the derivative of a metric tensor component with respect to a coordinate.
#
	#Parameters:
		#input_four_position (Tensor): The current four-position of the body.
		#mu, nu (String): Indices of the metric component to calculate the derivative for.
		#coordinate (String): The coordinate to differentiate with respect to (e.g., "r", "theta").
		#affecting_bodies_and_christoffels.keys() (Array): Array of bodies affecting the spacetime.
		#step_size (float): The small step value used for finite difference approximation.
#
	#Returns:
		#float: The approximate partial derivative of the metric component.
	#"""
	# Clone the input position to create perturbed positions
	var position_plus = input_four_position
	var position_minus = input_four_position
	
	# Adjust the specified coordinate by a small step forward and backward
	position_plus.additive_on_component([coordinate], step_size)
	position_minus.additive_on_component([coordinate], -step_size)
	
	# Calculate the metric tensor at the perturbed positions
	var metric_plus = calculate_metric_tensor(position_plus, input_affecting_bodies)
	var metric_minus = calculate_metric_tensor(position_minus, input_affecting_bodies)
	
	# Extract the specific component (mu, nu) of the metric tensor at both positions
	var metric_plus_component = metric_plus.get_components([mu, nu])
	var metric_minus_component = metric_minus.get_components([mu, nu])
	
	# Compute the finite difference approximation for the derivative
	return (metric_plus_component - metric_minus_component) / (2 * step_size)

func calculate_metric_derivatives(input_metric : Tensor, input_four_position : Tensor, input_affecting_bodies : Array, effective_r : float) -> Dictionary:
	# Compute partial derivatives of metric components (example for "t,t")
	var metric_derivatives : Tensor = Tensor.new()
	
	# Initialize metric_derivatives to zero for Minkowski spacetime
	metric_derivatives.set_multiple_components(Tensor.CHR_KEY_ARR, [0.0, 0.0, 0.0, 0.0])
	
	for mu in Tensor.FOUR_VEC_KEY_ARR:
		for nu in Tensor.FOUR_VEC_KEY_ARR:
			#Derivative with respect to "r"
			if nu == mu:
				var partial_r = (input_metric.get_components([mu[0], nu[0]]) - calculate_metric_change(input_four_position, mu[0], nu[0], "r", input_affecting_bodies)) / effective_r
				metric_derivatives.set_component([mu[0], nu[0], "r"], partial_r)
	
	return metric_derivatives.components



func integrate_geodesic(
	input_four_position: Tensor,
	input_four_velocity: Tensor,
	metric_tensor: Tensor,
	affecting_bodies: Array
) -> void:
	"""
	Integrates the geodesic equation using the RK4 method with curved spacetime.	
	Parameters:
	- input_four_position: Current four-position (Tensor).
	- input_four_velocity: Current four-velocity (Tensor).
	- step_size: Proper time step size (float).
	- metric_tensor: The metric tensor to calculate Christoffel symbols.
	- affecting_bodies: Array of nearby massive bodies influencing the motion.	
	Modifies:
	- input_four_position and input_four_velocity in-place.
	"""

	# RK4 Step Calculations
	var k1_velocity = calculate_effective_acceleration(
		input_four_position,
		input_four_velocity,
		metric_tensor,
		affecting_bodies
	)
	var k1_position = input_four_velocity.clone()	
	var temp_position = input_four_position.clone()
	var temp_velocity = input_four_velocity.clone()	
	var k2_position = Tensor.new()
	var k2_velocity = Tensor.new()

	var letters = ["t", "r", "phi", "theta"]

	for mu in letters:
		k2_position.set_component([mu], k1_position.get_components([mu]) + 0.5 * step_size * k1_velocity.get_components([mu]))
		temp_position.set_component([mu], input_four_position.get_components([mu]) + 0.5 * step_size * k1_position.get_components([mu]))
		temp_velocity.set_component([mu], input_four_velocity.get_components([mu]) + 0.5 * step_size * k1_velocity.get_components([mu]))

	k2_velocity = calculate_effective_acceleration(temp_position, temp_velocity, metric_tensor, affecting_bodies)	
	var k3_position = Tensor.new()
	var k3_velocity = Tensor.new()

	for mu in letters:
		k3_position.set_component([mu], k2_position.get_components([mu]) + 0.5 * step_size * k2_velocity.get_components([mu]))
		temp_position.set_component([mu], input_four_position.get_components([mu]) + 0.5 * step_size * k2_position.get_components([mu]))
		temp_velocity.set_component([mu], input_four_velocity.get_components([mu]) + 0.5 * step_size * k2_velocity.get_components([mu]))

	k3_velocity = calculate_effective_acceleration(temp_position, temp_velocity, metric_tensor, affecting_bodies)	
	var k4_position = Tensor.new()
	var k4_velocity = Tensor.new()

	for mu in letters:
		k4_position.set_component([mu], k3_position.get_components([mu]) + step_size * k3_velocity.get_components([mu]))
		temp_position.set_component([mu], input_four_position.get_components([mu]) + step_size * k3_position.get_components([mu]))
		temp_velocity.set_component([mu], input_four_velocity.get_components([mu]) + step_size * k3_velocity.get_components([mu]))

	k4_velocity = calculate_effective_acceleration(temp_position, temp_velocity, metric_tensor, affecting_bodies)

	# Final update to four-position and four-velocity
	for mu in letters:
		var position_update = (1.0 / 6.0) * step_size * (
			k1_position.get_components([mu]) + 2.0 * k2_position.get_components([mu]) +
			2.0 * k3_position.get_components([mu]) + k4_position.get_components([mu])
		)
		var velocity_update = (1.0 / 6.0) * step_size * (
			k1_velocity.get_components([mu]) + 2.0 * k2_velocity.get_components([mu]) +
			2.0 * k3_velocity.get_components([mu]) + k4_velocity.get_components([mu])
		)

		input_four_position.additive_on_component([mu], position_update)
		input_four_velocity.additive_on_component([mu], velocity_update)

	update_position_in_scene(input_four_position)
	update_velocity_in_scene(input_four_velocity)


# Helper function to compute acceleration from Christoffel symbols
func update_position_in_scene(input_four_position: Tensor) -> void:
	global_transform.origin = translate_spherical_to_cartesian(input_four_position)
	
	# Update the position of the Node3D

func translate_spherical_to_cartesian(input_four_position : Tensor) -> Vector3:
	var r = input_four_position.get_components(["r"])  # Radial distance
	var phi = input_four_position.get_components(["phi"])  # Azimuthal angle
	var theta = input_four_position.get_components(["theta"])  # Polar angle
	
	# Convert spherical coordinates to Cartesian coordinates for Godot
	# Spherical to Cartesian:
	# x = r * sin(theta) * cos(phi)
	# y = r * sin(theta) * sin(phi)
	# z = r * cos(theta)
	var x = r * sin(theta) * cos(phi)
	var y = 0
	var z = r * cos(theta)
	
	return Vector3(x, y, z)

func update_velocity_in_scene(input_four_velocity : Tensor):
	var v_r = input_four_velocity.get_components(["r"])
	var v_phi = input_four_velocity.get_components(["phi"])
	var v_theta = input_four_velocity.get_components(["theta"])
	
	# Convert spherical velocity components to Cartesian
	var v_x = v_r * sin(v_theta) * cos(v_phi) + v_r * cos(v_theta) * cos(v_phi) * v_theta - v_r * sin(v_theta) * sin(v_phi) * v_phi
	var v_y = 0#v_r * sin(v_theta) * sin(v_phi+ v_r * cos(v_theta) * sin(v_phi)) * v_theta + v_r * sin(v_theta) * cos(v_phi) * v_phi
	var v_z = v_r * cos(v_theta) - v_r * sin(v_theta) * v_theta
	
	
	velocity = Vector3(v_x, v_y, v_z)

func translate_to_four_position(reference_point: Vector3, object_position: Vector3) -> Tensor:
	# """
	# Translates the object's current position into a four-position Tensor
	# relative to a given reference point.
	#
	# Parameters:
	# - reference_point: Vector3 representing the origin/reference point in space.
	# - object_position: Vector3 representing the object's current position in Cartesian coordinates.
	#
	# Returns:
	# - four_position: Tensor with components [t, r, phi, theta].
	# """
	
	# Calculate displacement vector
	var displacement = object_position - reference_point
	
	# Convert Cartesian to spherical coordinates
	var r = displacement.length()  # Radial distance
	var theta = acos(displacement.z / r) if r > 0 else 0.0  # Polar angle
	var phi = atan2(displacement.y, displacement.x)  # Azimuthal angle
	
	# Create a new four-position Tensor
	var new_four_position = Tensor.new()
	
	# Set time component (initialize to 0; can adjust based on your needs)
	new_four_position.set_component(["t"], 0)
	
	# Set spatial components
	new_four_position.set_component(["r"], r)
	new_four_position.set_component(["phi"], phi)
	new_four_position.set_component(["theta"], theta)
	
	return four_position



func calculate_effective_acceleration(
	body_position: Tensor,
	body_velocity: Tensor,
	metric_tensor: Tensor,
	affecting_bodies_array : Array
) -> Tensor:
	var effective_r = 0.0
	
	if affecting_bodies_array.size() > 0:
		effective_r = Tensor.SMALL_VALUE
	else:
		effective_r = body_position.get_components(["r"])
	
	var effective_acceleration = Tensor.new()
	effective_acceleration.ZERO(Tensor.FOUR_VEC_KEY_ARR)	
	# 1. Add acceleration due to curvature via Christoffel symbols
	for mu in Tensor.FOUR_VEC_KEY_ARR:  # Iterate over spacetime dimensions
		var acc_component = 0.0
		for nu in Tensor.FOUR_VEC_KEY_ARR:
			for lambda in Tensor.FOUR_VEC_KEY_ARR:
				# Get the Christoffel symbol for this specific (mu, nu, lambda) set
				var christoffel = get_christoffel_symbol(mu[0], nu[0], lambda[0], body_position, affecting_bodies_array, effective_r)
				acc_component -= christoffel * body_velocity.get_components(nu) * body_velocity.get_components(lambda)
		effective_acceleration.set_component(mu, acc_component)	
	# 2. Add contributions from nearby bodies (optional, if dynamic gravitational sources exist)
	for other_body in affecting_bodies_and_christoffels:
		var other_position: Tensor = other_body.get("four_position")
		var other_mass: float = other_body.get("mass")	
		# Calculate spatial separation and distance
		var spatial_separation = Tensor.new()
		for mu in ["r", "phi", "theta"]:
			var separation_component = other_position.get_components([mu]) - body_position.get_components([mu])
			spatial_separation.set_component([mu], separation_component)
		var distance_squared = spatial_separation.norm_squared_spatial()
		if distance_squared < Tensor.SMALL_VALUE:
			continue  # Avoid singularities
		
		var distance = sqrt(distance_squared)	
		# Approximation: Add Newtonian term or use metric perturbations (if applicable)
		for mu in ["r", "phi", "theta"]:
			var acc_component = - (Tensor.G_CONST * other_mass / distance_squared) * \
								(spatial_separation.get_components([mu]) / distance)
			effective_acceleration.additive_on_component([mu], acc_component)	
	return effective_acceleration

# Helper Method 1: Get Christoffel Symbol
func get_christoffel_symbol(mu: String, nu: String, lambda: String, input_four_position: Tensor, affecting_bodies_array : Array, effective_r : float) -> float:
	# Retrieve the metric tensor at the given input_four_position
	var metric = calculate_metric_tensor(input_four_position, affecting_bodies_array)
	
	# Calculate the inverse of the metric tensor
	var inverse_metric = get_inverse_metric(metric)
	
	# Compute the partial derivatives of the metric tensor components with respect to the coordinates
	var metric_derivatives = calculate_metric_derivatives(metric, input_four_position, affecting_bodies_array, effective_r)
	
	# Initialize the Christoffel symbol
	var christoffel_symbol = 0.0
	
	# Sum over rho for the Christoffel symbol formula
	for sigma in Tensor.FOUR_VEC_KEY_ARR:
		# Calculate Christoffel symbols using the metric tensor and its derivatives
		var term_1 = metric_derivatives.get([nu, sigma[0], mu], 0.0)
		if term_1 != 0:
			print("TERM 1 :", term_1)
		var term_2 = metric_derivatives.get([sigma[0], mu, nu], 0.0)
		if term_2 != 0:
			print("TERM 2 :",term_2)
		var term_3 = metric_derivatives.get([mu, nu, sigma[0]], 0.0)
		if term_3 != 0:
			print("TERM 3 :",term_3)
		
		if lambda == sigma[0]:
			christoffel_symbol += inverse_metric.get_components([lambda, sigma[0]]) * (term_1 + term_2 - term_3)
		
		if christoffel_symbol != 0:
			#print(metric_tensor.get_components([mu[0], mu[0]]))
			print("Metric Inverse @ [", mu, ", ", mu, "] : ", inverse_metric)
			#print("Christoffel Symbol - ", [mu, nu, sigma], " : ", christoffel_value)
	  
	christoffel_symbol *= 0.5
	
	return christoffel_symbol

# Helper Method 3: Get Inverse Metric (Placeholder)
func get_inverse_metric(metric: Tensor) -> Tensor:
	var inverse_metric = metric.gaussian_elimination_invert()
	# Compute the inverse of the metric tensor here
	return inverse_metric

# Helper Method 4: Get Metric Derivatives (Placeholder)
func get_metric_derivatives(input_four_position: Tensor, metric_tensor : Tensor, affecting_bodies_array : Array, effective_r : float) -> Dictionary:
	#"""
	#Returns the partial derivatives of the metric tensor components at a given 	
	#Parameters:
	#- input_four_position: The four-input_four_position tensor representing the current location in 	
	#Returns:
	#- Dictionary where each key is a tuple of (mu, nu, lambda) and the value is the derivative of g_{\mu\nu}.
	#"""
	var derivatives = calculate_metric_derivatives(metric_tensor, input_four_position, affecting_bodies_array, effective_r)
	# Compute the partial derivatives of the metric tensor here
	return derivatives




func _on_cell_updater_body_checker_body_entered(body):
	if body.is_in_group("Body"):
		affecting_bodies_and_christoffels[body] = body.christoffel
		affecting_bodies_and_christoffels.keys().push_front(body)


func _on_cell_updater_body_checker_body_exited(body):
	if body.is_in_group("Body"):
		affecting_bodies_and_christoffels[body] = null
		affecting_bodies_and_christoffels.erase(body)
