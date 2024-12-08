extends Resource

class_name Tensor

var components : Dictionary

func set_components(indices : Array, value : float):
	components[indices] = value

func get_components(indices : Array) -> float:
	return components[indices]

func add(tensor_a: Tensor, tensor_b: Tensor) -> Tensor:
	var result = Tensor.new()
	for key in tensor_a.components.keys():
		if tensor_b.components.has(key):
			result.set_components([key], tensor_a.get_components([key]) + tensor_b.get_components([key]))
	return result

func subtract(tensor_a: Tensor, tensor_b: Tensor) -> Tensor:
	var result = Tensor.new()
	for key in tensor_a.components.keys():
		if tensor_b.components.has(key):
			result.set_components([key], tensor_a.get_components([key]) - tensor_b.get_components([key]))
	return result

func scale(tensor: Tensor, scalar: float) -> Tensor:
	var result = Tensor.new()
	for key in tensor.components.keys():
		result.set_components([key], tensor.get_components([key]) * scalar)
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
	var sum = 0.0
	for key in components.keys():
		if metric.components.has(key):
			sum += metric.components[key] * pow(components[key], 2)
	return sqrt(sum)

func metric_based_distance(other: Tensor, metric: Tensor) -> float:
	# Ensure both tensors have the same keys
	var keys = components.keys()
	
	if other.components.keys() != keys:
		push_error("Tensors must have the same keys to compute distance.")
		return -1.0  # Return an error value
	
	var sum = 0.0
	
	# If a metric is provided, use it for distance calculation
	if metric.components.size() > 0:
		for mu in keys:
			for nu in keys:
				# Compute the metric-weighted distance
				var delta_i = components[mu] - other.components[mu]
				var delta_j = components[nu] - other.components[nu]
				var g_ij 
				
				if mu[0] == nu[0]:
					g_ij = metric.get_components([mu[0], nu[0]])
				else:
					g_ij = 0.0
				
				sum += g_ij * delta_i * delta_j
	else:
	   # Euclidean distance if no metric is provided
		for key in keys:
			var delta = components[key] - other.components[key]
			sum += pow(delta, 2)
	
	
	return sqrt(sum)
