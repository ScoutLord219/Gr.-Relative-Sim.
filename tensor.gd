extends Resource

class_name Tensor

const EPI : float = 1e-6
const G_CONST : float = 6.6743e-11
const S_O_L : float = 299792458e-11
const ARR_IND : Array = ["t" , "r", "theta", "phi"]

var components : Dictionary = {}


#operational functions : facilitate operations perfomed between tensors or on tensor components
func set_value(mu : String, value : float, nu : String = "") -> void:
	if nu == "":
		components[mu] = value
	else:
		components[mu + nu] = value

func get_value(mu : String, nu : String = "") -> float:
	if nu == "":
		return components[mu]
	else:
		return components[mu + nu]

func self_additive(mu : String, nu : String, value : float) -> void:
	set_value(mu, get_value(mu, nu) + value, nu)

func self_multiply(mu : String, nu : String, value : float) -> void:
	set_value(mu, get_value(mu, nu) * value, nu)

#Utility functions : help in the ease and efficiency of creating new functions 
func duplicate_components() -> Dictionary:
	var tensor_duplicate : Dictionary = {}
	for mu in ARR_IND:
		for nu in ARR_IND:
			tensor_duplicate[mu + nu] = get_value(mu, nu)
	
	return tensor_duplicate

# Computes the finite difference of a parameter using the specified method
# `values` is an Array of the function's values at discrete points
# `t_values` is an Array of the corresponding independent variable values (time in this case)
# `method` specifies the finite difference method: "forward", "backward", or "central"
func finite_difference(values: Array, t_values: Array, index: int, method: String = "central") -> float:
	# Ensure valid inputs
	if index < 0 or index >= values.size():
		push_error("Index out of bounds for finite difference computation.")
		return 0.0
	# Handle forward difference (f(x+h) - f(x)) / h
	if method == "forward":
		if index + 1 >= values.size():
			push_error("Forward difference not possible at the end of the array.")
			return 0.0
		var h = t_values[index + 1] - t_values[index]
		return (values[index + 1] - values[index]) / h
		
	# Handle backward difference (f(x) - f(x-h)) / h
	elif method == "backward":
		if index - 1 < 0:
			push_error("Backward difference not possible at the start of the array.")
			return 0.0
		var h = t_values[index] - t_values[index - 1]
		return (values[index] - values[index - 1]) / h
		
	# Handle central difference (f(x+h) - f(x-h)) / (2h)
	elif method == "central":
		if index - 1 < 0 or index + 1 >= values.size():
			push_error("Central difference requires points on both sides.")
			return 0.0
		var h = t_values[index + 1] - t_values[index - 1]
		return (values[index + 1] - values[index - 1]) / h
		
	else:
		push_error("Invalid method for finite difference computation. Use 'forward', 'backward', or 'central'.")
		return 0.0
