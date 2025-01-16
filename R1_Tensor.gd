extends Tensor

class_name R1_Tensor


func R1_ZERO() -> void:
	for mu in ARR_IND:
		set_value(mu, 0.0)

func vec_3_to_R1_tensor(vector : Vector3, time_value : float = 0) -> R1_Tensor:
	var new_tensor : R1_Tensor
	new_tensor.set_value("t", time_value)
	new_tensor.set_value("r", vector.x)
	new_tensor.set_value("theta", vector.y)
	new_tensor.set_value("phi", vector.z)
	
	return new_tensor

func R1_tensor_to_vec_3() -> Vector3:
	var vector_version : Vector3
	vector_version = Vector3(get_value("r"), get_value("theta"), get_value("phi"))
	
	return vector_version

func spatial_seperation(position_1 : R1_Tensor, position_2 : R1_Tensor) -> Vector3:
	# Compute the spatial separation vector
	var separation_tensor = Tensor.new()
	var separation_vector : Vector3
	for mu in ["r", "theta", "phi"]:
		var separation_component = position_1.get_value(mu) - position_2.get_value(mu)
		separation_tensor.set_value(mu, separation_component)
		
	separation_vector = separation_tensor.R1_tensor_to_vec_3()
	return separation_vector

func norm_squared_spatial() -> float:
	return R1_tensor_to_vec_3().x ** 2 + R1_tensor_to_vec_3().y ** 2 + R1_tensor_to_vec_3().z ** 2
