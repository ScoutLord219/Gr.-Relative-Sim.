extends Tensor

class_name R4_Tensor


func R4_ZERO() -> void:
	for mu in ARR_IND:
		for nu in ARR_IND:
			for lambda in ARR_IND:
				for sigma in ARR_IND:
					set_R4_value(mu, nu, lambda, sigma, 0.0)

func get_R4_value(rho : String, sigma : String, mu : String, nu : String) -> float:
	return components[sigma + rho + mu + nu]

func set_R4_value(rho : String, sigma : String, mu : String, nu : String, value : float):
	components[sigma + rho + mu + nu] = value

func contract_R4_tensor(tensor : R4_Tensor) -> R2_Tensor:
	var result : R2_Tensor = R2_Tensor.new()
	
	for mu in ARR_IND:
		for nu in ARR_IND:
			var sum_terms = 0.0
			
			# Iterate over the index to be contracted
			for lambda in ARR_IND:
				sum_terms += tensor.get_R4_value(mu, lambda, lambda, nu)
			
			result.set_value(mu, sum_terms, nu)
	
	return result



