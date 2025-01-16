extends Tensor

class_name R3_Tensor


func R3_ZERO() -> void:
	for mu in ARR_IND:
		for nu in ARR_IND:
			for lambda in ARR_IND:
				set_R3_value(mu, nu, lambda, 0.0)


func get_R3_value(lambda : String, mu : String, nu : String) -> float:
	return components[lambda + mu + nu]

func set_R3_value(lambda : String, mu : String, nu : String, value : float):
	components[lambda + mu + nu] = value

#Utility functions : help in the ease and efficiency of creating new functions 
func duplicate_R3_tensor() -> R3_Tensor:
	var tensor_duplicate : R3_Tensor = R3_Tensor.new()
	
	for mu in ARR_IND:
		for nu in ARR_IND:
			for sigma in ARR_IND:
				tensor_duplicate.set_R3_value(mu, nu, sigma, get_R3_value(mu, nu, sigma))
	
	return tensor_duplicate



# RIEMANNNN

func compute_riemann_tensor(chr: R3_Tensor) -> R4_Tensor:
	# Riemann tensor to be returned
	var riemann : R4_Tensor = R4_Tensor.new()
	
	for rho in ARR_IND:
		for sigma in ARR_IND:
			for mu in ARR_IND:
				for nu in ARR_IND:
					# Partial derivatives of Christoffel symbols
					var d_mu_chr = R3_partial_derivative(chr, rho, nu, sigma)
					var d_nu_chr = R3_partial_derivative(chr, rho, mu, sigma)
					
					# Christoffel product terms
					var gamma_mu_term = 0.0
					var gamma_nu_term = 0.0
					
					for lambda in ARR_IND:
						gamma_mu_term += chr.get_R3_value(rho, mu, lambda) * chr.get_R3_value(lambda, nu, sigma)
						gamma_nu_term += chr.get_R3_value(rho, nu, lambda) * chr.get_R3_value(lambda, mu, sigma)
						
					# Compute Riemann tensor component
					riemann.set_R4_value(rho, sigma, mu, nu, d_mu_chr - d_nu_chr + gamma_mu_term - gamma_nu_term)
	return riemann

func R3_partial_derivative(chr: R3_Tensor, rho: String, mu: String, nu: String) -> float:
	var h = EPI  # Small step for finite differences
	var incremented = chr.duplicate_R3_tensor()
	var decremented = chr.duplicate_R3_tensor()
	
	# Adjust the Christoffel symbols by a small step in the coordinate direction
	incremented.set_R3_value(rho, mu, nu, chr.get_R3_value(rho, mu, nu) + h)
	decremented.set_R3_value(rho, mu, nu, chr.get_R3_value(rho, mu, nu) - h)
	
	# Compute central difference approximation
	var value_incremented = incremented.get_R3_value(rho, mu, nu)
	var value_decremented = decremented.get_R3_value(rho, mu, nu)
	return (value_incremented - value_decremented) / (2 * h)

