extends Area3D

class_name Cell

var four_position
var metric_tensor : Tensor = Tensor.new().metric_to_minkowski()
var metric_derivatives : Tensor = Tensor.new().ZERO(Tensor.MET_KEY_ARR)

func init(cell_dimensions : Vector2):
	pass
