extends Node3D

class_name SPACETIME

const CELL = preload("res://scenes/cell.tscn")

@export var cell_dimensions : Vector2 
@export var grid_size : int

var grid = {Vector3(cell_dimensions.x, 0, cell_dimensions.y): CELL}
var cell_array = [Cell]
var next_grid_position : Vector3
var current_cell

var observer_time : float

var initialized_cells : int = 0

var all_cells_successfully_loaded : bool = false

func _ready():
	for grid_z_position in grid_size:
		for grid_x_position in grid_size:
			current_cell = CELL.instantiate()
			self.add_child(current_cell)
			
			current_cell.global_position = next_grid_position
			grid[current_cell.global_position] = current_cell
			next_grid_position.x += cell_dimensions.x
			
			current_cell.init(cell_dimensions)
			
			cell_array.append(current_cell)
		next_grid_position.x = 0
		next_grid_position.z -= cell_dimensions.y
		grid_z_position += 1
	global_position = Vector3(-grid_size/2, 0, grid_size/2)

