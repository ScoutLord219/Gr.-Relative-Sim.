extends Node3D

class_name SPACETIME

const CELL = preload("res://scenes/cell.tscn")

const G_CONST : float = 6.6743e-11 #i had to remove the "e-11" from the g_const because godot wouldn't register such a small value
#in turn, since I basically multiplied the g_const by e11, i did the same with s_o_l. pretty sure thats how math works anyway
const S_O_L : float = 299792458

@export var cell_dimensions : Vector2 
@export var grid_size : int

signal all_cells_loaded

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
			
			current_cell.initializing_finished.connect(on_cell_initialized)
			current_cell.init_started(cell_dimensions, self)
			
			cell_array.append(current_cell)
		next_grid_position.x = 0
		next_grid_position.z -= cell_dimensions.y
		grid_z_position += 1
	global_position = Vector3(-grid_size/2, 0, grid_size/2)

func _process(delta):
	if all_cells_successfully_loaded:
		observer_time += delta

func on_cells_loaded():
	all_cells_successfully_loaded = true

func on_cell_initialized():
	initialized_cells += 1
	
	if initialized_cells == cell_array.size():
		all_cells_loaded.emit()
		on_cells_loaded()
	else:
		print("Failure to Load All Cells")
