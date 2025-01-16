extends Node3D

class_name SPACETIME

const CELL = preload("res://scenes/cell.tscn")

@export var grid_spacing : float 
@export var grid_size : Vector2

signal progress_simulation
signal evolve_simulation

var cell_array = [Node3D]
var current_cell

var metric_dictionary = { } 
# Vector3 : R2_Tensor

var observer_time : float = 0.0

var cells_computed : int = 0
var all_cells_computed : bool = false

func _ready():
	create_2d_grid(Vector3.ZERO, grid_size, grid_spacing)

func _input(event):
	if event.is_action_pressed("progress"):
		progress_simulation.emit()
		cells_computed = 0
	elif event.is_action_pressed("evolve") && all_cells_computed:
		evolve_simulation.emit()



# Function to create a grid centered on the Cartesian origin
func create_2d_grid(center: Vector3, dimensions: Vector2, cell_size: float) -> Array:
	#Creates a 2D grid of points in the x-z plane, centered at the specified origin.
	#
	#Parameters:
	#- center: Vector3 - The center of the grid (typically (0, 0, 0)).
	#- dimensions: Vector2 - Number of cells in the x and z directions.
	#- cell_size: float - Distance between points in the grid.
	#
	#Returns:
	#- Array of grid points as Vector3 positions in the x-z plane.
	
	var grid_points = []
	
	# Calculate half-dimensions for centering
	var half_dims = dimensions * cell_size * 0.5
	
	var cell_count : int = 0 
	
	# Iterate over x and z dimensions
	for x in range(dimensions.x):
		for z in range(dimensions.y):  # Using dimensions.y for z size
			# Compute the position for each grid point
			var _position = Vector3(
				center.x + (x * cell_size) - half_dims.x + (cell_size * 0.5),
				center.y,  # Always at the y-coordinate of the center
				center.z + (z * cell_size) - half_dims.y + (cell_size * 0.5)
			)
			
			current_cell = CELL.instantiate()
			self.add_child(current_cell)
			
			current_cell.global_transform.origin = _position
			
			current_cell.metric_computed.connect(check_all_cells_computed)
			current_cell.init(cell_size)
			cell_count += 1
			
			print("Initialized Cells : ", cell_count)
			
			grid_points.append(position)
	
	return grid_points

func check_all_cells_computed():
	cells_computed += 1
	print("Cells Computed : ", cells_computed)
	if cells_computed == grid_size.x * grid_size.y:
		all_cells_computed = true
