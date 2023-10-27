import numpy as np
import copy


class Cell:
	"""A (almost) general cell class. Keeps track on the cell's mother and children,
	its expression levels from birth to division and division time. Can also include
	methylation, but then the class is not completely general and must be adapted a 
	bit depending on model."""
	def __init__(self, ID, mother, t_div, x_0=[], y_0=[]):
		# If not the original cell, i.e. if there exists a mother.
		# The cell's initial values are the mother's final. 
		if(ID > 0):
			self.mother = mother.ID				# The mother's identification number
			self.gen = mother.gen + 1			# Generation number
			self.t_0 = mother.t_div				# Birth time of cell = division time of mother
			self.x_0 = mother.x_final			# Initial expression levels
			self.t_div = t_div + mother.t_div	# Time when cell divides
			
			# Methylation inheritage:
			if(len(mother.y_final) == 2):
				self.y_0 = mother.y_final		# Initial methylation levels
			elif(len(mother.y_final) == 3):
				# Only apply special division rules if X is absent.
				# C -> I, I -> 50% I + 50% O, O -> O
				if(self.x_0[4] <= 0):
					self.y_0 = np.zeros(3)
					self.y_0[1] = mother.y_final[0]
					self.y_0[2] = mother.y_final[2]
					random = np.random.rand(int(mother.y_final[1]))
					for r in random:
						if(r < 0.5):
							self.y_0[1] += 1
						else:
							self.y_0[2] += 1
					
				else:
					self.y_0 = mother.y_final	# Initial methylation levels
			else:
					self.y_0 = mother.y_final	# Initial methylation levels
					print(len(mother.y_final))

		# If the original cell.
		else:
			self.mother = -1
			self.gen = 0
			self.t_0 = 0
			self.x_0 = x_0
			self.y_0 = y_0
			self.t_div = t_div
			self.gen_i = 0	# The cell's indexing number within its generation
		
		self.ID = ID				# Cell's identification number
		self.t = np.array([self.t_0])			# Time series
		self.x = np.array([np.array(self.x_0)])		# Expression level at every timestep
		self.x_final = self.x_0		# The final expression levels
		
		self.t_meth = np.array([self.t_0])	# Time series methylation
		self.y = np.array([np.array(self.y_0)])		# Methylation level at every timestep
		self.y_final = self.y		# The final methylation levels
		
		self.children = [-1, -1]
		
	
	def set_division_time(self, t):
		self.t_div = t
	
	def set_final_expression_levels(self, x):
		self.x_final = x
	
	def set_final_methylation_levels(self, y):
		self.y_final = y
	
	def set_children(self, child_1, child_2):
		self.children = [child_1, child_2]	# List of the cell's children's ID number.

	def set_generation_index(self, gen_i):
		self.gen_i = gen_i
	
	def get_life_span(self):
		return [self.t_0, self.t_div]
	

	

