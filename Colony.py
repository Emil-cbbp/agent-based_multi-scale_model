import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
import seaborn as sns
import Cell
import OdeSolver as ODE
import Gillespie_methyl_collab as GILcoll
import copy
from itertools import cycle, islice
import timer
import types
from ete3 import Tree, TreeStyle, NodeStyle, TextFace, CircleFace
import Colourmap as CM



class Colony:
	"""A general class which keeps track on colonies of Cell-objects."""

	def __init__(self, x_0, division_scheme, rate_equations, params_x, gene_names, PATH='', day_length=75, methyl='', **kwargs):
		"""Initialises the Colony object. 
		Positional arguments:
			x_0 -- list of initial expression level of transcription factors
			division_scheme -- method which specifies when cells should divide
			rate_equations -- method which contains the rate equations for x on gillespie-format
			params_x -- list of parameters required for the rate equations
			gene_names -- dictionare of names of the genes with corresponding indices
		Keyword arguments:
			PATH -- string path where to store results
			day_length -- number to divide a day into
			methyl='' -- string which specifies what type of methylation should be used. If empty, no methylation is used.
		
		Possible **kwargs :
		related to methylation
			y_0 -- list of initial methylation values.
			rates_epi -- method which contains the rate equations for y on gillespie-format.
			params_epi -- list of parameters required for rates_epi
			N_CpG -- int with the number of of CpG sites used
			mediator_range -- int with the range of which mediators reach
			gillespie -- class in which collaborative gillespie solving methods is
			knockdown -- dictionary containing the names of all genes and how the desired knockdown-level (0,1) and an argument 'time_on' which is the time when the knockdown is initiated.
		"""
		self.methyl = methyl	# String which decides which methylation model to use
		if(methyl):
			if('rates_epi' in kwargs):
				if(isinstance(kwargs['rates_epi'], (types.FunctionType, types.MethodType))):
					self.rates_epi = kwargs['rates_epi']
				else:
					print('Error: rates_epi must be a function.')
			else:
				print('Error: rate equations for epigenetics not provided.')
			
			if('params_epi' in kwargs):
				self.params_epi = kwargs['params_epi']
			else:
				print('Error: parameters for epigenetic rate equations not provided.')
			if('y_0' in kwargs):
				y_0 = kwargs['y_0']
			else:
				y_0 = []
			if('N_CpG' in kwargs):
				self.N_CpG = kwargs['N_CpG']
			if('mediator_range' in kwargs):
				self.mediator_range = kwargs['mediator_range']
			else:
				self.mediator_range = False
			if('gillespie' in kwargs):
				self.gillespie = kwargs['gillespie']
			else:
				print('Error: Which Gillespie solver class to use not provided.')
		else:
			y_0 = []
		if('knockdown' in kwargs):
			self.knockdown = kwargs['knockdown']
			#example: knockdown = {'Runx1': 1.0, 'Tcf7': 0.2, 'Gata3': 1.0, 'Pu1': 1.0, 'time_on':200}
		else:
			print('Error: Knockdown dictionary not provided.')
		
		self.division_scheme = division_scheme
		self.rates = rate_equations
		self.params_x = params_x
		self.gene_names = list(gene_names.keys())
		self.gene_dic = gene_names
		self.day_length = day_length
		self.h_to_au = self.day_length/24
		self.PATH = PATH
		
		# If knockdown from time 0 adjust initial values
		x_0_kd = copy.deepcopy(x_0)
		if(self.knockdown['time_on'] == 0):
			tmp = [k for k in self.knockdown if k != 'time_on'] 
			for k in tmp:
				x_0_kd[self.gene_dic[k]] = round(self.knockdown[k] * x_0_kd[self.gene_dic[k]])
			

		"""Dictionary containing all Cell objects. {cell_ID:Cell_object}."""
		self.colony = {0: Cell.Cell(0, -1, self.division_scheme(0) * self.h_to_au, x_0=x_0_kd, y_0=y_0)}

		"""Dictionary containing the lineage, i.e. all the predecessors, of 
		each cell. {cell_ID:[ID of predecessors]}
		"""
		self.lineage = {0: []}
		
		"""Dictionary containing the total offsprings of 
		each cell. {cell_ID:[ID of offsprings]}
		"""
		self.offspring = {0: []}

		"""2-dimensional list with all the cell_IDs structured by generation [g][c]"""
		self.generations = [[0]]
		
		"""2-dimensional list with all the division_times structured by generation [g][c]"""
		self.division_times = [[self.colony[0].t_div]]
		
		
		"""Storage lists/dictionary for measurments"""
		self.measure_gens = []
		self.measure_x = {}
		self.measure_meth = []
		
			
		
	
	
	def get_division_times(self):
		"""Creates from scratch and returns a list with the division times of all cells."""
		times = copy.deepcopy(self.generations)
		for g in range(len(times)):
			for c in range(len(times[g])):
				times[g][c] = self.colony[self.generations[g][c]].t_div
		return times
	

	def update_division_times(self):
		"""Updates the list with division times."""
		if(len(self.division_times) == len(self.generations)):
			trueth = np.zeros(len(self.generations), dtype=bool)
			for i in range(len(self.division_times)):
				if(len(self.division_times[i]) == len(self.generations[i])):
					trueth[i] = True
			if(np.all(trueth)):
				return
			else:
				for g, b in enumerate(trueth):
					if(not b):
						tmp = [0 for x in range(len(self.generations[g]))]
						for c in range(len(self.generations[g])):
							tmp[c] = self.colony[self.generations[g][c]].t_div
						self.division_times[g] = tmp
		else:
			for g in range(len(self.division_times)):
				if(len(self.division_times[g]) == len(self.generations[g])):
					continue
				else:
					tmp = [0 for x in range(len(self.generations[g]))]
					for c in range(len(self.generations[g])):
						tmp[c] = self.colony[self.generations[g][c]].t_div
					self.division_times[g] = tmp
			for g in range(len(self.division_times), len(self.generations)):
				tmp = [0 for x in range(len(self.generations[g]))]
				for c in range(len(self.generations[g])):
					tmp[c] = self.colony[self.generations[g][c]].t_div
				self.division_times.append(tmp)
			

	def set_offspring(self, ID):
		"""For the given cell_ID, the cell is added as offspring to its lineage."""
		added_offspring = self.colony[ID].children
		self.offspring[ID] = added_offspring
		while(True):
			m_ID = self.colony[ID].mother
			if(m_ID < 0):
				break
			self.offspring[m_ID] += added_offspring
			ID = m_ID
	
	
	def set_lineage(self, ID):
		"""For the given cell_ID, the cell's lineage is created"""
		m_ID = self.colony[ID].mother
		self.lineage[ID] = [m_ID] + self.lineage[m_ID]
			
	
	def divide(self, ID):
		"""Carries out the cell division process for the given cell_ID.
		Creates two children and add them to the colony. Updates linage, 
		offspring, generations and division_times.
		"""

		cell = self.colony[ID]

		# Assign child-ID and add children to the colony.
		tmp = sum([2**g for g in range(cell.gen)])
		child_1 = tmp + 2**cell.gen + 2 * (cell.ID - tmp)
		child_2 = child_1 + 1
		t_1 = self.division_scheme(cell.gen + 1) * self.h_to_au
		t_2 = self.division_scheme(cell.gen + 1) * self.h_to_au
		self.colony[child_1] = Cell.Cell(child_1, cell, t_1)
		self.colony[child_2] = Cell.Cell(child_2, cell, t_2)
		cell.set_children(child_1, child_2)
		
		# Add to generations list.
		if(cell.gen + 1 > len(self.generations) - 1):
			self.generations.append([child_1])
		else:
			self.generations[cell.gen+1].append(child_1)
		self.colony[child_1].set_generation_index(len(self.generations[cell.gen+1]) - 1)
		self.generations[cell.gen+1].append(child_2)
		self.colony[child_2].set_generation_index(len(self.generations[cell.gen+1]) - 1)

		# Add to division_times list.
		if(cell.gen + 1 > len(self.division_times) - 1):
			self.division_times.append([self.colony[child_1].t_div])
		else:
			self.division_times[cell.gen+1].append(self.colony[child_1].t_div)
		self.division_times[cell.gen+1].append(self.colony[child_2].t_div)
		
		# Update offspring and lineage.
		self.set_offspring(cell.ID)
		self.set_lineage(child_1)
		self.set_lineage(child_2)
	
	
	
	def evolve_cell(self, ID):
		"""Evolces the gene expressions x for the given cell_ID by the rates 
		equations. If methylation is used in the colony, that is evolved as well.
		"""
		
		cell = self.colony[ID]
		
		# Evolve the expression levels.
		stoch_sol = ODE.StochOdeSolver(self.rates, self.params_x, cell.x_0, self.knockdown, 10000, 1)
		stoch_sol.gillespie_time(cell.t_div, cell.t_0)
		
		if(len(stoch_sol.t) > 0):
			cell.x = np.array(stoch_sol.x)
			cell.t = np.array(stoch_sol.t)
		cell.set_final_expression_levels(cell.x[-1])
		
		
		# Evolve the methylation.
		if(self.methyl):
			stoch_epi = self.gillespie(	self.rates_epi, cell.t_0, cell.y_0, self.params_epi,
										10000, cell.t, cell.x, self.mediator_range, self.knockdown)
			stoch_epi.gillespie_timeseries()
			
			if(len(stoch_epi.t) > 0):
				cell.y = np.array(stoch_epi.y)
				cell.t_meth = np.array(stoch_epi.t)
				cell.set_final_methylation_levels(cell.y[-1])
		
		
		
	
	def run_N_generations(self, N):
		""" Creates and evolve a colony for N generations."""
		g = 0	# Generation counter
		while(g < N):
			print('Generation', g)
			for i in self.generations[g]:
				cell = self.colony[i]
				self.evolve_cell(i)
				# Do not do cell division for the last generation.
				if(g < N - 1):
					self.divide(i)
			g += 1
				

	def run_max_time(self, T):
		"""Creates and evolve a colony until all the cells in a generation have 
		passed the given maximum time T.
		"""
		g = 0	# Generation counter
		t_min = min(self.division_times[g])
		while(True):
			print('Generation', g)
			for i in self.generations[g]:
				cell = self.colony[i]
				self.evolve_cell(i)
			self.update_division_times()
			t_min = min(self.division_times[g])
			# If at least 1 cell has not reached T, do one more division.
			if(t_min < T):
				for i in self.generations[g]:
					self.divide(i)
			else:
				break
			g += 1


	def measure(self, measure_points, measure_levels, measure_methylation=False, meth_frac=0.75):
		"""Peforms measurements at given timepoints. Measures generation distribution,
		expression levels of chosen genes, and optionally methylation.
		Arguments:
			measure_points -- list of time points at which measurments should be conducted
			measure_levels -- dictionary {gene_name: index_in_x} of genes to be measured
			measure_methylation -- bool if methylation should be measured of not
			meth_frac -- float with what ratio unmethylated/N_CpG should be considered as open
		"""
		self.measure_gens = [np.zeros(len(self.generations)) for i in range(len(measure_points))]
		self.measure_x = {g: [[] for i in range(len(measure_points))] for g in measure_levels}
		self.measure_meth = np.zeros(len(measure_points))
		
		self.meth_frac = meth_frac
		
		count = 0
		for i in self.colony:
			cell = self.colony[i]
			for k, m_p in enumerate(measure_points):
				T = cell.get_life_span()
				if(m_p > T[0] and m_p <= T[1]):
					count += 1
					self.measure_gens[k][cell.gen] += 1
					for g in measure_levels:
						self.measure_x[g][k].append(cell.x[:,measure_levels[g]][np.max(np.where(cell.t < m_p))])
						
					if(measure_methylation):
						if(cell.y[:,-1][np.max(np.where(cell.t_meth < m_p))] >= meth_frac * self.N_CpG):
							self.measure_meth[k] += 1
			
	
				
	
	def save_colony(self, sparse=False):
		""" Saves colony information in textfile colony.tsv
		"""
		print('Save')
		list_subfolders = [f.name for f in os.scandir(self.PATH) if f.is_dir()]
		name = 'gen_{:d}'.format(
					len(self.generations)-1)
		
		existing = [d for d in list_subfolders if name in d]
		name += '_run_{:d}'.format(len(existing) + 1)
		self.name = name
		os.mkdir(self.PATH + name)
		
		os.mkdir(self.PATH + name + '/Expression/')
		os.mkdir(self.PATH + name + '/Methylation/')
		
		info = open(self.PATH + name + '/info', 'w')
		info.write('x_0 = {:s}\n'.format(str(self.colony[0].x_0)))
		info.write('y_0 = {:s}\n'.format(str(self.colony[0].y_0)))
		var = self.__dict__
		for v in var:
			info.write('{:s} = {:s}\n'.format(v, str(var[v]).replace('\n', '')))
		info.close()
		
		col = open(self.PATH + name + '/colony.tsv', 'w')
		
		initial_x = '\t'.join([x + '_init' for x in self.gene_names])
		final_x = '\t'.join([x + '_final' for x in self.gene_names])
		
		if(len(self.colony[0].y_0) == 2):
			meth = 'C\tO'
			initial_meth = 'C_init\tO_init'
			final_meth = 'C_final\tO_final'
		elif(len(self.colony[0].y_0) == 3):
			meth = 'C\tI\tO'
			initial_meth = 'C_init\tI_init\tO_init'
			final_meth = 'C_final\tI_final\tO_final'
		col.write((	'cell_ID\tmother\tchild_1\tchild_2\tgen\tgen_ind\tt_birth\tt_div\t{:s}'
					+'\t{:s}\t{:s}\t{:s}\tBcl11b_frac\texpression_file\tmethylation_file\n'
					).format(initial_x, final_x, initial_meth, final_meth))
		
		for c in self.colony:
			cell = self.colony[c]
			out = '{:d}\t{:d}\t{:d}\t{:d}\t{:d}\t{:d}\t{:f}\t{:f}\t'.format(cell.ID, 
					cell.mother, cell.children[0], cell.children[1], cell.gen,
					cell.gen_i, cell.t_0, cell.t_div)
			for x in cell.x_0:
				out += str(x) + '\t'
			for x in cell.x_final:
				out += str(x) + '\t'
			for y in cell.y_0:
				out += str(y) + '\t'
			for y in cell.y_final:
				out += str(y) + '\t'
			out += str(cell.y_final[-1] / self.N_CpG) + '\t'
			out += 'Expression/cell_{:d}\t'.format(cell.ID)
			out += 'Methylation/cell_{:d}'.format(cell.ID)
			out += '\n'
			col.write(out)
			
			
			# Save every time step of the stochastic simulations for both 
			# transcriptional level and epigenetic level.
			if(not sparse):
				cell.t = cell.t.reshape((cell.t.shape[0],1))
				conc = np.concatenate((cell.t, cell.x), axis=1)
				np.savetxt(self.PATH + name + '/Expression/cell_{:d}.tsv'.format(cell.ID), conc, delimiter='\t', header='time\t' + '\t'.join(self.gene_names), comments='')
				
				cell.t_meth = cell.t_meth.reshape((cell.t_meth.shape[0],1))
				conc = np.concatenate((cell.t_meth, cell.y), axis=1)
				np.savetxt(self.PATH + name + '/Methylation/cell_{:d}.tsv'.format(cell.ID), conc, delimiter='\t', header='time\t' + meth, comments='')
		
		col.close()
		
		
	

	def create_tree(self, binary=False, save=False, node_labels=[], mark_X=False):
		"""Convert the colony to a ete3 tree-structure and plot it."""
		print('Create tree')
		max_val = np.max([cell.y_final[-1] / self.N_CpG for cell in list(self.colony.values())])
		cm = CM.ColourMap()
		
		gen_tree = copy.deepcopy(self.generations)

		gen_tree[0][0] = Tree(name='0', dist=self.colony[0].t_div)
		gen_tree[0][0].add_feature('cell', self.colony[0])
		if(len(node_labels)):
			gen_tree[0][0].add_face(TextFace('{:s}'.format(node_labels[self.generations[0][0]])), 0, position='branch-top')


		nstyle = NodeStyle()
		if(binary):
			if(mark_X):
				if(self.colony[0].y_final[-1]/self.N_CpG >= self.meth_frac):
					nstyle["fgcolor"] = 'white'
				else:
					if(self.colony[0].x_final[self.gene_dic['X']] <= 0):
						nstyle["fgcolor"] = 'red'
					else:
						nstyle["fgcolor"] = 'black'
			else:
				nstyle["fgcolor"] = 'white' if self.colony[0].y_final[-1]/self.N_CpG >= self.meth_frac else 'black'
		else:
			nstyle["fgcolor"] = cm.colour(self.colony[0].y_final[-1]/self.N_CpG)
		nstyle["bgcolor"] = "PaleGreen"
		nstyle["size"] = 15
		gen_tree[0][0].set_style(nstyle)
		
		
		ts = TreeStyle()
		ts.show_leaf_name = False
		ts.show_branch_length = False
		ts.mode = "c"
		ts.arc_start = 0
		ts.arc_span = 360
		
		
		for g in range(1, len(gen_tree)):
			for c in range(len(gen_tree[g])):
				cell = self.colony[self.generations[g][c]]
				
				if(g == 0):
					gen_tree[g][c] = t.add_child(name=str(cell.ID), dist=(cell.t_div-0))
				else:
					mother = self.colony[cell.mother]
					gen_tree[g][c] = gen_tree[g-1][mother.gen_i].add_child(name=str(cell.ID), dist=(cell.t_div - mother.t_div))
				gen_tree[g][c].add_feature('cell', cell)
				gen_tree[g][c].add_feature('age', cell.t_div)
				gen_tree[g][c].add_feature('Bcl11b', cell.y_final[-1]/self.N_CpG)
				
				if(len(node_labels) > 0):
					gen_tree[g][c].add_face(TextFace('{:s}'.format(node_labels[self.generations[g][c]])), 0, position='branch-top')
				
				nstyle = NodeStyle()
				if(binary):
					if(mark_X):
						if(cell.y_final[-1]/self.N_CpG >= self.meth_frac):
							nstyle["fgcolor"] = 'white'
						else:
							X_off = []
							for m in self.lineage[cell.ID]:
								X_off.append(True if 0 in self.colony[m].x[:,self.gene_dic['X']] else False)
							X_off.append(True if 0 in cell.x[:,self.gene_dic['X']] else False)
							if(cell.x_final[self.gene_dic['X']] <= 0):
								nstyle["fgcolor"] = 'red'
							else:
								nstyle["fgcolor"] = 'black'
					else:
						nstyle["fgcolor"] = 'white' if cell.y_final[-1]/self.N_CpG >= self.meth_frac else 'black'
				else:
					nstyle["fgcolor"] = cm.colour(cell.y_final[-1]/self.N_CpG)
				nstyle["size"] = 15
				gen_tree[g][c].set_style(nstyle)
		
		if(save):
			if(binary):
				if(mark_X):
					gen_tree[0][0].render(self.PATH + self.name + '/tree_binary_X.svg', tree_style=ts)
					gen_tree[0][0].render(self.PATH + self.name + '/tree_binary_X.png', tree_style=ts)
				else:
					gen_tree[0][0].render(self.PATH + self.name + '/tree_binary.svg', tree_style=ts)
					gen_tree[0][0].render(self.PATH + self.name + '/tree_binary.png', tree_style=ts)
			else:
				gen_tree[0][0].render(self.PATH + self.name + '/tree_colour.svg', tree_style=ts)
				gen_tree[0][0].render(self.PATH + self.name + '/tree_colour.png', tree_style=ts)
		return gen_tree[0][0], ts
		
	
	

	
	def plot_gen_dist(self, plot=True):
		N = len(self.colony)
		count = [c for day in self.measure_gens for c in day]
		generations = [x for x in range(len(self.generations))]*len(self.measure_gens)
		days = []
		for d in range(len(self.measure_gens)):
			days += ['Day {:d}'.format(d+2)]*len(self.generations)
		gens_dic = {'Generations': generations, 'Cells [%]': count, 'Days': days}
		
		self.gens_df = pd.DataFrame(gens_dic)
		if(plot):
			plot = sns.barplot(x='Generations', y='Cells [%]', hue='Days', data=self.gens_df)
		
		return self.gens_df
			
	
	
	def string_to_list(self, string, conv, delimiter=''):
		"""Possible conv args: 
			* 'str' -- convert/keep as string
			* 'np' -- convert to numpy array
			* 'list' -- convert to nested list
			* int -- convert to int
			* float -- convert to floar
		"""
		if(delimiter):
			if(delimiter == '],'):
				string = string.replace('],', ']],')
			elif(delimiter == '),'):
				string = string.replace('),', ')),')	
			tmp = string[1:-1].split(delimiter)
			tmp = [x.strip() for x in tmp]
		else:
			tmp = string[1:-1].split()
		if(conv == 'str'):
			tmp = [s.strip()[1:-1] for s in tmp]
		elif(conv == 'np'):
			int_or_float = tmp[0][7:-2].split()[0]
			if('.' in int_or_float):
				int_or_float = float
			else:
				int_or_float = int
			tmp = [np.array(self.string_to_list(x[6:-1], delimiter=',', conv=int_or_float)) for x in tmp]
		elif(conv == 'list'):
			int_or_float = tmp[0][1:-1].split()[0]
			if('.' in int_or_float):
				int_or_float = float
			else:
				int_or_float = int
			tmp = [self.string_to_list(x, delimiter=',', conv=int_or_float) for x in tmp]
		else:
			tmp = [s if not s.replace('.','',1).isdigit() else conv(s) for s in tmp]
			
		if(len(tmp) == 1):
			if(tmp[0] == ''):
				return []
		return tmp
	
	def string_to_dic(self, string, conv, delimiter=''):
		"""Possible conv args: 
			* 'str' -- convert/keep as string
#			* 'np' -- convert to numpy array
			* 'list' -- convert to nested list
			* int -- convert to int
			* float -- convert to floar
		"""
		if(delimiter):
			if(delimiter == '],'):
				string = string.replace('],', ']],')
			elif(delimiter == ']],'):
				string = string.replace(']],', ']]],')
			elif(delimiter == '),'):
				string = string.replace('),', ')),')	
			tmp = string[1:-1].split(delimiter)
			tmp = [x.strip() for x in tmp]
		else:
			tmp = string[1:-1].split()
		
		if(conv == 'str'):
			dic = {s if not s.split(':')[0].strip().isnumeric else int(s.split(':')[0].strip()): s.split(':')[1].strip() for s in tmp}
		elif(conv == 'list'):
			int_or_float = tmp[-1].split(':')[1]
			if('.' in int_or_float):
				int_or_float = float
			else:
				int_or_float = int
			if(delimiter == '],'):
				dic = {	s.split(':')[0].strip()[1:-1] if not s.split(':')[0].strip().isnumeric() else int(s.split(':')[0].strip()): 
						self.string_to_list(s.split(':')[1].strip(), conv=int_or_float, delimiter=',') for s in tmp}
			elif(delimiter == ']],'):
				dic = {	s.split(':')[0].strip()[1:-1] if not s.split(':')[0].strip().isnumeric() else int(s.split(':')[0].strip()): 
						self.string_to_list(s.split(':')[1].strip(), conv='list', delimiter='],') for s in tmp}
		else:
			dic = [s if not s.replace('.','',1).isdigit() else conv(s) for s in tmp]
	
		return dic

