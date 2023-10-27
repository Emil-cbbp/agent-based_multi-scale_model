import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import os
import timer
import Colony_plot_LCA_labels as col


def string_to_list(string, conv, delimiter=''):
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
		tmp = [np.array(string_to_list(x[6:-1], delimiter=',', conv=int_or_float)) for x in tmp]
	elif(conv == 'list'):
		int_or_float = tmp[0][1:-1].split()[0]
		if('.' in int_or_float):
			int_or_float = float
		else:
			int_or_float = int
		tmp = [string_to_list(x, delimiter=',', conv=int_or_float) for x in tmp]
	else:
		tmp = np.array([s if not s.replace('.','',1).isdigit() else conv(s) for s in tmp])
		
	if(len(tmp) == 1):
		if(tmp[0] == ''):
			return []
		
	
	return tmp

def string_to_dic(string, conv, delimiter=''):
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
					string_to_list(s.split(':')[1].strip(), conv=int_or_float, delimiter=',') for s in tmp}
		elif(delimiter == ']],'):
			dic = {	s.split(':')[0].strip()[1:-1] if not s.split(':')[0].strip().isnumeric() else int(s.split(':')[0].strip()): 
					string_to_list(s.split(':')[1].strip(), conv='list', delimiter='],') for s in tmp}
	else:
		dic = [s if not s.replace('.','',1).isdigit() else conv(s) for s in tmp]

	return dic


N_CpG = 500
meth_frac = 0.75



PATH_results = 'Results/'

folders = [f.name for f in os.scandir(PATH_results) if f.is_dir()]
PATHS = [PATH_results + f.name + '/' for f in os.scandir(PATH_results) if f.is_dir()]



statistics = pd.DataFrame({'Generation': [], 'State': [], 'Open': [], 'Closed': [], 'Same': [], 'Opposite': []})


tot_time_start = timer.start_time()
for PATH in PATHS:
	folder = [f.name for f in os.scandir(PATH) if f.is_dir() and f.name.startswith('gen')]
	for j, fol in enumerate(folder):
		print('Colony {:d}'.format(j))
		print(fol)
		colony = pd.read_csv(PATH + fol + '/colony.tsv', delimiter='\t')
		
		
		with open(PATH + fol + '/info', 'r') as info:
			info_list = info.readlines()
		info_dic = {x.split('=')[0].strip(): x.split('=')[1].strip() for x in info_list}
		
		lineage = string_to_dic(info_dic['lineage'], delimiter='],', conv='list')
		offspring = string_to_dic(info_dic['offspring'], delimiter='],', conv='list')
		generations = string_to_list(info_dic['generations'], delimiter='],', conv='list')
		
		opened = []
		closed = []
		gens = []
		for c in range(len(colony)):
			gens.append(colony.loc[c]['gen'])
			if(colony.loc[c]['Bcl11b_frac'] >= meth_frac):
				opened.append(c)
			else:
				closed.append(c)
			


		print('Start classification')

		labels = ['' for x in range(len(colony))]
		cells_to_label = [x for x in range(len(colony))]

		# If no cell is open
		if(len(closed) == len(colony)):
			labels = ['closed pre LCA' for c in range(len(colony))]
			cells_to_label = []

		# Assign open-label and remove children whose mothers are open.
		opened_copy = opened.copy()
		k = 0
		while(k < len(opened)):
			c = opened[k]
			labels[c] = 'open'
			cells_to_label.remove(c)
			if(colony.loc[c, 'child_1'] != -1):
				tmp = np.in1d(offspring[c], opened, assume_unique=True)
				open_offspring = offspring[c][tmp]
				for i in open_offspring:
					labels[i] = 'open'
					opened.remove(i)
					cells_to_label.remove(i)
				if(not np.all(tmp)):	#Should never happen
					closed_offspring = offspring[c][np.invert(tmp)]
					for i in closed_offspring:	
						labels[i] = 'closed offspring'
						cells_to_label.remove(i)
			k += 1
		opened_sparse = opened
		opened = opened_copy

		LCA = {}
		k = 0
		while(len(cells_to_label) > 0):
			c = cells_to_label[0]
			cell = colony.loc[c]
			if(c in opened):	#safety check, should have been labeled already
				labels[c] = 'open'
				cells_to_label.remove(c)
			if(cell['child_1'] != -1):
				if(cell['child_1'] in list(offspring.keys())):
					if(cell['child_1'] in opened):
						inter_1 = np.intersect1d(np.append(offspring[cell['child_1']], cell['child_1']), opened, assume_unique=True)
					else:
						inter_1 = np.intersect1d(offspring[cell['child_1']], opened, assume_unique=True)
					if(cell['child_2'] in opened):
						inter_2 = np.intersect1d(np.append(offspring[cell['child_2']], cell['child_2']), opened, assume_unique=True)
					else:
						inter_2 = np.intersect1d(offspring[cell['child_2']], opened, assume_unique=True)
				else:
					inter_1 = np.intersect1d([cell['child_1']], opened, assume_unique=True)
					inter_2 = np.intersect1d([cell['child_2']], opened, assume_unique=True)
				open_child_1_offspring = True if len(inter_1) > 0 else False
				open_child_2_offspring = True if len(inter_2) > 0 else False
				if(open_child_1_offspring and open_child_2_offspring):
					order = min(gens[np.min(inter_1)], gens[np.min(inter_2)]) - gens[c]
					LCA[c] = order
					labels[c] = 'LCA {:d}'.format(order)
					cells_to_label.remove(c)
				else:
					child_open = []
					child_closed = []
					if(open_child_1_offspring):
						child_open.append(1)
					else:
						child_closed.append(1)
					if(open_child_2_offspring):
						child_open.append(2)
					else:
						child_closed.append(2)
					LCA_lin = np.intersect1d(lineage[c], list(LCA.keys()), assume_unique=True)
					if(len(LCA_lin) == 0):
						labels[c] = 'closed pre LCA'
						cells_to_label.remove(c)
						for h in child_closed:
							labels[cell['child_{:d}'.format(h)]] = 'closed pre LCA'
							cells_to_label.remove(cell['child_{:d}'.format(h)])
							if(cell['child_1'] in list(offspring.keys())):
								for i in offspring[cell['child_{:d}'.format(h)]]:
									labels[i] = 'closed pre LCA'
									cells_to_label.remove(i)
					else:
						LCA_prev = np.max(LCA_lin)
						labels[c] = 'closed post LCA {:d}'.format(gens[c] - gens[LCA_prev])
						cells_to_label.remove(c)
						for h in child_closed:
							labels[cell['child_{:d}'.format(h)]] = 'closed post LCA {:d}'.format(gens[cell['child_{:d}'.format(h)]] - gens[LCA_prev])
							cells_to_label.remove(cell['child_{:d}'.format(h)])
							if(cell['child_1'] in list(offspring.keys())):
								for i in offspring[cell['child_{:d}'.format(h)]]:
									labels[i] = 'closed post LCA {:d}'.format(gens[i] - gens[LCA_prev])
									cells_to_label.remove(i)
			else:
				LCA_lin = np.intersect1d(lineage[c], list(LCA.keys()), assume_unique=True)
				if(len(LCA_lin) == 0):	#safety check, should have been labeled already
					labels[c] = 'closed pre LCA'
					cells_to_label.remove(c)
				else:
					LCA_prev = np.max(LCA_lin)
					labels[c] = 'closed post LCA {:d}'.format(gens[c] - gens[LCA_prev])
					cells_to_label.remove(c)


		colony['Labels'] = labels
		colony.to_csv(PATH + fol + '/colony_labels.tsv', '\t', index=False)
		timer.time_elapsed(tot_time_start, seconds=False)

























