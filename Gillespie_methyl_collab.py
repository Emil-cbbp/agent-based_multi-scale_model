import numpy as np
import copy


class Gillespie:
	def __init__(self, ode, t0, y0, params, iterations,
				 t_transcriptional, expr_transcriptional, mediator_range, knockdown):
		self.iterate = iterations
		self.t0 = t0
		self.y0 = y0
		self.ode = ode
		self.params = params
		self.t_trans = t_transcriptional
		self.x_trans = expr_transcriptional
		self.t = [self.t0]
		self.y = [np.array(copy.deepcopy(self.y0))]
		self.N_CpG = int(sum(self.y0))
		self.CpG = np.zeros((self.N_CpG))
		self.CpG[0:int(self.y0[1])] = 1		# Set Intermediate
		self.CpG[int(self.y0[1]):int(self.y0[1]) + int(self.y0[0])] = 2		# Set Closed
		self.mediator = mediator_range

	def gillespie_timeseries(self):
		standard_rates = True
		while(True):
			ind_trans = np.max(np.where(self.t_trans <= self.t[-1]))
			R = self.x_trans[ind_trans][0]
			X = self.x_trans[ind_trans][4]
			N = self.x_trans[ind_trans][5]
			x = [R, X, N]

			# compute reaction rates a_mu and sum a_0
			a_mu = self.ode(self.y[-1], x, self.t[-1], self.params)
			a_0 = np.sum(a_mu)

			# finish when a0 = 0
			if a_0 == 0:
				self.t = np.array(self.t)
				self.y = np.array(self.y)
				break

			# choose reaction
			a = 0
			idx = -1
			rand = np.random.rand(2)
			while(a < rand[1] * a_0):
				idx += 1
				a += a_mu[idx]

			# update time and y values
			self.t.append(self.t[-1] + 1./a_0 * np.log(1/rand[0]))
			self.y.append(np.array(self.y[-1]))
			rd_site = np.random.randint(0, len(self.CpG))
			if( idx == 0 and self.y[-1][2] != 0 and self.CpG[rd_site] == 0):
				self.y[-1][1] += 1
				self.y[-1][2] -= 1
				self.CpG[rd_site] = 1
			elif(idx == 1 and self.y[-1][1] != 0 and self.CpG[rd_site] == 1):
				self.y[-1][0] += 1
				self.y[-1][1] -= 1
				self.CpG[rd_site] = 2
			elif(idx == 2 and self.y[-1][0] != 0 and self.CpG[rd_site] == 2):
				self.y[-1][1] += 1
				self.y[-1][0] -= 1
				self.CpG[rd_site] = 1
			elif(idx == 3 and self.y[-1][1] != 0 and self.CpG[rd_site] == 1):
				self.y[-1][2] += 1
				self.y[-1][1] -= 1
				self.CpG[rd_site] = 0
				

			elif(idx == 4 and self.CpG[rd_site] == 0 and self.y[-1][2] != 0):
				rd_mediator = np.random.randint(max(0, rd_site - self.mediator),
												min(self.N_CpG, rd_site + self.mediator))
				if(self.CpG[rd_mediator] == 2):
					self.y[-1][1] += 1
					self.y[-1][2] -= 1
					self.CpG[rd_site] = 1
			elif(idx == 5 and self.CpG[rd_site] == 1 and self.y[-1][1] != 0):
				rd_mediator = np.random.randint(max(0, rd_site - self.mediator),
												min(self.N_CpG, rd_site + self.mediator))
				if(self.CpG[rd_mediator] == 2):
					self.y[-1][0] += 1
					self.y[-1][1] -= 1
					self.CpG[rd_site] = 2
			elif(idx == 6 and self.CpG[rd_site] == 1 and self.y[-1][1] != 0):
				rd_mediator = np.random.randint(max(0, rd_site - self.mediator),
												min(self.N_CpG, rd_site + self.mediator))
				if(self.CpG[rd_mediator] == 0):
					self.y[-1][2] += 1
					self.y[-1][1] -= 1
					self.CpG[rd_site] = 0
			elif(idx == 7 and self.CpG[rd_site] == 2 and self.y[-1][0] != 0):
				rd_mediator = np.random.randint(max(0, rd_site - self.mediator),
												min(self.N_CpG, rd_site + self.mediator))
				if(self.CpG[rd_mediator] == 0):
					self.y[-1][1] += 1
					self.y[-1][0] -= 1
					self.CpG[rd_site] = 1
			elif(idx == 8 and self.CpG[rd_site] == 1 and self.y[-1][1] != 0):
				rd_mediator = np.random.randint(max(0, rd_site - self.mediator),
												min(self.N_CpG, rd_site + self.mediator))
				if(self.CpG[rd_mediator] == 1):
					self.y[-1][0] += 1
					self.y[-1][1] -= 1
					self.CpG[rd_site] = 2

			# break when max time on transcriptional level is reached
			if self.t[-1] > self.t_trans[-1]:
				self.t = np.array(self.t)
				self.y = np.array(self.y)
				break




