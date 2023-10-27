import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import matplotlib
from matplotlib.collections import PatchCollection

class ColourMap:
	
	def __init__(self, max_val=1, min_val=0):
		self.max_val = max_val
		self.min_val = min_val
		self.norm = 0xff /(self.max_val - self.min_val)

	def colour(self, val, tup=False):
		""" Normalises the colour scale and yields a colour for a given values. """	
		c = int((val - self.min_val) * self.norm)	# Normalise value
		c = max(0, min(0xff,c))	# Limits value to 0 - 255
		c = self.mapColour(c, tup)	# applies mapping
		return c

	def mapColour(self, val, tup=False):
		"""	Imitates the default scale of pm3d in gnuplot. Similar to
		RGBColourScale, but goes through violet rather than blue.
		"""
		rgb = 0
		# RED
		if(val < 0x80):
			rgb |= (val * 2) << 16
		else:
			rgb |= 0xff << 16

		# GREEN
		if(val >= 0x80):
			rgb |= ((val - 0x80) * 2) << 8

		# BLUE
		if(val < 0x40):
			rgb |= val * 4
		elif(val < 0x80):
			rgb |= 0xff - (val - 0x40) * 4
		elif(val > 0xc0):
			rgb |= (val - 0xc0) * 4

		return rgb
	
	def convert_to_hex_string(self, num):
		h = '#{:0>6}'.format(hex(num)[2:])
		return(h)
		
	
	def colourbar(self, n):
		fig, ax = plt.subplots(1, figsize=(1,5))
		bar = [plt.Rectangle((0, x * 1/n), 0.2, 1/n, fill=True, edgecolor=None, color=self.convert_to_hex_string(self.colour(x * 1/n))) for x in range(n)]
		pc = PatchCollection(bar, match_original=True)
		ax.add_collection(pc)
		
		plt.xlim(0,0.2)
		plt.xticks([])
		plt.tight_layout()
		plt.show()
	


