import time

def start_time():
	return time.time()

# By default return in H M S. Can by option return only seconds
def time_elapsed(start_time, p_out=True, seconds=False, prec='', return_time=False, return_time_string=False, current_time=False):
		run_time = time.time() - start_time
		if(current_time):
			c_t = time.strftime("%H:%M:%S", time.localtime())
		if(seconds):
			t = ("--- Time elapsed: {:" + prec + "f} s ---").format(run_time)
			if(current_time):
				t +=  "Current time: {:s} ---".format(c_t)
			if(p_out):
				print(t)
			if(return_time and not return_time_string):
				return run_time
			elif(not return_time and return_time_string):
				return t
			elif(return_time and return_time_string):
				return run_time, t
		else:
			hours = run_time // 3600
			minutes = (run_time % 3600) // 60
			seconds = run_time - hours * 3600 - minutes * 60
			t = ("--- Time elapsed: {:d} h {:d} min and {:" + prec + "f} s ---").format(int(hours), int(minutes), seconds)
			if(current_time):
				t +=  " Current time: {:s} ---".format(c_t)
			if(p_out):
				print(t)
			if(return_time and not return_time_string):
				return run_time
			elif(not return_time and return_time_string):
				return t
			elif(return_time and return_time_string):
				return run_time, t
