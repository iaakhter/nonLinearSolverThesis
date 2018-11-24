# @author Mark Greenstreet
# Basic interval operations

import numpy as np
import functools

def my_reduce_last_dim_help(op, src, dst):
	if(src.ndim == 2):
		for i in range(len(src)):
			dst[i] = reduce(op, src[i])
	else:
		for i in range(len(src)):
			my_reduce_last_dim_help(op, src[i], dst[i])

def my_reduce_last_dim(op, x):
	if(not hasattr(x, 'ndim')):
		return x
	if(x.ndim == 1):
		return(np.array(functools.reduce(op, x)))
	dims = [];
	xx = x;
	for i in range(x.ndim - 1):
		dims.append(len(xx))
		xx = xx[0]
	result = np.zeros(dims)
	my_reduce_last_dim_help(op, x, result)
	return result

def my_min(x):
	return my_reduce_last_dim(lambda x, y: min(x,y), x)

def my_max(x):
	return my_reduce_last_dim(lambda x, y: max(x,y), x)

def interval_p(x):
	if(x is None): return(True)
	else: return hasattr(x, 'ndim') and (x.ndim == 1) and (len(x) == 2)

def tiny_p(x):
	if(interval_p(x)):
	  return(tiny_p(x[0]) and tiny_p(x[1]))
	else:
	  return(abs(x) < 1.0e-14)

def interval_fix(x):
	return (x if interval_p(x) else np.array([x, x]))

def interval_lo(x):
	if(x is None): return None
	if(interval_p(x)): return x[0]
	else: return x

def interval_mid(x):
	if(x is None): return None
	if(interval_p(x)): return (x[0] + x[1])/2.0
	else: return x

def interval_hi(x):
	if(x is None): return None
	if(interval_p(x)): return x[1]
	else: return x

def interval_r(x):
	if(x is None): return 0
	if(interval_p(x)): return max((x[1] - x[0])/2.0, 0)
	else: return 0

# Round interval x using rounded interval arithmetic
def interval_round(x):
	if interval_p(x):
		x[0] = np.nextafter(x[0], np.float("-inf"))
		x[1] = np.nextafter(x[1], np.float("inf"))
	return x

# Turn a floating point into an interval
def interval_applyRia(x):
	if interval_p(x):
		return x
	return np.array([np.nextafter(x, np.float("-inf")), np.nextafter(x, np.float("inf"))])

def interval_add(x, y):
	if((x is None) or (y is None)): return None
	if(interval_p(x) and interval_p(y)):
		return interval_round(np.array([x[0]+y[0], x[1]+y[1]]))
	elif(interval_p(x)):
		return interval_round(np.array([x[0]+y, x[1]+y]))
	elif(interval_p(y)):
		return interval_round(np.array([x+y[0], x+y[1]]))
	else: return(x+y)

def interval_neg(x):
	if(x is None): return None
	if(interval_p(x)):
		return(np.array([-x[1], -x[0]]))
	else: return(-x)

def interval_max(x,y):
	if not(interval_p(x) or interval_p(y)): return max(x,y)
	if interval_p(x) and interval_p(y): return np.array([max(x[0], y[0]), max(x[1], y[1])])
	if interval_p(x): return np.array([max(x[0], y), max(x[1], y)])
	if interval_p(y): return np.array([max(x, y[0]), max(x, y[1])])


def interval_sub(x, y):
	if((x is None) or (y is None)): return None
	return(interval_add(x, interval_neg(y)))

def interval_mult(x, y):
	if((x is None) or (y is None)): return None
	if(interval_p(x) and interval_p(y)):
		p = [xx*yy for xx in x for yy in y]
		return interval_round(np.array([min(p), max(p)]))
	elif(interval_p(x)):
		if(y >= 0):
			return interval_round(np.array([y*x[0], y*x[1]]))
		else:
			return interval_round(np.array([y*x[1], y*x[0]]))
	elif(interval_p(y)):
		return interval_mult(y,x)
	else: return(x*y)

def interval_div(x, y):
	if((x is None) or (y is None)): return None
	if((interval_p(y) and y[0]*y[1] <= 0) or tiny_p(y)):
		return interval_round(np.array([float('-inf'), float('+inf')]))
	elif(interval_p(y)):
		q = [xx/yy for xx in interval_fix(x) for yy in y]
		return interval_round(np.array([min(q), max(q)]))
	elif(interval_p(x)):
		if(y >= 0): 
			return interval_round(np.array([x[0]/y, x[1]/y]))
		else: 
			return interval_round(np.array[x[1]/y, x[0]/y])
	else: return((x+0.0)/y)

def interval_dotprod(x, y):
	try:
		if(len(x) != len(y)):
			raise ValueError('mismatched lengths for x and y: len(x) = ' + str(len(x)) + ', len(y) = ' + str(len(y)))
	except:
		raise ValueError('x and y must be vectors')
	return reduce(
		lambda p, q: interval_add(p, q),
		[ interval_mult(x[i], y[i]) for i in range(len(x))])

def interval_intersect(x, y):
	if((x is None) or (y is None)):
		r = None
	else:
		#print ("x", x)
		#print ("y", y)
		lo = max(x[0], y[0])
		hi = min(x[1], y[1])
		if(lo <= hi):
			r = np.array([lo, hi])
		else:
			r = None
	return r

def interval_union(x, y):
	if(x is None): return y
	elif(y is None): return x
	else: return np.array([min(x[0], y[0]), max(x[1], y[1])])
