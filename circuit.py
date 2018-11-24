# @author: Mark Greenstreet and Itrat Akhter
# # Implementation of long channel and short channel mosfet models
# Short channel model calls on a shared library built by C++ boost
# version which actually implements the model. 
# Short channel model is from https://nanohub.org/publications/15/4
# # Implementation of a circuit class that can combine the short
# channel and long channel transistors.
import math
import numpy as np
from cvxopt import matrix,solvers
from scipy.spatial import ConvexHull
from intervalBasics import *
from lpUtilsMark import *
# boost built shared library from C++
import stChannel_py



# @author Mark Greenstreet
class MosfetModel:
	def __init__(self, channelType, Vt = None, k = None, gds=0.0):
		self.channelType = channelType   # 'pfet' or 'nfet'
		self.Vt = Vt                     # threshold voltage
		self.k = k                       # carrier mobility
		if(gds == 'default'):
			self.gds = 1.0e-8  # for "leakage" -- help keep the Jacobians non-singular
		else: self.gds = gds

	def __str__(self):
		return "MosfetModel(" + str(self.channelType) + ", " + str(self.Vt) + ", " + str(self.k) + ", " + str(self.s) + ")"


# @author Itrat Akhter
# Short channel model
class ScMosfet:
	# s, g, d: Indices indicating source, gate and drain
	# model: MosfetModel
	# shape: shape
	def __init__(self, s, g, d, model, shape):
		self.s = s
		self.g = g
		self.d = d
		self.shape = shape
		self.model = model
		self.gradDict = {} # stores computations of gradients for specific Vs, Vg, Vd and channelType
		self.funDict = {} # stores computations of functions for specific Vs, Vg, Vd and channelType
		self.stMosfet = stChannel_py.StMosfet()

	# calculate ids from the given V
	def ids(self, V):
		model = self.model
		return(self.ids_help(V[self.s], V[self.g], V[self.d], model.channelType))

	# ids helper function. 
	# Vs, Vg, Vd can be interval or point values
	def ids_help(self, Vs, Vg, Vd, channelType):
		VbC = stChannel_py.MyList();
		if channelType == "nfet":
			fetFunc = self.stMosfet.mvs_idnMon
			VbC[:] = [0.0, 0.0]
		else:
			fetFunc = self.stMosfet.mvs_idpMon
			VbC[:] = [1.8, 1.8]

		VdC = stChannel_py.MyList()
		VgC = stChannel_py.MyList()
		VsC = stChannel_py.MyList()
		VsDict, VgDict, VdDict = None, None, None
		if interval_p(Vd):
			VdC[:] = [Vd[0], Vd[1]]
			VdDict = tuple(list(Vd))
		else:
			VdC[:] = [Vd, Vd]
			VdDict = Vd

		if interval_p(Vg):
			VgC[:] = [Vg[0], Vg[1]]
			VgDict = tuple(list(Vg))
		else:
			VgC[:] = [Vg, Vg]
			VgDict = Vg

		if interval_p(Vs):
			VsC[:] = [Vs[0], Vs[1]]
			VsDict = tuple(list(Vs))
		else:
			VsC[:] = [Vs, Vs]
			VsDict = Vs

		if (VsDict, VgDict, VdDict, channelType) in self.funDict:
			return self.funDict[(VsDict, VgDict, VdDict, channelType)]

		iVal = fetFunc(VdC, VgC, VsC, VbC)

		if(interval_p(Vs) or interval_p(Vg) or interval_p(Vd)):
			funVal =  self.shape*np.array([np.nextafter(iVal[0], float("-inf")), np.nextafter(iVal[1], float("inf"))])
		else:
			funVal = self.shape*iVal[0]

		self.funDict[(VsDict, VgDict, VdDict, channelType)] = funVal
		
		# Clean up funDict
		if len(self.funDict) >= 10000:
			self.funDict = {}

		return funVal


	# calculate gradient of ids with respect to Vs, Vg and Vd
	# from the given V
	def grad_ids(self, V):
		model = self.model
		return(self.grad_ids_help(V[self.s], V[self.g], V[self.d], model.channelType))

	# gradient helper function. 
	# Vs, Vg, Vd can be interval or point values
	def grad_ids_help(self, Vs, Vg, Vd, channelType):
		VbC = stChannel_py.MyList();
		if channelType == "nfet":
			fetGrad = self.stMosfet.mvs_idnGrad
			VbC[:] = [0.0, 0.0]
		else:
			fetGrad = self.stMosfet.mvs_idpGrad
			VbC[:] = [1.8, 1.8]

		VdC = stChannel_py.MyList()
		VgC = stChannel_py.MyList()
		VsC = stChannel_py.MyList()
		VsDict, VgDict, VdDict = None, None, None
		if interval_p(Vd):
			VdC[:] = [Vd[0], Vd[1]]
			VdDict = tuple(list(Vd))
		else:
			VdC[:] = [Vd, Vd]
			VdDict = Vd

		if interval_p(Vg):
			VgC[:] = [Vg[0], Vg[1]]
			VgDict = tuple(list(Vg))
		else:
			VgC[:] = [Vg, Vg]
			VgDict = Vg

		if interval_p(Vs):
			VsC[:] = [Vs[0], Vs[1]]
			VsDict = tuple(list(Vs))
		else:
			VsC[:] = [Vs, Vs]
			VsDict = Vs

		if (VsDict, VgDict, VdDict, channelType) in self.gradDict:
			return self.gradDict[(VsDict, VgDict, VdDict, channelType)]

		jac = fetGrad(VdC, VgC, VsC, VbC)

		if(interval_p(Vs) or interval_p(Vg) or interval_p(Vd)):
			grad = self.shape*np.array([np.array([np.nextafter(jac[2][0], float("-inf")), np.nextafter(jac[2][1], float("inf"))]),
							np.array([np.nextafter(jac[1][0], float("-inf")), np.nextafter(jac[1][1], float("inf"))]),
							np.array([np.nextafter(jac[0][0], float("-inf")), np.nextafter(jac[0][1], float("inf"))])])
		else:
			grad = self.shape*np.array([jac[2][0],jac[1][0],jac[0][0]])

		# Clean up gradDict
		self.gradDict[(VsDict, VgDict, VdDict, channelType)] = grad
		if len(self.gradDict) >= 10000:
			self.gradDict = {}
		return grad

	
	# Construct and return lp for short channel mosfet model
	# for interval Vs, Vg and Vd
	# In this case, this is just the hyperrectangle bounds.
	# TODO: Figure out the inflection points and construct tighter
	# linear bounds around short channel ids
	def lp_ids_help(self, Vs, Vg, Vd, channelType):
		
		Ids = self.ids_help(Vs, Vg, Vd, channelType)
		hyperLp = LP()
		if interval_p(Vs):
			hyperLp.ineq_constraint([-1.0, 0.0, 0.0, 0.0], -Vs[0])
			hyperLp.ineq_constraint([1.0, 0.0, 0.0, 0.0], Vs[1])
		else:
			hyperLp.ineq_constraint([-1.0, 0.0, 0.0, 0.0], -(Vs))
			hyperLp.ineq_constraint([1.0, 0.0, 0.0, 0.0], Vs)
		if interval_p(Vg):
			hyperLp.ineq_constraint([0.0, -1.0, 0.0, 0.0], -Vg[0])
			hyperLp.ineq_constraint([0.0, 1.0, 0.0, 0.0], Vg[1])
		else:
			hyperLp.ineq_constraint([0.0, -1.0, 0.0, 0.0], -(Vg))
			hyperLp.ineq_constraint([0.0, 1.0, 0.0, 0.0], Vg)
		if interval_p(Vd):
			hyperLp.ineq_constraint([0.0, 0.0, -1.0, 0.0], -Vd[0])
			hyperLp.ineq_constraint([0.0, 0.0, 1.0, 0.0], Vd[1])
		else:
			hyperLp.ineq_constraint([0.0, 0.0, -1.0, 0.0], -(Vd))
			hyperLp.ineq_constraint([0.0, 0.0, 1.0, 0.0], Vd)
		if interval_p(Ids):
			hyperLp.ineq_constraint([0.0, 0.0, 0.0, -1.0], -(Ids[0]))
			hyperLp.ineq_constraint([0.0, 0.0, 0.0, 1.0], (Ids[1]))
		else:
			hyperLp.ineq_constraint([0.0, 0.0, 0.0, -1.0], -(Ids))
			hyperLp.ineq_constraint([0.0, 0.0, 0.0, 1.0], (Ids))
		return hyperLp;

	# Return lp defining linear constraints containing ids function
	def lp_ids(self, V):
		model = self.model.channelType	
		return(self.lp_ids_help(V[self.s], V[self.g], V[self.d], model))


# @author Mark Greenstreet, Itrat Akhter
# Long channel mosfet model
# s, g, d: Indices indicating source, gate and drain
# model: MosfetModel
# shape: shape
class LcMosfet:
	def __init__(self, s, g, d, model, shape=3.0):
		self.s = s
		self.g = g
		self.d = d
		self.shape = shape
		self.model = model
		self.funDict = {}
		self.gradDict = {}

	# @author Itrat Akhter
	# ids helper function. 
	# Vs, Vg, Vd can be interval or point values
	def ids_help(self, Vs, Vg, Vd, channelType, Vt, ks):
		#print ("Vs", Vs, "Vg", Vg, "Vd", Vd)
		if interval_p(Vs) or interval_p(Vg) or interval_p(Vd):
			Vge = interval_sub(Vg, Vt)
			Vds = interval_sub(Vd, Vs)
			VgseTerm = interval_max(interval_sub(Vge, Vs), 0)
			VgdeTerm = interval_max(interval_sub(Vge, Vd), 0)
			i0 = interval_mult(ks/2.0, interval_sub(interval_mult(VgseTerm, VgseTerm),
												interval_mult(VgdeTerm, VgdeTerm)))
			if self.model.gds != 0.0:
				i_leak = interval_mult(Vds, self.model.gds)
				return interval_add(i0, i_leak)
			else:
				return i0
		else:
			Vge = Vg - Vt
			Vds = Vd - Vs
			VgseTerm = max(Vge - Vs, 0)
			VgdeTerm = max(Vge - Vd, 0)
			i0 = 0.5*ks*(VgseTerm*VgseTerm - VgdeTerm*VgdeTerm)
			i_leak = self.model.gds*Vds
			return i0 + i_leak

	# @author Itrat Akhter
	# ids helper function that takes advantage of monotonicity properties. 
	# Vs, Vg, Vd can be interval or point values
	def ids_help_mon(self, Vs, Vg, Vd, channelType, Vt, ks):
		if(channelType == 'pfet'):
			return interval_neg(self.ids_help(interval_neg(Vs), interval_neg(Vg), interval_neg(Vd), 'nfet', -Vt, -ks))

		elif interval_p(Vs) or interval_p(Vg) or interval_p(Vd):
			VsDict, VgDict, VdDict = None, None, None
			if interval_p(Vd):
				VdDict = tuple(list(Vd))
			else:
				VdDict = Vd

			if interval_p(Vg):
				VgDict = tuple(list(Vg))
			else:
				VgDict = Vg

			if interval_p(Vs):
				VsDict = tuple(list(Vs))
			else:
				VsDict = Vs

			if (VsDict, VgDict, VdDict, channelType) in self.funDict:
				return self.funDict[(VsDict, VgDict, VdDict, channelType)]

			Vs = interval_fix(Vs)
			Vg = interval_fix(Vg)
			Vd = interval_fix(Vd)
			if Vs[1] <= Vd[0]:
				currentVal = np.array([
				self.ids_help(interval_applyRia(np.amax(Vs)), interval_applyRia(np.amin(Vg)), interval_applyRia(np.amin(Vd)), channelType, Vt, ks)[0],
				self.ids_help(interval_applyRia(np.amin(Vs)), interval_applyRia(np.amax(Vg)), interval_applyRia(np.amax(Vd)), channelType, Vt, ks)[1]])
			elif Vd[1] < Vs[0]:
				currentVal = np.array([
				self.ids_help(interval_applyRia(np.amax(Vs)), interval_applyRia(np.amax(Vg)), interval_applyRia(np.amin(Vd)), channelType, Vt, ks)[0],
				self.ids_help(interval_applyRia(np.amin(Vs)), interval_applyRia(np.amin(Vg)), interval_applyRia(np.amax(Vd)), channelType, Vt, ks)[1]])
			else:
				currentVal = np.array([
				self.ids_help(interval_applyRia(np.amax(Vs)), interval_applyRia(np.amax(Vg)), interval_applyRia(np.amin(Vd)), channelType, Vt, ks)[0],
				self.ids_help(interval_applyRia(np.amin(Vs)), interval_applyRia(np.amax(Vg)), interval_applyRia(np.amax(Vd)), channelType, Vt, ks)[1]])

			self.funDict[(VsDict, VgDict, VdDict, channelType)] = currentVal
			# Clean up funDict
			if len(self.funDict) >= 10000:
				self.funDict = {}

		else: currentVal = self.ids_help(Vs, Vg, Vd, channelType, Vt, ks)
		return currentVal
	
	# @author Mark Greenstreet
	# calculate ids from the given V
	def ids(self, V):
		model = self.model
		#print ("V", V)
		#print ("self.s", self.s, "self.g", self.g, "self.d", self.d)
		return(self.ids_help_mon(V[self.s], V[self.g], V[self.d], model.channelType, model.Vt, model.k*self.shape))


	# @author Mark Greenstreet
	# grad_ids: compute the partials of ids wrt. Vs, Vg, and Vd
	#   This function is rather dense.  I would be happier if I could think of
	#    a way to make it more obvious.
	def dg_fun(self, Vs, Vg, Vd, Vt, ks):
		if(Vs[0] > Vd[1]): return None
		Vgse = interval_sub(interval_sub(Vg, np.array([Vs[0], min(Vs[1], Vd[1])])), Vt)
		Vgse[0] = max(Vgse[0], 0)
		Vgse[1] = max(Vgse[1], 0)
		Vds = interval_sub(Vd, Vs)
		Vds[0] = max(Vds[0], 0)
		Vx = np.array([Vg[0] - Vt[1] - Vd[1], Vg[1] - Vt[0] - max(Vs[0], Vd[0])])
		Vx[0] = max(Vx[0], 0)
		Vx[1] = max(Vx[1], 0)
		dg = interval_mult(ks, np.array([min(Vgse[0], Vds[0]), min(Vgse[1], Vds[1])]))
		dd = interval_add(interval_mult(ks, Vx), self.model.gds)
		# print "ks = " + str(ks) + ", gds = " + str(self.model.gds)
		# print "Vgse = " + str(Vgse) + ", Vds = " + str(Vds) + ", Vx = " + str(Vx)
		# print "dg = " + str(dg) + ", dd = " + str(dd)
		return np.array([interval_neg(interval_add(dg, dd)), dg, dd])

	# @author Mark Greenstreet
	def grad_ids_help(self, Vs, Vg, Vd, channelType, Vt, ks):
		if(channelType == 'pfet'):
			# self.grad_ids_help(-Vs, -Vg, -Vd, 'nfet', -Vt, -ks)
			# returns the partials of -Ids wrt. -Vs, -Vg, and -Vd,
			# e.g. (d -Ids)/(d -Vs).  The negations cancel out; so
			# we can just return that gradient.
			return self.grad_ids_help(interval_neg(Vs), interval_neg(Vg), interval_neg(Vd), 'nfet', -Vt, -ks)
		elif(interval_p(Vs) or interval_p(Vg) or interval_p(Vd)):
			VsDict, VgDict, VdDict = None, None, None
			if interval_p(Vd):
				VdDict = tuple(list(Vd))
			else:
				VdDict = Vd

			if interval_p(Vg):
				VgDict = tuple(list(Vg))
			else:
				VgDict = Vg

			if interval_p(Vs):
				VsDict = tuple(list(Vs))
			else:
				VsDict = Vs

			if (VsDict, VgDict, VdDict, channelType) in self.gradDict:
				return self.gradDict[(VsDict, VgDict, VdDict, channelType)]
			Vs = interval_applyRia(Vs)
			Vg = interval_applyRia(Vg)
			Vd = interval_applyRia(Vd)
			Vt = interval_applyRia(Vt)
			g0 = self.dg_fun(Vs, Vg, Vd, Vt, ks)
			g1x = self.dg_fun(Vd, Vg, Vs, Vt, ks)
			if(g1x is None): g1 = None
			else: g1 = np.array([interval_neg(g1x[2]), interval_neg(g1x[1]), interval_neg(g1x[0])])
			if g0 is None: grad = g1
			elif g1 is None: grad = g0
			else: grad = np.array([interval_union(g0[i], g1[i]) for i in range(len(g0))])
			self.gradDict[(VsDict, VgDict, VdDict, channelType)] = grad
			# Clean up funDict
			if len(self.gradDict) >= 10000:
				self.gradDict = {}
			return grad
		elif(Vd < Vs):
			gx = self.grad_ids_help(Vd, Vg, Vs, channelType, Vt, ks)
			return np.array([-gx[2], -gx[1], -gx[0]])
		Vgse = (Vg - Vs) - Vt
		Vds = Vd - Vs
		if(Vgse < 0):  # cut-off: Ids = 0
			return np.array([-self.model.gds, 0.0, self.model.gds])
		elif(Vgse < Vds): # saturation: Ids = (ks/2.0)*Vgse*Vgse
			return np.array([-ks*Vgse - self.model.gds, ks*Vgse, self.model.gds])
		else: # linear: ks*(Vgse - Vds/2.0)*Vds
			dg = ks*Vds
			dd = ks*(Vgse - Vds) + self.model.gds
			return np.array([-(dg + dd), dg, dd])

	# @author Mark Greenstreet
	def grad_ids(self, V):
		model = self.model
		return(self.grad_ids_help(V[self.s], V[self.g], V[self.d], model.channelType, model.Vt, model.k*self.shape))
			
	
	# @author Itrat Akhter
	# This function constructs linear program
	# in terms of src, gate, drain and Ids given the model
	# representing the cutoff, linear or saturation region.
	# Amatrix represents model for quadratic function
	# The quadratic function should be with respect
	# to at most 2 variables - for example, for linear region
	# the function is with respect to 2 variables and for saturation
	# or cutoff it is with respect to 1 variable.
	# For the mosfet case, if it is linear region, the variables
	# are Vg - Vs - Vt and Vd - Vs. For saturation, 
	# the variables are Vd - Vs
	def quad_lin_constraints(self, Amatrix, vertList, Vt):
		if Amatrix.shape[0] > 2:
			raise Exception("quad_lin_constraints: can only accept functions of at most 2 variables")

		# The costs used to find the extreme points
		# in terms of variables for Amatrix should be 
		# derived from the tangents. 
		costs = []
		ds = []
		for vert in vertList:
			grad = 2*np.dot(Amatrix, vert)
			costs.append(list(-grad) + [1])
			dVal = np.dot(np.transpose(vert), np.dot(Amatrix, vert)) - np.dot(grad, vert)
			ds.append(dVal)
		
		lp2d = LP()

		# handle the case where A is neither 
		# positive semidefinite or negative semidefinite
		if Amatrix.shape[0] > 1 and np.linalg.det(Amatrix) < 0:
			# We need to sandwitch the model by the tangents
			# on both sides so add negation of existing costs
			allCosts = [[grad for grad in cost] for cost in costs]
			for cost in costs:
				allCosts.append([-grad for grad in cost])

			for cost in allCosts:
				# Multiply A with coefficient for Ids from cost
				# In all cases, the cost for Ids is costs = [[0,0,1]]
				cA = cost[-1]*Amatrix

				eigVals, eigVectors = np.linalg.eig(cA)

				# sort the eigVals in descending order
				# order eigVectors in the same way
				sortedIndices = np.fliplr([np.argsort(eigVals)])[0]
				sortedEigVals = eigVals[sortedIndices]
				sortedEigVectors = eigVectors[:,sortedIndices]

				v0 = sortedEigVectors[:,0]
		
				# Find the intersection of eigen vector corresponding to positive
				# eigen value with the hyperrectangle
				intersectionPoints = self.intersection([v0[0], v0[1], 0.0], vertList)

				# Now test all the corners and all the intersection points to
				# check which one is the minimum. That is the constant for our cost
				pointsToTest = [point for point in vertList]
				pointsToTest += intersectionPoints

				valToCompare = float("inf")
				
				for point in pointsToTest:
					currentAtPoint = np.dot(np.transpose(point), np.dot(Amatrix, point))
					funVal = np.dot(cost, np.array([point[0], point[1], currentAtPoint]))			
					valToCompare = min(funVal, valToCompare)

				# transform expression in terms of vgse and vds, to Vs, Vg and Vd
				dConst =  -(valToCompare + cost[0]*Vt)

				lp2d.ineq_constraint([-cost[0], -cost[1], -cost[2]], -valToCompare)
		
		else:
			# handle the case where A is positive semi definite or negative semidefinite
			eigVals, eigVectors = np.linalg.eig(Amatrix)

			# positive semidefinite - convex
			cvx_flag = all([eigVals[ei] >= 0 for ei in range(len(eigVals))])

			# add the tangent constraints
			for ci in range(len(costs)):
				cost = costs[ci]
				if len(cost) == 3:
					gradgsdICons2d = [cost[0], cost[1], cost[2], ds[ci]]
				elif len(cost) == 2:
					gradgsdICons2d = [cost[0], 0.0, cost[1], ds[ci]]
				if(cvx_flag):
					gradgsdICons2d = [-grad for grad in gradgsdICons2d]
				
				lp2d.ineq_constraint(gradgsdICons2d[:-1], gradgsdICons2d[-1])


			# take average of cost this is needed for cap constraint
			avgCost = np.zeros((len(costs[0])))
			for cost in costs:
				avgCost += np.array(cost)
			avgCost = avgCost * (1.0/len(costs))

			# now find the additive constant for the cap constraint
			d = None
			for vert in vertList:
				IVal = np.dot(np.transpose(vert), np.dot(Amatrix, vert))
				bb = IVal - np.dot(-avgCost[:-1], vert)
				if(d is None): d = bb
				elif(cvx_flag): d = max(d, bb)
				else: d = min(d, bb)

			if len(cost) == 3:
				gradgsdICons2d = [avgCost[0], avgCost[1], avgCost[2], d]
			elif len(cost) == 2:
				gradgsdICons2d = [avgCost[0], 0.0, avgCost[1], d]

			if not(cvx_flag):
				gradgsdICons2d = [-grad for grad in gradgsdICons2d]

			lp2d.ineq_constraint(gradgsdICons2d[:-1], gradgsdICons2d[-1])


		# deal with leakage term
		if self.model.gds != 0.0:
			leakTermMat = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, -self.model.gds, 1.0]])
			newAList = list(np.dot(np.array(lp2d.A), leakTermMat))
			lp2d.A = newAList

		# Convert LP expressing Ids into a 3 variable LP (Vs, Vg, Vd) from a 2 variable 
		# LP (Vg - Vs - Vt, Vd - Vs)
		lp3d = LP()
		for i in range(len(lp2d.A)):
			lp3d.ineq_constraint([-lp2d.A[i][0] - lp2d.A[i][1], lp2d.A[i][0], lp2d.A[i][1], lp2d.A[i][2]], 
				lp2d.b[i] + lp2d.A[i][0]*Vt)
			
	
		return lp3d

	# @author Itrat Akhter
	# Construct linear constraints for the mosfet model
	# depending on whatever region they belong to - cutoff,
	# saturation or linear - Assume nfet at the point when
	# function is called
	# Vgse = Vg - Vs - Vt
	# Vds = Vd - Vs
	# ks = constant relevant to mosfet
	# hyperLp = hyperrectangle bounds
	def lp_grind(self, Vgse, Vds, Vt, ks, hyperLp):
		Vgse = interval_fix(Vgse)
		Vds = interval_fix(Vds)

		satA = (ks/2.0)*np.array([[1.0]])
		linA = (ks/2.0)*np.array([[0.0, 1.0], [1.0, -1.0]])

		if interval_hi(Vgse) <= 0.0: # cutoff everywhere in the hyperrectangle
			lp = LP()
			return lp.concat(hyperLp)

		elif(interval_lo(Vgse) >= 0 and interval_lo(Vds) >= interval_hi(Vgse)):  # saturation everywhere in the hyperrectangle
			vertList = [np.array([Vgse[0]]), \
					np.array([Vgse[1]])]
			return self.quad_lin_constraints(satA, vertList, Vt).concat(hyperLp)

		elif(interval_lo(Vgse) >= 0 and interval_hi(Vds) <= interval_lo(Vgse)):  # linear everywhere in the hyperrectangle
			vertList = [np.array([Vgse[0], Vds[0]]), \
						np.array([Vgse[1], Vds[0]]), \
						np.array([Vgse[1], Vds[1]]), \
						np.array([Vgse[0], Vds[1]])]
			return self.quad_lin_constraints(linA, vertList, Vt).concat(hyperLp)

		else: # When regions overlap
			# Taking the union of LPs slows down performance so we just return the hyper-rectangle
			# constraints
			return hyperLp

			# take the union of the separate region constraints
			# When taking the union, it is important to bound each variable
			# because of the way the union code works involves solving
			# linear programs and we cannot have degenerate columns in this case
			lp = LP()
			lp.concat(hyperLp)
			newVgse, newVds = Vgse, Vds
			
			# Check if there is a cutoff region. If there is add 
			# linear constraints for cutoff
			if Vgse[0] < 0.0 and Vgse[1] > 0.0:
				newVgse = np.array([Vgse[0], 0.0])
				regionLp = self.lp_grind(newVgse, newVds, Vt, ks, hyperLp)
				lp = lp.union(regionLp, (regionLp.num_constraints() - hyperLp.num_constraints(), regionLp.num_constraints()))
				newVgse = np.array([0.0, Vgse[1]])

			# Find intersection between Vgse = Vds line and hyperrectangle
			# To check if there is overlap between the saturation and linear
			# region
			vertList = [np.array([newVgse[0], newVds[0]]), \
						np.array([newVgse[1], newVds[0]]), \
						np.array([newVgse[1], newVds[1]]), \
						np.array([newVgse[0], newVds[1]])]
			intersectionPoints = self.intersection([-1, 1, 0.0], vertList)

			if len(intersectionPoints) == 0:
				# Deal with just saturation
				regionLp = self.lp_grind(newVgse, newVds, Vt, ks, hyperLp)
				lp = lp.union(regionLp, (regionLp.num_constraints() - hyperLp.num_constraints(), regionLp.num_constraints()))
			else:
				# deal with  ("saturation + linear")
				# find the two polygons caused by the line Vgse = Vds line
				# intersecting newVgse
				vertList1, vertList2 = [], []
				leftI, rightI = intersectionPoints[0], intersectionPoints[1]

				if leftI[1] > newVds[0]:
					vertList1 += [leftI, np.array([newVgse[0], newVds[0]])]
					vertList2 += [np.array([leftI[0]])]
				else:
					vertList1 += [leftI]
					vertList2 += [np.array([newVgse[0]]), np.array([leftI[0]])]

				vertList1.append(np.array([newVgse[1], newVds[0]]))

				if rightI[1] < newVds[1]:
					vertList1 += [rightI]
					vertList2 += [np.array([rightI[0]]), np.array([newVgse[1]])]
				else:
					vertList1 += [np.array([newVgse[1], newVds[1]]),rightI]
					vertList2 += [np.array([rightI[0]])]

				vertList2.append(np.array([newVgse[0]]))
				
				# linear region
				regionLp = self.quad_lin_constraints(linA, vertList1, Vt).concat(hyperLp)
				lp = lp.union(regionLp, (regionLp.num_constraints() - hyperLp.num_constraints(), regionLp.num_constraints()))

				# saturation region
				regionLp = self.quad_lin_constraints(satA, vertList2, Vt).concat(hyperLp)
				lp = lp.union(regionLp, (regionLp.num_constraints() - hyperLp.num_constraints(), regionLp.num_constraints()))

			return lp

	# @author Mark Greenstreet, Itrat Akhter
	# Return an lp formed around ids for intervals given by
	# Vs, Vg and Vd
	def lp_ids_help(self, Vs, Vg, Vd, channelType, Vt, ks):
		Vgs = interval_sub(Vg, Vs)
		Vgse = Vgs - Vt
		Vds = interval_sub(Vd, Vs)

		#print ("Vs", Vs, "Vg", Vg, "Vd", Vd)
		Ids = self.ids_help_mon(Vs, Vg, Vd, channelType, Vt, ks)
		#print ("Ids", Ids)
		# form the LP where the constraints for the Vg and Vd and Vs
		# bounds are added. This is important to make sure we don't get
		# unboundedness from our linear program
		hyperLp = LP()
		if interval_p(Vs):
			hyperLp.ineq_constraint([-1.0, 0.0, 0.0, 0.0], -Vs[0])
			hyperLp.ineq_constraint([1.0, 0.0, 0.0, 0.0], Vs[1])
		else:
			hyperLp.ineq_constraint([-1.0, 0.0, 0.0, 0.0], -(Vs))
			hyperLp.ineq_constraint([1.0, 0.0, 0.0, 0.0], Vs)
		if interval_p(Vg):
			hyperLp.ineq_constraint([0.0, -1.0, 0.0, 0.0], -Vg[0])
			hyperLp.ineq_constraint([0.0, 1.0, 0.0, 0.0], Vg[1])
		else:
			hyperLp.ineq_constraint([0.0, -1.0, 0.0, 0.0], -(Vg))
			hyperLp.ineq_constraint([0.0, 1.0, 0.0, 0.0], Vg)
		if interval_p(Vd):
			hyperLp.ineq_constraint([0.0, 0.0, -1.0, 0.0], -Vd[0])
			hyperLp.ineq_constraint([0.0, 0.0, 1.0, 0.0], Vd[1])
		else:
			hyperLp.ineq_constraint([0.0, 0.0, -1.0, 0.0], -(Vd))
			hyperLp.ineq_constraint([0.0, 0.0, 1.0, 0.0], Vd)
		if interval_p(Ids):
			hyperLp.ineq_constraint([0.0, 0.0, 0.0, -1.0], -(Ids[0]))
			hyperLp.ineq_constraint([0.0, 0.0, 0.0, 1.0], (Ids[1]))
		else:
			hyperLp.ineq_constraint([0.0, 0.0, 0.0, -1.0], -(Ids))
			hyperLp.ineq_constraint([0.0, 0.0, 0.0, 1.0], (Ids))

		#return hyperLp

		# If hyperrectangles are tiny do not try to construct 
		# linear convex regions and just return the lp representing
		# the bounds
		if interval_r(Vs) < 1e-9 or interval_r(Vd) < 1e-9:
			return hyperLp
		
		if(not(interval_p(Vs) or interval_p(Vg) or interval_p(Vd))):
			# Vs, Vg, and Vd are non-intervals -- they define a point
			# Add an inequality that our Ids is the value for this point
			return(LP(None, None, None, [[0,0,0,1.0]], self.ids_help_mon(Vs, Vg, Vd, channelType, Vt, ks)))
		elif(channelType == 'pfet'):
			LPnfet = self.lp_ids_help(interval_neg(Vs), interval_neg(Vg), interval_neg(Vd), 'nfet', -Vt, -ks)
			return LPnfet.neg_A()
		elif((interval_lo(Vs) < interval_hi(Vd)) and (interval_hi(Vs) > interval_lo(Vd))):
			# If the Vs interval overlaps the Vd interval, then Vds can change sign.
			# That means we've got a saddle.  We won't try to generated LP constraints
			# for the saddle.  So, we just return a hyperLp.
			return hyperLp
		elif(interval_lo(Vs) >= interval_hi(Vd)):
			LPswap = self.lp_ids_help(Vd, Vg, Vs, channelType, Vt, ks)
			A = []
			for i in range(len(LPswap.A)):
				row = LPswap.A[i]
				A.append([row[2], row[1], row[0], -row[3]])
			Aeq = []
			for i in range(len(LPswap.Aeq)):
				row = LPswap.Aeq[i]
				Aeq.append([row[2], row[1], row[0], -row[3]])
			return LP(LPswap.c, A, LPswap.b, Aeq, LPswap.beq)

		else:
			return self.lp_grind(Vgse, Vds, Vt, ks, hyperLp)

	# @author Itrat Akhter
	# Find the intersection points between the line
	# represented by [a, b, c] where a, b, and c represent
	# the line ax + by = c and polygon represented by vertList
	# a list of arrays where each array represents the vertex
	def intersection(self, line, vertList):
		intersectionPoints = []
		for vi in range(len(vertList)):
			vert1 = vertList[vi]
			vert2 = vertList[(vi + 1)%len(vertList)]
			minVert = np.minimum(vert1, vert2)
			maxVert = np.maximum(vert1, vert2)
			if vert2[0] - vert1[0] == 0:
				x0 = vert2[0]
				x1 = (line[2] - line[0]*x0)/line[1]
				if x1 >= minVert[1] and x1 <= maxVert[1]:
					intersectionPoints.append(np.array([x0, x1]))
			else:
				m = (vert2[1] - vert1[1])/(vert2[0] - vert1[0])
				c = vert1[1] - m*vert1[0]
				x0 = (line[2] - line[1]*c)/(line[0] + line[1]*m)
				x1 = m*x0 + c
				if x0 >= minVert[0] and x0 <= maxVert[0] and x1 >= minVert[1] and x1 <= maxVert[1]:
					intersectionPoints.append(np.array([x0, x1]))

		return intersectionPoints

	# @ author Mark Greenstreet
	# Return an lp
	def lp_ids(self, V, numDivisions = None, mainDict = None):
		model = self.model
		idsLp = self.lp_ids_help(V[self.s], V[self.g], V[self.d], model.channelType, model.Vt, model.k*self.shape)
		return idsLp


# @author Mark Greenstreet, Itrat Akhter
# Circuit combining different transistor models
class Circuit:
	def __init__(self, tr):
		self.tr = tr # list of transistor objects (long channel - Mosfet or short channel - StMosfet)

	# @author Mark Greenstreet
	# Return node currents given the voltages at nodes
	def f(self, V):
		intervalVal = any([interval_p(x) for x in V])
		if intervalVal:
			I_node = np.zeros((len(V),2))
		else:
			I_node = np.zeros(len(V))
		for i in range(len(self.tr)):
			tr = self.tr[i]
			Ids = tr.ids(V)
			#print "Circuit.f: i = " + str(i) + ", tr.s = " + str(tr.s) + "(" + str(V[tr.s]) + "), tr.g = " + str(tr.g) + "(" + str(V[tr.g]) + "), tr.d = " + str(tr.d) + "(" + str(V[tr.d]) + "), ids = " + str(Ids)
			I_node[tr.s] = interval_add(I_node[tr.s], Ids)
			I_node[tr.d] = interval_sub(I_node[tr.d], Ids)
		return I_node

	# @author Mark Greenstreet
	# Return jacobian of node currents with respect to node voltages
	def jacobian(self, V):
		#print ("Calculating jacobian for V")
		#print (V)
		intervalVal = any([interval_p(x) for x in V])	
		
		if intervalVal:
			J = np.zeros([len(V), len(V), 2])
		else:
			J = np.zeros([len(V), len(V)])
		for i in range(len(self.tr)):
			tr = self.tr[i]
			g = tr.grad_ids(V)
			#print 'i = ' + str(i) + ', tr.s = ' + str(tr.s) + ', tr.g = ' + str(tr.g) + ', tr.d = ' + str(tr.d) + ', g = ' + str(g)
			sgd = [tr.s, tr.g, tr.d]
			for i in range(len(sgd)):
				J[tr.s, sgd[i]] = interval_add(J[tr.s, sgd[i]], g[i])
				J[tr.d, sgd[i]] = interval_sub(J[tr.d, sgd[i]], g[i])
		return J

	

	# @author Mark Greenstreet, Itrat Akhter
	# The lp we return has one variable for each node and one variable for each transistor.
	# For 0 <= i < #nodes, variable[i] is the voltage on node i.
	# For 0 <= j < #transistors, variable[#nodes+j] is Ids for transistor j
	# This function collects all the linear constraints related to the 
	# transistors and set the node currents to 0
	# grndPower Index is the list of indices indicating
	# ground and power in V. We need this when we are setting sum of node currents to 0
	# because during the summation we want to avoid the ground and power nodes
	def lp(self, V, grndPowerIndex):
		lp = LP()
		n_nodes = len(V)
		n_tr = len(self.tr)
		nvars = len(V) + n_tr

		#print ("nvars", nvars)
		#print ("n_tr", n_tr)

		eqCoeffs = np.zeros((n_nodes, nvars))

		# keep an array of nvars with boolean values
		# if they haven't been dealt with yet, add the hyper constraint
		# for that specific variable
		hyperAdded = [False]*nvars
		for i in range(n_tr):
			#print ("transistor number", i)
			tr = self.tr[i]
			lptr = tr.lp_ids(V)

			# remove the hyper constraints that we have added 
			# (hyper constraints have been repeated as these are
			# the hyper constraints from each transistor. We had needed
			# to add these if we created unions of two lps and avoid unboundedness).
			# We need
			# to avoid duplicate constraints so that we do not slow down
			# the linear program solver. We will add the hyper constraints
			# for each variable exactly once in a bit. 
			# TODO: Construct a better way to do this
			modLptr = LP()
			for ai in range(6,len(lptr.A)):
				modLptr.ineq_constraint(lptr.A[ai], lptr.b[ai])
			
			lp.concat(modLptr.varmap(nvars, [tr.s, tr.g, tr.d, n_nodes+i]))
			
			# Add the hyper constraints for tr.s,tr.g and tr.d if they have not been
			# added yet. This makes sure that we have added hyper constraints for 
			# each variable exactly once
			for ind in [tr.s, tr.g, tr.d]:
				if not(hyperAdded[ind]):
					hyperAdded[ind] = True
					greatConstr = np.zeros((nvars))
					greatConstr[ind] = -1.0
					lessConstr = np.zeros((nvars))
					lessConstr[ind] = 1.0
					if interval_p(V[ind]):
						lp.ineq_constraint(greatConstr, -V[ind][0])
						lp.ineq_constraint(lessConstr, V[ind][1])
					else:
						lp.ineq_constraint(greatConstr, -V[ind])
						lp.ineq_constraint(lessConstr, V[ind])

			eqCoeffs[tr.s, n_nodes + i] += 1.0
			eqCoeffs[tr.d, n_nodes + i] += -1.0


		# need to add equality constraints that the sum of the currents into each node is zero
		for i in range(n_nodes):
			if all([i != gpi for gpi in grndPowerIndex]):
				lp.ineq_constraint(list(-eqCoeffs[i]), 0.0)
				lp.ineq_constraint(list(eqCoeffs[i]), 0.0)

		#print ("lp")
		#print (lp)
		return lp



	# @author Itrat Akhter
	# This function solves the linear program 
	# returns a tighter hyperrectangle if feasible
	# grndPower Index is the list of indices indicating
	# ground and power in V
	def linearConstraints(self, V, grndPowerIndex):
		lp = self.lp(V, grndPowerIndex)
		n_nodes = len(V)
		n_tr = len(self.tr)
		nvars = len(V) + n_tr
		#print ("nvars", nvars)
	
		tighterHyper = [x for x in V]
		feasible = True
		numSuccessLp, numUnsuccessLp, numTotalLp = 0, 0, 0
		for i in range(n_nodes):
			if interval_p(tighterHyper[i]):
				numTotalLp += 2
				cost = np.zeros((nvars))

				#minimize variable i
				cost[i] = 1.0
				lp.add_cost(list(cost))
				minSol = lp.solve()

				#maximize variable i
				cost[i] = -1.0
				lp.add_cost(list(cost))
				maxSol = lp.solve()

				if minSol is None or maxSol is None:
					numUnsuccessLp += 2
					continue

				if minSol["status"] == "primal infeasible" and maxSol["status"] == "primal infeasible":
					numSuccessLp += 2
					feasible = False
					break

				if minSol["status"] == "optimal":
					tighterHyper[i][0] = minSol['x'][i]
					numSuccessLp += 1
				else:
					numUnsuccessLp += 1
				if maxSol["status"] == "optimal":
					tighterHyper[i][1] = maxSol['x'][i]
					numSuccessLp += 1
				else:
					numUnsuccessLp += 1


		#print ("numTotalLp", numTotalLp, "numSuccessLp", numSuccessLp, "numUnsuccessLp", numUnsuccessLp)
		return [feasible, tighterHyper, numTotalLp, numSuccessLp, numUnsuccessLp]

