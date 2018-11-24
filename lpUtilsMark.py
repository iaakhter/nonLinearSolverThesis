# @author: Mark Greenstreet, Itrat Akhter
# Class to help construct, manipulate and solve linear programs

import numpy as np
from cvxopt import matrix,solvers
from scipy.spatial import ConvexHull
from intervalBasics import *


class LP:
	# We'll store the A and Aeq matrices by rows: each row corresponds
	#   to a single constraint.  Each matrix is a python list of rows,
	#   each row is a list of numbers.
	#   Note that cvxopt expects the matrices to stored by colums -- I
	#   guess there's some Fortran in the genealogy of cvxopt.  They also
	#   use scipy matrices.  We'll have to convert when calling the cvxopt
	#   routines.
	def __init__(self, c=None, A=None, b=None, Aeq=None, beq=None):
		if(c is None):
			self.c = None
		else:
			self.c = [ e for e in c ]
		if(A is None):
			self.A = []
		else:
			self.A = [ [ e for e in row ] for row in A]
		if(b is None):
			self.b = []
		else:
			self.b = [ e for e in b ]
		if(Aeq is None):
			self.Aeq = []
		else:
			self.Aeq = [ [ e for e in row ] for row in Aeq]
		if(beq is None):
			self.beq = []
		else:
			self.beq = [ e for e in beq ]

	def num_constraints(self):
		return len(self.A) + len(self.Aeq)

	def __str__(self):
		ineqAStr = ""
		for ineq in self.A:
			ineqAStr += str(np.array(ineq)) + "\n"
		return "ineqA\n" + ineqAStr + "\n" + "ineqB\n" + str(np.array(self.b)) + "\n" + \
				"eqA\n" + str(np.array(self.Aeq)) + "\n" + "eqB\n" + str(np.array(self.beq)) + "\n" + \
				"cost\n" + str(np.array(self.c)) + "\n"

	# @author: Mark Greenstreet
	# concat: add the constraints from LP2 to our constraints.
	# We keep our cost vector and ignore LP2.c.  We should probably
	# check to makes sure that LP2 has the same number of variables as
	# we do, but that will be a later enhancement.
	# concat modifies self.
	def concat(self, LP2):
		self.A = self.A + LP2.A
		self.b = self.b + LP2.b
		self.Aeq = self.Aeq + LP2.Aeq
		self.beq = self.beq + LP2.beq
		return(self)

	# @author: Mark Greenstreet
	# sometimes, we want to replace a variable with its negation (i.e.
	#   when using nfet code to model pfet currents).  That means we need
	#   to negate the elements of A and Aeq.
	def neg_A(self):
		nA = [[-e for e in row] for row in self.A]
		nAeq = [[-e for e in row] for row in self.Aeq]
		return LP(self.c, nA, self.b, nAeq, self.beq)

	# @author:Mark Greenstreet
	# add an equality constraint
	# eq_constraint modifies self.
	def eq_constraint(self, aeq, beq):
		self.Aeq.append(aeq)
		self.beq.append(beq)

	# @author:Mark Greenstreet
	# add an inequality constraint
	# ineq_constraint modifies self.
	def ineq_constraint(self, a, b):
		self.A.append(a)
		self.b.append(b)


	# @author: Itrat Akhter
	# Use inequality constraints of 
	# one LP (self) as costs (both min and max) of other LP
	# to construct constraints that satisfy both LP
	# return a new LP from the new constraints.
	# otherLpExcludingInds indicates the indices of the inequality
	# constraints in self that are not considered as cost for consideration
	def constraint_as_cost(self, otherLp, otherLpExcludingInds):
		newLp = LP()
		# Use the normals of inequalities of
		# self's as cost (minimize and maximize) 
		# for other LP, solve the other LP to
		# find the constraint (from self's, minimization and maximization)
		# that satisfies both the LPs
		validConstraints = []
		for i in range(len(self.A)):
			if otherLpExcludingInds is None or not(i>= otherLpExcludingInds[0] and i < otherLpExcludingInds[1]):
				minCost = [ val*1.0 for val in self.A[i]]
				maxCost = [-val*1.0 for val in self.A[i]]
				possibleValidConstraints = []
				possibleBs = []
				
				# list containing [x, y] where x and y
				# represent a constraint to be added to
				# A and b respectively
				possibleValidConstraints.append([[ val for val in self.A[i]], self.b[i]])
				possibleBs.append(self.b[i])
				

				# Maximize
				otherLp.add_cost(maxCost)
				maxSol = otherLp.solve()
				if maxSol is not None and maxSol["status"] == "optimal":
					maxB = np.dot(np.array(maxCost), np.array(maxSol['x']))[0]
					possibleValidConstraints.append([minCost, -maxB])
					possibleBs.append(-maxB)
				

				# Add the constraint with the highest constant to the newLp
				maxB = max(possibleBs)
				maxBindex = possibleBs.index(maxB)
				if maxSol is not None and maxSol["status"] == "optimal":
					newLp.ineq_constraint(possibleValidConstraints[maxBindex][0], possibleValidConstraints[maxBindex][1])

		
		return newLp


	
	# @author: Itrat Akhter
	# Create a union of self and another LP
	# The union of LPs should satisfy both
	# the LPs
	def union(self, otherLp, otherLpExcludingInds):
		newLp = LP()

		# Use inequality constraints of self as 
		# costs of otherLp to find new constraints
		# that satisfy self and otherLP
		newLp.concat(self.constraint_as_cost(otherLp, None))

		# Use inequality constraints of otherLp as 
		# costs of self to find new constraints
		# that satisfy self and otherLP
		newLp.concat(otherLp.constraint_as_cost(self, otherLpExcludingInds))

		# add the equality constraints
		for i in range(len(self.Aeq)):
			newLp.eq_constraint([x for x in self.Aeq[i]], self.beq[i])
		for i in range(len(otherLp.Aeq)):
			newLp.eq_constraint([x for x in otherLp.Aeq[i]], otherLp.beq[i])

		
		return newLp

	# @author: Itrat Akhter
	# This will replace the current cost
	# function. 
	def add_cost(self, c):
		self.c = c

	# @author: Itrat Akhter
	# return cvxopt solution after solving lp
	def solve(self):
		if self.num_constraints() == 0:
			return None
		cocantenatedA = [ e for e in self.A ]
		for eqConstraint in self.Aeq:
			cocantenatedA.append(eqConstraint)
			negConstraint = [ -x for x in eqConstraint]
			cocantenatedA.append(negConstraint)

		cocantenatedb = [ e for e in self.b ]
		for eqb in self.beq:
			cocantenatedb.append(eqb)
			cocantenatedb.append(-eqb)

		AMatrix = matrix(np.array(cocantenatedA))
		bMatrix = matrix(cocantenatedb)
		cMatrix = matrix(self.c)
		solvers.options["show_progress"] = False

		try:
			sol = solvers.lp(cMatrix, AMatrix, bMatrix)
			return sol
		except ValueError:
			return None

	# @author: Itrat Akhter
	# Calculate the slack for each inequality constraint
	def slack(self):
		sol = self.solve()
		if sol is None:
			raise Exception('LP:slack - could not solve LP')

		if sol["status"] == "primal infeasible":
			raise Exception('LP:slack - LP infeasible')

		calculatedB = np.dot(np.array(self.A), np.array(sol['x']))
		slack = np.array(self.b) - np.transpose(calculatedB)
		return slack


	# @author: Mark Greenstreet
	# Map thre current LP to a new LP where
	# the new Lp has nvar number of variables
	# vmap indicates to which index in the new LP
	# the variables of the current LP will correspond to
	def varmap(self, nvar, vmap):
		n_ineq = len(self.A)
		A = []
		for i in range(n_ineq):
			r = [0 for j in range(nvar)]
			for k in range(len(vmap)):
				r[vmap[k]] = self.A[i][k]
			A.append(r)
		n_eq = len(self.Aeq)
		Aeq = []
		for i in range(n_eq):
			r = [0 for j in range(nvar)]
			for k in range(len(vmap)):
				r[vmap[k]] = self.Aeq[i][k]
			Aeq.append(r)
		return LP(self.c, A, self.b, Aeq, self.beq)

