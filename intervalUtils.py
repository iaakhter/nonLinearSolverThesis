# Functions implementing interval verification algorithm - the Krawczyk
# operator and its helper functions
# @author Itrat Ahmed Akhter

import numpy as np
import random
import math
from intervalBasics import *


'''
Multiply 2 matrices and return the resulting matrix. 
Any matrix can be an interval matrix
'''
def multiplyMats(mat1, mat2):
	isInterval = len(mat1.shape) == 3 or len(mat2.shape) == 3
	if isInterval:
		result = np.zeros((mat1.shape[0],mat2.shape[1],2))
	else:
		result = np.zeros((mat1.shape[0],mat2.shape[1]))
	
	for i in range(mat1.shape[0]):
		for j in range(mat2.shape[1]):
			if isInterval:
				intervalVal = np.zeros(2)
			else:
				intervalVal = 0.0
			for k in range(mat2.shape[1]):
				intervalVal += interval_mult(mat1[i,k],mat2[k,j])
				intervalVal = interval_round(intervalVal)
			result[i,j] = interval_round(intervalVal)

	return result


'''
Subtract 2 matrices. Either of the matrix can be an interval matrix
and return the resulting matrix
'''
def subtractMats(mat1, mat2):
	isInterval = len(mat1.shape) == 3 or len(mat2.shape) == 3
	if isInterval:
		result = np.zeros((mat1.shape[0],mat2.shape[1],2))
	else:
		result = np.zeros((mat1.shape[0],mat2.shape[1]))

	for i in range(result.shape[0]):
		for j in range(result.shape[1]):
			result[i,j] = interval_sub(mat1[i,j], mat2[i,j])
	return result

'''
Multiply interval or regular matrix with interval or regular vector
and regular or interval vector depending on whether mat and vec 
are interval or regular
'''
def multiplyMatWithVec(mat,vec):
	isInterval = interval_p(mat[0,0]) or interval_p(vec[0])

	if isInterval:
		result = np.zeros((mat.shape[0],2))
	else:
		result = np.zeros((mat.shape[0],1))
	for i in range(mat.shape[0]):
		if isInterval:
			intervalVal = np.zeros((2))
		else:
			intervalVal = np.zeros((1))
		for j in range(mat.shape[1]):
			mult = interval_mult(mat[i,j],vec[j])
			intervalVal += mult
			intervalVal = interval_round(intervalVal)

		result[i,:] = interval_round(intervalVal)
	
	return result


'''
Turn a regular matrix into an interval matrix by 
subtracting ulp from each element to create a lower bound and 
adding ulp to each element to create an upper bound
'''
def turnRegMatToIntervalMat(mat):
	intervalMat = np.zeros((mat.shape[0], mat.shape[1], 2))
	for i in range(mat.shape[0]):
		for j in range(mat.shape[1]):
			intervalMat[i,j] = np.array([np.nextafter(mat[i,j], float("-inf")), np.nextafter(mat[i,j], float("inf"))])

	return intervalMat

'''
Turn a regular vector into an interval vector by 
subtracting ulp from each element to create a lower bound and 
adding ulp to each element to create an upper bound
'''
def turnRegVecToIntervalVec(vec):
	intervalVec = np.zeros((len(vec), 2))
	for i in range(len(vec)):
		intervalVec[i] = np.array([np.nextafter(vec[i], float("-inf")), np.nextafter(vec[i], float("inf"))])

	return intervalVec


'''
Return the volume of the hyperrectangle
'''
def volume(hyperRectangle):
	if hyperRectangle is None:
		return None
	vol = 1
	for i in range(hyperRectangle.shape[0]):
		vol *= (hyperRectangle[i,1] - hyperRectangle[i,0])
	return vol

'''
Apply newton's method to find a solution using
function defined by model
@param model defines the problem
@param soln the starting point for Newton's method
@return (False, soln) if Newton's method doesn't find a solution
					within the bounds defined by model
@return (True, soln) if Newton's method can find a solution within
					the bounds defined by model
'''
def newton(model,soln,normThresh=1e-8):
	h = soln
	count = 0
	maxIter = 100
	bounds = model.bounds
	lenV = len(soln)
	overallHyper = np.zeros((lenV,2))
	for i in range(lenV):
		overallHyper[i,0] = bounds[i][0]
		overallHyper[i,1] = bounds[i][1]
	while count < maxIter and (np.linalg.norm(h) > normThresh or count == 0):
		res = model.f(soln)
		res = -np.array(res)
		jac = model.jacobian(soln)
		try:
			h = np.linalg.solve(jac,res)
		except np.linalg.LinAlgError:
			h = np.linalg.lstsq(jac, res)[0]
		soln = soln + h
		if np.less(soln, overallHyper[:,0] - 0.001).any() or np.greater(soln, overallHyper[:,1]+0.001).any():
			return (False, soln)
		count+=1
	if count >= maxIter and np.linalg.norm(h) > normThresh:
		return(False, soln)
	return (True,soln)



'''
Do a krawczyk update on hyperrectangle defined by startBounds
@param startBounds hyperrectangle
@param jacInterval interval jacobian over startBounds
@param samplePoint mid point in startBounds
@param fSamplePoint function evaluation at samplePoint
@param jacSamplePoint jacobian at samplePoint
@return (True, refinedHyper) if hyperrectangle contains a unique solution.
		refinedHyper also contains the solution and might be smaller
		than hyperRectangle
@return (False, refinedHyper) if hyperrectangle may contain more
		than one solution. refinedHyper also contains all the solutions
		that hyperRectangle might contain and might be smaller than
		hyperRectangle
@return (False, None) if hyperrectangle contains no solution
'''
def krawczykHelp(startBounds, jacInterval, samplePoint, fSamplePoint, jacSamplePoint):
	numV = startBounds.shape[0]
	I = np.identity(numV)
	

	'''print ("startBounds")
	printHyper(startBounds)
	print ("samplePoint", samplePoint)
	print ("fSamplePoint", fSamplePoint)
	print ("jacInterval", jacInterval)
	print ("jacSamplePoint", jacSamplePoint)'''


	try:
		C = np.linalg.inv(jacSamplePoint)
	except:
		# In case jacSamplePoint is singular
		C = np.linalg.pinv(jacSamplePoint)

	#print ("C", C)
	#print ("fSamplePoint", fSamplePoint)
	C_fSamplePoint = multiplyMatWithVec(C,fSamplePoint)
	#print ("C_fSamplePoint", C_fSamplePoint)

	C_jacInterval = multiplyMats(C,jacInterval)

	#print ("C_jacInterval", C_jacInterval)

	I_minus_C_jacInterval = subtractMats(I,C_jacInterval)

	#print ("I_minus_C_jacInterval", I_minus_C_jacInterval)
	

	xi_minus_samplePoint = np.zeros((numV, 2))
	for i in range(numV):
		xi_minus_samplePoint[i] = interval_sub(startBounds[i], samplePoint[i])

	#print ("xi_minus_samplePoint", xi_minus_samplePoint)
	lastTerm = multiplyMatWithVec(I_minus_C_jacInterval, xi_minus_samplePoint)

	#print ("lastTerm", lastTerm)
	
	kInterval = np.zeros((numV,2))

	for i in range(numV):
		kInterval[i,:] = interval_add(interval_sub(samplePoint[i], C_fSamplePoint[i]), lastTerm[i])

	#print ("startBounds")
	#printHyper(startBounds)
	#print ("kInterval")
	#printHyper(kInterval)
	# if kInterval is in the interior of startBounds, found a unique solution
	if(all([ (interval_lo(kInterval[i]) > interval_lo(startBounds[i])) and (interval_hi(kInterval[i]) < interval_hi(startBounds[i]))
		 for i in range(numV) ])):
		return (True, kInterval)

	
	intersect = np.zeros((numV,2))
	for i in range(numV):
		#print ("i", i)
		#print ("startBounds", startBounds[i][0], startBounds[i][1])
		#print ("kInterval", kInterval[i][0], kInterval[i][1])
		intersectVar = interval_intersect(kInterval[i], startBounds[i])
		if intersectVar is not None:
			intersect[i] = intersectVar
		else:
			# no solution
			return (False, None)

	return (False, intersect)



'''Print hyperrectangle hyper'''
def printHyper(hyper):
	for i in range(hyper.shape[0]):
		print (hyper[i,0], hyper[i,1])


'''
Find a Newton's solution in hyper. If there exists a solution
inflate current hyper so that solution is in centre and
check the inflated hyper with Krawczyk
@param model defines the problem
@param hyper hyperRectangle
@param epsilonInflation indicates the proportion of hyper-rectangle distance by which the 
 	hyper-rectangle needs to be inflated before the Krawczyk operator is applied
@return Krawczyk result of inflated hyper
'''
def checkInflatedHyper(model, hyper, epsilonBounds):
	startBounds = np.copy(hyper)
	prevIntersect, intersect = newtonInflation(model, startBounds, epsilonBounds)

	if prevIntersect is None:
		return [False, hyper]

	startBounds = np.copy(intersect)
	samplePoint = (startBounds[:,0] + startBounds[:,1])/2.0
	fSamplePoint = np.array(model.f(samplePoint))
	jacSamplePoint = model.jacobian(samplePoint)
	jacInterval = model.jacobian(startBounds)
	
	# Krawczyk update
	kHelpResult = krawczykHelp(startBounds, jacInterval, samplePoint, fSamplePoint, jacSamplePoint)

	if kHelpResult[0]:
		return [True, kHelpResult[1]]
	else:
		return [False, hyper]


'''
Check whether hyperrectangle hyperRectangle contains
a unique solution, no solution or maybe more than one solution
to function identified by the model with Krawczyk operator. 
Use rounded interval arithmetic for every operation in the 
Krawczyk update
@param model defines the problem
@param hyperRectangle the hyperRectangle
@param alpha indicates how many times the Krawczyk operator is 
		used to refine hyperRectangle before the function returns
		If the reduction in volume is below alpha, then we are done
@param epsilonInflation the amount by which either side of the hyper-rectangle
		is inflated before applying the Krawczyk method. This allows for quicker
		convergence to a unique solution if one exists
@return (True, refinedHyper) if hyperrectangle contains a unique solution.
		refinedHyper also contains the solution and might be smaller
		than hyperRectangle
@return (False, refinedHyper) if hyperrectangle may contain more
		than one solution. refinedHyper also contains all the solutions
		that hyperRectangle might contain and might be smaller than
		hyperRectangle
@return (False, None) if hyperrectangle contains no solution
'''
def checkExistenceOfSolution(model,hyperRectangle, alpha = 1.0, epsilonInflation=0.01):
	epsilonBounds = 1e-12
	numV = len(hyperRectangle[0])

	startBounds = np.zeros((numV,2))
	startBounds[:,0] = hyperRectangle[0,:]
	startBounds[:,1] = hyperRectangle[1,:]

	# First do an interval arithmetic test
	# Calculate the interval evaluation of the function
	# for hyperrectangle. If any component of the result
	# does not contain zero then the hyperrectangle does not
	# contain any solution
	if hasattr(model, 'f'):
		#print ("startBounds", startBounds)
		funVal = model.f(startBounds)
		#print ("funVal", funVal)
		if(any([np.nextafter(funVal[i,0], np.float("-inf"))*np.nextafter(funVal[i,1], np.float("inf")) > np.nextafter(0.0, np.float("inf")) 
				for i in range(numV)])):
			return [False, None]

	# Start the Krawczyk update
	constructBiggerHyper = False
	iteration = 0
	prevIntersect = None
	while True:
		oldVolume = volume(startBounds)
		#print ("startBounds before")
		#printHyper(startBounds)
		dist = startBounds[:,1] - startBounds[:,0]
		startBounds[:,0] = startBounds[:,0] - (epsilonInflation*dist + epsilonBounds)
		startBounds[:,1] = startBounds[:,1] + (epsilonInflation*dist + epsilonBounds)
	
		#print ("startBounds after")
		#printHyper(startBounds)
		samplePointSing = (startBounds[:,0] + startBounds[:,1])/2.0
		samplePoint = turnRegVecToIntervalVec(samplePointSing)
		fSamplePoint = np.array(model.f(samplePoint))
		jacSamplePoint = model.jacobian(samplePointSing)
		jacInterval = model.jacobian(startBounds)
		
		# Krawczyk update
		kHelpResult = krawczykHelp(startBounds, jacInterval, samplePoint, fSamplePoint, jacSamplePoint)

		if kHelpResult[0] or kHelpResult[1] is None:
			return kHelpResult
		
		intersect = kHelpResult[1]
		#print("intersect")
		#printHyper(intersect)

		newVolume = volume(intersect)
		volReduc = (oldVolume - newVolume)/(oldVolume*1.0)
		#print ("volReduc", volReduc)

		# If the reduction of volume is less than equal to alpha
		# then do no more Krawczyk updates. We are done
		if (math.isnan(volReduc) or volReduc <= alpha):
			intersect[:,0] = np.maximum(hyperRectangle[0,:], intersect[:,0])
			intersect[:,1] = np.minimum(hyperRectangle[1,:], intersect[:,1])
			return [False,intersect]
		else:
			startBounds = intersect

		iteration += 1




'''
Use newton's method to find a solution in hyperrectangle hyper.
If a solution exists then construct a hyperrectnagle with the 
solution in center and contains hyper. Return the original hyper
and the inflated hyper.
@param model defines the problem
@param hyper hyperRectangle
@epsilonBounds defines how far away from the hyperrectangle hyper
we can allow the newton solution to be to start the hyperrectangle 
inflation process
@return (hyper, inflatedHyper) if there is a Newton solution within
epsilon bounds of hyper. Otherwise return (None, hyper)
'''
def newtonInflation(model, hyper, epsilonBounds):
	numV = hyper.shape[0]
	exampleVolt = (hyper[:,0] + hyper[:,1])/2.0
	soln = newton(model,exampleVolt)
	prevHyper = None
	solnInHyper = True
	if soln[0]:
		# the new hyper must contain the solution in the middle and enclose old hyper
		#print ("soln ", soln[1][0], soln[1][1], soln[1][2])
		for sl in range(numV):
			if (soln[1][sl] < hyper[sl][0] - epsilonBounds or
				soln[1][sl] > hyper[sl][1] + epsilonBounds):
				solnInHyper = False
				break
	if soln[0] and solnInHyper:
		prevHyper = np.copy(hyper)
		for si in range(numV):
			maxDiff = max(abs(hyper[si,1] - soln[1][si]), abs(soln[1][si] - hyper[si,0]))
			#print ("si", si, "maxDiff", maxDiff)
			if maxDiff < epsilonBounds:
				maxDiff = epsilonBounds
			hyper[si,0] = soln[1][si] - maxDiff
			hyper[si,1] = soln[1][si] + maxDiff

	return (prevHyper, hyper)


