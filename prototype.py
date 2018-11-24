# @author Itrat Ahmed Akhter
# Main file containing implementation of our solver

import numpy as np
import time
import intervalUtils
from intervalBasics import *
from circuitModels import RambusTanh, RambusMosfet
from circuitModels import SchmittMosfet
from circuitModels import InverterTanh, InverterMosfet
from circuitModels import InverterLoopTanh, InverterLoopMosfet
import dcUtils
import random
import math
import circuit


# A mechanism to bisect a hyperrectangle
# Find the dimension in hyper with the highest
# length and bisect the hyperrectangle at that dimension
def bisectMax(hyper, options=None):
	intervalLength = hyper[:,1] - hyper[:,0]
	bisectIndex = np.argmax(intervalLength)
	#print ("bisectIndex", bisectIndex)
	lHyper = np.copy(hyper)
	rHyper = np.copy(hyper)
	midVal = (hyper[bisectIndex][0] + hyper[bisectIndex][1])/2.0
	lHyper[bisectIndex][1] = midVal
	rHyper[bisectIndex][0] = midVal
	return [lHyper, rHyper]

# Bisect guided by the Newton's method
# Try to find a Newton's solution in hyper.
# If it exists bisect so that one half contains
# the Newton's solution at a significant distance from the border
# and the other half does not. If it does not perform bisectMax
def bisectNewton(hyper, model):
	# Figure out the dimension to bisect over and the value at which
	# to visect over
	bisectIndex, cutoffVal = findBisectingIndexProportion(model, hyper)
	if bisectIndex is None:
		lHyp, rHyp = bisectMax(hyper)
	else:
		lHyp, rHyp = bisectAtIndex(hyper, [bisectIndex, cutoffVal])
	return [lHyp, rHyp]

# options = [bisectIndex, cutoffVal]
# Bisect the hyper-rectangle at the dimension indicated by
# bisectIndex and at the curoffVal
def bisectAtIndex(hyper, options):
	bisectIndex, cutoffVal = options
	lHyper = np.copy(hyper)
	rHyper = np.copy(hyper)
	lHyper[bisectIndex][1] = cutoffVal
	rHyper[bisectIndex][0] = cutoffVal
	return [lHyper, rHyper]

# Figure out the dimension to bisect over and the value at which
# to bisect over
def findBisectingIndexProportion(model, hyper):
	bisectIndex, cutoffVal = None, None
	# Find Newton's solution with the mid point as the initial solution
	hyperDist = hyper[:,1] - hyper[:,0]
	trialSoln = hyper[:,0] + 0.5*hyperDist
	finalSoln = intervalUtils.newton(model, trialSoln)
	if finalSoln[0]:
		# If the Newton's solution exists find the dimension for which
		# Newton's solution is farthest from the border.
		# Choose the cutoff val to be 0.3 times the distance between the Newton's solution
		# and the farthest border
		if np.all(finalSoln[1] >= hyper[:,0]) and np.all(finalSoln[1] <= hyper[:,1]):
			distFromLow = finalSoln[1] - hyper[:,0]
			distFromHigh = hyper[:,1] - finalSoln[1]
			maxDistIndexFromLo = np.argmax(distFromLow)
			maxDistIndexFromHi = np.argmax(distFromHigh)
			if distFromLow[maxDistIndexFromLo] > distFromHigh[maxDistIndexFromHi]:
				bisectIndex = maxDistIndexFromLo
				cutoffVal = finalSoln[1][bisectIndex] - 0.3*(finalSoln[1][bisectIndex] - hyper[bisectIndex][0])
			else:
				bisectIndex = maxDistIndexFromHi		
				cutoffVal = finalSoln[1][bisectIndex] + 0.3*(hyper[bisectIndex][1] - finalSoln[1][bisectIndex])

	return bisectIndex, cutoffVal			

# solver's main loop that uses LP
# @param uniqueHypers is a list of hyperrectangle containing unique solutions
#	found by solverLoop
# @param model indicates the problem we are trying to solve rambus/schmitt/metitarski
# @param statVars holds statistical information about the operations performed by the solver.
#	For example, number of bisections, number of Lp's performed
# @param volRedThreshold is indicates the stopping criterion for the loop of
#	Krawczyk and LP (implemented by the function ifFeasibleHyper) is applied
# @param bisectFun is a function that takes in a hyperrectangle and employes some mechanism
#	to bisect it	
# @param numSolutions indicates the number of solutions wanted by the user
# @param kAlpha is the threshold which indicates the stopping criterion for the Krawczyk loop
# @param epsilonInflation indicates the proportion of hyper-rectangle distance by which the 
# 	hyper-rectangle needs to be inflated before the Krawczyk operator is applied
# @param hyperRectangle the initial hyperrectangle over which the search for solutions
#	is done by solverLoop. If this argument is None then the hyperrectangle defined
#	by the bounds of the model is used
def solverLoop(uniqueHypers, model, statVars=None, volRedThreshold=1.0, bisectFun=bisectNewton, numSolutions="all", kAlpha=1.0, epsilonInflation=0.01, hyperRectangle = None):
	if not(hasattr(model, 'linearConstraints')):
		raise Exception("model has no instance of linearConstraints. Define a method called linearConstraints in the model class to use linear programming feature. Or use the solver without the linear programming feature\n")
	if statVars is None:
		statVars = {}
		statVars.update({'numBisection':0, 'numLp':0, 'numK':0, 'numSingleKill':0, 'numDoubleKill':0,
					'totalKTime':0, 'totalLPTime':0, 'avgKTime':0, 'avgLPTime':0, 'stringHyperList':[],
					'numLpCalls':0, 'numSuccessLpCalls':0, 'numUnsuccessLpCalls':0})
	lenV = len(model.bounds)
	
	if hyperRectangle is None:
		hyperRectangle = np.zeros((lenV,2))

		for i in range(lenV):
			hyperRectangle[i,0] = model.bounds[i][0]
			hyperRectangle[i,1] = model.bounds[i][1]

	
	statVars['stringHyperList'].append(("i", intervalUtils.volume(hyperRectangle)))
	
	start = time.time()
	feas = intervalUtils.checkExistenceOfSolution(model,hyperRectangle.transpose(),kAlpha,epsilonInflation)
	end = time.time()
	statVars['totalKTime'] += end - start
	statVars['numK'] += 1
	
	statVars['stringHyperList'].append(("g", intervalUtils.volume(feas[1])))


	#stack containing hyperrectangles about which any decision
	#has not been made - about whether they contain unique solution
	#or no solution
	stackList = []
	if feas[1] is not None:
		stackList.append(feas[1])

	while len(stackList) > 0:
		#pop the hyperrectangle
		#print ("len(stackList)", len(stackList))
		hyperPopped = stackList.pop(-1)
		#print ("solver loop hyperPopped")
		#intervalUtils.printHyper(hyperPopped)
		
		#if the popped hyperrectangle is contained in a hyperrectangle
		#that is already known to contain a unique solution, then do not
		#consider this hyperrectangle for the next steps
		hyperAlreadyConsidered = False
		for hyper in uniqueHypers:
			if np.greater_equal(hyperPopped[:,0], hyper[:,0]).all() and np.less_equal(hyperPopped[:,1], hyper[:,1]).all():
				hyperAlreadyConsidered = True
				break
		if hyperAlreadyConsidered:
			continue

		#Apply the Krawczyk + Lp loop
		feasibility = ifFeasibleHyper(hyperPopped, statVars, volRedThreshold, model, kAlpha, epsilonInflation)
		
		#print ("feasibility", feasibility)
		if feasibility[0]:
			#If the Krawczyk + Lp loop indicate uniqueness, then add the hyperrectangle
			#to our list
			if numSolutions == "all" or len(uniqueHypers) < numSolutions:
				addToSolutions(model, uniqueHypers, feasibility[1], kAlpha, epsilonInflation)

		elif feasibility[0] == False and feasibility[1] is not None:
			#If the Krawczyk + Lp loop cannot make a decision about
			#the hyperrectangle, the do the bisect and kill loop - keep
			#bisecting as long atleast one half either contains a unique
			#solution or no solution. Otherwise, add the two halves to
			#the stackList to be processed again.
			hypForBisection = feasibility[1]
			while hypForBisection is not None:
				#print ("hypForBisection")
				#intervalUtils.printHyper(hypForBisection)
				lHyp, rHyp = bisectFun(hypForBisection, model)
				statVars['numBisection'] += 1
				statVars['stringHyperList'].append(("b", [intervalUtils.volume(lHyp), intervalUtils.volume(rHyp)]))
				#print ("lHyp")
				#intervalUtils.printHyper(lHyp)
				start = time.time()
				lFeas = intervalUtils.checkExistenceOfSolution(model, lHyp.transpose(), kAlpha, epsilonInflation=epsilonInflation)
				end = time.time()
				#print ("lFeas", lFeas)
				statVars['totalKTime'] += end - start
				statVars['numK'] += 1
				statVars['stringHyperList'].append(("g", intervalUtils.volume(lFeas[1])))
				#print ("rHyp")
				#intervalUtils.printHyper(rHyp)
				start = time.time()
				rFeas = intervalUtils.checkExistenceOfSolution(model, rHyp.transpose(), kAlpha, epsilonInflation=epsilonInflation)
				end = time.time()
				#print ("rFeas", rFeas)
				statVars['totalKTime'] += end - start

				statVars['numK'] += 1
				statVars['stringHyperList'].append(("g", intervalUtils.volume(rFeas[1])))
				if lFeas[0] or rFeas[0] or (lFeas[0] == False and lFeas[1] is None) or (rFeas[0] == False and rFeas[1] is None):
					if lFeas[0] and rFeas[0]:
						statVars['numDoubleKill'] += 1
					elif lFeas[0] == False and lFeas[1] is None and rFeas[0] == False and rFeas[1] is None:
						statVars['numDoubleKill'] += 1
					else:
						statVars['numSingleKill'] += 1
					
					if lFeas[0]:
						if numSolutions == "all" or len(uniqueHypers) < numSolutions:
							addToSolutions(model, uniqueHypers, lFeas[1], kAlpha, epsilonInflation)
					if rFeas[0]:
						if numSolutions == "all" or len(uniqueHypers) < numSolutions:
							addToSolutions(model, uniqueHypers, rFeas[1], kAlpha, epsilonInflation)

					if lFeas[0] == False and lFeas[1] is not None:
						hypForBisection = lFeas[1]
					elif rFeas[0] == False and rFeas[1] is not None:
						hypForBisection = rFeas[1]
					else:
						hypForBisection = None

				
				else:
					stackList.append(lFeas[1])
					stackList.append(rFeas[1])
					hypForBisection = None


# solver's main loop that doesn't use LP
# @param uniqueHypers is a list of hyperrectangle containing unique solutions
#	found by solverLoop
# @param model indicates the problem we are trying to solve rambus/schmitt/metitarski
# @param statVars holds statistical information about the operations performed by the solver.
#	For example, number of bisections, number of Lp's performed
# @param bisectFun is a function that takes in a hyperrectangle and employes some mechanism
#	to bisect it	
# @param numSolutions indicates the number of solutions wanted by the user
# @param kAlpha is the threshold which indicates the stopping criterion for the Krawczyk loop
# @param epsilonInflation indicates the proportion of hyper-rectangle distance by which the 
# 	hyper-rectangle needs to be inflated before the Krawczyk operator is applied
# @param hyperRectangle the initial hyperrectangle over which the search for solutions
#	is done by solverLoop. If this argument is None then the hyperrectangle defined
#	by the bounds of the model is used
def solverLoopNoLp(uniqueHypers, model, statVars=None, bisectFun=bisectMax, numSolutions="all", kAlpha=1.0, epsilonInflation=0.001, hyperRectangle = None):
	if statVars is None:
		statVars = {}
		statVars.update({'numBisection':0, 'numLp':0, 'numK':0, 'numSingleKill':0, 'numDoubleKill':0,
					'totalKTime':0, 'totalLPTime':0, 'avgKTime':0, 'avgLPTime':0, 'stringHyperList':[],
					'numLpCalls':0, 'numSuccessLpCalls':0, 'numUnsuccessLpCalls':0})
	lenV = len(model.bounds)
	
	if hyperRectangle is None:
		hyperRectangle = np.zeros((lenV,2))

		for i in range(lenV):
			hyperRectangle[i,0] = model.bounds[i][0]
			hyperRectangle[i,1] = model.bounds[i][1]

	
	statVars['stringHyperList'].append(("i", intervalUtils.volume(hyperRectangle)))
	
	#stack containing hyperrectangles about which any decision
	#has not been made - about whether they contain unique solution
	#or no solution
	stackList = [hyperRectangle]

	while len(stackList) > 0:
		#pop the hyperrectangle
		#print ("len(stackList)", len(stackList))
		hyperPopped = stackList.pop(-1)
		#print ("solver loop hyperPopped")
		#intervalUtils.printHyper(hyperPopped)
		
		#if the popped hyperrectangle is contained in a hyperrectangle
		#that is already known to contain a unique solution, then do not
		#consider this hyperrectangle for the next steps
		hyperAlreadyConsidered = False
		for hyper in uniqueHypers:
			if np.greater_equal(hyperPopped[:,0], hyper[:,0]).all() and np.less_equal(hyperPopped[:,1], hyper[:,1]).all():
				hyperAlreadyConsidered = True
				break
		if hyperAlreadyConsidered:
			continue

		start = time.time()
		feasibility = intervalUtils.checkExistenceOfSolution(model, hyperPopped.transpose(), kAlpha, epsilonInflation=epsilonInflation)
		end = time.time()
		statVars['totalKTime'] += end - start
		statVars['numK'] += 1
		statVars['stringHyperList'].append(("g", intervalUtils.volume(feasibility[1])))
		
		#print ("feasibility", feasibility)
		if feasibility[0]:
			#If the Krawczyk loop indicate uniqueness, then add the hyperrectangle
			#to our list
			if numSolutions == "all" or len(uniqueHypers) < numSolutions:
				#print ("solution found")
				#print ("hyper")
				#intervalUtils.printHyper(hyperPopped)
				#print ("feas")
				#intervalUtils.printHyper(feasibility[1])
				addToSolutions(model, uniqueHypers, feasibility[1], kAlpha, epsilonInflation)

		elif feasibility[0] == False and feasibility[1] is not None:
			#If the Krawczyk loop cannot make a decision about
			#the hyperrectangle, bisect and add the two halves to
			#the stackList to be processed again.
			hypForBisection = feasibility[1]
			lHyp, rHyp = bisectFun(hypForBisection, model)
			statVars['numBisection'] += 1
			statVars['stringHyperList'].append(("b", [intervalUtils.volume(lHyp), intervalUtils.volume(rHyp)]))
			stackList.append(lHyp)
			stackList.append(rHyp)


# Apply Krawczyk and linear programming to refine the hyperrectangle
# @param hyperRectangle 
# @param statVars statVars holds statistical information about the operations performed by the solver.
#	For example, number of Lp's performed
# @param volRedTheshold indicates the stopping criterion for the loop of
#	Krawczyk and LP (implemented by this function). Basically keep repeating the loop
#	as long as the percentage of volume reduction of the hyperrectangle is atleast volRedThreshold. 
# @param model indicates the problem we are trying to solve
# @param kAlpha is the threshold which indicates the stopping criterion for the Krawczyk loop
# @param epsilonInflation indicates the proportion of hyper-rectangle distance by which the 
# 	hyper-rectangle needs to be inflated before the Krawczyk operator is applied
# @return (True, hyper) if hyperRectangle contains a unique
# 	solution and hyper maybe smaller than hyperRectangle containing the solution
# @return (False, None) if hyperRectangle contains no solution
# @return (False, hyper) if hyperRectangle may contain more
# 	than 1 solution and hyper maybe smaller than hyperRectangle containing the solutions
def ifFeasibleHyper(hyperRectangle, statVars, volRedThreshold, model, kAlpha,epsilonInflation):
	lenV = hyperRectangle.shape[0]
	iterNum = 0
	while True:
		newHyperRectangle = np.copy(hyperRectangle)

		# Apply linear programming step
		start = time.time()
		feasible, newHyperRectangle, numTotalLp, numSuccessLp, numUnsuccessLp = model.linearConstraints(newHyperRectangle)
		end = time.time()
		statVars['totalLPTime'] += end - start
		statVars['numLpCalls'] += numTotalLp
		statVars['numSuccessLpCalls'] += numSuccessLp
		statVars['numUnsuccessLpCalls'] += numUnsuccessLp
		statVars['numLp'] += 1
		if feasible:
			vol = intervalUtils.volume(newHyperRectangle)
		else:
			vol = None
		statVars['stringHyperList'].append(("l", vol))
		#print ("newHyperRectangle", newHyperRectangle)
		#intervalUtils.printHyper(newHyperRectangle)
		if feasible == False:
			return (False, None)

		for i in range(lenV):
			if newHyperRectangle[i,0] < hyperRectangle[i,0]:
				newHyperRectangle[i,0] = hyperRectangle[i,0]
			if newHyperRectangle[i,1] > hyperRectangle[i,1]:
				newHyperRectangle[i,1] = hyperRectangle[i,1]


		start = time.time()
		
		#Apply Krawczyk
		kResult = intervalUtils.checkExistenceOfSolution(model, newHyperRectangle.transpose(), kAlpha, epsilonInflation=epsilonInflation)
		end = time.time()
		statVars['totalKTime'] += (end - start)
		statVars['numK'] += 1
		statVars['stringHyperList'].append(("g", intervalUtils.volume(kResult[1])))
		
		#Unique solution or no solution
		if kResult[0] or kResult[1] is None:
			#print ("uniqueHyper", hyperRectangle)
			return kResult

		newHyperRectangle = kResult[1]			 

		hyperVol = intervalUtils.volume(hyperRectangle)

		newHyperVol = intervalUtils.volume(newHyperRectangle)

		propReduc = (hyperVol - newHyperVol)/hyperVol
		
		# If the proportion of volume reduction is not atleast
		# volRedThreshold then return
		if math.isnan(propReduc) or propReduc <= volRedThreshold:
			return (False, newHyperRectangle)
		hyperRectangle = newHyperRectangle
		iterNum+=1
	


# A function that adds a new solution to the list of existing hyperrectangles
# containing unique solutions if the new hyperrectangle does not contain
# the same solution as the solution in any of the existing hyperrectangles
# @param model indicates the problem we are trying to solve rambus/schmitt/metitarski
# @param allHypers list of hyperrectangles containing unique solutions
# @param solHyper new hyperrectangle containing unique solution
# @param kAlpha is the threshold which indicates the stopping criterion for the Krawczyk loop
# @param epsilonInflation indicates the proportion of hyper-rectangle distance by which the 
# 	hyper-rectangle needs to be inflated before the Krawczyk operator is applied
def addToSolutions(model, allHypers, solHyper, kAlpha,epsilonInflation):
	epsilon = 1e-12
	lenV = len(model.bounds)
	foundOverlap = False
	exampleVolt = (solHyper[:,0] + solHyper[:,1])/2.0
	soln = intervalUtils.newton(model,exampleVolt)
	
	if not(soln[0]):
		raise Exception("prototype.py addToSolutions: Something went wrong. Should contain a unique solution")
	
	for hi in range(len(allHypers)):
		oldHyper = allHypers[hi]

		'''for i in range(lenV):
			if interval_intersect(solHyper[i], oldHyper[i]) is None:
				print ("not intersecting", solHyper[i], oldHyper[i])'''
		#Check if solHyper overlaps with oldHyper
		if all(interval_intersect(solHyper[i], oldHyper[i]) is not None for i in range(lenV)):
			intersectHyper = np.zeros((lenV,2))
			for ui in range(lenV):
				intersectHyper[ui,:] = interval_intersect(solHyper[ui], oldHyper[ui])

			#print ("intersectHyper", intersectHyper)
			if np.all(soln[1] >= intersectHyper[:,0]) and np.all(soln[1] <= intersectHyper[:,1]):
				hyperAroundNewton = np.zeros((lenV, 2))
				for si in range(lenV):
					minDiff = min(abs(intersectHyper[si,1] - soln[1][si]), abs(soln[1][si] - intersectHyper[si,0]))
					hyperAroundNewton[si,0] = soln[1][si] - minDiff
					hyperAroundNewton[si,1] = soln[1][si] + minDiff
				feasibility = intervalUtils.checkExistenceOfSolution(model, hyperAroundNewton.transpose(), alpha = kAlpha, epsilonInflation=epsilonInflation)
				if feasibility[0]:
					foundOverlap = True
					break

	if not(foundOverlap):
		allHypers.append(solHyper)
		return True
	else:
		return False




# Find the dc equilibrium points for the schmitt trigger for a specific
# input voltage
# @param modelType indicates the type of transistor model used for the schmitt
#	trigger. If modelType == "lcMosfet", use the long channel mosfet model.
#	If modelType == "scMosfet", use the short channel mosfet model
# @param inputVoltage the value of the specific input voltage for which 
# 	the dc equilibrium points are found
# @param statVars dictionary to hold statistics like number of bisections, number of Krawczyk calls
# @param kAlpha is the threshold which indicates the stopping criterion for the Krawczyk loop
# @param epsilonInflation indicates the proportion of hyper-rectangle distance by which the 
# 	hyper-rectangle needs to be inflated before the Krawczyk operator is applied
# @param bisectType indicates the type of bisection used in the solver
# 	A bisectType of "bisectMax" indicates the solver that it should use bisectMax
# 	method. A bisecType of "bisectNewton" indicates the solver that it should use
#	bisectNewton method
# @param numSolutions, number of dc equilibrium points we are looking for
# @param useLp flag to decide whether to use linear programming in our method or not
# @return a list of hyperrectangles containing unique dc equilibrium points
def schmittTrigger(modelType, inputVoltage, statVars, kAlpha = 1.0, epsilonInflation=0.001, bisectType="bisectMax", numSolutions = "all", useLp = False):
	statVars.update({'numBisection':0, 'numLp':0, 'numK':0, 'numSingleKill':0, 'numDoubleKill':0,
					'totalKTime':0, 'totalLPTime':0, 'avgKTime':0, 'avgLPTime':0, 'stringHyperList':[],
					'numLpCalls':0, 'numSuccessLpCalls':0, 'numUnsuccessLpCalls':0})

	#load the schmitt trigger model
	if modelType == "lcMosfet":
		#modelParam = [Vtp, Vtn, Vdd, Kn, Kp, Sn]
		modelParam = [-0.4, 0.4, 1.8, 270*1e-6, -90*1e-6, 8/3.0]
		model = SchmittMosfet(modelType = modelType, modelParam = modelParam, inputVoltage = inputVoltage)
	elif modelType == "scMosfet":
		modelParam = [1.0] #Vdd
		model = SchmittMosfet(modelType = modelType, modelParam = modelParam, inputVoltage = inputVoltage)

	startExp = time.time()

	allHypers = []
	#print ("model val", model.f(np.array([1.7, 0.1, 0.1])))

	if bisectType == "bisectMax":
		bisectFun = bisectMax
	if bisectType == "bisectNewton":
		bisectFun = bisectNewton
	if useLp:
		volRedThreshold = 1.0
		solverLoop(uniqueHypers=allHypers, model=model, statVars=statVars, volRedThreshold=volRedThreshold, bisectFun=bisectFun, numSolutions=numSolutions, kAlpha=kAlpha, epsilonInflation=epsilonInflation)
	else:
		solverLoopNoLp(uniqueHypers=allHypers, model=model, statVars=statVars, bisectFun=bisectFun, numSolutions=numSolutions, kAlpha=kAlpha, epsilonInflation=epsilonInflation)

	#print ("allHypers")
	#print (allHypers)
	#print ("numSolutions", len(allHypers))

	
	#dcUtils.printSol(allHypers, model)
	endExp = time.time()
	#print ("TOTAL TIME ", endExp - startExp)
	if statVars['numLp'] != 0:
		statVars['avgLPTime'] = (statVars['totalLPTime']*1.0)/statVars['numLp']
	if statVars['numK'] != 0:
		statVars['avgKTime'] = (statVars['totalKTime']*1.0)/statVars['numK']

	#print ("numBisection", statVars['numBisection'], "numLp", statVars['numLp'], "numK", statVars['numK'],
	#	"numSingleKill", statVars['numSingleKill'], "numDoubleKill", statVars['numDoubleKill'])
	#print ("totalKTime", statVars['totalKTime'], "totalLPTime", statVars['totalLPTime'], "avgKTime", 
	#	statVars['avgKTime'], "avgLPTime", statVars['avgLPTime'])
	#print ("numLpCalls", statVars['numLpCalls'], "numSuccessLpCalls", statVars['numSuccessLpCalls'], "numUnsuccessLpCalls", statVars['numUnsuccessLpCalls'])
	return allHypers

# Find the dc equilibrium points for an inverter for a specific
# input voltage
# @param modelType indicates the type of transistor model used for the inverter. 
# If modelType == "lcMosfet", use the long channel mosfet model.
#	If modelType == "scMosfet", use the short channel mosfet model
# @param inputVoltage the value of the specific input voltage for which 
# 	the dc equilibrium points are found
# @param statVars dictionary to hold statistics like number of bisections, number of Krawczyk calls
# @param kAlpha is the threshold which indicates the stopping criterion for the Krawczyk loop
# @param epsilonInflation indicates the proportion of hyper-rectangle distance by which the 
# 	hyper-rectangle needs to be inflated before the Krawczyk operator is applied
# @param bisectType indicates the type of bisection used in the solver
# 	A bisectType of "bisectMax" indicates the solver that it should use bisectMax
# 	method. A bisecType of "bisectNewton" indicates the solver that it should use
#	bisectNewton method
# @param numSolutions, number of dc equilibrium points we are looking for
# @param useLp flag to decide whether to use linear programming in our method or not
# @return a list of hyperrectangles containing unique dc equilibrium points
def inverter(modelType, inputVoltage, statVars, kAlpha=1.0, epsilonInflation=0.001, bisectType="bisectMax", numSolutions="all" , useLp=False):
	statVars.update({'numBisection':0, 'numLp':0, 'numK':0, 'numSingleKill':0, 'numDoubleKill':0,
					'totalKTime':0, 'totalLPTime':0, 'avgKTime':0, 'avgLPTime':0, 'stringHyperList':[],
					'numLpCalls':0, 'numSuccessLpCalls':0, 'numUnsuccessLpCalls':0})
	

	#load the inverter model
	if modelType == "tanh":
		modelParam = [-5.0, 0.0] # y = tanh(modelParam[0]*x + modelParam[1])
		model = InverterTanh(modelParam = modelParam, inputVoltage = inputVoltage)
	if modelType == "lcMosfet":
		#modelParam = [Vtp, Vtn, Vdd, Kn, Kp, Sn]
		modelParam = [-0.4, 0.4, 1.8, 270*1e-6, -90*1e-6, 8/3.0]
		model = InverterMosfet(modelType = modelType, modelParam = modelParam, inputVoltage = inputVoltage)
	if modelType == "scMosfet":
		modelParam = [1.0] #Vdd
		model = InverterMosfet(modelType = modelType, modelParam = modelParam, inputVoltage = inputVoltage)

	startExp = time.time()
	
	allHypers = []
	if bisectType == "bisectMax":
		bisectFun = bisectMax
	if bisectType == "bisectNewton":
		bisectFun = bisectNewton
	if useLp:
		volRedThreshold = 1.0
		solverLoop(uniqueHypers=allHypers, model=model, statVars=statVars, volRedThreshold=volRedThreshold, bisectFun=bisectFun, numSolutions=numSolutions, kAlpha=kAlpha, epsilonInflation=epsilonInflation)
	else:
		solverLoopNoLp(uniqueHypers=allHypers, model=model, statVars=statVars, bisectFun=bisectFun, numSolutions=numSolutions, kAlpha=kAlpha, epsilonInflation=epsilonInflation)
	
	#print ("allHypers")
	#print (allHypers)
	#print ("numSolutions", len(allHypers))
	
	endExp = time.time()
	#print ("TOTAL TIME ", endExp - startExp)

	if statVars['numLp'] != 0:
		statVars['avgLPTime'] = (statVars['totalLPTime']*1.0)/statVars['numLp']
	if statVars['numK'] != 0:
		statVars['avgKTime'] = (statVars['totalKTime']*1.0)/statVars['numK']
	
	#print ("numBisection", statVars['numBisection'], "numLp", statVars['numLp'], "numK", statVars['numK'],
	#	"numSingleKill", statVars['numSingleKill'], "numDoubleKill", statVars['numDoubleKill'])
	#print ("totalKTime", statVars['totalKTime'], "totalLPTime", statVars['totalLPTime'], "avgKTime", 
	#	statVars['avgKTime'], "avgLPTime", statVars['avgLPTime'])
	#print ("numLpCalls", statVars['numLpCalls'], "numSuccessLpCalls", statVars['numSuccessLpCalls'], "numUnsuccessLpCalls", statVars['numUnsuccessLpCalls'])
	return allHypers


# Find the dc equilibrium points for an inverter loop 
# @param modelType indicates the type of transistor model used for the inverter. 
# If modelType == "lcMosfet", use the long channel mosfet model.
#	If modelType == "scMosfet", use the short channel mosfet model
# @param numInverters the number of inverters in the inverter loop
# @param statVars dictionary to hold statistics like number of bisections, number of Krawczyk calls
# @param kAlpha is the threshold which indicates the stopping criterion for the Krawczyk loop
# @param epsilonInflation indicates the proportion of hyper-rectangle distance by which the 
# 	hyper-rectangle needs to be inflated before the Krawczyk operator is applied
# @param bisectType indicates the type of bisection used in the solver
# 	A bisectType of "bisectMax" indicates the solver that it should use bisectMax
# 	method. A bisecType of "bisectNewton" indicates the solver that it should use
#	bisectNewton method
# @param numSolutions, number of dc equilibrium points we are looking for
# @param useLp flag to decide whether to use linear programming in our method or not
# @return a list of hyperrectangles containing unique dc equilibrium points
def inverterLoop(modelType, numInverters, statVars, kAlpha=1.0, epsilonInflation=0.001, bisectType="bisectMax", numSolutions="all" , useLp=False):
	statVars.update({'numBisection':0, 'numLp':0, 'numK':0, 'numSingleKill':0, 'numDoubleKill':0,
					'totalKTime':0, 'totalLPTime':0, 'avgKTime':0, 'avgLPTime':0, 'stringHyperList':[],
					'numLpCalls':0, 'numSuccessLpCalls':0, 'numUnsuccessLpCalls':0})
	

	#load the inverter model
	if modelType == "tanh":
		modelParam = [-5.0, 0.0] # y = tanh(modelParam[0]*x + modelParam[1])
		model = InverterLoopTanh(modelParam = modelParam, numInverters = numInverters)
	if modelType == "lcMosfet":
		#modelParam = [Vtp, Vtn, Vdd, Kn, Kp, Sn]
		modelParam = [-0.4, 0.4, 1.8, 270*1e-6, -90*1e-6, 8/3.0]
		model = InverterLoopMosfet(modelType = modelType, modelParam = modelParam, numInverters = numInverters)
	if modelType == "scMosfet":
		modelParam = [1.0] #Vdd
		model = InverterLoopMosfet(modelType = modelType, modelParam = modelParam, numInverters = numInverters)

	startExp = time.time()
	
	allHypers = []
	if bisectType == "bisectMax":
		bisectFun = bisectMax
	if bisectType == "bisectNewton":
		bisectFun = bisectNewton
	if useLp:
		volRedThreshold = 1.0
		solverLoop(uniqueHypers=allHypers, model=model, statVars=statVars, volRedThreshold=volRedThreshold, bisectFun=bisectFun, numSolutions=numSolutions, kAlpha=kAlpha, epsilonInflation=epsilonInflation)
	else:
		solverLoopNoLp(uniqueHypers=allHypers, model=model, statVars=statVars, bisectFun=bisectFun, numSolutions=numSolutions, kAlpha=kAlpha, epsilonInflation=epsilonInflation)
	
	#print ("allHypers")
	#print (allHypers)
	#print ("numSolutions", len(allHypers))
	
	endExp = time.time()
	#print ("TOTAL TIME ", endExp - startExp)

	if statVars['numLp'] != 0:
		statVars['avgLPTime'] = (statVars['totalLPTime']*1.0)/statVars['numLp']
	if statVars['numK'] != 0:
		statVars['avgKTime'] = (statVars['totalKTime']*1.0)/statVars['numK']
	
	#print ("numBisection", statVars['numBisection'], "numLp", statVars['numLp'], "numK", statVars['numK'],
	#	"numSingleKill", statVars['numSingleKill'], "numDoubleKill", statVars['numDoubleKill'])
	#print ("totalKTime", statVars['totalKTime'], "totalLPTime", statVars['totalLPTime'], "avgKTime", 
	#	statVars['avgKTime'], "avgLPTime", statVars['avgLPTime'])
	#print ("numLpCalls", statVars['numLpCalls'], "numSuccessLpCalls", statVars['numSuccessLpCalls'], "numUnsuccessLpCalls", statVars['numUnsuccessLpCalls'])
	return allHypers



# Find the dc equilibrium points for a rambus ring oscillator
# @param modelType indicates the type of inverter used in the rambus oscillator
# 	If modelType == "tanh", use the tanh inverter model.
#	If modelType == "lcMosfet", use transistor with two long channel mosfet models.
#	If modelType == "scMosfet", use transistor with two short channel mosfet models.
# @param numStages the number of stages in the rambus ring oscillator
# @param g_cc strength of the cross coupled inverter (as compared to that of the forward)
# @param statVars dictionary to hold statistics like number of bisections, number of Krawczyk calls
# @param kAlpha is the threshold which indicates the stopping criterion for the Krawczyk loop
# @param epsilonInflation indicates the proportion of hyper-rectangle distance by which the 
# 	hyper-rectangle needs to be inflated before the Krawczyk operator is applied
# @param bisectType indicates the type of bisection used in the solver
# 	A bisectType of "bisectMax" indicates the solver that it should use bisectMax
# 	method. A bisecType of "bisectNewton" indicates the solver that it should use
#	bisectNewton method
# @param numSolutions, number of dc equilibrium points we are looking for
# @param useLp flag to decide whether to use linear programming in our method or not
# @return a list of hyperrectangles containing unique dc equilibrium points
def rambusOscillator(modelType, numStages, g_cc, statVars, kAlpha=1.0, epsilonInflation=0.01, bisectType="bisectMax", numSolutions="all", useLp=False):
	statVars.update({'numBisection':0, 'numLp':0, 'numK':0, 'numSingleKill':0, 'numDoubleKill':0,
					'totalKTime':0, 'totalLPTime':0, 'avgKTime':0, 'avgLPTime':0, 'stringHyperList':[],
					'numLpCalls':0, 'numSuccessLpCalls':0, 'numUnsuccessLpCalls':0})
	
	if modelType == "tanh":
		modelParam = [-5.0, 0.0] # y = tanh(modelParam[0]*x + modelParam[1])
		model = RambusTanh(modelParam = modelParam, g_cc = g_cc, g_fwd = 1.0, numStages=numStages)
	elif modelType == "lcMosfet":
		#modelParam = [Vtp, Vtn, Vdd, Kn, Kp, Sn]
		#modelParam = [-0.25, 0.25, 1.0, 1.0, -0.5, 1.0]
		#modelParam = [-0.4, 0.4, 1.8, 1.5, -0.5, 8/3.0]
		modelParam = [-0.4, 0.4, 1.8, 270*1e-6, -90*1e-6, 8/3.0]
		model = RambusMosfet(modelType = modelType, modelParam = modelParam, g_cc = g_cc, g_fwd = 1.0, numStages = numStages)	
	elif modelType == "scMosfet":
		modelParam = [1.0] #Vdd
		model = RambusMosfet(modelType = modelType, modelParam = modelParam, g_cc = g_cc, g_fwd = 1.0, numStages = numStages)

	startExp = time.time()
	
	allHypers = []

	'''hyper1 = np.array([[1.35, 1.8],
						[0.9, 1.8],
						[0.0, 0.9],
						[0.0, 0.9]])'''
	hyper1 = None
	#feasibility = ifFeasibleHyper(hyperRectangle=hyper1, statVars=statVars, volRedThreshold=1.0, model=model, kAlpha=kAlpha,epsilonInflation=epsilonInflation)
	#print ("feasibility", feasibility)
	'''hyper1 = np.array([[1.7999998100803811, 1.8],
						[0.0, 4.4122565118110095e-12],
						[0.5204758706791751, 0.5539991578829824],
						[1.795431201546126, 1.8],
						[0.0, 8.105592737403754e-08],
						[1.7999999999986622, 1.8],
						[1.0453925769488364, 1.0685168999355248],
						[0.00396157418470353, 0.009124314562001691]])
	fVal = model.f(hyper1)
	print ("fVal")
	intervalUtils.printHyper(fVal)
	feas = intervalUtils.checkExistenceOfSolution(model,hyper1.transpose(), alpha = kAlpha, epsilonInflation=epsilonInflation)
	print ("feasibility")
	print (feas)'''


	'''hyper1 = np.array([[1.6278178919881217, 1.6278178942899837],
						[1.55407062766072, 1.5540706309561325],
						[6.523159433547119e-17, 1.2288759830968326e-17],
						[0.16870818405562973, 0.16870818565396228],
						[0.1076139571020803, 0.10761395883468822],
						[0.12766813879925398, 0.12766814077766875],
						[1.7999999999999992, 1.800000000000001],
						[1.452495225275838, 1.4524952267582616]])

	hyper2 = np.array([[1.6278178920290314, 1.6278178942563253],
						[1.5540706276747582, 1.5540706309385526],
						[5.106161194109459e-17, 1.456900429145809e-17],
						[0.16870818412176625, 0.16870818560404693],
						[0.10761395711218251, 0.10761395882091715],
						[0.1276681388034702, 0.1276681407746626],
						[1.7999999999999992, 1.800000000000001],
						[1.4524952253462955, 1.452495226650087]])
	allHypers = [hyper1]
	addToSolutions(model, allHypers, hyper2, kAlpha,epsilonInflation)
	print ("allHypers")
	print (allHypers)'''

	if bisectType == "bisectMax":
		bisectFun = bisectMax
	if bisectType == "bisectNewton":
		bisectFun = bisectNewton
	if useLp:
		volRedThreshold = 1.0
		solverLoop(uniqueHypers=allHypers, model=model, statVars=statVars, volRedThreshold=volRedThreshold, bisectFun=bisectFun, numSolutions=numSolutions, kAlpha=kAlpha, epsilonInflation=epsilonInflation, hyperRectangle=hyper1)
	else:
		solverLoopNoLp(uniqueHypers=allHypers, model=model, statVars=statVars, bisectFun=bisectFun, numSolutions=numSolutions, kAlpha=kAlpha, epsilonInflation=epsilonInflation, hyperRectangle = hyper1)
	
	#print ("allHypers")
	#print (allHypers)
	#dcUtils.printSol(allHypers, model)
	#print ("numSolutions", len(allHypers))

	endExp = time.time()
	#print ("TOTAL TIME ", endExp - startExp)

	if statVars['numLp'] != 0:
		statVars['avgLPTime'] = (statVars['totalLPTime']*1.0)/statVars['numLp']
	if statVars['numK'] != 0:
		statVars['avgKTime'] = (statVars['totalKTime']*1.0)/statVars['numK']
	
	#print ("numBisection", statVars['numBisection'], "numLp", statVars['numLp'], "numK", statVars['numK'],
	#	"numSingleKill", statVars['numSingleKill'], "numDoubleKill", statVars['numDoubleKill'])
	#print ("totalKTime", statVars['totalKTime'], "totalLPTime", statVars['totalLPTime'], "avgKTime", 
	#	statVars['avgKTime'], "avgLPTime", statVars['avgLPTime'])
	#print ("numLpCalls", statVars['numLpCalls'], "numSuccessLpCalls", statVars['numSuccessLpCalls'], "numUnsuccessLpCalls", statVars['numUnsuccessLpCalls'])
	return allHypers


if __name__ == "__main__":
	statVars = {}
	start = time.time()
	#allHypers = schmittTrigger(modelType="lcMosfet", inputVoltage = 1.8, statVars=statVars, numSolutions = "all", useLp = False)
	#allHypers = inverter(modelType="lcMosfet", inputVoltage=0.9, statVars=statVars, numSolutions="all")
	#allHypers = inverterLoop(modelType="tanh", numInverters=4, statVars=statVars, numSolutions="all", useLp = False)
	allHypers = rambusOscillator(modelType="scMosfet", numStages=4, g_cc=4.0, statVars=statVars, kAlpha = 1.0, epsilonInflation=0.001, numSolutions="all", bisectType="bisectMax", useLp = False)
	end = time.time()
	#print ("allHypers")
	#for hyper in allHypers:
	#	print ("hyper")
	#	intervalUtils.printHyper(hyper)
	print ("numSolutions", len(allHypers))
	print ("time taken", end - start)
