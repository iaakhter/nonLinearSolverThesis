# @author: Itrat Ahmed Akhter
# Find dc equilibrium points for rambus oscillator with the inverter
# modeled as a tanh function with Z3

from z3 import *
import numpy as np
import time



# Return list of constraints approximating tanh with piecewise polynomial function
# a <= 0 -- gain for tanh
# VinVar is a Z3 variable representing input voltage
# VoutVar is a Z3 variable representing output voltage
def tanhCurrent(a, VinVar, VoutVar):
	constraints = []
	constraints.append(Or(And(VinVar > -2.0/a, VoutVar == -1-(a*VinVar)**(-5)),
							And(VinVar <= -2.0/a, VinVar >= 2.0/a, VoutVar == -(13/256.0)*(a*VinVar)**3 + (11/16.0)*(a*VinVar)),
							And(VinVar < 2.0/a, VoutVar == 1-(a*VinVar)**(-5))))
	return constraints

# Try and find dc equilibrium points for rambus oscillator with the inverter
# being modeled by tanh function
# numStages indicates the number of stages in the rambus oscillator
# g_cc is the strength of the cross coupled inverter as compared to that of forward (g_fwd = 1.0)
# numSolutions indicates the number of solutions we want Z3 to find
# a is the gain of the inverter. Should be negative
def rambusOscillatorTanh(numStages, g_cc = 0.5, numSolutions = "all", a = -5.0):
	start = time.time()
	g_fwd = 1.0
	lenV = numStages*2

	set_option(rational_to_decimal=True)

	vs = RealVector("v", lenV)
	vfwds = RealVector("vfwd", lenV)
	vccs = RealVector("vcc", lenV)

	allConstraints = []
		
	# Store rambus oscillator constraints
	for i in range(lenV):
		allConstraints.append(vs[i] >= -1)
		allConstraints.append(vs[i] <= 1)
		fwdInd = (i-1)%lenV
		ccInd = (i+lenV//2)%lenV
		allConstraints += tanhCurrent(a, vs[fwdInd], vfwds[i])
		allConstraints += tanhCurrent(a, vs[ccInd], vccs[i])
		allConstraints.append(g_fwd*vfwds[i] + (-g_fwd-g_cc)*vs[i] + g_cc*vccs[i] == 0)

	s = Solver()
	allSolutions = []
	while True:
		if numSolutions != "all" and len(allSolutions) == numSolutions:
			break

		# Store constraints pruning search space so that
		# old solutions are not considered
		excludingConstraints = []
		for solution in allSolutions:
			singleExcludingConstraints = []
			for i in range(lenV):
				singleExcludingConstraints.append(vs[i] != solution[i])
			excludingConstraints.append(singleExcludingConstraints)
		
		#print ("numConstraints", len(allConstraints))
		# Add all the rambus oscillator constraints
		f_sat = And(*allConstraints)
		# Add constraints so that old solutions are not considered
		if len(excludingConstraints) > 0:
			for constraints in excludingConstraints:
				f_sat = And(f_sat, Or(*constraints))
		
		# Add the constraints to Z3 with a push and pop operation
		#print ("f_sat")
		#print (f_sat)
		s.push()
		s.add(f_sat)
		#print ("s")
		#print (s)
		result = s.check()
		#print (result)
		if result != sat:
			break
		
		m = s.model()
		sol = [None]*lenV
		for d in m.decls():
			dName = str(d.name())
			firstLetter = dName[0]
			if (dName[0] == "v" and dName[1] == "_"):
				index = int(dName[len(dName) - 1])
				
				sol[index] = m[d]


		print ("sol", sol)
		s.pop()

		allSolutions.append(sol)
		
		print ("num solutions found", len(allSolutions))
		
		
	'''print ("all solutions")
	for solution in allSolutions:
		print (solution)
	
	print ("num solutions", len(allSolutions))'''
	end = time.time()
	print ("time taken", end - start)
	return allSolutions


# Try and find dc equilibrium points for inverter modeled by a tanh function
# for a specific input voltage
# a indicates the gain by the tanh function
# numSolutions indicates the number of solutions we want Z3 to find
def inverterTanh(inputVoltage, a=-5.0, numSolutions = "all"):
	start = time.time()

	set_option(rational_to_decimal=True)

	outputVolt = Real("outputVolt")
	vRes = Real("vRes")

	allConstraints = []	
	allConstraints.append(outputVolt >= -1.0)
	allConstraints.append(outputVolt <= 1.0)
	tanhVal = np.tanh(a*inputVoltage)
	allConstraints.append(vRes == tanhVal)
	allConstraints.append(vRes - outputVolt == 0.0)
	
	s = Solver()
	allSolutions = []
	while True:
		if numSolutions != "all" and len(allSolutions) == numSolutions:
			break
		
		# Store constraints pruning search space so that
		# old solutions are not considered
		excludingConstraints = []
		for solution in allSolutions:
			singleExcludingConstraints = []
			singleExcludingConstraints.append(outputVolt != solution)
			excludingConstraints.append(singleExcludingConstraints)
		
		#print ("allConstraints")
		#print (allConstraints)
		#print ("numConstraints", len(allConstraints))
		f_sat = And(*allConstraints)
		if len(excludingConstraints) > 0:
			for constraints in excludingConstraints:
				f_sat = And(f_sat, Or(*constraints))
		
		# Add the constraints to Z3 with a push and pop operation
		#print ("f_sat")
		#print (f_sat)
		s.push()
		s.add(f_sat)
		#print ("s")
		#print (s)
		result = s.check()
		#print (result)
		if result != sat:
			break
		
		m = s.model()
		#print ("m", m)

		sol = None
		for d in m.decls():
			dName = str(d.name())
			if dName == "outputVolt":
				sol = m[d]
		
		print ("sol", sol)
		s.pop()
		
		allSolutions.append(sol)

		print ("num solutions found", len(allSolutions))

	
	'''print ("all solutions")
	for solution in allSolutions:
		print (solution)'''

	end = time.time()
	print ("time taken", end - start)
	return allSolutions


def inverterLoopTanh(numInverters, numSolutions = "all", a= -5.0):
	start = time.time()
	
	set_option(rational_to_decimal=True)

	vs = RealVector("v", numInverters)
	vRes = RealVector("vRes", numInverters)

	allConstraints = []
	s = Solver()
		
	# Store rambus oscillator constraints
	for i in range(numInverters):
		allConstraints.append(vs[i] >= -1)
		allConstraints.append(vs[i] <= 1)
		inputInd = i
		outputInd = (i+1)%numInverters
		allConstraints += tanhCurrent(a, vs[inputInd], vRes[i])
		allConstraints.append(vRes[i] - vs[outputInd] == 0.0)

	allSolutions = []
	while True:
		if numSolutions != "all" and len(allSolutions) == numSolutions:
			break

		# Store constraints pruning search space so that
		# old hyperrectangles are not considered
		excludingConstraints = []
		for solution in allSolutions:
			singleExcludingConstraints = []
			for i in range(numInverters):
				singleExcludingConstraints.append(vs[i] != solution[i])
			excludingConstraints.append(singleExcludingConstraints)
		
		# Add all the rambus oscillator constraints
		f_sat = And(*allConstraints)
		# Add constraints so that old hyperrectangles are not considered
		if len(excludingConstraints) > 0:
			for constraints in excludingConstraints:
				f_sat = And(f_sat, Or(*constraints))
		
		# Add the constraints to Z3 with a push and pop operation
		#print ("f_sat")
		#print (f_sat)
		s.push()
		s.add(f_sat)
		#print ("s")
		#print (s)
		result = s.check()
		#print (result)
		if result != sat:
			break
		
		m = s.model()
		sol = [None]*numInverters
		for d in m.decls():
			dName = str(d.name())
			firstLetter = dName[0]
			if (dName[0] == "v" and dName[1] == "_"):
				index = int(dName[len(dName) - 1])		
				sol[index] = m[d]


		print ("sol", sol)

		s.pop()
		allSolutions.append(sol)
		
		print ("num solutions found", len(allSolutions))
		

	end = time.time()
	print ("time taken", end - start)
	return allSolutions


if __name__ == "__main__":
	#rambusOscillatorTanh(a = -5.0, numStages = 2, numSolutions = "all", g_cc = 0.5)
	#allSolutions = inverterTanh(-1.0, -5.0)
	allSolutions = inverterLoopTanh(1)
	print ("allSolutions")
	print (allSolutions)
