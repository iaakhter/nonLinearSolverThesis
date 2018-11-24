# @author Itrat Ahmed Akhter
# Find dc equilibrium points for circuits involving long channel Mosfet models
# using z3

from z3 import *
import time
import numpy as np

# Constraints for nfet with leakage current
def nFetLeak(Vtn, Vdd, Kn, Sn, src, gate, drain, tI):
	constraints = []
	gds = 1e-8
	constraints.append(Or(And(src > drain, gate - drain <= Vtn, tI == -(src - drain)*gds),
							And(src > drain, gate - drain >= Vtn, src - drain >= gate - drain - Vtn, tI == -0.5*Sn*Kn*(gate - drain - Vtn)*(gate - drain - Vtn) -(src - drain)*gds),
							And(src > drain, gate - drain >= Vtn, src - drain <= gate - drain - Vtn, tI == -Sn*Kn*(gate - drain - Vtn - (src - drain)/2.0)*(src - drain) -(src - drain)*gds),
							And(src <= drain, gate - src <= Vtn, tI == (drain - src)*gds),
							And(src <= drain, gate - src >= Vtn, drain - src >= gate - src - Vtn, tI == 0.5*Sn*Kn*(gate - src - Vtn)*(gate - src - Vtn) + (drain - src)*gds),
							And(src <= drain, gate - src >= Vtn, drain - src <= gate - src - Vtn, tI == Sn*Kn*(gate - src - Vtn - (drain - src)/2.0)*(drain - src) + (drain - src)*gds)))
	
	return constraints

# Constraints for nfet without leakage current
def nFet(Vtn, Vdd, Kn, Sn, src, gate, drain, tI):
	constraints = []
	constraints.append(Or(And(src > drain, gate - drain <= Vtn, tI == 0.0),
							And(src > drain, gate - drain >= Vtn, src - drain >= gate - drain - Vtn, tI == -0.5*Sn*Kn*(gate - drain - Vtn)*(gate - drain - Vtn)),
							And(src > drain, gate - drain >= Vtn, src - drain <= gate - drain - Vtn, tI == -Sn*Kn*(gate - drain - Vtn - (src - drain)/2.0)*(src - drain)),
							And(src <= drain, gate - src <= Vtn, tI == 0.0),
							And(src <= drain, gate - src >= Vtn, drain - src >= gate - src - Vtn, tI == 0.5*Sn*Kn*(gate - src - Vtn)*(gate - src - Vtn)),
							And(src <= drain, gate - src >= Vtn, drain - src <= gate - src - Vtn, tI == Sn*Kn*(gate - src - Vtn - (drain - src)/2.0)*(drain - src))))
	
	return constraints

# Constraints for pfet with leakage current
def pFetLeak(Vtp, Vdd, Kp, Sp, src, gate, drain, tI):
	constraints = []
	gds = 1e-8
	constraints.append(Or(And(src < drain, gate - drain >= Vtp, tI == -(src - drain)*gds),
							And(src < drain, gate - drain <= Vtp, src - drain <= gate - drain - Vtp, tI == -0.5*Sp*Kp*(gate - drain - Vtp)*(gate - drain - Vtp) -(src - drain)*gds),
							And(src < drain, gate - drain <= Vtp, src - drain >= gate - drain - Vtp, tI == -Sp*Kp*(gate - drain - Vtp - (src - drain)/2.0)*(src - drain) -(src - drain)*gds),
							And(src >= drain, gate - src >= Vtp, tI == (drain - src)*gds),
							And(src >= drain, gate - src <= Vtp, drain - src <= gate - src - Vtp, tI == 0.5*Sp*Kp*(gate - src - Vtp)*(gate - src - Vtp) + (drain - src)*gds),
							And(src >= drain, gate - src <= Vtp, drain - src >= gate - src - Vtp, tI == Sp*Kp*(gate - src - Vtp - (drain - src)/2.0)*(drain - src) + (drain - src)*gds)))
	
	return constraints

# Constraints for pfet without leakage current
def pFet(Vtp, Vdd, Kp, Sp, src, gate, drain, tI):
	constraints = []
	constraints.append(Or(And(src < drain, gate - drain >= Vtp, tI == 0.0),
							And(src < drain, gate - drain <= Vtp, src - drain <= gate - drain - Vtp, tI == -0.5*Sp*Kp*(gate - drain - Vtp)*(gate - drain - Vtp)),
							And(src < drain, gate - drain <= Vtp, src - drain >= gate - drain - Vtp, tI == -Sp*Kp*(gate - drain - Vtp - (src - drain)/2.0)*(src - drain)),
							And(src >= drain, gate - src >= Vtp, tI == 0.0),
							And(src >= drain, gate - src <= Vtp, drain - src <= gate - src - Vtp, tI == 0.5*Sp*Kp*(gate - src - Vtp)*(gate - src - Vtp)),
							And(src >= drain, gate - src <= Vtp, drain - src >= gate - src - Vtp, tI == Sp*Kp*(gate - src - Vtp - (drain - src)/2.0)*(drain - src))))
	
	return constraints


# Try and find dc equilibrium points for rambus oscillator with long channel mosfet model
# numStages indicates the number of stages in the rambus oscillator
# numSolutions indicates the number of solutions we want Z3 to find
# g_cc is the strength of the cross coupled inverter as compared to that of forward (g_fwd = 1.0)
# Vtp, Vtn, Vdd, Kn, Kp, Sn are parameters for the long channel model
def rambusOscillatorLcMosfet(numStages, numSolutions = "all", g_cc = 0.5, Vtp = -0.4, Vtn = 0.4, Vdd = 1.8, Kn = 270*1e-6, Kp = -90*1e-6, Sn = 3.0):
	start = time.time()
	#print ("Vtp", Vtp, "Vtn", Vtn, "Vdd", Vdd, "Kn", Kn, "Kp", Kp, "Sn", Sn)
	g_fwd = 1.0
	Sp = Sn *2.0
	lenV = numStages*2

	set_option(rational_to_decimal=True)

	vs = RealVector('v', lenV)
	ifwdNs = RealVector('ifwdN', lenV)
	ifwdPs = RealVector('ifwdP', lenV)
	iccNs = RealVector('iccN', lenV)
	iccPs = RealVector('iccP', lenV)

	allConstraints = []	
	for i in range(lenV):
		allConstraints.append(vs[i] >= 0.0)
		allConstraints.append(vs[i] <= Vdd)
		allConstraints.append(g_fwd*(-ifwdNs[i]-ifwdPs[i]) + g_cc*(-iccNs[i]-iccPs[i]) == 0)
		fwdInd = (i-1)%lenV
		ccInd = (i+lenV//2)%lenV
		fwdConstraints = nFet(Vtn, Vdd, Kn, Sn, 0.0, vs[fwdInd], vs[i], ifwdNs[i])
		fwdConstraints += pFet(Vtp, Vdd, Kp, Sp, Vdd, vs[fwdInd], vs[i], ifwdPs[i])
		ccConstraints = nFet(Vtn, Vdd, Kn, Sn, 0.0, vs[ccInd], vs[i], iccNs[i])
		ccConstraints += pFet(Vtp, Vdd, Kp, Sp, Vdd, vs[ccInd], vs[i], iccPs[i])
		allConstraints += fwdConstraints + ccConstraints

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
		
		#print ("allConstraints")
		#print (allConstraints)
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
		hyper = np.zeros((lenV, 2))
		sol = np.zeros((lenV))
		for d in m.decls():
			dName = str(d.name())
			firstLetter = dName[0]
			if (dName[0] == "v" and dName[1] == "_"):
				index = int(dName[len(dName) - 1])
				if str(m[d])[-1] == "?":
					val = float(str(m[d])[:-1])
				else:
					val = float(Fraction(str(m[d])))
				sol[index] = val

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
		print (solution)'''

	end = time.time()
	print ("time taken", end - start)
	return allSolutions

# Try and find dc equilibrium points for schmitt trigger with long channel mosfet model
# for a specific input voltage
# Vtp, Vtn, Vdd, Kn, Kp, Sn are parameters for the long channel model
# numSolutions indicates the number of solutions we want Z3 to find
def schmittTriggerLcMosfet(inputVoltage, Vtp = -0.4, Vtn = 0.4, Vdd = 1.8, Kn = 270*1e-6, Kp = -90*1e-6, Sn = 3.0, numSolutions = "all"):
	start = time.time()
	#print ("Vtp", Vtp, "Vtn", Vtn, "Vdd", Vdd, "Kn", Kn, "Kp", Kp, "Sn", Sn)
	Sp = Sn *2.0

	set_option(rational_to_decimal=True)

	vs = RealVector('v', 3)
	tIs = RealVector('tI', 6)
	nIs = RealVector('nI',3)
	lenV = 3
	s = Solver()

	allConstraints = []
	for i in range(lenV):
		allConstraints.append(vs[i] >= 0.0)
		allConstraints.append(vs[i] <= Vdd)
		allConstraints.append(nIs[i] == 0.0)
	allConstraints += nFetLeak(Vtn, Vdd, Kn, Sn, 0.0, inputVoltage, vs[1], tIs[0])
	allConstraints += nFetLeak(Vtn, Vdd, Kn, Sn, vs[1], inputVoltage, vs[0], tIs[1])
	allConstraints += nFetLeak(Vtn, Vdd, Kn, Sn, vs[1], vs[0], Vdd, tIs[2])
	allConstraints += pFetLeak(Vtp, Vdd, Kp, Sp, Vdd, inputVoltage, vs[2], tIs[3])
	allConstraints += pFetLeak(Vtp, Vdd, Kp, Sp, vs[2], inputVoltage, vs[0], tIs[4])
	allConstraints += pFetLeak(Vtp, Vdd, Kp, Sp, vs[2], vs[0], 0.0, tIs[5])
	allConstraints.append(nIs[0] == -tIs[4] - tIs[1])
	allConstraints.append(nIs[1] == -tIs[0] + tIs[1] + tIs[2])
	allConstraints.append(nIs[2] == -tIs[3] + tIs[5] + tIs[4])

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
		print (solution)'''

	end = time.time()
	print ("time taken", end - start)
	return allSolutions

# Try and find dc equilibrium points for inverter with long channel mosfet model
# for a specific input voltage
# Vtp, Vtn, Vdd, Kn, Kp, Sn are parameters for the long channel model
# numSolutions indicates the number of solutions we want Z3 to find
def inverterLcMosfet(inputVoltage, Vtp = -0.4, Vtn = 0.4, Vdd = 1.8, Kn = 270*1e-6, Kp = -90*1e-6, Sn = 3.0, numSolutions = "all"):
	start = time.time()
	Sp = Sn*3.0
	set_option(rational_to_decimal=True)
	outputVolt = Real("outputVolt")
	iP = Real("iP")
	iN = Real("iN")

	allConstraints = []	
	allConstraints.append(outputVolt >= 0.0)
	allConstraints.append(outputVolt <= Vdd)
	allConstraints.append(-iP-iN == 0)
	allConstraints += nFet(Vtn, Vdd, Kn, Sn, 0.0, inputVoltage, outputVolt, iN)
	allConstraints += pFet(Vtp, Vdd, Kp, Sp, Vdd, inputVoltage, outputVolt, iP)
	
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


def inverterLoopLcMosfet(numInverters, numSolutions = "all", Vtp = -0.4, Vtn = 0.4, Vdd = 1.8, Kn = 270*1e-6, Kp = -90*1e-6, Sn = 3.0):
	start = time.time()
	#print ("Vtp", Vtp, "Vtn", Vtn, "Vdd", Vdd, "Kn", Kn, "Kp", Kp, "Sn", Sn)
	g_fwd = 1.0
	Sp = Sn *2.0

	set_option(rational_to_decimal=True)

	vs = RealVector('v', numInverters)
	iNs = RealVector('iN', numInverters)
	iPs = RealVector('iP', numInverters)

	allConstraints = []	
	for i in range(numInverters):
		allConstraints.append(vs[i] >= 0.0)
		allConstraints.append(vs[i] <= Vdd)
		inputInd = i
		outputInd = (i+1)%numInverters
		allConstraints += nFet(Vtn, Vdd, Kn, Sn, 0.0, vs[inputInd], vs[outputInd], iNs[i])
		allConstraints += pFet(Vtp, Vdd, Kp, Sp, Vdd, vs[inputInd], vs[outputInd], iPs[i])
		allConstraints.append(-iNs[i]-iPs[i] == 0.0)

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
			for i in range(numInverters):
				singleExcludingConstraints.append(vs[i] != solution[i])
			excludingConstraints.append(singleExcludingConstraints)
		
		#print ("allConstraints")
		#print (allConstraints)
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

	
	'''print ("all solutions")
	for solution in allSolutions:
		print (solution)'''

	end = time.time()
	print ("time taken", end - start)
	return allSolutions

if __name__ == "__main__":
	#allSolutions = inverterLcMosfet(0.0)
	allSolutions = inverterLoopLcMosfet(3)
	print ("allSolutions")
	print (allSolutions)



