# @author Itrat Ahmed Akhter
# Find dc equilibrium points for circuits involving long channel Mosfet models
# using dReal

from dreal.symbolic import Variable, logical_and, logical_or, tanh
from dreal.symbolic import logical_not
from dreal.symbolic import forall
from dreal.api import CheckSatisfiability, Minimize
import time
import numpy as np


# Constraints for nfet with leakage current
def nFetLeak(Vtn, Vdd, Kn, Sn, src, gate, drain, tI):
	constraints = []
	gds = 1e-8
	if type(src) == Variable or type(gate) == Variable:
		constraints.append(logical_or(logical_and(src > drain, gate - drain <= Vtn, tI == -(src - drain)*gds),
										logical_and(src > drain, gate - drain >= Vtn, src - drain >= gate - drain - Vtn, tI == -0.5*Sn*Kn*(gate - drain - Vtn)*(gate - drain - Vtn) -(src - drain)*gds),
										logical_and(src > drain, gate - drain >= Vtn, src - drain <= gate - drain - Vtn, tI == -Sn*Kn*(gate - drain - Vtn - (src - drain)/2.0)*(src - drain) -(src - drain)*gds),
										logical_and(src <= drain, gate - src <= Vtn, tI == (drain - src)*gds),
										logical_and(src <= drain, gate - src >= Vtn, drain - src >= gate - src - Vtn, tI == 0.5*Sn*Kn*(gate - src - Vtn)*(gate - src - Vtn) + (drain - src)*gds),
										logical_and(src <= drain, gate - src >= Vtn, drain - src <= gate - src - Vtn, tI == Sn*Kn*(gate - src - Vtn - (drain - src)/2.0)*(drain - src) + (drain - src)*gds)))
	else:
		if gate - src <= Vtn:
			constraints.append(logical_or(logical_and(src > drain, gate - drain <= Vtn, tI == -(src - drain)*gds),
											logical_and(src > drain, gate - drain >= Vtn, src - drain >= gate - drain - Vtn, tI == -0.5*Sn*Kn*(gate - drain - Vtn)*(gate - drain - Vtn) -(src - drain)*gds),
											logical_and(src > drain, gate - drain >= Vtn, src - drain <= gate - drain - Vtn, tI == -Sn*Kn*(gate - drain - Vtn - (src - drain)/2.0)*(src - drain) -(src - drain)*gds),
											logical_and(src <= drain, tI == (drain - src)*gds)))
		else:
			constraints.append(logical_or(logical_and(src > drain, gate - drain <= Vtn, tI == -(src - drain)*gds),
											logical_and(src > drain, gate - drain >= Vtn, src - drain >= gate - drain - Vtn, tI == -0.5*Sn*Kn*(gate - drain - Vtn)*(gate - drain - Vtn) -(src - drain)*gds),
											logical_and(src > drain, gate - drain >= Vtn, src - drain <= gate - drain - Vtn, tI == -Sn*Kn*(gate - drain - Vtn - (src - drain)/2.0)*(src - drain) -(src - drain)*gds),
											logical_and(src <= drain, drain - src >= gate - src - Vtn, tI == 0.5*Sn*Kn*(gate - src - Vtn)*(gate - src - Vtn) + (drain - src)*gds),
											logical_and(src <= drain, drain - src <= gate - src - Vtn, tI == Sn*Kn*(gate - src - Vtn - (drain - src)/2.0)*(drain - src) + (drain - src)*gds)))
	
	return constraints

# Constraints for nfet without leakage current
def nFet(Vtn, Vdd, Kn, Sn, src, gate, drain, tI):
	constraints = []
	if type(src) == Variable or type(gate) == Variable:
		constraints.append(logical_or(logical_and(src > drain, gate - drain <= Vtn, tI == 0.0),
										logical_and(src > drain, gate - drain >= Vtn, src - drain >= gate - drain - Vtn, tI == -0.5*Sn*Kn*(gate - drain - Vtn)*(gate - drain - Vtn)),
										logical_and(src > drain, gate - drain >= Vtn, src - drain <= gate - drain - Vtn, tI == -Sn*Kn*(gate - drain - Vtn - (src - drain)/2.0)*(src - drain)),
										logical_and(src <= drain, gate - src <= Vtn, tI == 0.0),
										logical_and(src <= drain, gate - src >= Vtn, drain - src >= gate - src - Vtn, tI == 0.5*Sn*Kn*(gate - src - Vtn)*(gate - src - Vtn)),
										logical_and(src <= drain, gate - src >= Vtn, drain - src <= gate - src - Vtn, tI == Sn*Kn*(gate - src - Vtn - (drain - src)/2.0)*(drain - src))))
	else:
		if gate - src <= Vtn:
			constraints.append(logical_or(logical_and(src > drain, gate - drain <= Vtn, tI == 0.0),
											logical_and(src > drain, gate - drain >= Vtn, src - drain >= gate - drain - Vtn, tI == -0.5*Sn*Kn*(gate - drain - Vtn)*(gate - drain - Vtn)),
											logical_and(src > drain, gate - drain >= Vtn, src - drain <= gate - drain - Vtn, tI == -Sn*Kn*(gate - drain - Vtn - (src - drain)/2.0)*(src - drain)),
											logical_and(src <= drain, tI == 0.0)))
		else:
			constraints.append(logical_or(logical_and(src > drain, gate - drain <= Vtn, tI == 0.0),
											logical_and(src > drain, gate - drain >= Vtn, src - drain >= gate - drain - Vtn, tI == -0.5*Sn*Kn*(gate - drain - Vtn)*(gate - drain - Vtn)),
											logical_and(src > drain, gate - drain >= Vtn, src - drain <= gate - drain - Vtn, tI == -Sn*Kn*(gate - drain - Vtn - (src - drain)/2.0)*(src - drain)),
											logical_and(src <= drain, drain - src >= gate - src - Vtn, tI == 0.5*Sn*Kn*(gate - src - Vtn)*(gate - src - Vtn)),
											logical_and(src <= drain, drain - src <= gate - src - Vtn, tI == Sn*Kn*(gate - src - Vtn - (drain - src)/2.0)*(drain - src))))
	
	return constraints

# Constraints for pfet with leakage current
def pFetLeak(Vtp, Vdd, Kp, Sp, src, gate, drain, tI):
	constraints = []
	gds = 1e-8
	if type(src) == Variable or type(gate) == Variable:
		constraints.append(logical_or(logical_and(src < drain, gate - drain >= Vtp, tI == -(src - drain)*gds),
										logical_and(src < drain, gate - drain <= Vtp, src - drain <= gate - drain - Vtp, tI == -0.5*Sp*Kp*(gate - drain - Vtp)*(gate - drain - Vtp) -(src - drain)*gds),
										logical_and(src < drain, gate - drain <= Vtp, src - drain >= gate - drain - Vtp, tI == -Sp*Kp*(gate - drain - Vtp - (src - drain)/2.0)*(src - drain) -(src - drain)*gds),
										logical_and(src >= drain, gate - src >= Vtp, tI == (drain - src)*gds),
										logical_and(src >= drain, gate - src <= Vtp, drain - src <= gate - src - Vtp, tI == 0.5*Sp*Kp*(gate - src - Vtp)*(gate - src - Vtp) + (drain - src)*gds),
										logical_and(src >= drain, gate - src <= Vtp, drain - src >= gate - src - Vtp, tI == Sp*Kp*(gate - src - Vtp - (drain - src)/2.0)*(drain - src) + (drain - src)*gds)))
	else:
		if gate - src >= Vtp:
			constraints.append(logical_or(logical_and(src < drain, gate - drain >= Vtp, tI == -(src - drain)*gds),
											logical_and(src < drain, gate - drain <= Vtp, src - drain <= gate - drain - Vtp, tI == -0.5*Sp*Kp*(gate - drain - Vtp)*(gate - drain - Vtp) -(src - drain)*gds),
											logical_and(src < drain, gate - drain <= Vtp, src - drain >= gate - drain - Vtp, tI == -Sp*Kp*(gate - drain - Vtp - (src - drain)/2.0)*(src - drain) -(src - drain)*gds),
											logical_and(src >= drain, tI == (drain - src)*gds)))
		else:
			constraints.append(logical_or(logical_and(src < drain, gate - drain >= Vtp, tI == -(src - drain)*gds),
											logical_and(src < drain, gate - drain <= Vtp, src - drain <= gate - drain - Vtp, tI == -0.5*Sp*Kp*(gate - drain - Vtp)*(gate - drain - Vtp) -(src - drain)*gds),
											logical_and(src < drain, gate - drain <= Vtp, src - drain >= gate - drain - Vtp, tI == -Sp*Kp*(gate - drain - Vtp - (src - drain)/2.0)*(src - drain) -(src - drain)*gds),
											logical_and(src >= drain, drain - src <= gate - src - Vtp, tI == 0.5*Sp*Kp*(gate - src - Vtp)*(gate - src - Vtp) + (drain - src)*gds),
											logical_and(src >= drain, drain - src >= gate - src - Vtp, tI == Sp*Kp*(gate - src - Vtp - (drain - src)/2.0)*(drain - src) + (drain - src)*gds)))
	
	return constraints

# Constraints for pfet without leakage current
def pFet(Vtp, Vdd, Kp, Sp, src, gate, drain, tI):
	constraints = []
	if type(src) == Variable or type(gate) == Variable:
		constraints.append(logical_or(logical_and(src < drain, gate - drain >= Vtp, tI == 0.0),
										logical_and(src < drain, gate - drain <= Vtp, src - drain <= gate - drain - Vtp, tI == -0.5*Sp*Kp*(gate - drain - Vtp)*(gate - drain - Vtp)),
										logical_and(src < drain, gate - drain <= Vtp, src - drain >= gate - drain - Vtp, tI == -Sp*Kp*(gate - drain - Vtp - (src - drain)/2.0)*(src - drain)),
										logical_and(src >= drain, gate - src >= Vtp, tI == 0.0),
										logical_and(src >= drain, gate - src <= Vtp, drain - src <= gate - src - Vtp, tI == 0.5*Sp*Kp*(gate - src - Vtp)*(gate - src - Vtp)),
										logical_and(src >= drain, gate - src <= Vtp, drain - src >= gate - src - Vtp, tI == Sp*Kp*(gate - src - Vtp - (drain - src)/2.0)*(drain - src))))
	else:
		if gate - src >= Vtp:
			constraints.append(logical_or(logical_and(src < drain, gate - drain >= Vtp, tI == 0.0),
											logical_and(src < drain, gate - drain <= Vtp, src - drain <= gate - drain - Vtp, tI == -0.5*Sp*Kp*(gate - drain - Vtp)*(gate - drain - Vtp)),
											logical_and(src < drain, gate - drain <= Vtp, src - drain >= gate - drain - Vtp, tI == -Sp*Kp*(gate - drain - Vtp - (src - drain)/2.0)*(src - drain)),
											logical_and(src >= drain, tI == 0.0)))
		else:
			constraints.append(logical_or(logical_and(src < drain, gate - drain >= Vtp, tI == 0.0),
											logical_and(src < drain, gate - drain <= Vtp, src - drain <= gate - drain - Vtp, tI == -0.5*Sp*Kp*(gate - drain - Vtp)*(gate - drain - Vtp)),
											logical_and(src < drain, gate - drain <= Vtp, src - drain >= gate - drain - Vtp, tI == -Sp*Kp*(gate - drain - Vtp - (src - drain)/2.0)*(src - drain)),
											logical_and(src >= drain, drain - src <= gate - src - Vtp, tI == 0.5*Sp*Kp*(gate - src - Vtp)*(gate - src - Vtp)),
											logical_and(src >= drain, drain - src >= gate - src - Vtp, tI == Sp*Kp*(gate - src - Vtp - (drain - src)/2.0)*(drain - src))))
	
	return constraints

# Try and find dc equilibrium points for rambus oscillator with long channel mosfet model
# numStages indicates the number of stages in the rambus oscillator
# numSolutions indicates the number of solutions we want dReal to find
# g_cc is the strength of the cross coupled inverter as compared to that of forward (g_fwd = 1.0)
# Vtp, Vtn, Vdd, Kn, Kp, Sn are parameters for the long channel model
def rambusOscillatorLcMosfet(numStages, numSolutions = "all", g_cc = 0.5, Vtp = -0.4, Vtn = 0.4, Vdd = 1.8, Kn = 270*1e-6, Kp = -90*1e-6, Sn = 3.0):
	epsilon = 1e-14
	start = time.time()
	#print ("Vtp", Vtp, "Vtn", Vtn, "Vdd", Vdd, "Kn", Kn, "Kp", Kp, "Sn", Sn)
	g_fwd = 1.0
	lenV = numStages*2
	Sp = 2*Sn

	vs = []
	ifwdNs = []
	ifwdPs = []
	iccNs = []
	iccPs = []
	for i in range(lenV):
		vs.append(Variable("v" + str(i)))
		ifwdNs.append(Variable("ifwdN" + str(i)))
		ifwdPs.append(Variable("ifwdP" + str(i)))
		iccNs.append(Variable("iccN" + str(i)))
		iccPs.append(Variable("iccP" + str(i)))

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

	allSolutions = []
	while True:
		if numSolutions != "all" and len(allSolutions) == numSolutions:
			break
		
		# Store constraints pruning search space so that
		# old hyperrectangles are not considered
		excludingConstraints = []
		for solution in allSolutions:
			singleExcludingConstraints = []
			for i in range(lenV):
				singleExcludingConstraints.append(vs[i] < solution[i][0])
				singleExcludingConstraints.append(vs[i] > solution[i][1])
			excludingConstraints.append(singleExcludingConstraints)
		
		#print ("allConstraints")
		#print (allConstraints)
		f_sat = logical_and(*allConstraints)
		if len(excludingConstraints) > 0:
			for constraints in excludingConstraints:
				f_sat = logical_and(f_sat, logical_or(*constraints))
		
		#print ("f_sat")
		#print (f_sat)
		result = CheckSatisfiability(f_sat, epsilon)
		#print (result)
		if result is None:
			break
		hyper = np.zeros((lenV,2))
		for i in range(lenV):
			hyper[i,:] = [result[vs[i]].lb() - 1000*epsilon, result[vs[i]].ub() + 1000*epsilon]

		#print ("hyper", hyper)
		allSolutions.append(hyper)

		print ("num solutions found", len(allSolutions))


	end = time.time()
	print ("time taken", end - start)
	return allSolutions


# Try and find dc equilibrium points for schmitt trigger with long channel mosfet model
# for a specific input voltage
# Vtp, Vtn, Vdd, Kn, Kp, Sn are parameters for the long channel model
# numSolutions indicates the number of solutions we want dReal to find
def schmittTriggerLcMosfet(inputVoltage, Vtp = -0.4, Vtn = 0.4, Vdd = 1.8, Kn = 270*1e-6, Kp = -90*1e-6, Sn = 3.0, numSolutions = "all"):
	epsilon = 1e-14
	start = time.time()
	#print ("Vtp", Vtp, "Vtn", Vtn, "Vdd", Vdd, "Kn", Kn, "Kp", Kp, "Sn", Sn)
	Sp = Sn *2.0

	lenV = 3

	vs = []
	tIs = []
	nIs = []

	for i in range(lenV):
		vs.append(Variable("v" + str(i)))
		nIs.append(Variable("nI" + str(i)))
	
	for i in range(lenV*2):
		tIs.append(Variable("tI" + str(i)))

	allConstraints = []	
	for i in range(lenV):
		allConstraints.append(vs[i] >= 0.0)
		allConstraints.append(vs[i] <= 1.8)
	allConstraints += nFetLeak(Vtn, Vdd, Kn, Sn, 0.0, inputVoltage, vs[1], tIs[0])
	allConstraints += nFetLeak(Vtn, Vdd, Kn, Sn, vs[1], inputVoltage, vs[0], tIs[1])
	allConstraints += nFetLeak(Vtn, Vdd, Kn, Sn, vs[1], vs[0], Vdd, tIs[2])
	allConstraints += pFetLeak(Vtp, Vdd, Kp, Sp, Vdd, inputVoltage, vs[2], tIs[3])
	allConstraints += pFetLeak(Vtp, Vdd, Kp, Sp, vs[2], inputVoltage, vs[0], tIs[4])
	allConstraints += pFetLeak(Vtp, Vdd, Kp, Sp, vs[2], vs[0], 0.0, tIs[5])
	allConstraints.append(nIs[0] == -tIs[4] - tIs[1])
	allConstraints.append(nIs[1] == -tIs[0] + tIs[1] + tIs[2])
	allConstraints.append(nIs[2] == -tIs[3] + tIs[5] + tIs[4])
	for i in range(lenV):
		allConstraints.append(nIs[i] == 0.0)

	allSolutions = []
	while True:
		if numSolutions != "all" and len(allSolutions) == numSolutions:
			break
		
		# Store constraints pruning search space so that
		# old hyperrectangles are not considered
		excludingConstraints = []
		for solution in allSolutions:
			singleExcludingConstraints = []
			for i in range(lenV):
				singleExcludingConstraints.append(vs[i] < solution[i][0])
				singleExcludingConstraints.append(vs[i] > solution[i][1])
			excludingConstraints.append(singleExcludingConstraints)
		
		#print ("allConstraints")
		#print (allConstraints)
		#print ("numConstraints", len(allConstraints))

		f_sat = logical_and(*allConstraints)
		if len(excludingConstraints) > 0:
			for constraints in excludingConstraints:
				f_sat = logical_and(f_sat, logical_or(*constraints))
		
		#print ("f_sat")
		#print (f_sat)
		result = CheckSatisfiability(f_sat, epsilon)
		#print (result)
		if result is None:
			break
		hyper = np.zeros((lenV,2))
		for i in range(lenV):
			hyper[i,:] = [result[vs[i]].lb() - 1000*epsilon, result[vs[i]].ub() + 1000*epsilon]

		#print ("hyper", hyper)
		allSolutions.append(hyper)

		print ("num solutions found", len(allSolutions))


	end = time.time()
	print ("time taken", end - start)
	return allSolutions


# Try and find dc equilibrium points for inverter with long channel mosfet model
# for a specific input voltage
# Vtp, Vtn, Vdd, Kn, Kp, Sn are parameters for the long channel model
# numSolutions indicates the number of solutions we want dReal to find
def inverterLcMosfet(inputVoltage, Vtp = -0.4, Vtn = 0.4, Vdd = 1.8, Kn = 270*1e-6, Kp = -90*1e-6, Sn = 3.0, numSolutions = "all"):
	epsilon = 1e-14
	start = time.time()
	Sp = Sn*3.0

	outputVolt = Variable("outputVolt")
	iP = Variable("iP")
	iN = Variable("iN")
	
	allConstraints = []	
	allConstraints.append(outputVolt >= 0.0)
	allConstraints.append(outputVolt <= Vdd)
	allConstraints.append(-iP-iN == 0)
	allConstraints += nFet(Vtn, Vdd, Kn, Sn, 0.0, inputVoltage, outputVolt, iN)
	allConstraints += pFet(Vtp, Vdd, Kp, Sp, Vdd, inputVoltage, outputVolt, iP)

	allSolutions = []
	while True:
		if numSolutions != "all" and len(allSolutions) == numSolutions:
			break
		
		# Store constraints pruning search space so that
		# old hyperrectangles are not considered
		excludingConstraints = []
		for solution in allSolutions:
			singleExcludingConstraints = []
			singleExcludingConstraints.append(outputVolt < solution[0][0])
			singleExcludingConstraints.append(outputVolt > solution[0][1])
			excludingConstraints.append(singleExcludingConstraints)
		
		#print ("allConstraints")
		#print (allConstraints)
		f_sat = logical_and(*allConstraints)
		if len(excludingConstraints) > 0:
			for constraints in excludingConstraints:
				f_sat = logical_and(f_sat, logical_or(*constraints))
		
		#print ("f_sat")
		#print (f_sat)
		result = CheckSatisfiability(f_sat, epsilon)
		#print (result)
		if result is None:
			break
		hyper = np.zeros((1,2))
		hyper[0,:] = [result[outputVolt].lb() - 1000*epsilon, result[outputVolt].ub() + 1000*epsilon]

		#print ("hyper", hyper)
		allSolutions.append(hyper)

		print ("num solutions found", len(allSolutions))

	
	end = time.time()
	print ("time taken", end - start)
	return allSolutions


def inverterLoopLcMosfet(numInverters, numSolutions = "all", Vtp = -0.4, Vtn = 0.4, Vdd = 1.8, Kn = 270*1e-6, Kp = -90*1e-6, Sn = 3.0):
	epsilon = 1e-14
	start = time.time()
	#print ("Vtp", Vtp, "Vtn", Vtn, "Vdd", Vdd, "Kn", Kn, "Kp", Kp, "Sn", Sn)
	Sp = 2*Sn

	vs = []
	iNs = []
	iPs = []

	for i in range(numInverters):
		vs.append(Variable("v" + str(i)))
		iNs.append(Variable("iN" + str(i)))
		iPs.append(Variable("iP" + str(i)))

	allConstraints = []	
	for i in range(numInverters):
		allConstraints.append(vs[i] >= 0.0)
		allConstraints.append(vs[i] <= Vdd)
		allConstraints.append(-iNs[i]-iPs[i] == 0)
		inputInd = i
		outputInd = (i+1)%numInverters
		allConstraints += nFet(Vtn, Vdd, Kn, Sn, 0.0, vs[inputInd], vs[outputInd], iNs[i])
		allConstraints += pFet(Vtp, Vdd, Kp, Sp, Vdd, vs[inputInd], vs[outputInd], iPs[i])

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
				singleExcludingConstraints.append(vs[i] < solution[i][0])
				singleExcludingConstraints.append(vs[i] > solution[i][1])
			excludingConstraints.append(singleExcludingConstraints)
		
		#print ("allConstraints")
		#print (allConstraints)
		f_sat = logical_and(*allConstraints)
		if len(excludingConstraints) > 0:
			for constraints in excludingConstraints:
				f_sat = logical_and(f_sat, logical_or(*constraints))
		
		#print ("f_sat")
		#print (f_sat)
		result = CheckSatisfiability(f_sat, epsilon)
		#print (result)
		if result is None:
			break
		hyper = np.zeros((numInverters,2))
		for i in range(numInverters):
			hyper[i,:] = [result[vs[i]].lb() - 1000*epsilon, result[vs[i]].ub() + 1000*epsilon]

		#print ("hyper", hyper)
		allSolutions.append(hyper)

		print ("num solutions found", len(allSolutions))


	end = time.time()
	print ("time taken", end - start)
	return allSolutions

if __name__ == "__main__":
	allSolutions = inverterLoopLcMosfet(1)
	print ("allSolutions")
	for solution in allSolutions:
		print ("solution")
		for i in range(solution.shape[0]):
			print (solution[i,0], solution[i,1])

