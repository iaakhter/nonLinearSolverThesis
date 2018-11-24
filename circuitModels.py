# @author Itrat Ahmed Akhter
# Classes defining different models used of our circuits
# Each class has an f method that gives a point and interval evaluation
# of function depending on the type of function argument (point or interval)
# Each class has a jacobian method that gives a point and interval
# jacobian of function depending on the type of function argument
# (point or interval)
# Some classes also have a linearConstraints method that constructs
# linear bounds for the model, solves optimization problems to minimize
# and maximize each variable to give a tighter hyperrectangle

import numpy as np
import lpUtils
from cvxopt import matrix,solvers
from scipy.spatial import ConvexHull
import circuit
import funCompUtils as fcUtils
from intervalBasics import *

'''
Rambus ring oscillator with tanh as the inverter model
@param modelParam is constants in tanh where y = tanh(modelParam[0]*x + modelParam[1])
@param g_cc strength of cross coupled inverters
@param g_fwd strength of forward inverters
@param numStages number of stages in the rambus oscillator
'''
class RambusTanh:
	def __init__(self, modelParam, g_cc, g_fwd, numStages):
		self.modelParam = modelParam # y = tanh(modelParam[0]*x + modelParam[1])
		self.g_cc = g_cc
		self.g_fwd = g_fwd
		self.numStages = numStages
		self.xs = []
		self.ys = []
		self.zs = []
		for i in range(numStages*2):
			self.xs.append("x" + str(i))
			self.ys.append("y" + str(i))
			self.zs.append("z" + str(i))

		self.bounds = []
		for i in range(numStages*2):
			self.bounds.append([-1.0,1.0])

	def f(self, V):
		intervalVal = any([interval_p(x) for x in V])
		lenV = len(V)
		
		if intervalVal:
			fVal = np.zeros((lenV,2))
		else:
			fVal = np.zeros((lenV))

		for i in range(lenV):
			fwdInd = (i-1)%lenV
			ccInd = (i + lenV//2) % lenV
			tanhFwd = fcUtils.tanhFun(V[fwdInd], self.modelParam[0], self.modelParam[1])
			tanhCc = fcUtils.tanhFun(V[ccInd], self.modelParam[0], self.modelParam[1])
			fwdTerm = interval_sub(tanhFwd, V[i])
			ccTerm = interval_sub(tanhCc, V[i])
			fVal[i] = interval_add(interval_mult(self.g_fwd, fwdTerm),interval_mult(self.g_cc, ccTerm))

		return fVal


	'''Get jacobian of rambus oscillator at V
	'''
	def jacobian(self,V):
		lenV = len(V)
		intervalVal = any([interval_p(x) for x in V])

		if intervalVal:
			jac = np.zeros((lenV, lenV, 2))
		else:
			jac = np.zeros((lenV, lenV))
		for i in range(lenV):
			if intervalVal:
				jac[i,i] = interval_fix(-(self.g_fwd + self.g_cc))
			else:
				jac[i,i] = -(self.g_fwd + self.g_cc)
			jac[i,(i-1)%lenV] = interval_mult(self.g_fwd, fcUtils.tanhFunder(V[(i-1)%lenV], self.modelParam[0], self.modelParam[1]))
			jac[i,(i + lenV//2) % lenV] = interval_mult(self.g_cc, fcUtils.tanhFunder(V[(i + lenV//2) % lenV], self.modelParam[0], self.modelParam[1]))

		return jac



	def linearConstraints(self, hyperRectangle):
		solvers.options["show_progress"] = False
		allConstraints = ""
		lenV = self.numStages*2

		allConstraints = ""
		for i in range(lenV):
			fwdInd = (i-1)%lenV
			ccInd = (i+lenV//2)%lenV
			#print ("fwdInd ", fwdInd, " ccInd ", ccInd)
			#print "hyperRectangle[fwdInd][0]", hyperRectangle[fwdInd][0], "hyperRectangle[fwdInd][1]", hyperRectangle[fwdInd][1]
			
			triangleClaimFwd = fcUtils.tanhLinearConstraints(self.modelParam[0], self.modelParam[1], self.xs[fwdInd], self.ys[i], hyperRectangle[fwdInd,0],hyperRectangle[fwdInd,1])
			#print ("triangleClaimFwd", triangleClaimFwd)
			allConstraints += triangleClaimFwd

			triangleClaimCc = fcUtils.tanhLinearConstraints(self.modelParam[0], self.modelParam[1], self.xs[ccInd], self.zs[i], hyperRectangle[ccInd,0],hyperRectangle[ccInd,1])
			#print ("triangleClaimCc", triangleClaimCc)
			allConstraints += triangleClaimCc
				
			allConstraints += str(self.g_fwd) + " " + self.ys[i] + " + " + str(-self.g_fwd-self.g_cc) + \
			" " + self.xs[i] + " + " + str(self.g_cc) + " "  + self.zs[i] + " >= 0.0\n"
			allConstraints += str(self.g_fwd) + " " + self.ys[i] + " + " + str(-self.g_fwd-self.g_cc) + \
			" " + self.xs[i] + " + " + str(self.g_cc) + " "  + self.zs[i] + " <= 0.0\n"
		

		#print "allConstraints"
		#print allConstraints
		variableDict, A, B = lpUtils.constructCoeffMatrices(allConstraints)
		newHyperRectangle = np.copy(hyperRectangle)
		feasible = True
		numSuccessLp, numUnsuccessLp, numTotalLp = 0, 0, 0
		for i in range(lenV):
			numTotalLp += 2
			#print "min max ", i
			minObjConstraint = "min 1 " + self.xs[i]
			maxObjConstraint = "max 1 " + self.xs[i]
			Cmin = lpUtils.constructObjMatrix(minObjConstraint,variableDict)
			Cmax = lpUtils.constructObjMatrix(maxObjConstraint,variableDict)
			minSol = solvers.lp(Cmin,A,B)
			maxSol = solvers.lp(Cmax,A,B)
			if minSol["status"] == "primal infeasible" and maxSol["status"] == "primal infeasible":
				feasible = False
				numSuccessLp += 2
				break
			else:
				if minSol["status"] == "optimal":
					newHyperRectangle[i,0] = minSol['x'][variableDict[self.xs[i]]] - 1e-6
					numSuccessLp += 1
				else:
					numUnsuccessLp += 1
					#print ("min lp not optimal", minSol["status"])
					#print (allConstraints)
				if maxSol["status"] == "optimal":
					newHyperRectangle[i,1] = maxSol['x'][variableDict[self.xs[i]]] + 1e-6
					numSuccessLp += 1
				else:
					numUnsuccessLp += 1

		return [feasible, newHyperRectangle, numTotalLp, numSuccessLp, numUnsuccessLp]


'''
Inverter with tanh as the inverter model
@param modelParam is constants in tanh where y = tanh(modelParam[0]*x + modelParam[1])
@param inputVoltage specific inputVoltage for which we are looking for 
				DC equilibrium points
'''
class InverterTanh:
	def __init__(self, modelParam, inputVoltage):
		self.modelParam = modelParam # y = tanh(modelParam[0]*x + modelParam[1])
		self.bounds = []
		self.bounds.append([-1.0,1.0])
		self.inputVoltage = inputVoltage

	def f(self, V):
		return interval_sub(fcUtils.tanhFun(self.inputVoltage, self.modelParam[0], self.modelParam[1]), V)

	def jacobian(self,V):
		lenV = len(V)
		intervalVal = any([interval_p(x) for x in V])

		if intervalVal:
			jac = np.zeros((1, 1, 2))
			jac[0,0] = [-1.0, -1.0]
		else:
			jac = np.zeros((1, 1))
			jac[0,0] = -1.0
		return jac


'''
Inverter loops with tanh as the inverter model
@param modelParam is constants in tanh where y = tanh(modelParam[0]*x + modelParam[1])
@param numInverters number of inverters in the inverter loop
'''
class InverterLoopTanh:
	def __init__(self, modelParam, numInverters):
		self.modelParam = modelParam # y = tanh(modelParam[0]*x + modelParam[1])
		self.numInverters = numInverters
		self.bounds = []
		for i in range(numInverters):
			self.bounds.append([-1.0, 1.0])

	def f(self, V):
		intervalVal = any([interval_p(x) for x in V])
		lenV = len(V)
		
		if intervalVal:
			fVal = np.zeros((lenV,2))
		else:
			fVal = np.zeros((lenV))

		for i in range(lenV):
			inputInd = i
			outputInd = (i+1)%lenV
			fVal[i] = interval_sub(fcUtils.tanhFun(V[inputInd], self.modelParam[0], self.modelParam[1]),
									V[outputInd])

		return fVal

	def jacobian(self,V):
		lenV = len(V)
		intervalVal = any([interval_p(x) for x in V])

		if intervalVal:
			jac = np.zeros((lenV, lenV, 2))
		else:
			jac = np.zeros((lenV, lenV))
		for i in range(lenV):
			inputInd = i
			outputInd = (i+1)%lenV
			jac[i,inputInd] = fcUtils.tanhFunder(V[inputInd], self.modelParam[0], self.modelParam[1])
			if intervalVal:
				jac[i,outputInd] = interval_fix(-1.0)
			else:
				jac[i,outputInd] = -1.0

		return jac

'''
Rambus ring oscillator with 2 CMOS transistors having
	Mosfet models used as an inverter
@param modelType "lcMosfet" if transistors are long channel and
				"scMosfet" if transistors are short channel
@param modelParam parameters of model depending on whether modelType
				is long channel or short channel
@param g_cc strength of cross coupled inverters
@param g_fwd strength of forward inverters
@param numStages number of stages in the rambus oscillator
'''
class RambusMosfet:
	def __init__(self, modelType, modelParam, g_cc, g_fwd, numStages):
		s0 = 3.0
		self.g_cc = g_cc
		self.g_fwd = g_fwd
		self.numStages = numStages
		lenV = numStages*2

		if modelType == "lcMosfet":
			model = circuit.LcMosfet
			self.Vtp = modelParam[0]
			self.Vtn = modelParam[1]
			self.Vdd = modelParam[2]
			self.Kn = modelParam[3]
			self.Kp = modelParam[4]
			nfet = circuit.MosfetModel('nfet', self.Vtn, self.Kn)
			pfet = circuit.MosfetModel('pfet', self.Vtp, self.Kp)

		elif modelType == "scMosfet":
			model = circuit.ScMosfet
			self.Vdd = modelParam[0]
			nfet = circuit.MosfetModel('nfet')
			pfet = circuit.MosfetModel('pfet')

		# for 4 stage oscillator for example
		# V is like this for Mark's model V = [V0, V1, V2, V3, V4, V5, V6, V7, src, Vdd]

		transistorList = []
		for i in range(lenV):
			fwdInd = (i-1)%(lenV)
			ccInd = (i+lenV//2)%(lenV)
			transistorList.append(model(lenV, fwdInd, i, nfet, s0*g_fwd))
			transistorList.append(model(lenV+1, fwdInd, i, pfet, 2.0*s0*g_fwd))
			transistorList.append(model(lenV, ccInd, i, nfet, s0*g_cc))
			transistorList.append(model(lenV+1, ccInd, i, pfet, 2.0*s0*g_cc))

		self.c = circuit.Circuit(transistorList)

		self.bounds = []
		for i in range(numStages*2):
			self.bounds.append([0.0, self.Vdd])


	def f(self,V):
		myV = [x for x in V] + [0.0, self.Vdd]
		# print 'rambusMosfetMark.f: myV = ' + str(myV)
		funcVal = self.c.f(myV)
		# print 'rambusMosfetMark.f: funcVal = ' + str(funcVal)
		return funcVal[:-2]

	def jacobian(self,V):
		myV = [x for x in V] + [0.0, self.Vdd]
		myJac = self.c.jacobian(myV)
		return np.array(myJac[:-2,:-2])	

	def linearConstraints(self, hyperRectangle):
		lenV = len(hyperRectangle)
		cHyper = [x for x in hyperRectangle] + [0.0, self.Vdd]
		[feasible, newHyper, numTotalLp, numSuccessLp, numUnsuccessLp] = self.c.linearConstraints(cHyper, [lenV, lenV + 1])
		newHyper = newHyper[:-2]
		newHyperMat = np.zeros((lenV,2))
		for i in range(lenV):
			newHyperMat[i,:] = [newHyper[i][0], newHyper[i][1]]
		return [feasible, newHyperMat, numTotalLp, numSuccessLp, numUnsuccessLp]



'''
An inverter modeled by 2 CMOS transistors having long channel or 
short channel Mosfet models
@param modelType "lcMosfet" if transistors are long channel and
				"scMosfet" if transistors are short channel
@param modelParam parameters of model depending on whether modelType
				is long channel or short channel
@param inputVoltage specific inputVoltage for which we are looking for 
				DC equilibrium points
'''
class InverterMosfet:
	def __init__(self, modelType, modelParam, inputVoltage):
		self.inputVoltage = inputVoltage
		s0 = 3.0
		if modelType == "lcMosfet":
			model = circuit.LcMosfet
			self.Vtp = modelParam[0]
			self.Vtn = modelParam[1]
			self.Vdd = modelParam[2]
			self.Kn = modelParam[3]
			self.Kp = modelParam[4]
			nfet = circuit.MosfetModel('nfet', self.Vtn, self.Kn)
			pfet = circuit.MosfetModel('pfet', self.Vtp, self.Kp)

		elif modelType == "scMosfet":
			model = circuit.ScMosfet
			self.Vdd = modelParam[0]
			nfet = circuit.MosfetModel('nfet')
			pfet = circuit.MosfetModel('pfet')


		# V = [outputVoltage, inputVoltage, grnd, Vdd]	
		transistorList = []
		transistorList.append(model(2, 1, 0, nfet, s0))
		transistorList.append(model(3, 1, 0, pfet, s0*2))

		self.c = circuit.Circuit(transistorList)

		self.bounds = []
		# output voltage
		self.bounds.append([0.0, self.Vdd])



	def f(self,V):
		myV = [x for x in V] + [self.inputVoltage, 0.0, self.Vdd]
		funcVal = self.c.f(myV)
		return funcVal[:-3]

	def jacobian(self,V):
		myV = [x for x in V] + [self.inputVoltage, 0.0, self.Vdd]
		myJac = self.c.jacobian(myV)
		return np.array(myJac[:-3,:-3])	

	def linearConstraints(self, hyperRectangle):
		lenV = len(hyperRectangle)
		cHyper = [x for x in hyperRectangle] + [self.inputVoltage, 0.0, self.Vdd]
		[feasible, newHyper, numTotalLp, numSuccessLp, numUnsuccessLp] = self.c.linearConstraints(cHyper, [lenV, lenV + 1])
		newHyper = newHyper[:-3]
		newHyperMat = np.zeros((lenV,2))
		for i in range(lenV):
			newHyperMat[i,:] = [newHyper[i][0], newHyper[i][1]]
		return [feasible, newHyperMat, numTotalLp, numSuccessLp, numUnsuccessLp]

'''
An inverter loop where each inverter is modeled by 2 CMOS transistors 
having long channel or short channel Mosfet models
@param modelType "lcMosfet" if transistors are long channel and
				"scMosfet" if transistors are short channel
@param modelParam parameters of model depending on whether modelType
				is long channel or short channel
@param numInverters number of inverters in the inverter loop
'''
class InverterLoopMosfet:
	def __init__(self, modelType, modelParam, numInverters):
		s0 = 3.0
		if modelType == "lcMosfet":
			model = circuit.LcMosfet
			self.Vtp = modelParam[0]
			self.Vtn = modelParam[1]
			self.Vdd = modelParam[2]
			self.Kn = modelParam[3]
			self.Kp = modelParam[4]
			nfet = circuit.MosfetModel('nfet', self.Vtn, self.Kn)
			pfet = circuit.MosfetModel('pfet', self.Vtp, self.Kp)

		elif modelType == "scMosfet":
			model = circuit.ScMosfet
			self.Vdd = modelParam[0]
			nfet = circuit.MosfetModel('nfet')
			pfet = circuit.MosfetModel('pfet')

		self.numInverters = numInverters

		# V = [V0, V1, V2, ... grnd, Vdd]
		transistorList = []
		for i in range(numInverters):
			transistorList.append(model(numInverters, i, (i+1)%numInverters, nfet, s0))
			transistorList.append(model(numInverters+1, i, (i+1)%numInverters, pfet, s0*2))

		self.c = circuit.Circuit(transistorList)

		self.bounds = []
		for i in range(numInverters):
			self.bounds.append([0.0, self.Vdd])



	def f(self,V):
		myV = [x for x in V] + [0.0, self.Vdd]
		funcVal = self.c.f(myV)
		return funcVal[:-2]

	def jacobian(self,V):
		myV = [x for x in V] + [0.0, self.Vdd]
		myJac = self.c.jacobian(myV)
		return np.array(myJac[:-2,:-2])	

	def linearConstraints(self, hyperRectangle):
		lenV = len(hyperRectangle)
		cHyper = [x for x in hyperRectangle] + [0.0, self.Vdd]
		[feasible, newHyper, numTotalLp, numSuccessLp, numUnsuccessLp] = self.c.linearConstraints(cHyper, [lenV, lenV + 1])
		newHyper = newHyper[:-2]
		newHyperMat = np.zeros((lenV,2))
		for i in range(lenV):
			newHyperMat[i,:] = [newHyper[i][0], newHyper[i][1]]
		return [feasible, newHyperMat, numTotalLp, numSuccessLp, numUnsuccessLp]


'''
Schmitt Trigger where each transistor is modeled by a CMOS transistor
with long channel and short channel Mosfet models
@param modelType "lcMosfet" if transistors are long channel and
				"scMosfet" if transistors are short channel
@param modelParam parameters of model depending on whether modelType
				is long channel or short channel
@param inputVoltage specific inputVoltage for which we are looking for 
				DC equilibrium points
'''
class SchmittMosfet:
	def __init__(self, modelType, modelParam, inputVoltage):
		s0 = 3.0
		self.inputVoltage = inputVoltage

		if modelType == "lcMosfet":
			model = circuit.LcMosfet
			self.Vtp = modelParam[0]
			self.Vtn = modelParam[1]
			self.Vdd = modelParam[2]
			self.Kn = modelParam[3]
			self.Kp = modelParam[4]
			nfet = circuit.MosfetModel('nfet', self.Vtn, self.Kn, gds = "default")
			pfet = circuit.MosfetModel('pfet', self.Vtp, self.Kp, gds = "default")

		elif modelType == "scMosfet":
			model = circuit.ScMosfet
			self.Vdd = modelParam[0]
			nfet = circuit.MosfetModel('nfet')
			pfet = circuit.MosfetModel('pfet')


		# with the voltage array containing [grnd, Vdd, input, X[0], X[1], X[2]]
		# where X[0] is the output voltage and X[1] is the voltage at node with 
		# nfets and X[2] is the voltage at node with pfets

		# src, gate, drain = grnd, input, X[1]
		m0 = model(0, 2, 4, nfet, s0)

		# src, gate, drain = X[1], input, X[0]
		m1 = model(4, 2, 3, nfet, s0)

		# src, gate, drain = X[1], X[0], Vdd
		m2 = model(4, 3, 1, nfet, s0)

		# src, gate, drain = Vdd, input, X[2]
		m3 = model(1, 2, 5, pfet, s0*2.0)

		# src, gate, drain = X[2], in, X[0]
		m4 = model(5, 2, 3, pfet, s0*2.0)

		# src, gate, drain = X[2], X[0], grnd
		m5 = model(5, 3, 0, pfet, s0*2.0)

		self.c = circuit.Circuit([m0, m1, m2, m3, m4, m5])

		self.bounds = []
		for i in range(3):
			self.bounds.append([0.0,self.Vdd])


	def f(self,V):
		myV = [0.0, self.Vdd, self.inputVoltage] + [x for x in V]
		#print ("myV", myV)
		funcVal = self.c.f(myV)
		#print ("myV", myV)
		#print ("funcVal",funcVal)
		return funcVal[3:]


	def jacobian(self,V):
		#print ("calculating jacobian")
		myV = [0.0, self.Vdd, self.inputVoltage] + [x for x in V]
		#print ("V", V)
		jac = self.c.jacobian(myV)
		#print ("jac before")
		#print (jac)
		return jac[3:,3:]


	def linearConstraints(self, hyperRectangle):
		lenV = len(hyperRectangle)
		cHyper = [0.0, self.Vdd, self.inputVoltage] + [x for x in hyperRectangle]
		[feasible, newHyper, numTotalLp, numSuccessLp, numUnsuccessLp]  = self.c.linearConstraints(cHyper, [0, 1])
		newHyper = newHyper[3:]
		newHyperMat = np.zeros((lenV,2))
		for i in range(lenV):
			newHyperMat[i,:] = [newHyper[i][0], newHyper[i][1]]
		return [feasible, newHyperMat, numTotalLp, numSuccessLp, numUnsuccessLp]


