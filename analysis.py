import numpy as np
import matplotlib.pyplot as plt
from prototype import *
from intervalUtils import *
from circuitModels import InverterLoopMosfet

def analyzeSchmittTrigger():
	inputVoltages = np.linspace(0.35, 0.55, 50)
	posEigenValues = []
	for i in range(len(inputVoltages)):
		statVars = {}
		inputVoltage = inputVoltages[i]
		modelParam = [1.0] #Vdd
		model = SchmittMosfet(modelType = "scMosfet", modelParam = modelParam, inputVoltage = inputVoltage)
		allHypers = []
		solverLoopNoLp(allHypers, model)
		if len(allHypers) != 3:
			raise Exception("analyzeSchmittTrigger inputVoltage "+str(inputVoltage)+
				": there should have been 3 solutions but got " + str(len(allHypers)) + " solutions ")
		
		foundPosEig = False
		for hyper in allHypers:
			exampleVolt = (hyper[:,0] + hyper[:,1])/2.0
			soln = newton(model,exampleVolt)
			if not(soln[0]):
				raise Exception("analyzeSchmittTrigger inputVoltage "+str(inputVoltage)+
				": newton's method should have found solution in  " + str(allHypers))
			soln = soln[1]
			jacAtSoln = model.jacobian(soln)
			eigVals,_ = np.linalg.eig(jacAtSoln)
			maxEig = np.amax(eigVals.real)
			if maxEig > 0:
				foundPosEig = True
				posEigenValues.append(maxEig)

		if not(foundPosEig):
			raise Exception("analyzeSchmittTrigger inputVoltage "+str(inputVoltage)+
				": there should have been one equilibrium point with pos eig value  ")

	posEigenValues = np.array(posEigenValues)
	print ("Vin")
	print (inputVoltages)
	print ("Pos Eig")
	print (posEigenValues)
	plt.plot(inputVoltages, posEigenValues, "*")
	plt.xlabel("Vin")
	plt.ylabel("Real Part of Eigen Values")
	plt.show()


def analyzeInverterLoop():
	posEigenValues = []
	statVars = {}
	modelParam = [1.0] #Vdd
	model = InverterLoopMosfet(modelType = "scMosfet", modelParam = modelParam)
	allHypers = []
	solverLoopNoLp(allHypers, model)
	if len(allHypers) != 3:
		raise Exception("analyzeInverterLoop inputVoltage: there should have been 3 solutions but got " + str(len(allHypers)) + " solutions ")
	
	print ("len(allHypers)", len(allHypers))
	foundPosEig = False
	for hyper in allHypers:
		exampleVolt = (hyper[:,0] + hyper[:,1])/2.0
		soln = newton(model,exampleVolt)
		if not(soln[0]):
			raise Exception("analyzeInverterLoop: newton's method should have found solution in  " + str(allHypers))
		soln = soln[1]
		#print ("soln", soln)
		#print ("current", model.f(soln))
		jacAtSoln = model.jacobian(soln)
		eigVals,_ = np.linalg.eig(jacAtSoln)
		maxEig = np.amax(eigVals.real)
		if maxEig > 0:
			print ("pos eig value at", hyper)
			foundPosEig = True
			posEigenValues.append(maxEig)

	if not(foundPosEig):
		raise Exception("analyzeInverter inputVoltage : there should have been one equilibrium point with pos eig value  ")

	posEigenValues = np.array(posEigenValues)
	print ("Pos Eig")
	print (posEigenValues)


#analyzeSchmittTrigger()
analyzeInverterLoop()
