# @author Itrat Ahmed Akhter
# Functions to analyze dc equlibrium points
import numpy as np
import intervalUtils

def printSol(allHypers, model):
	# categorize solutions found
	sampleSols, rotatedSols, stableSols, unstableSols = categorizeSolutions(allHypers,model)

	for hi in range(len(sampleSols)):
		print ("equivalence class# ", hi)
		print ("main member ", sampleSols[hi])
		print ("check current ", model.f(sampleSols[hi]))
		print ("number of other members ", len(rotatedSols[hi]))
		print ("other member rotationIndices: ")
		for mi in range(len(rotatedSols[hi])):
			print (rotatedSols[hi][mi])
		print ("")

	'''for hi in range(len(sampleSols)):
		if len(rotatedSols[hi]) > lenV - 1 or (len(rotatedSols[hi]) >= 1 and rotatedSols[hi][0] == 0):
			print ("problem equivalence class# ", hi)
			print ("main member ", sampleSols[hi])
			print ("num other Solutions ", len(rotatedSols[hi]))

	print ("")
	print ("numSolutions, ", len(allHypers))
	print ("num stable solutions ", len(stableSols))'''
	'''for si in range(len(stableSols)):
		print stableSols[si]'''
	#print "num unstable solutions ", len(unstableSols)
	'''for si in range(len(unstableSols)):
		print unstableSols[si]'''

'''
Return true if stable and false otherwise
'''
def determineStability(equilibrium,model):
	jac = model.jacobian(equilibrium)
	#print ("equilibrium", equilibrium)
	#print ("jac", jac)
	eigVals,_ = np.linalg.eig(jac)
	maxEig = np.amax(eigVals.real)
	if maxEig > 0:
		return False
	return True


'''
Divide solutions according to equivalence classes.
Each equivalence class consists of a set of solution 
that are array rotations of each other.
Equivalence classes are represented by a sample solution
and a list of rotation indices. The list of rotation indices
hold the indices over which the sample solution must be rotated
to get the rest of the solutions in the equivalence class.
This equivalence class is more to help debug the rambus oscillator
problem.
Categorize solutions as stable or unstable
'''
def categorizeSolutions(allHypers,model):
	sampleSols = []
	rotatedSols = {}
	stableSols = []
	unstableSols = []
	allSols = []
	for hyper in allHypers:
		exampleSoln = (hyper[:,0] + hyper[:,1])/2.0
		lenV = len(exampleSoln)
		finalSoln = intervalUtils.newton(model,exampleSoln)
		#print "exampleSoln ", exampleSoln
		#print "finalSoln ", finalSoln
		stable = determineStability(finalSoln[1],model)
		if stable:
			stableSols.append(finalSoln[1])
		else:
			unstableSols.append(finalSoln[1])
		allSols.append(finalSoln[1])
		
		# Classify the solutions into equivalence classes
		if len(sampleSols) == 0:
			sampleSols.append(finalSoln[1])
			rotatedSols[0] = []
		else:
			foundSample = False
			for si in range(len(sampleSols)):
				sample = sampleSols[si]
				for ii in range(lenV):
					if abs(finalSoln[1][0] - sample[ii]) < 1e-8:
						rotatedSample = np.zeros_like(finalSoln[1])
						for ri in range(lenV):
							rotatedSample[ri] = sample[(ii+ri)%lenV]
						if np.less_equal(np.absolute(rotatedSample - finalSoln[1]), np.ones((lenV))*1e-8 ).all():
							foundSample = True
							rotatedSols[si].append(ii)
							break
				if foundSample:
					break

			if foundSample == False:
				sampleSols.append(finalSoln[1])
				rotatedSols[len(sampleSols)-1] = []


	return sampleSols, rotatedSols, stableSols, unstableSols

