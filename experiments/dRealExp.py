import sys
import multiprocessing

sys.path.append('../')
from prototype import *
sys.path.append('../dRealCode')
from dRealLcMosfet import *
from dRealScMosfet import *
from dRealTanh import *

def runRambusExperiment(modelType, numStages, g_cc):
	print ("rambus modelType", modelType, "numStages", numStages, "g_cc", g_cc)
	if modelType == "tanh":
		allHypers = rambusOscillatorTanh(numStages = numStages, g_cc = g_cc)
	if modelType == "lcMosfet":
		allHyper = rambusOscillatorLcMosfet(numStages = numStages, g_cc = g_cc)
	if modelType == "scMosfet":
		allHypers = rambusOscillatorScMosfet(numStages = numStages, g_cc = g_cc)


def runSchmittExperiment(modelType, inputVoltage):
	print ("schmitt modelType", modelType, "inputVoltage", inputVoltage)
	if modelType == "lcMosfet":
		allHypers = schmittTriggerLcMosfet(inputVoltage = inputVoltage)
	if modelType == "scMosfet":
		allHypers = schmittTriggerScMosfet(inputVoltage = inputVoltage)

def runInverterExperiment(modelType, inputVoltage):
	print ("inverter modelType", modelType, "inputVoltage", inputVoltage)
	if modelType == "tanh":
		allHypers = inverterTanh(inputVoltage = inputVoltage)
	if modelType == "lcMosfet":
		allHypers = inverterLcMosfet(inputVoltage = inputVoltage)
	if modelType == "scMosfet":
		allHypers = inverterScMosfet(inputVoltage = inputVoltage)

def runInverterLoopExperiment(modelType, numInverters):
	print ("inverterLoop modelType", modelType, "numInverters", numInverters)
	if modelType == "tanh":
		allHypers = inverterLoopTanh(numInverters = numInverters)
	if modelType == "lcMosfet":
		allHypers = inverterLoopLcMosfet(numInverters = numInverters)
	if modelType == "scMosfet":
		allHypers = inverterLoopScMosfet(numInverters = numInverters)



if __name__ == '__main__':
	timeout = 36000 # in seconds (10 hours)

	# Run rambus experiments
	'''modelTypesL = ["tanh", "lcMosfet", "scMosfet"]
	#modelTypesL = ["scMosfet"]
	numStagesL = [2, 4, 6]
	#numStagesL = [2]
	gccL = [0.5, 4.0]
	#gccL = [4.0]
	for modelType in modelTypesL:
		for numStages in numStagesL:
			for gcc in gccL:	
				p = multiprocessing.Process(target=runRambusExperiment, name="RunRambusExperiment", args=(modelType, numStages, gcc))
				
				p.start()

				# Wait a maximum of timeout seconds for process
				# Usage: join([timeout in seconds])
				p.join(timeout)

				# If thread is active
				if p.is_alive():
					print "After 10 hours... let's kill it..."

					# Terminate process
					p.terminate()
					p.join()

	# Run schmitt trigger experiments
	modelTypesL = ["lcMosfet", "scMosfet"]
	inputVoltages = []
	for modelType in modelTypesL:
		if modelType == "lcMosfet":
			inputVoltages = [0.0, 0.9, 1.8]
		elif modelType == "scMosfet":
			inputVoltages = [0.0, 0.5, 1.0]
		for inputVoltage in inputVoltages:
			p = multiprocessing.Process(target=runSchmittExperiment, name="RunSchmittExperiment", args=(modelType, inputVoltage))
			p.start()

			# Wait a maximum of timeout seconds for process
			# Usage: join([timeout in seconds])
			p.join(timeout)

			# If thread is active
			if p.is_alive():
				print "After 10 hours... let's kill it..."

				# Terminate process
				p.terminate()
				p.join()'''

	# Run inverter experiments
	'''modelTypesL = ["tanh","lcMosfet", "scMosfet"]
	inputVoltages = []
	for modelType in modelTypesL:
		if modelType == "tanh":
			inputVoltages = [-1.0, 1.0]
		elif modelType == "lcMosfet":
			inputVoltages = [0.0, 1.8]
		elif modelType == "scMosfet":
			inputVoltages = [0.0, 1.0]
		for inputVoltage in inputVoltages:
			p = multiprocessing.Process(target=runInverterExperiment, name="RunInverterExperiment", args=(modelType, inputVoltage))
			p.start()

			# Wait a maximum of timeout seconds for process
			# Usage: join([timeout in seconds])
			p.join(timeout)

			# If thread is active
			if p.is_alive():
				print "After 10 hours... let's kill it..."

				# Terminate process
				p.terminate()
				p.join()'''


	# Run inverter loop experiments
	modelTypesL = ["tanh", "lcMosfet", "scMosfet"]
	numInvertersL = [1, 2, 3, 4]
	for modelType in modelTypesL:
		for numInverters in numInvertersL:
			p = multiprocessing.Process(target=runInverterLoopExperiment, name="RunInverterLoopExperiment", args=(modelType, numInverters))
			p.start()


			# Wait a maximum of timeout seconds for process
			# Usage: join([timeout in seconds])
			p.join(timeout)

			# If thread is active
			if p.is_alive():
				print "After 10 hours... let's kill it..."

				# Terminate process
				p.terminate()
				p.join()

