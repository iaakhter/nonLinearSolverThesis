import sys
import csv
import multiprocessing

sys.path.append('../')
from prototype import *

def runRambusExperiment(timingFilename, modelType, numStages, gcc):
	useLp = False
	numTrials = 1

	headers = ["Problem", "bisectType", "kAlpha", "NumSolutions", "numBisection", "numLP", "numK", "numSingleKill", "numDoubleKill", "totalLPTime", "totalKTime", "avgLPTime", "avgKTime"]
	for i in range(numTrials):
		headers.append("Run_"+str(i))
	headers.append("UseLp")
	csvTimeData = [headers]

	#kAlphaArray = [0.2, 0.4, 0.6, 0.8, 1.0]
	kAlphaArray = [1.0]
	#bisectionTypes = ["bisectMax", "bisectNewton"]
	bisectionTypes = ["bisectNewton"]

	for bisectType in bisectionTypes:
		for kAlpha in kAlphaArray:	
			print ("rambus modelType", modelType, "numStages", numStages, "gcc", gcc)
			print ("bisectType", bisectType, "kAlpha", kAlpha)
			problemName = "rambus_"+modelType+"_stgs_"+str(numStages)+"_gcc_"+str(gcc)

			timeData = [problemName]
			if not(useLp):
				problemName += "_useLP_false"

			for n in range(numTrials):
				start = time.time()
				statVars = {}
				allHypers = rambusOscillator(modelType=modelType, numStages=numStages, g_cc=gcc, statVars=statVars, kAlpha=kAlpha, bisectType=bisectType, numSolutions="all", useLp=useLp)
				end = time.time()
				numSolutions = len(allHypers)
				timeTaken = end - start
				print ("run ", n, "timeTaken", timeTaken)

				if n == 0:
					timeData.append(bisectType)
					timeData.append(kAlpha)
					timeData.append(numSolutions)
					timeData.append(statVars["numBisection"])
					timeData.append(statVars["numLp"])
					timeData.append(statVars["numK"])
					timeData.append(statVars["numSingleKill"])
					timeData.append(statVars["numDoubleKill"])
					timeData.append(statVars["totalLPTime"])
					timeData.append(statVars["totalKTime"])
					timeData.append(statVars["avgLPTime"])
					timeData.append(statVars["avgKTime"])

				timeData.append(timeTaken)

			if not(useLp):
				timeData.append("False")
			else:
				timeData.append("True")

			csvTimeData.append(timeData)
					
	timeFile = open(timingFilename, 'a')
	with timeFile:
		writer = csv.writer(timeFile)
		writer.writerows(csvTimeData)

def runSchmittExperiment(timingFilename, modelType, inputVoltage, bisectType):
	numTrials = 1
	useLp = False
	
	headers = ["Problem", "bisectType","NumSolutions", "numBisection", "numLP", "numK", "numSingleKill", "numDoubleKill", "totalLPTime", "totalKTime", "avgLPTime", "avgKTime"]
	for i in range(numTrials):
		headers.append("Run_"+str(i))
	headers.append("UseLp")
	csvTimeData = [headers]

	print ("schmitt modelType", modelType, "inputVoltage", inputVoltage)
	problemName = "schmitt_"+modelType+"_input_voltage_"+str(inputVoltage)
	timeData = [problemName]

	for n in range(numTrials):
		statVars = {}
		start = time.time()
		allHypers = schmittTrigger(modelType=modelType, inputVoltage = inputVoltage, bisectType=bisectType, statVars=statVars, numSolutions = "all", useLp = useLp)
		end = time.time()
		numSolutions = len(allHypers)
		timeTaken = end - start
		print ("run ", n, "timeTaken", timeTaken)

		if n == 0:
			timeData.append(bisectType)
			timeData.append(numSolutions)
			timeData.append(statVars["numBisection"])
			timeData.append(statVars["numLp"])
			timeData.append(statVars["numK"])
			timeData.append(statVars["numSingleKill"])
			timeData.append(statVars["numDoubleKill"])
			timeData.append(statVars["totalLPTime"])
			timeData.append(statVars["totalKTime"])
			timeData.append(statVars["avgLPTime"])
			timeData.append(statVars["avgKTime"])

		timeData.append(timeTaken)
	
	if not(useLp):
		timeData.append("False")
	else:
		timeData.append("True")
	csvTimeData.append(timeData)
					
	timeFile = open(timingFilename, 'a')
	with timeFile:
		writer = csv.writer(timeFile)
		writer.writerows(csvTimeData)

def runInverterLoopExperiment(timingFilename, modelType, numInverters, bisectType):
	numTrials = 1
	useLp = False
	
	headers = ["Problem", "bisectType","NumSolutions", "numBisection", "numLP", "numK", "numSingleKill", "numDoubleKill", "totalLPTime", "totalKTime", "avgLPTime", "avgKTime"]
	for i in range(numTrials):
		headers.append("Run_"+str(i))
	headers.append("UseLp")
	csvTimeData = [headers]

	print ("inverterLoop modelType", modelType, "numInverters", numInverters)
	problemName = "inverterLoop_"+modelType+"_numInverters_"+str(numInverters)
	timeData = [problemName]

	for n in range(numTrials):
		statVars = {}
		start = time.time()
		allHypers = inverterLoop(modelType=modelType, numInverters=numInverters, bisectType=bisectType, statVars=statVars, numSolutions="all" , useLp=useLp)
		end = time.time()
		numSolutions = len(allHypers)
		timeTaken = end - start
		print ("run ", n, "timeTaken", timeTaken)

		if n == 0:
			timeData.append(bisectType)
			timeData.append(numSolutions)
			timeData.append(statVars["numBisection"])
			timeData.append(statVars["numLp"])
			timeData.append(statVars["numK"])
			timeData.append(statVars["numSingleKill"])
			timeData.append(statVars["numDoubleKill"])
			timeData.append(statVars["totalLPTime"])
			timeData.append(statVars["totalKTime"])
			timeData.append(statVars["avgLPTime"])
			timeData.append(statVars["avgKTime"])

		timeData.append(timeTaken)
	
	if not(useLp):
		timeData.append("False")
	else:
		timeData.append("True")
	csvTimeData.append(timeData)
					
	timeFile = open(timingFilename, 'a')
	with timeFile:
		writer = csv.writer(timeFile)
		writer.writerows(csvTimeData)


def runInverterExperiment(timingFilename, modelType, inputVoltage, bisectType):
	numTrials = 1
	useLp = False
	
	headers = ["Problem", "bisectType", "NumSolutions", "numBisection", "numLP", "numK", "numSingleKill", "numDoubleKill", "totalLPTime", "totalKTime", "avgLPTime", "avgKTime"]
	for i in range(numTrials):
		headers.append("Run_"+str(i))
	headers.append("UseLp")
	csvTimeData = [headers]

	print ("inverter modelType", modelType, "inputVoltage", inputVoltage)
	problemName = "inverter_"+modelType+"_input_voltage_"+str(inputVoltage)
	timeData = [problemName]

	for n in range(numTrials):
		statVars = {}
		start = time.time()
		allHypers = inverter(modelType=modelType, inputVoltage=inputVoltage, bisectType=bisectType, statVars=statVars, numSolutions="all" , useLp=useLp)
		end = time.time()
		numSolutions = len(allHypers)
		timeTaken = end - start
		print ("run ", n, "timeTaken", timeTaken)

		if n == 0:
			timeData.append(bisectType)
			timeData.append(numSolutions)
			timeData.append(statVars["numBisection"])
			timeData.append(statVars["numLp"])
			timeData.append(statVars["numK"])
			timeData.append(statVars["numSingleKill"])
			timeData.append(statVars["numDoubleKill"])
			timeData.append(statVars["totalLPTime"])
			timeData.append(statVars["totalKTime"])
			timeData.append(statVars["avgLPTime"])
			timeData.append(statVars["avgKTime"])

		timeData.append(timeTaken)
	
	if not(useLp):
		timeData.append("False")
	else:
		timeData.append("True")
	csvTimeData.append(timeData)
					
	timeFile = open(timingFilename, 'a')
	with timeFile:
		writer = csv.writer(timeFile)
		writer.writerows(csvTimeData)



if __name__ == '__main__':
	timeout = 36000 # in seconds (10 hours)
	timingFilename = "../data/solver_time_main.csv"
	#timingFilename = "../data/solver_time_mainLp.csv"
	#timingFilename = "../data/solver_time_bisectType.csv"
	#timingFilename = "../data/solver_time_largeEpsilon.csv"
	#timingFilename = "../data/solver_time_kAlpha.csv"

	# Run rambus experiments
	modelTypesL = ["tanh"]
	#modelTypesL = ["tanh"]
	numStagesL = [6]
	#numStagesL = [2]
	gccL = [4.0]
	for modelType in modelTypesL:
		for numStages in numStagesL:
			for gcc in gccL:	
				p = multiprocessing.Process(target=runRambusExperiment, name="RunRambusExperiment", args=(timingFilename, modelType, numStages, gcc))
				
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
	'''modelTypesL = ["lcMosfet", "scMosfet"]
	inputVoltages = []
	for modelType in modelTypesL:
		if modelType == "lcMosfet":
			inputVoltages = [0.0, 0.9, 1.8]
			bisectType = "bisectNewton"
		elif modelType == "scMosfet":
			inputVoltages = [0.0, 0.5, 1.0]
			bisectType = "bisectMax"
		for inputVoltage in inputVoltages:
			p = multiprocessing.Process(target=runSchmittExperiment, name="RunSchmittExperiment", args=(timingFilename, modelType, inputVoltage, bisectType))
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
	'''modelTypesL = ["tanh", "lcMosfet", "scMosfet"]
	inputVoltages = []
	for modelType in modelTypesL:
		if modelType == "tanh":
			inputVoltages = [-1.0, 1.0]
			bisectType = "bisectNewton"
		elif modelType == "lcMosfet":
			inputVoltages = [0.0, 1.8]
			bisectType = "bisectNewton"
		elif modelType == "scMosfet":
			inputVoltages = [0.0, 1.0]
			bisectType = "bisectMax"
		for inputVoltage in inputVoltages:
			p = multiprocessing.Process(target=runInverterExperiment, name="RunInverterExperiment", args=(timingFilename, modelType, inputVoltage, bisectType))
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
	'''modelTypesL = ["tanh", "lcMosfet", "scMosfet"]
	numInvertersL = [1, 2, 3, 4]
	for modelType in modelTypesL:
		if modelType == "tanh":
			bisectType = "bisectNewton"
		elif modelType == "lcMosfet":
			bisectType = "bisectNewton"
		elif modelType == "scMosfet":
			bisectType = "bisectMax"
		for numInverters in numInvertersL:
			p = multiprocessing.Process(target=runInverterLoopExperiment, name="RunInverterLoopExperiment", args=(timingFilename, modelType, numInverters, bisectType))
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


	


