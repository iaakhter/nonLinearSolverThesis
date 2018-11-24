from circuit import *
import numpy as np
from cvxopt import matrix,solvers
import random
import funCompUtils as fcUtils

def testLPUnionTanh(low, high):
	function = fcUtils.tanhFun
	functionDer = fcUtils.tanhFunder
	constant = -5
	# test with linear constraints for tanh function
	
	# construct LP for bounds from low to 0
	lowLp = LP()
	inputLow, inputHigh = low, 0.0
	funLow = function(inputLow, constant)
	dLow = functionDer(inputLow, constant)
	funHigh = function(inputHigh, constant)
	dHigh = functionDer(inputHigh, constant)
	
	cLow = funLow - dLow*inputLow
	cHigh = funHigh - dHigh*inputHigh

	diff = inputHigh - inputLow
	if(diff == 0):
		diff = 1e-10
	dThird = (funHigh - funLow)/diff
	cThird = funLow - dThird*inputLow
	# secant
	lowLp.ineq_constraint([-1, dThird], -cThird)
	#tangent constraints
	lowLp.ineq_constraint([1, -dLow], cLow)
	lowLp.ineq_constraint([1, -dHigh], cHigh)
	lowLp.ineq_constraint([0, -1], -inputLow)
	lowLp.ineq_constraint([0, 1], inputHigh)
	lowLp.add_cost([1.0,0.0])
	lowSol = lowLp.solve()
	if lowSol["status"] == "primal_infeasible":
		print ("oops lowLp is infeasible")
		return


	# construct LP for bounds from 0 to high
	highLp = LP()
	inputLow, inputHigh = 0.0, high
	funLow = function(inputLow, constant)
	dLow = functionDer(inputLow, constant)
	funHigh = function(inputHigh, constant)
	dHigh = functionDer(inputHigh, constant)
	
	cLow = funLow - dLow*inputLow
	cHigh = funHigh - dHigh*inputHigh

	diff = inputHigh - inputLow
	if(diff == 0):
		diff = 1e-10
	dThird = (funHigh - funLow)/diff
	cThird = funLow - dThird*inputLow
	# secant
	highLp.ineq_constraint([1, -dThird], cThird)
	#tangent constraints
	highLp.ineq_constraint([-1, dLow], -cLow)
	highLp.ineq_constraint([-1, dHigh], -cHigh)
	highLp.ineq_constraint([0, -1], -inputLow)
	highLp.ineq_constraint([0, 1], inputHigh)
	highLp.add_cost([1.0,0.0])
	highSol = highLp.solve()
	if highSol["status"] == "primal_infeasible":
		print ("oops highLp is infeasible")
		return


	# construct union of lowLp and highLp
	unionLp = highLp.union(lowLp)

	'''print ("lowLp") 
	print (lowLp)
	print ("highLp")
	print (highLp)
	print ("unionLp")
	print (unionLp)'''

	# check if lowLp and unionLp are feasible
	lpCheck1 = LP()
	lpCheck1 = lpCheck1.concat(lowLp)
	lpCheck1 = lpCheck1.concat(unionLp)
	lpCheck1.add_cost([1.0,0])
	lpSol1 = lpCheck1.solve()
	if lpSol1["status"] == "primal infeasible":
		print ("oops lowLp and unionLp together are infeasible")
		print ("lowLp")
		print (lowLp)
		print ("highLp")
		print (highLp)
		print ("unionLp")
		print (unionLp)
		return

	# check if highLp and unionLp are feasible
	lpCheck2 = LP()
	lpCheck2 = lpCheck2.concat(highLp)
	lpCheck2 = lpCheck2.concat(unionLp)
	lpCheck2.add_cost([1.0,0])
	lpSol2 = lpCheck2.solve()
	if lpSol2["status"] == "primal infeasible":
		print ("oops highLp and unionLp together are infeasible")
		print ("lowLp")
		print (lowLp)
		print ("highLp")
		print (highLp)
		print ("unionLp")
		print (unionLp)
		return


def testIdsTransistorHelp(m1, Vs, Vg, Vd, numDivisions = None, mainDict = None):
	#print ("testing", m1.model.channelType, "Vs", Vs, "Vg", Vg, "Vd", Vd)
	if m1.model.channelType == "nfet":
		v = np.array([Vs, Vg, Vd, np.array([1.8])])
	elif m1.model.channelType == "pfet":
		v = np.array([np.array([0.0]), Vg, Vd, Vs])
	IDS = m1.ids(v)

	sampleDelta = 0.0001
	vsSamples = np.linspace(Vs[0]+sampleDelta, Vs[1]-sampleDelta, 10)
	vgSamples = np.linspace(Vg[0]+sampleDelta, Vg[1]-sampleDelta, 10)
	vdSamples = np.linspace(Vd[0]+sampleDelta, Vd[1]-sampleDelta, 10)
	totalSamples = 0
	totalLps = 0
	for vs in vsSamples:
		for vg in vgSamples:
			for vd in vdSamples:
				totalSamples += 1
				if m1.model.channelType == "nfet":
					ids = m1.ids(np.array([vs, vg, vd, 1.8]))
				elif m1.model.channelType == "pfet":
					ids = m1.ids(np.array([0.0, vg, vd, vs]))
		
				if(((ids < IDS[0]) or (ids > IDS[1])) and not (tiny_p(ids - IDS[0]) or tiny_p(ids - IDS[1]))):
					print "oops ids not in interval"
					print "intervalVs", Vs, "intervalVg", Vg, "intervalVd", Vd
					print "sampleVs", vs, "sampleVg", vg, "sampleVd", vd
					print "IDS", IDS
					print "actualIds", ids
					exit()

def testIdsGradTransistorHelp(m1, Vs, Vg, Vd, numDivisions = None, mainDict = None):
	#print ("testing", m1.model.channelType, "Vs", Vs, "Vg", Vg, "Vd", Vd)
	if m1.model.channelType == "nfet":
		v = np.array([Vs, Vg, Vd, np.array([1.8])])
	elif m1.model.channelType == "pfet":
		v = np.array([np.array([0.0]), Vg, Vd, Vs])
	GRAD = m1.grad_ids(v)

	sampleDelta = 0.0001
	vsSamples = np.linspace(Vs[0]+sampleDelta, Vs[1]-sampleDelta, 10)
	vgSamples = np.linspace(Vg[0]+sampleDelta, Vg[1]-sampleDelta, 10)
	vdSamples = np.linspace(Vd[0]+sampleDelta, Vd[1]-sampleDelta, 10)
	totalSamples = 0
	totalLps = 0
	for vs in vsSamples:
		for vg in vgSamples:
			for vd in vdSamples:
				totalSamples += 1
				if m1.model.channelType == "nfet":
					grad = m1.grad_ids(np.array([vs, vg, vd, 1.8]))
				elif m1.model.channelType == "pfet":
					grad = m1.grad_ids(np.array([0.0, vg, vd, vs]))
		
				delV = 0.0;
				if (grad[0] < GRAD[0][0] - delV or grad[0] > GRAD[0][1] + delV or
					grad[1] < GRAD[1][0] - delV or grad[1] > GRAD[1][1] + delV or
					grad[2] < GRAD[2][0] - delV or grad[2] > GRAD[2][1] + delV):
					print "oops grad not in interval"
					print "intervalVs", Vs, "intervalVg", Vg, "intervalVd", Vd
					print "sampleVs", vs, "sampleVg", vg, "sampleVd", vd
					print "GRAD", GRAD
					print "actualgrad", grad
					exit()

def testLpTransistorHelp(m1, Vs, Vg, Vd, numDivisions = None, mainDict = None):
	#print ("testing", m1.model.channelType, "Vs", Vs, "Vg", Vg, "Vd", Vd)
	if m1.model.channelType == "nfet":
		v = np.array([Vs, Vg, Vd, np.array([1.8])])
	elif m1.model.channelType == "pfet":
		v = np.array([np.array([0.0]), Vg, Vd, Vs])
	IDS = m1.ids(v)
	lp = LP()
	lp = lp.concat(m1.lp_ids(v))


	lp.add_cost([0.0, 0.0, 0.0, 1.0])

	minsol = lp.solve()
	#print ("minsol")
	#print (minsol['x'])

	
	if minsol is None:
		#print ("non minsol")
		return

	lp.add_cost([0.0, 0.0, 0.0, -1.0])
	maxsol = lp.solve()
	#print ("maxsol")
	#print (maxsol['x'])
	
	if minsol["status"] == "primal infeasible" or maxsol["status"] == "primal infeasible":
		print ("oops infeasible lp")
		print ("intervalVs", Vs, "intervalVg", Vg, "intervalVd", Vd)
		print ("lp", str(lp))
		return

	#print ("lp")
	#print (lp)

	#print ("minsol['status']", minsol["status"])
	#print ("maxsol['status']", maxsol["status"])

	#print ("maxsol slack")
	#print (lp.slack())

	minIds = minsol['x'][3] - 1e-6
	maxIds = maxsol['x'][3] + 1e-6

	#minIds = minsol['x'][3]
	#maxIds = maxsol['x'][3]



	sampleDelta = 0.0001
	vsSamples = np.linspace(Vs[0]+sampleDelta, Vs[1]-sampleDelta, 10)
	vgSamples = np.linspace(Vg[0]+sampleDelta, Vg[1]-sampleDelta, 10)
	vdSamples = np.linspace(Vd[0]+sampleDelta, Vd[1]-sampleDelta, 10)
	totalSamples = 0
	totalLps = 0
	for vs in vsSamples:
		for vg in vgSamples:
			for vd in vdSamples:
				totalSamples += 1
				if m1.model.channelType == "nfet":
					ids = m1.ids(np.array([vs, vg, vd, 1.8]))
				elif m1.model.channelType == "pfet":
					ids = m1.ids(np.array([0.0, vg, vd, vs]))
		
				if(((ids < minIds) or (ids > maxIds)) and not (tiny_p(ids - minIds) or tiny_p(ids - maxIds))):
					print "oops ids not in lp interval"
					print "intervalVs", Vs, "intervalVg", Vg, "intervalVd", Vd
					print "sampleVs", vs, "sampleVg", vg, "sampleVd", vd
					print "minIds", minIds + 1e-4, "maxIds", maxIds - 1e-4
					#print "IDS", IDS
					print "actualIds", ids
					print "lp", str(lp)
					exit()


def testNfetIds():
	print ("testing nfet ids")
	nfet = MosfetModel(channelType = 'nfet', Vt =0.4, k = 270.0e-6)
	m1 = LcMosfet(0, 1, 2, nfet)
	Vss = list(np.linspace(0.0, 1.5, 10))
	Vgs = list(np.linspace(0.0, 1.5, 10))
	Vds = list(np.linspace(0.0, 1.5, 10))
	for vs in Vss:
		for vg in Vgs:
			for vd in Vds:
				Vs = np.array([vs, vs+0.3])
				Vd = np.array([vd, vd+0.3])
				Vg = np.array([vg, vg+0.3])
				testIdsTransistorHelp(m1, Vs, Vg, Vd)


def testPfetIds():
	print ("testing pfet ids")
	pfet = MosfetModel(channelType = 'pfet', Vt = -0.4, k = -90.0e-6)
	m1 = LcMosfet(3, 1, 2, pfet)

	Vss = list(np.linspace(0.0, 1.5, 10))
	Vgs = list(np.linspace(0.0, 1.5, 10))
	Vds = list(np.linspace(0.0, 1.5, 10))
	for vs in Vss:
		for vg in Vgs:
			for vd in Vds:
				Vs = np.array([vs, vs+0.3])
				Vd = np.array([vd, vd+0.3])
				Vg = np.array([vg, vg+0.3])
				testIdsTransistorHelp(m1, Vs, Vg, Vd)

def testNfetLp():
	print ("testing nfet lp")
	nfet = MosfetModel(channelType = 'nfet', Vt =0.4, k = 270.0e-6)
	m1 = LcMosfet(0, 1, 2, nfet)
	Vss = list(np.linspace(0.0, 1.5, 10))
	Vgs = list(np.linspace(0.0, 1.5, 10))
	Vds = list(np.linspace(0.0, 1.5, 10))
	for vs in Vss:
		for vg in Vgs:
			for vd in Vds:
				Vs = np.array([vs, vs+0.3])
				Vd = np.array([vd, vd+0.3])
				Vg = np.array([vg, vg+0.3])
				testLpTransistorHelp(m1, Vs, Vg, Vd)


def testPfetLp():
	print ("testing pfet lp")
	pfet = MosfetModel(channelType = 'pfet', Vt = -0.4, k = -90.0e-6)
	m1 = LcMosfet(3, 1, 2, pfet)

	Vss = list(np.linspace(0.0, 1.5, 10))
	Vgs = list(np.linspace(0.0, 1.5, 10))
	Vds = list(np.linspace(0.0, 1.5, 10))
	for vs in Vss:
		for vg in Vgs:
			for vd in Vds:
				Vs = np.array([vs, vs+0.3])
				Vd = np.array([vd, vd+0.3])
				Vg = np.array([vg, vg+0.3])
				testLpTransistorHelp(m1, Vs, Vg, Vd)

def testNfetGrads():
	print ("testing nfet grads")
	nfet = MosfetModel(channelType = 'nfet', Vt =0.4, k = 270.0e-6)
	m1 = LcMosfet(0, 1, 2, nfet)
	Vss = list(np.linspace(0.0, 1.5, 10))
	Vgs = list(np.linspace(0.0, 1.5, 10))
	Vds = list(np.linspace(0.0, 1.5, 10))
	for vs in Vss:
		for vg in Vgs:
			for vd in Vds:
				Vs = np.array([vs, vs+0.3])
				Vd = np.array([vd, vd+0.3])
				Vg = np.array([vg, vg+0.3])
				testIdsGradTransistorHelp(m1, Vs, Vg, Vd)


def testPfetGrads():
	print ("testing pfet grads")
	pfet = MosfetModel(channelType = 'pfet', Vt = -0.4, k = -90.0e-6)
	m1 = LcMosfet(3, 1, 2, pfet)

	Vss = list(np.linspace(0.0, 1.5, 10))
	Vgs = list(np.linspace(0.0, 1.5, 10))
	Vds = list(np.linspace(0.0, 1.5, 10))
	for vs in Vss:
		for vg in Vgs:
			for vd in Vds:
				Vs = np.array([vs, vs+0.3])
				Vd = np.array([vd, vd+0.3])
				Vg = np.array([vg, vg+0.3])
				testIdsGradTransistorHelp(m1, Vs, Vg, Vd)

def testCircuit():
	# set up a two stage oscillator and check if the lp constraints
	# original hyper to that containing solution
	g_fwd, g_cc = 1.0, 0.5
	Vtn, Vtp = 0.4, -0.4
	Kn, Kp = 270*1e-6, -90*1e-6
	s0 = 3.0
	grnd, Vdd = 0.0, 1.8

	nfet = MosfetModel('nfet', Vtn, Kn)
	pfet = MosfetModel('pfet', Vtp, Kp)

	numStages = 2
	lenV = numStages*2

	transistorList = []
	for i in range(lenV):
		fwdInd = (i-1)%(lenV)
		ccInd = (i+lenV//2)%(lenV)
		transistorList.append(Mosfet(lenV, fwdInd, i, nfet, s0*g_fwd))
		transistorList.append(Mosfet(lenV+1, fwdInd, i, pfet, 2.0*s0*g_fwd))
		transistorList.append(Mosfet(lenV, ccInd, i, nfet, s0*g_cc))
		transistorList.append(Mosfet(lenV+1, ccInd, i, pfet, 2.0*s0*g_cc))

	c = Circuit(transistorList)

	origHyper = np.array([[1.3, 1.6], [0.1, 0.4], [1.3, 1.6], [0.1, 0.4]])
	#origHyper = np.array([[0.1, 0.4], [1.3, 1.6], [0.1, 0.4], [1.3, 1.6]])
	#origHyper = np.array([[0.7, 0.9], [0.7, 0.9], [0.7, 0.9], [0.7, 0.9]])
	#origHyper = np.array([ 0.84948974,  0.84948974,  0.84948974,  0.84948974])
	cHyper = [x for x in origHyper] + [0.0, Vdd]
	#print ("current in origHyper", c.f(cHyper))
	newHyper = c.linearConstraints(cHyper, [4, 5])
	print (newHyper)

def testCircuitSchmitt():
	Vtn, Vtp = 0.4, -0.4
	Kn, Kp = 270*1e-6, -90*1e-6
	s0 = 3.0
	grnd, Vdd = 0.0, 1.8
	inputVoltage = 1.0


	nfet = MosfetModel('nfet', Vtn, Kn, "default")
	pfet = MosfetModel('pfet', Vtp, Kp, "default")

	# with the voltage array containing [grnd, Vdd, input, X[0], X[1], X[2]]
	# where X[0] is the output voltage and X[1] is the voltage at node with 
	# nfets and X[2] is the voltage at node with pfets

	# src, gate, drain = grnd, input, X[1]
	m0 = Mosfet(0, 2, 4, nfet, s0)

	# src, gate, drain = X[1], input, X[0]
	m1 = Mosfet(4, 2, 3, nfet, s0)

	# src, gate, drain = X[1], X[0], Vdd
	m2 = Mosfet(4, 3, 1, nfet, s0)

	# src, gate, drain = Vdd, input, X[2]
	m3 = Mosfet(1, 2, 5, pfet, s0*2.0)

	# src, gate, drain = X[2], in, X[0]
	m4 =Mosfet(5, 2, 3, pfet, s0*2.0)

	# src, gate, drain = X[2], X[0], grnd
	m5 = Mosfet(5, 3, 0, pfet, s0*2.0)

	c = Circuit([m0, m1, m2, m3, m4, m5])
	origHyper = np.array([  6.99592428e-05,   5.34975542e-05,   8.00042178e-01])
	cHyper = [0.0, Vdd, inputVoltage] + [x for x in origHyper]
	print ("current in origHyper", c.f(cHyper))



if __name__ == "__main__":
	testNfetIds()
	testPfetIds()
	#testNfetGrads()
	#testPfetGrads()
	testNfetLp()
	testPfetLp()
	#testCircuit()
	#testCircuitSchmitt()
	#testLPUnionTanh(-1.0, 1.0)
	#testLPUnionTanh(-1.0, 0.1)
	#testLPUnionTanh(-0.1, 1.0)
	#testLPUnionTanh(-0.1, 0.1)