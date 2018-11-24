# @author Itrat Ahmed Akhter
# Find dc equilibrium points for circuits involving short channel Mosfet models
# using dReal taken from https://nanohub.org/publications/15/4


from dreal.symbolic import Variable, logical_and, logical_or
from dreal.symbolic import abs as dabs
from dreal.symbolic import sqrt, exp, log
from dreal.symbolic import logical_not
from dreal.symbolic import forall
from dreal.api import CheckSatisfiability, Minimize
import time
import numpy as np

def model_params(mType):
		input_params = {}
		if mType == 'n':
			input_params = {'version':1.10, 'mType':1,
							'W':450e-7, 'Lgdr':45e-7,
							'dLg':3.75e-7, 'Cg':2.837e-6,
							'etov':1.67e-7, 'delta':0.1370, 
							'n0':1.4075, 'Rs0':100.5558,
							'Rd0':100.5558, 'Cif':1.9856e-12,
							'Cof':1.9856e-12, 'vxo':1.2277,
							'mu':641.7222, 'beta':1.8,
							'Tjun':298, 'phib':1.2,
							'gamma':0.2, 'Vt0':0.5496,
							'alpha': 3.5, 'mc':0.2,
							'CTM_select':1, 'CC':0,
							'nd':3.053e-14, 'zeta':1.0}

		elif mType == 'p':
			input_params = {'version':1.10, 'mType':-1,
							'W':450e-7, 'Lgdr':45e-7,
							'dLg':3.75e-7, 'Cg':2.837e-6,
							'etov':1.67e-7, 'delta':0.1637,
							'n0':1.4001, 'Rs0':153.8674,
							'Rd0':153.8674, 'Cif':1.9856e-12,
							'Cof':1.9856e-12, 'vxo':1.0751,
							'mu':285.6056, 'beta':1.6,
							'Tjun':298, 'phib':1.2,
							'gamma':0.2, 'Vt0':0.6252,
							'alpha':3.5, 'mc':0.2,
							'CTM_select':1, 'CC':0,
							'nd':0.0235, 'zeta':1.0}

		else:
			input_params = {'TNOM': 300, 'Wint':5e-7,
							  'Lgdr':45e-7, 'dLg':3.75e-7,
							  'Xl':-20e-9, 'Xw':0,
							  'Xj':1.4e-6, 'CjdTnom':0.0005,
							  'Tjun':298, 'Tcj':0.001,
							  'PbdTnom':1.0, 'Tpb':0.005,
							  'CjswdTnom':5e-10, 'Tcjsw':0.001,
							  'PbswdTnom':1.0, 'Tpbsw':0.005,
							  'CjswgdTnom':5e-10, 'Tcjswg':0.001,
							  'PbswgdTnom':1.0, 'Tpbswg':0.005,
							  'Mjd':0.5, 'Mjswd':0.33,
							  'Mjswgd':0.33, 'W':450e-7}

		return input_params

# dReal implemented short channel formula
# any of Vs, Vg, Vd can be dReal Variables
# I is a variable
# returns constraints relating I to Vs, Vg, Vd and
# other intermediate variables
def mvs_id(fetType, Vs, Vg, Vd, I, shape):
	params = model_params(fetType)
	version = params['version'];
	mType = params['mType'];
	W = params['W'];
	Lgdr = params['Lgdr'];
	dLg = params['dLg'];
	Cg = params['Cg'];
	etov = params['etov'];
	delta = params['delta'];
	n0 = params['n0'];
	Rs0 = params['Rs0'];
	Rd0 = params['Rd0'];
	Cif = params['Cif'];
	Cof = params['Cof'];
	vxo = params['vxo']*1e7;
	mu = params['mu'];
	beta = params['beta'];
	Tjun = params['Tjun'];
	phib = params['phib'];
	gamma = params['gamma'];
	Vt0 = params['Vt0'];
	alpha = params['alpha'];
	mc = params['mc'];
	CTM_select = params['CTM_select'];
	CC = params['CC'];
	nd = params['nd'];
	zeta = params['zeta'];

	# SMALL_VALUE
	SMALL_VALUE = 1e-10;
	# LARGE_VALUE
	LARGE_VALUE = 40;
	
	if mType == 1.0:
		Vb = 0.0
		Vdsi = Variable("Vdsin")
		Vgsi = Variable("Vgsin")
		Vbsi = Variable("Vbsin")
		n = Variable("nn")
		nphit = Variable("nphitn")
		phibVbs = Variable("phibVbsn")
		Vtpcorr = Variable("Vtpcorrn")
		eVgpre = Variable("eVgpren")
		FFpre = Variable("FFpren")
		ab = Variable("abn")
		Vcorr = Variable("Vcorrn")
		Vgscorr = Variable("Vgscorrn")
		Vbscorr = Variable("Vbscorrn")
		phibVbscorr = Variable("phibVbscorrn")
		Vt0bs = Variable("Vt0bsn")
		phibVbsi = Variable("phibVbsin")
		Vt0bs0 = Variable("Vt0bs0n")
		Vtp = Variable("Vtpn")
		Vtp0 = Variable("Vtp0n")
		eVg = Variable("eVgn")
		FF = Variable("FFn")
		eVg0 = Variable("eVg0n")
		FF0 = Variable("FF0n")
		Qref = Variable("Qrefn")
		eta = Variable("etan")
		Qinv_corr = Variable("Qinv_corrn")
		Vdsat = Variable("Vdsatn")
		VdsiVdsat = Variable("VdsiVdsatn")
		powVal = Variable("powValn")
		Fsat = Variable("Fsatn")

	else:
		Vb = 1.8
		Vdsi = Variable("Vdsip")
		Vgsi = Variable("Vgsip")
		Vbsi = Variable("Vbsip")
		n = Variable("np")
		nphit = Variable("nphitp")
		phibVbs = Variable("phibVbsp")
		Vtpcorr = Variable("Vtpcorrp")
		eVgpre = Variable("eVgprep")
		FFpre = Variable("FFprep")
		ab = Variable("abp")
		Vcorr = Variable("Vcorrp")
		Vgscorr = Variable("Vgscorrp")
		Vbscorr = Variable("Vbscorrp")
		phibVbscorr = Variable("phibVbscorrp")
		Vt0bs = Variable("Vt0bsp")
		phibVbsi = Variable("phibVbsip")
		Vt0bs0 = Variable("Vt0bs0p")
		Vtp = Variable("Vtpp")
		Vtp0 = Variable("Vtp0p")
		eVg = Variable("eVgp")
		FF = Variable("FFp")
		eVg0 = Variable("eVg0p")
		FF0 = Variable("FF0p")
		Qref = Variable("Qrefp")
		eta = Variable("etap")
		Qinv_corr = Variable("Qinv_corrp")
		Vdsat = Variable("Vdsatp")
		VdsiVdsat = Variable("VdsiVdsatp")
		powVal = Variable("powValp")
		Fsat = Variable("Fsatp")

	constraints = []

	if mType == 1:
		constraints.append(logical_or(logical_and(Vs <= Vd, Vdsi == mType*(Vd - Vs)),
										logical_and(Vs > Vd, Vdsi == mType*(Vs - Vd))))
		constraints.append(logical_or(logical_and(Vs <= Vd, Vgsi == mType*(Vg - Vs)),
										logical_and(Vs > Vd, Vgsi == mType*(Vg - Vd))))
		constraints.append(logical_or(logical_and(Vs <= Vd, Vbsi == mType*(Vb - Vs)),
										logical_and(Vs > Vd, Vbsi == mType*(Vb - Vd))))

	else:
		constraints.append(logical_or(logical_and(Vd <= Vs, Vdsi == mType*(Vd - Vs)),
										logical_and(Vd > Vs, Vdsi == mType*(Vs - Vd))))
		constraints.append(logical_or(logical_and(Vd <= Vs, Vgsi == mType*(Vg - Vs)),
										logical_and(Vd > Vs, Vgsi == mType*(Vg - Vd))))
		constraints.append(logical_or(logical_and(Vd <= Vs, Vbsi == mType*(Vb - Vs)),
										logical_and(Vd > Vs, Vbsi == mType*(Vb - Vd))))

	Cofs = 0*( 0.345e-12/ etov ) * dLg/ 2.0 + Cof;  # s-terminal outer fringing cap [F/cm]
	Cofd = 0*( 0.345e-12/ etov ) * dLg/ 2.0 + Cof;  # d-terminal outer fringing cap [F/cm]
	Leff = Lgdr - dLg;                            # Effective channel length [cm]. After subtracting overlap lengths on s and d side
	kB = 8.617e-5;                                # Boltzmann constant [eV/K]
	phit = kB*Tjun;                               # Thermal voltage, kT/q [V]
	me = (9.1e-31) * mc;                          # Carrier mass [Kg]

	constraints.append(n == n0 + nd *Vdsi)
	constraints.append(nphit == n*phit)
	aphit = alpha*phit

	constraints.append(phibVbs == dabs(phib - Vbsi))
	constraints.append(Vtpcorr == Vt0 + gamma * (sqrt(phibVbs)- sqrt(phib))- Vdsi * delta)
	constraints.append(eVgpre == exp(( Vgsi - Vtpcorr )/ ( aphit * 1.5 )))
	constraints.append(FFpre == 1.0/( 1.0 + eVgpre ))
	constraints.append(ab == 2 * ( 1 - 0.99 * FFpre ) * phit)
	constraints.append(Vcorr == ( 1.0 + 2.0 * delta ) * ( ab/ 2.0 ) * ( exp( -Vdsi/ ab )))
	constraints.append(Vgscorr == Vgsi + Vcorr)
	constraints.append(Vbscorr == Vbsi + Vcorr)
	constraints.append(phibVbscorr == dabs(phib - Vbscorr))
	constraints.append(Vt0bs == Vt0 + gamma * (sqrt(phibVbscorr) - sqrt( phib )))
	constraints.append(phibVbsi == dabs(phib - Vbsi))
	constraints.append(Vt0bs0 == Vt0 + gamma * (sqrt( phibVbsi) - sqrt( phib )))
	constraints.append(Vtp == Vt0bs - Vdsi * delta - 0.5 * aphit)
	constraints.append(Vtp0 == Vt0bs0 - Vdsi * delta - 0.5 * aphit)
	constraints.append(eVg == exp(( Vgscorr - Vtp )/ ( aphit )))
	constraints.append(FF == 1.0/ ( 1.0 + eVg ))
	constraints.append(eVg0 == exp(( Vgsi - Vtp0 )/ ( aphit )))
	constraints.append(Qref == Cg * nphit)
	constraints.append(eta == ( Vgscorr - ( Vt0bs - Vdsi * delta - FF * aphit ))/ ( nphit ))
	constraints.append(Qinv_corr == Qref*log(1.0 + exp(eta)))
	vx0 = vxo
	Vdsats = vx0*Leff/mu
	constraints.append(Vdsat == Vdsats * ( 1.0 - FF ) + phit * FF)
	constraints.append(VdsiVdsat ==  Vdsi/Vdsat)
	constraints.append(powVal == pow(VdsiVdsat, beta))
	constraints.append(Fsat==VdsiVdsat/(pow((1+powVal), (1.0/beta))))
	if mType == 1:
		constraints.append(logical_or(logical_and(Vs <= Vd, I == Qinv_corr * vx0 * Fsat * W * mType * shape),
										logical_and(Vs > Vd, I == Qinv_corr * vx0 * Fsat * W * mType * -1.0 * shape)))
	else:
		constraints.append(logical_or(logical_and(Vd <= Vs, I == Qinv_corr * vx0 * Fsat * W * mType * shape),
										logical_and(Vd > Vs, I == Qinv_corr * vx0 * Fsat * W * mType * -1.0 * shape)))

	return constraints

# Try and find dc equilibrium points for inverter with short channel mosfet model
# for a specific input voltage
# numSolutions indicates the number of solutions we want dReal to find
def inverterScMosfet(inputVoltage, numSolutions = "all"):
	epsilon = 1e-14
	start = time.time()
	Vdd = 1.0
	sn = 3
	sp = 2*sn

	outputVolt = Variable("outputVolt")
	iP = Variable("iP")
	iN = Variable("iN")

	allConstraints = []	
	allConstraints.append(outputVolt >= 0.0)
	allConstraints.append(outputVolt <= Vdd)
	allConstraints.append(-iP-iN == 0)
	allConstraints += mvs_id("n", 0.0, inputVoltage, outputVolt, iN, sn)
	allConstraints += mvs_id("p", Vdd, inputVoltage, outputVolt, iP, sp)
	

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
		hyper[0,:] = [result[outputVolt].lb() - 2*epsilon, result[outputVolt].ub() + 2*epsilon]

		#print ("hyper", hyper)
		allSolutions.append(hyper)

		print ("num solutions found", len(allSolutions))


	end = time.time()
	print ("time taken", end - start)
	return allSolutions

def inverterLoopScMosfet(numInverters, numSolutions = "all"):
	epsilon = 1e-14
	start = time.time()
	#print ("Vtp", Vtp, "Vtn", Vtn, "Vdd", Vdd, "Kn", Kn, "Kp", Kp, "Sn", Sn)
	Vdd = 1.0
	sn = 3
	sp = 2*sn

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
		allConstraints += mvs_id("n", 0.0, vs[inputInd], vs[outputInd], iNs[i], sn)
		allConstraints += mvs_id("p", Vdd, vs[inputInd], vs[outputInd], iPs[i], sp)

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
			hyper[i,:] = [result[vs[i]].lb() - 2*epsilon, result[vs[i]].ub() + 2*epsilon]

		print ("hyper", hyper)
		allSolutions.append(hyper)

		print ("num solutions found", len(allSolutions))


	end = time.time()
	print ("time taken", end - start)
	return allSolutions


# Try and find dc equilibrium points for rambus oscillator with short channel mosfet model
# numStages indicates the number of stages in the rambus oscillator
# g_cc is the strength of the cross coupled inverter as compared to that of forward (g_fwd = 1.0)
# numSolutions indicates the number of solutions we want dReal to find
def rambusOscillatorScMosfet(numStages, g_cc = 0.5, numSolutions = "all"):
	epsilon = 1e-14
	start = time.time()
	g_fwd = 1.0
	lenV = numStages*2
	Vdd = 1.0
	sn = 3
	sp = 2*sn

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
		fwdConstraints = mvs_id("n", 0.0, vs[fwdInd], vs[i], ifwdNs[i], sn)
		fwdConstraints += mvs_id("p", Vdd, vs[fwdInd], vs[i], ifwdPs[i], sp)
		ccConstraints = mvs_id("n", 0.0, vs[ccInd], vs[i], iccNs[i], sn)
		ccConstraints += mvs_id("p", Vdd, vs[ccInd], vs[i], iccPs[i], sp)
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
			hyper[i,:] = [result[vs[i]].lb() - 2*epsilon, result[vs[i]].ub() + 2*epsilon]

		print ("hyper", hyper)
		allSolutions.append(hyper)

		print ("num solutions found", len(allSolutions))

	end = time.time()
	print ("time taken", end - start)
	return allSolutions

# Try and find dc equilibrium points for schmitt trigger with short channel mosfet model
# for a specific input voltage
# numSolutions indicates the number of solutions we want dReal to find
def schmittTriggerScMosfet(inputVoltage, numSolutions = "all"):
	epsilon = 1e-14
	start = time.time()

	lenV = 3
	Vdd = 1.0
	sn = 3
	sp = 2*sn

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
		allConstraints.append(vs[i] <= Vdd)
	allConstraints += mvs_id('n', 0.0, inputVoltage, vs[1], tIs[0], sn)
	allConstraints += mvs_id('n', vs[1], inputVoltage, vs[0], tIs[1], sn)
	allConstraints += mvs_id('n', vs[1], vs[0], Vdd, tIs[2], sn)
	allConstraints += mvs_id('p', Vdd, inputVoltage, vs[2], tIs[3], sp)
	allConstraints += mvs_id('p', vs[2], inputVoltage, vs[0], tIs[4], sp)
	allConstraints += mvs_id('p', vs[2], vs[0], 0.0, tIs[5], sp)
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
			hyper[i,:] = [result[vs[i]].lb() - 2*epsilon, result[vs[i]].ub() + 2*epsilon]

		#print ("hyper", hyper)
		allSolutions.append(hyper)

		print ("num solutions found", len(allSolutions))
		

	end = time.time()
	print ("time taken", end - start)
	return allSolutions

if __name__ == "__main__":
	allSolutions = inverterLoopScMosfet(1)
	print ("allSolutions")
	for solution in allSolutions:
		print ("solution")
		for i in range(solution.shape[0]):
			print (solution[i,0], solution[i,1])

