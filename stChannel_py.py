import stChannel_py
import numpy as np

def compareSample(fetType, VdI, VgI, VsI, vdSam, vgSam, vsSam):
	stMosfet = stChannel_py.StMosfet();
	Vb = stChannel_py.MyList();
	if fetType == "n":
		fetFunc = stMosfet.mvs_idn;
		Vb[:] = [0.0, 0.0];
	else:
		fetFunc = stMosfet.mvs_idp;
		Vb[:] = [1.8, 1.8];

	idI = fetFunc(VdI, VgI, VsI, Vb);
	idSam = fetFunc(vdSam, vgSam, vsSam, Vb);
	if (idSam[0] < idI[0] or idSam[1] > idI[1]):
		print ("vd", vdSam[0], " in ", VdI[0], VdI[1])
		print ("vg", vgSam[0], " in ", VgI[0], VgI[1])
		print ("vs", vsSam[0], " in ", VsI[0], VsI[1])
		print ("gives i ", idSam[0], idSam[1], "not in i interval")
		print (idI[0], idI[1])
		print ("")


def testFetCurrent(fetType):
	stMosfet = stChannel_py.StMosfet();
	vds = list(np.linspace(0.0, 1.5, 10));
	vgs = list(np.linspace(0.0, 1.5, 10));
	vss = list(np.linspace(0.0, 1.5, 10));
	Vb = stChannel_py.MyList();
	if fetType == "n":
		fetFunc = stMosfet.mvs_idnMon
		Vb.append(0.0);
		Vb.append(0.0);
	else:
		fetFunc = stMosfet.mvs_idpMon
		Vb.append(1.8);
		Vb.append(1.8);
	for vd in vds:
		for vg in vgs:
			for vs in vss:
				Vd = stChannel_py.MyList()
				Vg = stChannel_py.MyList()
				Vs = stChannel_py.MyList()
				Vd.append(vd)
				Vd.append(vd + 0.3)
				Vg.append(vg)
				Vg.append(vg + 0.3)
				Vs.append(vs)
				Vs.append(vs + 0.3)

				'''print ("CALCULATING INTERVAL FOR ")
				print ("Vd", Vd[0], Vd[1])
				print ("Vg", Vg[0], Vg[1])
				print ("Vs", Vs[0], Vs[1])'''

				idn = fetFunc(Vd, Vg, Vs, Vb)
				
				epsilon = 0.0
				vdSample = list(np.linspace(Vd[0] + epsilon, Vd[1] - epsilon, 10))
				vgSample = list(np.linspace(Vg[0] + epsilon, Vg[1] - epsilon, 10))
				vsSample = list(np.linspace(Vs[0] + epsilon, Vs[1] - epsilon, 10))

				for vdSam in vdSample:
					for vgSam in vgSample:
						for vsSam in vsSample:
							VdSam = stChannel_py.MyList()
							VgSam = stChannel_py.MyList()
							VsSam = stChannel_py.MyList()
							VdSam.append(vdSam)
							VdSam.append(vdSam)
							VgSam.append(vgSam)
							VgSam.append(vgSam)
							VsSam.append(vsSam)
							VsSam.append(vsSam)
							iVal = fetFunc(VdSam, VgSam, VsSam, Vb)
							epsilon = 0
							if (iVal[0] < idn[0] - epsilon or iVal[1] > idn[1] + epsilon):
								print ("vd", vdSam, " in ", Vd[0], Vd[1])
								print ("vg", vgSam, " in ", Vg[0], Vg[1])
								print ("vs", vsSam, " in ", Vs[0], Vs[1])
								print ("gives i ", iVal[0], iVal[1], "not in i interval")
								print (idn[0], idn[1])
								print ("")
								return


def testFetJac(fetType):
	delta = 0.001;
	Vb = stChannel_py.MyList();
	if fetType == "n":
		fetFun = stMosfet.mvs_idn
		fetJac = stMosfet.mvs_idnGrad
		Vb.append(0.0);
		Vb.append(0.0);
	else:
		fetFun = stMosfet.mvs_idp
		fetJac = stMosfet.mvs_idpGrad
		Vb.append(1.8);
		Vb.append(1.8);

	vds = list(np.linspace(0.0, 1.8 - delta, 100));
	vgs = list(np.linspace(0.0, 1.8 - delta, 100));
	vss = list(np.linspace(0.0, 1.8 - delta, 100));
	for vd in vds:
		for vg in vgs:
			for vs in vss:
				Vd = stChannel_py.MyList()
				Vg = stChannel_py.MyList()
				Vs = stChannel_py.MyList()
				VdDelta = stChannel_py.MyList()
				VgDelta = stChannel_py.MyList()
				VsDelta = stChannel_py.MyList()
				Vd.append(vd)
				Vd.append(vd)
				Vg.append(vg)
				Vg.append(vg)
				Vs.append(vs)
				Vs.append(vs)
				VdDelta.append(vd + delta)
				VdDelta.append(vd + delta)
				VgDelta.append(vg + delta)
				VgDelta.append(vg + delta)
				VsDelta.append(vs + delta)
				VsDelta.append(vs + delta)

				idJac = fetJac(Vd, Vg, Vs, Vb)

				idVd0 = fetFun(Vd, Vg, Vs, Vb);
				idVd1 = fetFun(VdDelta, Vg, Vs, Vb);

				idVg0 = fetFun(Vd, Vg, Vs, Vb);
				idVg1 = fetFun(Vd, VgDelta, Vs, Vb);

				idVs0 = fetFun(Vd, Vg, Vs, Vb);
				idVs1 = fetFun(Vd, Vg, VsDelta, Vb);


				compJac0 = (idVd1[0] - idVd0[0])/delta
				compJac1 = (idVg1[0] - idVg0[0])/delta
				compJac2 = (idVs1[0] - idVs0[0])/delta


				if (abs(idJac[0][0] - compJac0) > 1e-3 or abs(idJac[1][0] - compJac1) > 1e-3 or
					abs(idJac[2][0] - compJac2) > 1e-3):
					print ("oops vd", vd, "vg", vg, "vs", vs)
					print ("jacobian ", idJac[0][0], idJac[1][0], idJac[2][0], "not close to ")
					print ("numeric jac ", compJac0, compJac1, compJac2)
					return;


def testFetJacInterval(fetType):
	vds = list(np.linspace(0.0, 1.5, 10));
	vgs = list(np.linspace(0.0, 1.5, 10));
	vss = list(np.linspace(0.0, 1.5, 10));
	Vb = stChannel_py.MyList();
	if fetType == "n":
		fetFun = stMosfet.mvs_idn
		fetJac = stMosfet.mvs_idnGrad
		Vb.append(0.0);
		Vb.append(0.0);
	else:
		fetFun = stMosfet.mvs_idp
		fetJac = stMosfet.mvs_idpGrad
		Vb.append(1.8);
		Vb.append(1.8);
	
	for vd in vds:
		for vg in vgs:
			for vs in vss:
				Vd = stChannel_py.MyList()
				Vg = stChannel_py.MyList()
				Vs = stChannel_py.MyList()
				Vd.append(vd)
				Vd.append(vd + 0.3)
				Vg.append(vg)
				Vg.append(vg + 0.3)
				Vs.append(vs)
				Vs.append(vs + 0.3)

				jac = fetJac(Vd, Vg, Vs, Vb)
				
				vdSample = list(np.linspace(Vd[0], Vd[1], 10))
				vgSample = list(np.linspace(Vg[0], Vg[1], 10))
				vsSample = list(np.linspace(Vs[0], Vs[1], 10))

				for vdSam in vdSample:
					for vgSam in vgSample:
						for vsSam in vsSample:
							VdSam = stChannel_py.MyList()
							VgSam = stChannel_py.MyList()
							VsSam = stChannel_py.MyList()
							VdSam.append(vdSam)
							VdSam.append(vdSam)
							VgSam.append(vgSam)
							VgSam.append(vgSam)
							VsSam.append(vsSam)
							VsSam.append(vsSam)
							jacSam = fetJac(VdSam, VgSam, VsSam, Vb)

							delV = 0.0;
							if (jacSam[0][0] < jac[0][0] - delV or jacSam[0][1] > jac[0][1] + delV or
								jacSam[1][0] < jac[1][0] - delV or jacSam[1][1] > jac[1][1] + delV or
								jacSam[2][0] < jac[2][0] - delV or jacSam[2][1] > jac[2][1] + delV or
								jacSam[3][0] < jac[3][0] - delV or jacSam[3][1] > jac[3][1] + delV):
								print ("vd", vd, "in", Vd[0], Vd[1])
								print ("vg", vg, "in", Vg[0], Vg[1])
								print ("vs", vs, "in", Vs[0], Vs[1])
								print ("jacSam[0]", jacSam[0][0], jacSam[1][0], "jacSam[1]", jacSam[1][0], jacSam[1][1], 
									"jacSam[2]", jacSam[2][0], jacSam[2][0], "jacSam[3]", jacSam[3][0], jacSam[3][1])
								print ("not in")
								print ("jac[0]", jac[0][0], jac[1][0], "jac[1]", jac[1][0], jac[1][1], 
									"jac[2]", jac[2][0], jac[2][0], "jac[3]", jac[3][0], jac[3][1])
								return;



#ids = StMosfet.mvs_idn([0.7, 0.7], [0.6, 0.6], [0.9, 0.9], [0.0, 0.0])
stMosfet = stChannel_py.StMosfet()
Vd = stChannel_py.MyList()
Vg = stChannel_py.MyList()
Vs = stChannel_py.MyList()
Vb = stChannel_py.MyList()
Vd[:] = [0.94, 0.94]
Vg[:] = [0.04, 0.04]
Vs[:] = [0.86, 0.86]
Vb[:] = [0.0, 0.0]
'''idn = stMosfet.mvs_idnMon(Vd, Vg, Vs, Vb)
print ("idn", idn[0], idn[1])
idnJac = stMosfet.mvs_idnJac(Vd, Vg, Vs, Vb)
print ("idn", idn[0], idn[1])
print ("jac der0", idnJac[0][0], idnJac[0][1],
	"der1", idnJac[1][0], idnJac[1][1],
	"der2", idnJac[2][0], idnJac[2][1],
	"der3", idnJac[3][0], idnJac[3][1])

Vd[:] = [1.799993521250602, 1.799993521250602]
Vg[:] = [0.0, 0.0]
Vs[:] = [1.8, 1.8] 
Vbp = stChannel_py.MyList()
Vbp.append(1.8)
Vbp.append(1.8)
idp = stMosfet.mvs_idp(Vd, Vg, Vs, Vbp, 0.0)
idpJac = stMosfet.mvs_idpJac(Vd, Vg, Vs, Vbp, 0.0)
print ("idp", idp[0], idp[1])
print ("jac der0", idpJac[0][0], idpJac[0][1],
	"der1", idpJac[1][0], idpJac[1][1],
	"der2", idpJac[2][0], idpJac[2][1],
	"der3", idpJac[3][0], idpJac[3][1])'''

print ("testing nfet")
testFetCurrent('n')
print ("testing pfet")
testFetCurrent('p')
print ("testing nfetJac")
testFetJac('n')
print ("testing pfetJac")
testFetJac('p')
print ("testing nfetJacInterval")
testFetJacInterval('n')
print ("testing pfetJacInterval")
testFetJacInterval('p')

'''VdI, VgI, VsI = stChannel_py.MyList(), stChannel_py.MyList(), stChannel_py.MyList();
VdI[:] = [0.3333333333333333, 0.6333333333333333]
VgI[:] = [0.0, 0.3]
VsI[:] = [1.5, 1.8]

vdSam, vgSam, vsSam = stChannel_py.MyList(), stChannel_py.MyList(), stChannel_py.MyList();
vdSam[:] = [0.33333333433333334, 0.33333333433333334]
vgSam[:] = [0.29999999899999996, 0.29999999899999996]
vsSam[:] = [1.799999999, 1.799999999]
compareSample("n", VdI, VgI, VsI, vdSam, vgSam, vsSam)'''

