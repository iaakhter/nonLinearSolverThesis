function vDot = schmittDot(t, V, timeScale)
	if(nargin < 3) timeScale = 1.0; end
	V;
	Vtn = 0.4;
	Vtp = -0.4;
	KSn = 270*1e-6 * 3.0;
	KSp = -90*1e-6 * 3*2.0;
	%Cg = 180e-9*1800e-9*(KSn-KSp);
	Cg = 1e-3;
	CI = 2*1e-4;
	gnd = 0.0;
	Vdd = 1.8;
	vDot = zeros(3,1);
	tIs = zeros(6,1);
	tIs(1) = lcMosfet(gnd, V(4), V(2), 1, Vtn, KSn);
	tIs(2) = lcMosfet(V(2), V(4), V(1), 1, Vtn, KSn);
	tIs(3) = lcMosfet(V(2), V(1), Vdd, 1, Vtn, KSn);
	tIs(4) = lcMosfet(Vdd, V(4), V(3), -1, Vtp, KSp);
	tIs(5) = lcMosfet(V(3), V(4), V(1), -1, Vtp, KSp);
	tIs(6) = lcMosfet(V(3), V(1), gnd, -1, Vtp, KSp);
	tIs;
	vDot(1) = (-tIs(5) -tIs(2))*(1/Cg);
	vDot(2) = (-tIs(1) + tIs(2) + tIs(3))*(1/CI);
	vDot(3) = (-tIs(4) + tIs(6) + tIs(5))*(1/CI);

	vDot = vDot*timeScale;

	%vDot = (1/Cg)*vDot;