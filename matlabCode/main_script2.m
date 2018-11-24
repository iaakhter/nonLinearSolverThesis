currentFun = @(x,y)z_foo(x,y);
%currentFun = @(x,y)z_default(x,y);
hyperIn = [-0.8; 1.2; 1.2; -0.8];
hyperOut = [-0.3; -0.3; 1.0; 1.0];

Vin = min(hyperIn):0.01:max(hyperIn);
Vout = min(hyperOut):0.01:max(hyperOut);

I = zeros(length(Vout), length(Vin));
firDerIn = zeros(length(Vout), length(Vin));
firDerOut = zeros(length(Vout), length(Vin));
secDerIn = zeros(length(Vout), length(Vin));
secDerOut = zeros(length(Vout), length(Vin));
secDerOutIn = zeros(length(Vout), length(Vin));
regTest = zeros(length(Vout), length(Vin));
for i = 1:length(Vin)
	for j = 1:length(Vout)
		[I(j,i), firDerIn(j,i), firDerOut(j,i), secDerIn(j,i), secDerOut(j,i), secDerOutIn(j,i)] = currentFun(Vin(i), Vout(j));
		if secDerIn(j,i) > 0.0 && secDerOut(j,i) > 0.0 && Vin(i) <= 0.75 && Vout(j) >= 0.13761195 && -Vin(i) + Vout(j) <= -0.25
			regTest(j,i) = 1;
		end;
	end;
end;

[posSecDerInR, posSecDerInC] = find(secDerIn > 0);
[negSecDerInR, negSecDerInC] = find(secDerIn < 0);

[posSecDerOutR, posSecDerOutC] = find(secDerOut > 0);
[negSecDerOutR, negSecDerOutC] = find(secDerOut < 0);

[posSecDerOutInR, posSecDerOutInC] = find(secDerOutIn > 0);
[negSecDerOutInR, negSecDerOutInC] = find(secDerOutIn < 0);

CO(:,:,1) = zeros(length(Vout), length(Vin));
CO(:,:,2) = zeros(length(Vout), length(Vin));
CO(:,:,3) = zeros(length(Vout), length(Vin));

for i = 1:length(posSecDerInR)
	%CO(posSecDerInR(i),posSecDerInC(i),1) = 0.0; 
	CO(posSecDerInR(i),posSecDerInC(i),2) = 1.0; 
	%CO(posSecDerInR(i),posSecDerInC(i),3) = 0.0; 
end;

for i = 1:length(negSecDerInR)
	%CO(negSecDerInR(i),negSecDerInC(i),1) = 0.0; 
	CO(negSecDerInR(i),negSecDerInC(i),2) = 1.0; 
	%CO(negSecDerInR(i),negSecDerInC(i),3) = 0.0; 
end;

for i = 1:length(posSecDerOutR)
	%CO(negSecDerOutR(i),negSecDerOutC(i),1) = 0.0; 
	%CO(negSecDerOutR(i),negSecDerOutC(i),2) = 0.0; 
	CO(posSecDerOutR(i),posSecDerOutC(i),3) = 1.0;  
end;

for i = 1:length(negSecDerOutR)
	%CO(negSecDerOutR(i),negSecDerOutC(i),1) = 0.0; 
	%CO(negSecDerOutR(i),negSecDerOutC(i),2) = 0.0; 
	CO(negSecDerOutR(i),negSecDerOutC(i),3) = 1.0;  
end;

currentsAtHypers = zeros(length(hyperIn));
firDerInAtHypers = zeros(length(hyperIn));
firDerOutAtHypers = zeros(length(hyperIn));
for i = 1:length(hyperIn)
	[currentAtHypers(i),firDerInAtHypers(i),firDerOutAtHypers(i),blah3,blah4] = currentFun(hyperIn(i), hyperOut(i));
end;

planeVin = min(hyperIn):0.01:max(hyperIn);
planeVout = min(hyperOut):0.01:max(hyperOut);

% tangent planes
tanPlanePoint1 = [hyperIn(1), hyperOut(1), currentAtHypers(1)];
tanPlanePoint2 = [hyperIn(2), hyperOut(2), currentAtHypers(2)];
tanPlanePoint3 = [hyperIn(3), hyperOut(3), currentAtHypers(3)];
tanPlanePoint4 = [hyperIn(4), hyperOut(4), currentAtHypers(4)];
tanPlaneCurrent1 = zeros(length(planeVout), length(planeVin));
tanPlaneCurrent2 = zeros(length(planeVout), length(planeVin));
tanPlaneCurrent3 = zeros(length(planeVout), length(planeVin));
tanPlaneCurrent4 = zeros(length(planeVout), length(planeVin));

% secant planes
plane1Point1 = [hyperIn(1), hyperOut(1), currentAtHypers(1)];
plane1Point2 = [hyperIn(2), hyperOut(2), currentAtHypers(2)];
plane1Point3 = [hyperIn(3), hyperOut(3), currentAtHypers(3)];
plane1Current = zeros(length(planeVout), length(planeVin));
%plane1Current = zeros(length(Vout), length(Vin));

plane2Point1 = [hyperIn(1), hyperOut(1), currentAtHypers(1)];
plane2Point2 = [hyperIn(2), hyperOut(2), currentAtHypers(2)];
plane2Point3 = [hyperIn(4), hyperOut(4), currentAtHypers(4)];
plane2Current = zeros(length(planeVout), length(planeVin));

plane3Point1 = [hyperIn(2), hyperOut(2), currentAtHypers(2)];
plane3Point2 = [hyperIn(3), hyperOut(3), currentAtHypers(3)];
plane3Point3 = [hyperIn(4), hyperOut(4), currentAtHypers(4)];
plane3Current = zeros(length(planeVout), length(planeVin));

plane4Point1 = [hyperIn(1), hyperOut(1), currentAtHypers(1)];
plane4Point2 = [hyperIn(3), hyperOut(3), currentAtHypers(3)];
plane4Point3 = [hyperIn(4), hyperOut(4), currentAtHypers(4)];
plane4Current = zeros(length(planeVout), length(planeVin));

for i = 1:length(planeVin)
	for j = 1:length(planeVout)
		%plane1Current(j,i) = planeEqn(plane1Point1, plane1Point2, plane1Point3, Vin(i), Vout(j));
		plane1Current(j,i) = planeEqn(plane1Point1, plane1Point2, plane1Point3, planeVin(i), planeVout(j));
		plane2Current(j,i) = planeEqn(plane2Point1, plane2Point2, plane2Point3, planeVin(i), planeVout(j));
		plane3Current(j,i) = planeEqn(plane3Point1, plane3Point2, plane3Point3, planeVin(i), planeVout(j));
		plane4Current(j,i) = planeEqn(plane4Point1, plane4Point2, plane4Point3, planeVin(i), planeVout(j));
		%convexHullPlane1(j,i) = planeEqn(cPlane1Pt1, cPlane1Pt2, cPlane1Pt3, planeVin(i), planeVout(j));
		%convexHullPlane2(j,i) = planeEqn(cPlane2Pt1, cPlane2Pt2, cPlane2Pt3, planeVin(i), planeVout(j));
		tanPlaneCurrent1(j,i) = tanEqn(tanPlanePoint1, firDerInAtHypers(1), firDerOutAtHypers(1), planeVin(i), planeVout(j));
		tanPlaneCurrent2(j,i) = tanEqn(tanPlanePoint2, firDerInAtHypers(2), firDerOutAtHypers(2), planeVin(i), planeVout(j));
		tanPlaneCurrent3(j,i) = tanEqn(tanPlanePoint3, firDerInAtHypers(3), firDerOutAtHypers(3), planeVin(i), planeVout(j));
		tanPlaneCurrent4(j,i) = tanEqn(tanPlanePoint4, firDerInAtHypers(4), firDerOutAtHypers(4), planeVin(i), planeVout(j));
	end;
end;

[surfVin,surfVout] = meshgrid(Vin,Vout);
[surfPlaneVin,surfPlaneVout] = meshgrid(planeVin,planeVout);

I(find(I - (tanPlaneCurrent1 - 1e-12) < 0))
planeCurrent = plane1Current + 1e-15;
tanPlaneCurrent = tanPlaneCurrent4 - 1e-15;

if all(I >= tanPlaneCurrent)
	disp('tanyay')
end;

if all(I <= planeCurrent)
	disp('secyay')
end;

figure;
s1 = surf(surfVin, surfVout, I);
shading interp;
%hold on;
%tanPlane = surf(surfPlaneVin,surfPlaneVout,tanPlaneCurrent);
%hold on;
%tanPlane = surf(surfPlaneVin,surfPlaneVout,tanPlaneCurrent1 - 0.1);
%hold on;
%tanPlane = surf(surfPlaneVin,surfPlaneVout,tanPlaneCurrent3 - 0.1);
%hold on;
%secPlane = surf(surfPlaneVin,surfPlaneVout,plane1Current + 0.1);
%hold on;
%secPlane = surf(surfPlaneVin,surfPlaneVout,plane4Current + 0.1);

