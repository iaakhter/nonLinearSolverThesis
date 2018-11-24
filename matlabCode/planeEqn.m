function I = planeEqn(point1, point2, point3, Vin, Vout)
	normal = cross(point1 - point2, point1 - point3);
	%if normal(3) < 0
	%	normal = -normal;
	%end;
	I = (-normal(1)*(Vin - point1(1)) - normal(2)*(Vout - point1(2)))/normal(3) + point1(3);
end