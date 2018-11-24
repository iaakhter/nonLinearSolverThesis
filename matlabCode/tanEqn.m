function I = tanEqn(point, derIn, derOut, Vin, Vout)
	I = derIn*(Vin - point(1)) + derOut*(Vout - point(2)) + point(3);
end