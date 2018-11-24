function Vin = vin(t)

	if (t <= 1.8)
		Vin = t;
	else
		Vin = -t+3.6;
	end
end