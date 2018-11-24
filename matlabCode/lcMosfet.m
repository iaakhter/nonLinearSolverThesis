% when channelType = 1, nfet
% when channelType = -1, pfet
function ids = lcMosfet(Vs, Vg, Vd, channelType, Vt, ks)
	ids = 0.0;
	if (channelType == -1)
		ids = -lcMosfet(-Vs, -Vg, -Vd, 1, -Vt, -ks);
	elseif (Vd < Vs)
		ids = -lcMosfet(Vd, Vg, Vs, channelType, Vt, ks);
	else
		Vgse = (Vg - Vs) - Vt;
		Vds = Vd - Vs;
		if(Vgse < 0)  % cut-off
			ids = 0;
		elseif(Vgse < Vds) % saturation
			ids = (ks/2.0)*Vgse*Vgse;
		else % linear
			ids = ks*(Vgse - Vds/2.0)*Vds;
		end
	end

end