% Itrat Ahmed Akhter
% CPSC 538G Proposal
% oscl.m

% Rambus ring oscillator model
% v - input voltage array
% a - input for inverter model
% inverterFunc - inverter model
% g_cc - strength of the cross-coupled inverters
% g_fwd - strength of the forward inverters
% shift - input for inverter model


function I = oscITransistor(v, g_cc, g_fwd, Vtp, Vtn, Vdd, Kn, Sn)
  if(nargin < 7) Sn = 8/3.0; end;
  if(nargin < 7) Kn = 1.5; end;
  if(nargin < 6) Vdd = 1.8; end;
  if(nargin < 5) Vtn = 0.4; end;
  if(nargin < 4) Vtp = -0.4; end;
  if(nargin < 3) g_fwd = 1;   end;
  if(nargin < 2) g_cc  = 0.5; end;
  I = zeros(size(v)); % the current
  vin = [v(2,end); v(1,end)]; % wrap-around with a twist
  for m = 1:size(v,2) % for each stage
    for n = 1:2 % top and bottom halves
      [I_fwd, secDerFwdIn, secDerFwdOut] = currentMosfet(vin(n), v(n,m), Vtp, Vtn, Vdd, Kn);
      [I_cc, secDerCcIn, secDerCcOut] = currentMosfet(v(3-n, m), v(n,m), Vtp, Vtn, Vdd, Kn);
      I(n,m) = g_fwd*I_fwd + g_cc*I_cc;
    end % for n
    vin = v(:,m);
  end % for m
end % oscI

