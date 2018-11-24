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


function I = oscI(v, a, inverterFunc, g_cc, g_fwd, shift)
  if(nargin < 6) shift = 0; end;
  if(nargin < 5) g_fwd = 1;   end;
  if(nargin < 4) g_cc  = 0.5; end;
  if(nargin < 3) inverterFunc = @inverter; end;
  if(nargin < 2) a = 5; end;
  I = zeros(size(v)); % the current
  vin = [v(2,end); v(1,end)]; % wrap-around with a twist
  for m = 1:size(v,2) % for each stage
    for n = 1:2 % top and bottom halves
      I_fwd = (inverterFunc(vin(n), 1, a) - v(n,m))*g_fwd;
      I_cc  = (inverterFunc(v(3-n, m), 1, a) - v(n,m))*g_cc; % beware -- tricky indexing with "3-n"
      I(n,m) = I_fwd + I_cc;
    end % for n
    vin = v(:,m);
  end % for m
end % oscI

