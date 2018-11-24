% Itrat Ahmed Akhter
% CPSC 538G Proposal
% inverterTanh.m

% Hyperbolic tangent function model for inverter
% in - input voltage
% n - the number of cycles through an inverter
% a - slope of the middle portion of the function.
% shift - the amount by which the inverter function gets
%         translated by so that nonlinear models can be
%         bound by non-linear model


function out = inverterTanh(in, n, a)
  if(nargin < 3) a = 5; end
  if(nargin < 2) n = 1; end
  out = in;
  smallA = a*0.25;
  for i = 1:n
    out = -tanh(a*out);
  end
end % inverter
