% Itrat Ahmed Akhter
% CPSC 538G Proposal
% inverter.m

% 3 piece linear piecewise function model for inverter
% in - input voltage
% n - the number of cycles through an inverter
% a - slope of the middle portion of the function.
% shift - the amount by which the inverter function gets
%         translated by so that nonlinear models can be
%         bound by non-linear model
function out = inverter(in, n, a, shift)
  if(nargin < 4) shift = 0; end
  if(nargin < 3) a = 5; end
  if(nargin < 2) n = 1; end
  out = in;
  for i = 1:n
    out = min(max(-a*(out-shift), -1), +1) + shift;
  end
end % inverter
