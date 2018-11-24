% Itrat Ahmed Akhter
% CPSC 538G Proposal
% oscDot.m

% Helper function to simulate the ring oscillator over time
% using ode solvers
% t - time for simulation
% y - input voltage
% a - input for inverter model
% inverterFunc - inverter model
% g_cc - strength of the cross-coupled inverters
% g_fwd - strength of the forward inverters
% shift - input for inverter model

function vdot = oscDot(t, y,a,inverterFunc,g_cc,shift)
  g_fwd = 1;
  v = reshape(y, [2, size(y,1)/2]);
  I = oscI(v, a, inverterFunc,g_cc,g_fwd,shift);
  vdot = reshape(I, size(y));
% oscDot
