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

function vdot = oscDotTransistor(t,y,g_cc)
  g_fwd = 1;
  v = reshape(y, [2, size(y,1)/2]);
  I = oscITransistor(v,g_cc,g_fwd);
  vdot = reshape(I, size(y));
% oscDot
