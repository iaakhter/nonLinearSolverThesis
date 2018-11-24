% Itrat Ahmed Akhter
% CPSC 538G Proposal
% main_script.m
% figures for project report



%Figure 4 to show convergence around unstable equilibrium for 3 pc
t = 0:0.1:10;
y = [0.34; -0.34; 0.34; -0.34];
g_cc = 0.5;
inverterFunc = @inverterTanh;
a = 5.0;
[t,y] = ode45(@(t,y)oscDot(t,y,a,inverterFunc,g_cc,0),t,y);
slice = 100;
y(slice,:)
figure;
plot(t(1:slice),y(1:slice, :));
xlim([0,10]);
ylim([-1.0,1.0]);
xlabel('Time','FontSize', 15);
ylabel('Voltage','FontSize', 15);
lgd = legend('V(0,0)', 'V(0,1)', 'V(1,0)', 'V(1,1)');
lgd.FontSize = 15;
%function schmittPlot(timeScale, LineStyle)
%	if(nargin < 1) timeScale = 200.0; end
%	if(nargin < 2) LineStyle = 'b-'; end
%	t = [0 3.6];
%	V = [1.8; 0.1; 1.7];

%	[t,V] = ode45(@(t,V)schmittDot(t,[V(1),V(2),V(3),vin(t)], timeScale),t,V);
%	if(nargin < 2) figure; end;
%	plot(t, V(:,1));
%	hold on;

%	in = zeros(1,length(t))
%	for i = 1:length(t)
%		in(i) = vin(t(i));
%	end
%	plot(t,in);
%	xlabel('Time', 'FontSize', 15);
%	ylabel('Voltage','FontSize', 15);
%	lgd = legend('Vout', 'Vin');
%	lgd.FontSize = 15;
%end % schmittPlot
