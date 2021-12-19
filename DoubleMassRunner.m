function [std_ref] = DoubleMassRunner(x, varargin)
%DoubleMassRunner Run the Simulation DoubleMassSimplified 
% 
%   DoubleMassRunner([1 2 3]) Passing the constants as parameter to the
%   simulink model and run the simulation, returning the standard deviation
%   in comparison to the reference signal
%
%   DoubleMassRunner([1 2 3], S) where S is boolean (true/false) value. If
%   true then the plot diagram comparing the simulation and the reference
%   signal. S is by default false.
%
%   DoubleMassRunner([1 2 3], 'SimTime', time) adding option of the
%   duration of simulation time of the Simulink Model DoubleMassSimplified
%

%% Using input parser for optional parameter;
p = inputParser;
validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
p.CaseSensitive = false;
p.addOptional('verbose',false,@islogical);
p.addParameter('SimTime', 10, validScalarPosNum);
p.parse(varargin{:});
SimTime = p.Results.SimTime;
show_plot = p.Results.verbose;

%% Load reference file once
persistent ref;
if isempty(ref)
    ref = load('var4.mat');
end

%% Specifying the constants for the Simulation
Ts = ref.Ts;
k1 = x(1);
k2 = x(2);
C12 = x(3);

%% Run simulink simulation
Out = sim('DoubleMassSimplified.slx', 'SrcWorkspace', 'current');
Yout_sim = Out.simout;
t_sim = Out.tout;

%% show the plot
if show_plot
    plot(ref.t, ref.Yout, 'b'); grid on; hold on;
    plot(t_sim, Yout_sim, 'r--'); grid on; hold off;
    legend('Reference', 'Simulation');
    xlabel('Time (s)');
    ylabel('Amplitude');
    pause(0.001);
end

%% calculate the standard deviation simulation with the reference
std_ref = std(ref.Yout-Yout_sim);

end

