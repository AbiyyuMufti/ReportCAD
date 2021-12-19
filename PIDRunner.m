function [std_ref] = PIDRunner(X, varargin)
%PIDCONTROLRUNNER Summary of this function goes here
%   Detailed explanation goes here

%% Using input parser for optional parameter;
p = inputParser;
validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
p.CaseSensitive = false;
p.addParameter('Omega', 10, validScalarPosNum);
p.addParameter('SimTime', 10, validScalarPosNum);
p.addParameter('t', 'PID', @(x) any(validatestring(x,{'PID', 'Filtered'})));
p.addParameter('StepTime', 1e-3, validScalarPosNum);
p.addParameter('Object', [0.2244 0.5636 16.6338], @ismatrix);
p.addOptional('verbose', false, @islogical);

p.parse(varargin{:});
SimTime = p.Results.SimTime;
show_plot = p.Results.verbose;
Ts = p.Results.StepTime;
type = p.Results.t;
OParam = p.Results.Object;
omega = p.Results.Omega;


%% Run simulink simulation
k1 = OParam(1);
k2 = OParam(2);
C12 = OParam(3);

Kp = X(1);
Ki = X(2);
Kd = X(3);

% assignin('base','k1',k1)
% assignin('base','k2',k2)
% assignin('base','C12',C12)
% assignin('base','SimTime',SimTime)
% assignin('base','Ts',Ts)

As_et = poly(-omega*ones(1,3));
Bs_et = As_et(end);
% [Bs_et, As_et] = linmod('DoubleMassSimplified');
tf_ref = tf(Bs_et, As_et);
% step(tf_ref, 0:Ts:SimTime-Ts)
[ref.Yout, ref.t] = step(tf_ref, 0:Ts:SimTime-Ts);

switch type
    case 'Filtered'
        Simulation = 'PIDControlDMSFiltered';
        PID_Lgnd = 'PID with Filter';
    case 'PID'
        Simulation = 'PIDControlDMS';
        PID_Lgnd = 'PID';
end

Out = sim(Simulation, 'SrcWorkspace', 'current');
Yout_sim = Out.simout;
t_sim = Out.tout;
% deltaY = out.simout(1:end,2) - out.simout(1:end,1);
% std_ref = std(deltaY);
std_ref = std(ref.Yout-Yout_sim);

%% Show plot
% if show_plot
%     plot(out.tout, out.simout(1:end,2)); grid on; hold on;
%     plot(out.tout, out.simout(1:end,1)); hold off;
%     legend('Reference', PID_Lgnd);
%     pause(0.001);
% end

if show_plot
    plot(ref.t, ref.Yout, 'b'); grid on; hold on;
    plot(t_sim, Yout_sim, 'r'); grid on; hold off;
    legend('Reference', PID_Lgnd);
    xlabel('Time (s)');
    ylabel('Amplitude');
    pause(0.001);
end

end

