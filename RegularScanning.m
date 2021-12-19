function [out, varargout] = RegularScanning(func, xInit, varargin)
%REGULARSCANNING Identification method using regular scanning approach
%
%   [out] = RegularScanning(@func, xInit) searching the unknown constants
%   that compatible to the simulation that being run in 'func' using
%   regular scanning approach. 'xInit' will be used as as starting value 
%   that will be generated as searching area.
%
%   [out, std] = RegularScanning(@func, xInit) also returning the standard
%   deviation of the end result.
%
%   [out, std, iter] = RegularScanning(@func, xInit) also returning the
%   total number of iterations that being processed during the computation.
%
%   [out, ...] = RegularScanning(@func, xInit, 'improved') using improved 
%   method of regular scanning that will iterate the standard regular 
%   scanning and optimize the searching area. Best to use if the searching 
%   area is unknown.
% 
%   [out, ...] = RegularScanning(___, 'dim', Integer) Option to specifying 
%   the dimension each searching area. The input must be an Integer. By 
%   default dim = 5;
%
%   [out, ...] = RegularScanning(___, 'v', bool) Option to see the 
%   comparison diagram of the simulated signal and the reference signal 
%   in every single loop. By default v = false;
%
%   [out, ...] = RegularScanning(___, 'fSize', Numeric) Option to change 
%   the field size of the generated searching area. By default fSize = 1;
%
%   [out, ...] = RegularScanning(___, 'improved', 'std', Numeric) Option to 
%   specified the desired standard deviation of improved regular scanning 
%   algorithm. By default std = 0.01;
%
%   [out, ...] = RegularScanning(___, 'Params', Cell) Option to add
%   parameters that will be passed to the 'func'. This must be in cell
%   format
%

%% Input Parser
p = inputParser;
p.CaseSensitive = false;
p.addParameter('dim', 5, @(x) isnumeric(x) && (x > 1));
p.addParameter('fSize', 1, @(x) isnumeric(x));
p.addParameter('v', false, @islogical);
p.addParameter('std', 0.01, @isnumeric);
p.addOptional('type', 'standard', @(x) any(validatestring(x,{'standard', 'improved'})));
p.addParameter('Params', {}, @iscell);
p.parse(varargin{:});
show_plot = p.Results.v;
type = p.Results.type;
dim = p.Results.dim;
fieldSize = p.Results.fSize;
desired_std = p.Results.std;
func_param = p.Results.Params;

%% Running Regular Scanning
if ~show_plot
    disp("Busy...");
end
if type == "improved"
    [out, std_signal, iter] = ImprovedRegularScanning(func, xInit, fieldSize, dim, show_plot, desired_std, func_param);
else
    [out, std_signal, iter] = StandardRegularScanning(func, xInit, fieldSize, dim, show_plot, func_param);    
end
if ~show_plot
    disp("Done");
end
func(out, true, func_param{:});
varargout{1} = std_signal;
varargout{2} = iter;
end


%% Standard approach of Regular Scanning
function [out, varargout] = StandardRegularScanning(func, xInit, fieldSize, dim, show_plot, func_param)
minimum_value = 10^10;
iter = 0;
xMin = xInit*(1-fieldSize);
xMax = xInit*(1+fieldSize);
if all(xMin == 0)
    xMin = 0.1*ones(1,length(xInit));
end
X = zeros(dim, length(xInit));
for i = 1:length(xInit)
    X(:,i) = linspace(xMin(i), xMax(i), dim);
end
for i=1:dim
    for j=1:dim
        for k=1:dim 
            iter = iter + 1;
            x = [X(i,1), X(j,2), X(k,3)];
            std_signal = func(x, show_plot, func_param{:});
            if(std_signal < minimum_value)
                minimum_value = std_signal;
                x_opt = x;
            end
        end
    end
end
varargout{1} = minimum_value;
varargout{2} = iter;
out = x_opt;
end


%% Improved Regular Scanning with searching algorithm
function [out, varargout] = ImprovedRegularScanning(func, xInit, fieldSize, dim, show_plot, desired_std, func_param)

[X, minimum_value, iter] = StandardRegularScanning(func, xInit, fieldSize, dim, show_plot, func_param);
fieldSize = 0.5;
while fieldSize >= 0.1
    [sig, std_signal, lp] = StandardRegularScanning(func, X, fieldSize, dim, show_plot, func_param);
    iter = iter + lp;
    if(std_signal < minimum_value)
        minimum_value = std_signal;
        X = sig;
        fieldSize = fieldSize - 0.1;
    end
    if(std_signal <= desired_std)
        minimum_value = std_signal;
        X = sig;
        break;
    end
end
varargout{1} = minimum_value;
varargout{2} = iter;
out = X;
end