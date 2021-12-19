function [out, varargout] = MonteCarloSearch(func, xInit, varargin)
%MONTECARLOSEARCH Identification method using monte-carlo approach
%
%   [out] = MonteCarloSearch(@func, xInit) searching the unknown constants
%   that compatible to the simulation that being run in 'func' using
%   monte-carlo approach. 'xInit' will be used as as starting value 
%   that will be generated as searching area in which the random value will
%   be generated.
%
%   [out, std] = MonteCarloSearch(@func, xInit) also returning the standard
%   deviation of the end result.
%
%   [out, std, iter] = MonteCarloSearch(@func, xInit) also returning the
%   total number of iterations that being processed during the computation.
%
%   [out, ...] = MonteCarloSearch(@func, xInit, 'type', Type) specifying
%   the type of Monte-Carlo approach. The possible types are the following:
%   'standard': Standard Monte-Carlo approach (will be used by default)
%   'LasVegas': Using forever loop until desired std achieved.
%   'improved': Optimized Identification using combination of Monte-Carlo 
%   and Las Vegas approach.
% 
%   [out, ...] = MonteCarloSearch(___, 'dim', Integer) Option to specifying 
%   the number of iterations for standard Monte Carlod. By default dim = 100;
%
%   [out, ...] = MonteCarloSearch(___, 'v', bool) Option to see the 
%   comparison diagram of the simulated signal and the reference signal 
%   each time random value generated. By default v = false;
%
%   [out, ...] = MonteCarloSearch(___, 'fSize', Numeric) Option to change 
%   the field size of the generated searching area. By default fSize = 1;
%
%   [out, ...] = MonteCarloSearch(___, 'std', Numeric) Option to 
%   specified the desired standard deviation for 'LasVegas' or 'improved' 
%   algorithm. By default std = 0.01;
%
%   [out, ...] = MonteCarloSearch(___, 'timeout', Numeric) Option to 
%   specified the time out of 'LasVegas' or 'improved' algortihm, so the 
%   process can be ended when timeout is reached even the desired std not
%   yet achieved.
%
%   [out, ...] = MonteCarloSearch(___, 'Params', Cell) Option to add
%   parameters that will be passed to the 'func'. This must be in cell
%   format
%

%% Input Parser
p = inputParser;
p.CaseSensitive = false;
p.addParameter('dim', 100, @(x) isnumeric(x) && (x > 1));
p.addParameter('fSize', 1, @(x) isnumeric(x));
p.addParameter('v', false, @islogical);
p.addParameter('std', 0.01, @isnumeric);
p.addParameter('timeout', 5000, @isnumeric);
p.addParameter('P', {}, @iscell);
types = {'standard', 'LasVegas', 'improved'};
p.addOptional('type', 'standard',...
    @(x) any(validatestring(x, types)));
p.parse(varargin{:});
show_plot = p.Results.v;
type = p.Results.type;
dim = p.Results.dim;
fieldSize = p.Results.fSize;
desired_std = p.Results.std;
timeout = p.Results.timeout;
func_param = p.Results.P;

%% Run Identification using MonteCarlo approach
if ~show_plot
    disp("Busy..");
end
switch type
    case 'LasVegas'
        [res, std_signal, iter] = StandardLasVegas(func, xInit, fieldSize, desired_std, show_plot, timeout, 'Params', func_param);
    case 'improved'
        if ~any(strcmp(varargin,'dim'))
            dim = 10;
        end
        [res, std_signal, iter] = ImprovedMonteCarlo(func, xInit, fieldSize, dim, desired_std, show_plot, timeout, 'Params', func_param);
    otherwise
        [res, std_signal, iter] = StandardMonteCarlo(func, xInit, fieldSize, dim, show_plot, 'Params', func_param);
end
if ~show_plot
    disp("Done..");
end
out = res;
func(out, true, func_param{:});
varargout{1} = std_signal;
varargout{2} = iter;
end

%% Standard Monte Carlo
function [out, varargout] = StandardMonteCarlo(func, xInit, fieldSize, dim, show, varargin)
xMin = xInit*(1-fieldSize);
xMax = xInit*(1+fieldSize);

p = inputParser;
p.CaseSensitive = false;
p.addOptional('minimum_value', 10^10, @(x) isnumeric(x));
p.addParameter('Params', {}, @iscell);
p.parse(varargin{:});
minimum_value = p.Results.minimum_value;
func_param = p.Results.Params;

out = xInit;

for i=1:dim
    x = xMin + (xMax-xMin).*rand(size(xInit));
    std_signal = func(x, show, func_param{:});
    if (std_signal < minimum_value)
        minimum_value = std_signal;
        out = x;
    end
end

std_signal = minimum_value;
varargout{1} = std_signal;
varargout{2} = i;

end

%% Las Vegas
function [out, varargout] = StandardLasVegas(func, xInit, fieldSize, desired_std, show, timeout, varargin)
p = inputParser;
p.CaseSensitive = false;
p.addParameter('Params', {}, @iscell);
p.parse(varargin{:});
func_param = p.Results.Params;

xMin = xInit*(1-fieldSize);
xMax = xInit*(1+fieldSize);
minimum_value = 10^10;
iter = 0;
while minimum_value > desired_std
    iter = iter + 1;
    x = xMin + (xMax-xMin).*rand(size(xInit));
    
    try
        std_signal = func(x, show, func_param{:});    
    catch
        std_signal = 10^10;
    end
    
    if (std_signal < minimum_value)
        minimum_value = std_signal;
        sig = x;
    end
    if iter > timeout
        break;
    end
end
varargout{1} = minimum_value;
varargout{2} = iter; 
out = sig;
end

%% Improved Monte Carlo
function [out, varargout] = ImprovedMonteCarlo(func, xInit, fieldSize, dim, desired_std, show, timeout, varargin)
p = inputParser;
p.CaseSensitive = false;
p.addParameter('Params', {}, @iscell);
p.parse(varargin{:});
func_param = p.Results.Params;

[X, minimum_value, iter] = StandardMonteCarlo(func, xInit, fieldSize, dim, show, 'Params', func_param);
min_prev = minimum_value;
fieldSize = fieldSize/2;
while minimum_value > desired_std
    [X, minimum_value, lp] = StandardMonteCarlo(func, X, fieldSize, dim, show, minimum_value, 'Params', func_param);
    iter = iter + lp;
    if round(log10(min_prev/minimum_value)) > 0
        fieldSize = fieldSize/2;
        min_prev = minimum_value;
    end
    if iter>timeout
        break;
    end
end
out = X;
varargout{1} = minimum_value;
varargout{2} = iter; 
end
