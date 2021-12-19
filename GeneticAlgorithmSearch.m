function [out, varargout] = GeneticAlgorithmSearch(func, xInit, varargin)
%GENETICALGORITHMSEARCH Identification method using genetic algorithm approach
%
%   [out] = GeneticAlgorithmSearch(@func, xInit) searching the unknown 
%   constants that compatible to the simulation that being run in 'func' 
%   using genetic algorithm approach. 'xInit' will be used as as starting 
%   value that will be generated as searching area in which the first
%   population will be generated
%
%   [out, std] = GeneticAlgorithmSearch(@func, xInit) also returning the 
%   standard deviation of the end result.
%
%   [out, std, iter] = GeneticAlgorithmSearch(@func, xInit) also returning 
%   the total number of iterations that being processed during the 
%   computation.
% 
%   [out, ...] = GeneticAlgorithmSearch(___, 'n_pop', Numeric) Option to 
%   specifying the number of unit in a population. By default n_pop = 25;
%
%   [out, ...] = GeneticAlgorithmSearch(___, 'n_gen', Numeric) Option to 
%   specifying the maximum number of generations. By default n_gen = 10;
%
%   [out, ...] = GeneticAlgorithmSearch(___, 'n_best', Numeric) Option to 
%   specifying the number of best units to be selected. By default n_best = 10;
%
%   [out, ...] = GeneticAlgorithmSearch(___, 'mutation', Numeric) Option to 
%   specifying the maximum number of generations. By default n_gen = 10;
%
%   [out, ...] = GeneticAlgorithmSearch(___, 'v', bool) Option to see the 
%   comparison diagram of the simulated signal and the reference signal 
%   each time random value generated. By default v = false;
%
%   [out, ...] = GeneticAlgorithmSearch(___, 'fSize', Numeric) Option to change 
%   the field size of the generated searching area. By default fSize = 1;
%
%   [out, ...] = GeneticAlgorithmSearch(___, 'std', Numeric) Option to 
%   specified the desired standard deviation. If the desired std has been
%   reached then the program will ended.
%
%   [out, ...] = GeneticAlgorithmSearch(___, 'Params', Cell) Option to add
%   parameters that will be passed to the 'func'. This must be in cell
%   format
%
%% Input Parser
p = inputParser;
p.CaseSensitive = false;
p.addParameter('n_pop', 25, @(x) isnumeric(x) && (x > 1));
p.addParameter('n_best', 10, @(x) isnumeric(x) && (x > 1));
p.addParameter('mut', 3, @(x) isnumeric(x) && (x > 1));
p.addParameter('fSize', 1, @(x) isnumeric(x));
p.addParameter('v', false, @islogical);
p.addParameter('std', 0.01, @isnumeric);
p.addParameter('n_gen', 10, @isnumeric);
p.addParameter('P', {}, @iscell);
p.parse(varargin{:});
show_plot = p.Results.v;
n_best = p.Results.n_best;
number_pop = p.Results.n_pop;
fieldSize = p.Results.fSize;
desired_std = p.Results.std;
n_gen = p.Results.n_gen;
mut_occ = p.Results.mut;
func_param = p.Results.P;

%% Run the identification process
[out, std_signal, iter] = StandardGeneticAlgorithm(func, xInit, fieldSize,...
    number_pop, n_best, mut_occ, show_plot, n_gen, desired_std, 'Params', func_param);
func(out, true, func_param{:});
varargout{1} = std_signal;
varargout{2} = iter;

end

%% Improved Algortihm
function [out, varargout] = StandardGeneticAlgorithm(func, xInit, fieldSize, n_pop, n_best, mut_occ, show, max_gen, des_std, varargin)
p = inputParser;
p.CaseSensitive = false;
p.addParameter('Params', {}, @iscell);
p.parse(varargin{:});
func_param = p.Results.Params;

xMin = xInit*(1-fieldSize);
xMax = xInit*(1+fieldSize);
iter = 0;
n_genetic = length(xInit);
X = zeros(n_pop, n_genetic + 1);
minimum_value = 10^10;
min_prev = minimum_value;

% Populate first generation
for i=1:n_pop
    x_gen = xMin + (xMax-xMin).*(rand(size(xMin)));
    std_res = func(x_gen, show, func_param{:});
    X(i,:) = [x_gen, std_res];
end

% Sort the matrix X and slice the best
Xbest = sort_slice(X, n_best);
for j=1:max_gen
    iter = iter + 1;
    Child_prev = zeros(1, n_genetic);
    std_prev = 0;
    for i=1:n_pop
        Father = Xbest(randperm(length(Xbest),1), 1:end-1);
        Mother = Xbest(randperm(length(Xbest),1), 1:end-1);
        sel = randperm(2,1);
        sep = randperm(n_genetic-1, 1);
        if sel == 1
            Child = [Father(1:sep), Mother(sep+1:end)];
        else 
            Child = [Mother(1:sep), Father(sep+1:end)];
        end
        
        if mod(iter, mut_occ) == 0
            mutation_value = xMin + (xMax-xMin).*(rand(size(xMin)));
            mut_pos = randperm(n_genetic, 1);
            Child(mut_pos) = mutation_value(mut_pos);
        end
        
        % Run simulink only if new gen passed, else copy
        if all(Child_prev == Child)
            std_res = std_prev;
        else
            std_res = func(Child, show, func_param{:});    
        end
        
        X(i,:) = [Child, std_res];
        Child_prev = Child;
        std_prev = std_res;
    end
    Xbest = sort_slice(X, n_best);
    
    if minimum_value > Xbest(1,end)
        x_opt = Xbest(1,1:end-1);
        minimum_value = Xbest(1,end);
    end

    if round(log10(min_prev/minimum_value)) > 0
        fieldSize = fieldSize/2;
        min_prev = minimum_value;
    end

    xMin = x_opt*(1-fieldSize);
    xMax = x_opt*(1+fieldSize);
    
    if des_std > minimum_value
        break;
    end

end
out = x_opt;
std_signal = minimum_value;
varargout{1} = std_signal;
varargout{2} = iter*n_pop;
end

%% sort and slice
function [out] = sort_slice(X, n_best)
    [~, I] = sort(X(:,end));
    Xsorted = X(I,:);
    out = Xsorted(1:n_best,:);
end