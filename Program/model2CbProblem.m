function problem = model2CbProblem(varargin)
%% problem = model2CbProblem(varargin)
% This function constructs an optimization problem structure
% for a given solver.
% Input:
%   struct model:       metabolic model
%   char solver:        solver for which the problem should be constructed
% Output:
%   struct problem:     optimization problem

%% Parse input
p = inputParser;

addRequired(p,'model',@isstruct);
addRequired(p,'solver',@ischar)

parse(p,varargin{:});

model = p.Results.model;
solver = p.Results.solver;

%% Construct constraint-based problem
problem = struct;

switch solver
    
    case 'gurobi'
        % define objective
        problem.obj = model.c;
        % objective sense
        problem.modelsense = model.osenseStr;
        % equality LHS
        problem.A = sparse(model.S);
        % equality RHS
        problem.rhs = model.b;
        % constraint sense
        problem.sense = translateSenseCobra2Gurobi(model.csense);
        % lower bounds
        problem.lb = model.lb;
        % upper bounds
        problem.ub = model.ub;

    case 'ibm_cplex'
        disp('Not implemented')
end

end