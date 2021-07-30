function I = fcs(A,negativeWeightsOpt)
%fcs computes the functional connectivity strength for nodes of the adjacency matrix A.
%   FCS is computed as the mean weight between a node and all other nodes
%%% INPUT %%%
% A: m x m adjacency matrix
% negativeWeightsOpt: Specify how to treat negative weights. Options:
% 'standard' (default): No treatment of negative weights
% 'absolute': A is treated as abs(A)
% 'zero': Negative weights are set to 0
%%% OUTPUT %%%
% I: m element vector of FCS for each node in A
if nargin == 1
    negativeWeightsOpt = 'standard';
end
if strcmp(negativeWeightsOpt,'absolute')
    A = abs(A)
elseif strcmp(negativeWeightsOpt,'zero')
    A(A<0) = 0;
elseif ~strcmp(negativeWeightsOpt,'standard')
    error('Undefined 2nd argument');
end
m = size(A,1); %no. of nodes
I = (sum(A,2)-1)/(m-1); %Average of the off-diagonal elements (all 1s)
end
