function I = fcs(W,negativeWeightsOpt)
%fcs computes the functional connectivity strength for nodes of the weight matrix W.
%   FCS is computed as the mean weight between a node and all other nodes
%%% INPUT %%%
% W: m x m weight matrix
% negativeWeightsOpt: Specify how to treat negative weights. Options:
% 'standard' (default): No treatment of negative weights
% 'absolute': W is treated as abs(A)
% 'zero': Negative weights are set to 0
%%% OUTPUT %%%
% I: m element vector of FCS for each node in W
if nargin == 1
    negativeWeightsOpt = 'standard';
end
if strcmp(negativeWeightsOpt,'absolute')
    W = abs(W);
elseif strcmp(negativeWeightsOpt,'zero')
    W(W<0) = 0;
elseif ~strcmp(negativeWeightsOpt,'standard')
    error('Undefined 2nd argument');
end
m = size(W,1); %no. of nodes
I = (sum(W,2)-1)/(m-1); %Average of the off-diagonal elements (all 1s)
end

