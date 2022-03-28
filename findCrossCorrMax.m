function [corrMax, lag] = findCrossCorrMax(x1,x2,range)
%FINDCROSSCORRMAX computes the maximum cross-correlation value beween two signals in a specified
%range.
%   ARGUMENTS
%   x1, x2: signal vectors
%   range: maximum lag (from to -range to +range)considered for the cross-correlation, in number of samples.
%   OUTPUTS
%   corrMax: max cross correlation value
%   lag: lag corresponding to corrMax
[r,lags] = xcorr(x1,x2,range);
[~,ind] = max(abs(r));
corrMax = r(ind);
lag = lags(ind);
end

