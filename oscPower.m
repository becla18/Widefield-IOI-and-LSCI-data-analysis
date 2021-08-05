function P = oscPower(X,band,fs)
%oscAmplitude calculates the average power of a signal in a given frequency band by performing a
%fast fourier transform
%%% INPUTS %%%
%X: input signal, specified as a vector or matrix. If X is a matrix, then the function calculates
%   the average power for each column.
%band: frequency band [Hz] in which to average power, specified as a two-element vector
%      [flow,fhigh].
%fs: sampling frequency of X, specified in Hz.
%%% OUTPUT %%%
%P: Average power of X, returned as a scalar if X is a vector and a vector if X is a matrix.

%no. of time domain samples
if isvector(X)
    n = length(X);
    if size(X,1) == 1 %Make sure time is in row dimension
        X = X';
    end
else
    n = size(X,1);
end
Y = fft(X);
%frequency vector for power spectrum
fmax = fs/2;
df = fs/n;
f=0:fs/n:fs/2;
%Power spectrum
S = abs(Y).^2/n;
%Select first half of DFT
if mod(n,2) %if impair number of point
    S = S(1:(n-1)/2+1,:);
else
    S = S(1:n/2+1,:);
end
%Find indexes of values at band limits
[~,i1] = min(abs(f-band(1)));
[~,i2] = min(abs(f-band(2)));
P = mean(S(i1:i2,:));
end

