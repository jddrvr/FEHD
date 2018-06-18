function [trans_func,spectm,freq] = AR_calc_spectrum(A,R,cparams,lagList,fRange)
% Jon Drover
% Version complete October 19 2015.
%
% Computes the spectral matrices using the auto-regressive model.
%
% Inputs:
%   A - AR lag matrices
%   R - AR residuals
%   cparams - sampling rate (cparams.Fs) and frequency range
%   (cparams.fpass).
%   lagList - lags in the AR model, eg [1:5].
%   fRange - Holds frequency range. Redundant, but necessary to enter. 

[N M O] = size(A);
% N - channels
% M = N
% O - number of lags.

% For each frequency, set up a new matrix layer

numFreqs = 11;

f = linspace(fRange(1),fRange(2),numFreqs);

f_length = length(f);

rr=1/cparams.Fs;
spectm = zeros(N,N,f_length);

A = permute(A,[2 1 3]);

for i=1:f_length
    spectm(:,:,i) = eye(N);
    for k=1:O % Goes through the lags
        spectm(:,:,i) = spectm(:,:,i)-A(:,:,k)*exp(-2*pi*1i*lagList(k)*rr*f(i));
    end
end

for i=1:f_length
    spectm(:,:,i) = inv(spectm(:,:,i));
    trans_func(:,:,i) = spectm(:,:,i);
end

% Want to put the scaling in here...

Rmat = zeros(M);

Rmat = (R'*R)/size(R,1);

%for i=1:M
%    Rmat(i,i) = var(R(:,i));
%    Rmat(i,i) = 1
%end

for i=1:f_length
    spectm(:,:,i) = spectm(:,:,i)*Rmat*spectm(:,:,i)';
    freq(i)=f(i);
end


spectm = spectm*rr;

spectm = permute(spectm,[2 1 3]);

end
