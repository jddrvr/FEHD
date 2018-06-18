function [CI,CIofF,f] = grangerInt(A,R,k,lagList,cparams)
% Jon Drover
% Version Complete October 19 2015.
%
% Calculates the integral of the causal index from f1 to f2
%
% Inputs:
%   Same as previously.
% Outputs:
%   CI - the value of the Granger integral.
%   CIofF - Granger causality at each frequency.
%   f - frequnecy array.
%

[M,~,L] = size(A);

% Calculate the AR estimate spectrum

[H,S,f] = AR_calc_spectrum(A,R,cparams,lagList,cparams.fpass);

Rmat = (R'*R)/size(R,1)*1.0/cparams.Fs;

Hxx = H(1:k,1:k,:);

sigmaX = Rmat(1:k,1:k);

Sfull = S(1:k,1:k,:);

Sself = zeros(k,k,size(Sfull,3));

CIofF = zeros(1,size(Sfull,3));

for i=1:size(H,3)
    Sself(:,:,i) = Hxx(:,:,i)*sigmaX*Hxx(:,:,i)';
    
    CIofF(i) = real(log(det((Sfull(:,:,i)))/det((Sself(:,:,i)))));
    
end
% changed things a little here - /length(CIofF)
CI = (f(2)-f(1))*sum(CIofF);

end

