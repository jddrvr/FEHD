function [minSoS,x,Q] = grangerBlockMin(A,R,p,q,theta,k,lagList,cparams)
% Jon Drover
% Version complete October 19 2015.
%
% This calculates the Granger causality for each of the rotations, the
% minimum of which is determined and returned.
%
% Inputs: 
%   A - AR lag matrices
%   R - Residuals
%   p,q - size of caused and casual groups. For FEHD, q=1.
%   theta - list of angles of rotation.
%   k - Size of caused group (yes, it is redundant).
%   lagList - same as in FEHD.
%   cparams - same as in FEHD.
%
% Outputs:
%   minSOS - the minimal Granger Causality integral.
%   x - the value of theta that minimized the Granger causality.
%   Q - The rotation matrix corrsepnding to theta=x.
%   

[M,~,L] = size(A);

AHD = zeros(M,M,L);

thetaResults = zeros(length(theta),1);

for i = 1:length(theta)
    
    thVal = theta(i);
    
    cc = cos(thVal);
    ss = sin(thVal);

    Qrot = eye(M);

    Qrot(p,p) = cc;
    Qrot(p,q) = -ss;
    Qrot(q,q) = cc;
    Qrot(q,p) = ss;

    for j=1:L
        AHD(:,:,j) = Qrot*A(:,:,j)*Qrot';
    end

    thetaResults(i) = grangerInt(AHD,R,k,lagList,cparams);
    
end

[minSoS,ix] = min(thetaResults);

x = theta(ix);

Q = eye(M);

Q(p,p) = cos(x);
Q(p,q) = -sin(x);
Q(q,q) = cos(x);
Q(q,p) = sin(x);

end

