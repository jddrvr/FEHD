function [A_HM,Q_total,SOS] = minGrangerInt(A,R,k,tol,lagList,cparams)
% Jon Drover
% Version complete October 19 2015
%
% This function determines the rotation that minimizes the Granger
% causality over the given frequency range.
%
% Inputs:
%   A - the AR lag matrices (comps,comps,lags)
%
%   R - the residuals 
%
%   k - the size of the "caused" group.
%
%   tol - a tolerance parameter that determines how hard the minimizer has
%   to work.
%
%   lagList - same as in FEHD.
%   cparams - same as in FEHD.

% Storage
A_tmp = A;
[M,N,O] = size(A);

% Create a matrix to be used for rotations.
Q_total = eye(M);
% A counter
iterations = 0;

dif = 1;
minVal=Inf;

while dif>tol
  
    iterations = iterations + 1;    
   
    % Breakout
    if iterations > 2000
        keep_looking = 0;
        disp(['too many iterations. Breaking.']);
        A_HM = eye(M);
        Q = eye(M);
        break;
    end
    
    dif = 0.0;
    
    for p = 1:k
      for q = k+1:M
  
        upper = pi/2;
        lower = -pi/2;
        
    	%divisor = 125;
        divisor = 60;
    	step = (upper-lower)/divisor;
    	fevalspots = lower:step:upper-step;
    
    	while upper-lower > 0
            
            [SOS,thVal,Q] = grangerBlockMin(A_tmp,R,p,q,fevalspots,k,lagList,cparams);
            
            lower = thVal-step;
            upper = thVal+step;
            step = (upper-lower)/divisor;
            fevalspots = lower:step:upper-step;                        
        end
        
        for l=1:O
            A_tmp(:,:,l) = Q*A_tmp(:,:,l)*Q';
        end 
        
        Q_total = Q*Q_total;
        
      end
    end
    
    dif = minVal-SOS;
    if(dif < 0)
        disp(['negative where it should not be']);
    end
    
    minVal = SOS;
        
end

A_HM = A_tmp;
end

        