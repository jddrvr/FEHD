function [HD,SW] = FEHD(dataSegments,lagList,cparams,maxPerms)
% Jon Drover
% This Version complete October 19 2015
%
% FEHD algorithm
% Inputs:
%   dataSegments - data as a list of segments, d{1},d{2},d{3}, etc. Must
%   use the curly brackets, even if there is only one segment.
%
%   lagList - The lags to use for the auto-regressive model. For example,
%   lagList=[1:5] would use the 1st through 5th time lags. The functions
%   will tolerate non-consecutive lags.
%
%   cparams - CHRONUX style parameters - needs cparams.Fs = sampling rate
%   and cparams.fpass = [flow fhigh] to execute FEHD on the bandwidth
%   f=[flow,fhigh]. 
%   
%   maxPerms - to avoid getting trapped in a non-optimal local minimum,
%   FEHD permutes the data and repeats multiple times. maxPerms sets the
%   maximum number of permutations to consider. If the number of possible
%   permutations (comps!) is less that maxPerms, all of them are used. 
%
% Outputs:
%   HD - Transformed data as a list of segments.
%   SW - the transformation matrix - for each i, HD{i}=SW*dataSegments{i}.
%
% This routine can be tested by running FEHDsample at the prompt.

%
%
%

% Determine the size of the data (Epochs, channels, timePoints).

numEpochs = length(dataSegments);
[numComps,~] = size(dataSegments{1});

% Initialize the transformation matrices.

SWstep = eye(numComps);
SW = zeros(numComps);

% A switch
doRandom = 0;

% Iterate through the number of components
for comps=numComps:-1:2
    oldMin = Inf; % Initialize minimal value
    nPerms = factorial(comps); % number of permutations for number of comps.
    
    % Set up the permutation pool
    
    if(nPerms > maxPerms)
        
        numPerms = maxPerms;
        doRandom = 1;
        
    else
        numPerms = nPerms;
        V = [1:comps];
        possPerms = perms(V);
        doRandom = 0;
    end
    
    % Loop through the permutations
    
    for perm = 1:numPerms
        
        % Obtain a permutation matrix
        
        P = zeros(comps);
        
        if(doRandom==1)
    
            perm_vals = randperm(comps);

            for i=1:comps
                P(i,perm_vals(i)) = 1;
            end
    
        else
            
            perm_vals = possPerms(perm,:);
            
            for i=1:comps
                P(i,perm_vals(i)) = 1;
            end
            
        end
        
        % Construct AR model
        
        [R,A] = mkAR(dataSegments,lagList);
    
        A = convertA(A);
    
        % Ortho-normalize the residuals and permute the components.
        
        clear rdecor
    
        [rdecor,~,~,decor] = pca(R');
    
        clear Adecor
    
        for i=1:length(lagList)
            Adecor(:,:,i) = inv(decor)'*A(:,:,i)*decor';
            Adecor(:,:,i) = P*Adecor(:,:,i)*P';
        end
        
        % Find the best rotation, Q.
        
        [AHM,Qtmp,SOS] = minGrangerInt(Adecor,rdecor',comps-1,1e-6,lagList,cparams);

        if(SOS < oldMin)
            
            Q{comps} = Qtmp*P;
            oldMin = SOS;
            
        end
    end
        
    SWstep = Q{comps}*decor*SWstep(1:comps,:);
    
    SW(comps,:) = SWstep(comps,:);
        
    for i=1:numEpochs
        rawd{i} = Q{comps}*decor*dataSegments{i};
        HD{i}(comps,:) = rawd{i}(comps,:);
        dataSegments{i} = rawd{i}(1:comps-1,:);
    end
    
    disp(['Completed a component'])
    
end

for i=1:numEpochs
    HD{i}(1,:) = rawd{i}(1,:);
end

SW(1,:) = SWstep(1,:);

end