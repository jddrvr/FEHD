function [Xintegral,XofF,f] = PGC(ETS,cparams,lagList)

numEpochs = length(ETS); 
[M,N] = size(ETS{1}); % Assumes that the epochs all have the same length.
numComps = M;

X = zeros(numComps,numComps);

for i=1:numEpochs
    dset{i} = zeros(2,N); % Initialize
end

for i=1:numComps
    for j=1:numComps
        
        if(j~=i)
        
            for k=1:numEpochs
                dset{k}(1,:) = ETS{k}(i,:);
                dset{k}(2,:) = ETS{k}(j,:);
            end
        
            [R,A] = mkAR(dset,lagList);
            
            A = convertA(A);
            
            Rmat = R'*R/size(R,1);
            
            C = Rmat(1,2);
            sig2 = Rmat(1,1);
            
            P = [1 0;-C'*inv(sig2) 1];
            
            A2 = permute(A,[2 1 3]);
            
            for count=lagList
                
                A3(:,:,count) = P*A2(:,:,count)*inv(P);
                
            end
            
            A = permute(A3,[2 1 3]);
            
            R = (P*R')';
            
            % R'*R/size(R,1)
            
            [CI,CIofF,f] = grangerInt(A,R,1,lagList,cparams);
        
            X(i,j) = CI;
            XI{i,j} = CIofF;
        end
        
    end
end
%f = linspace(cparams.fpass(1),cparams.fpass(2),21);
Xintegral = X;
XofF = XI;
