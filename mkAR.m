function [R,A] = mkAR(data,lagList)

numLags = size(lagList,2);

numEpochs = length(data);

numChannels = size(data{1},1);

npts = size(data{1},2);

for i=1:numEpochs
  if(size(data{i},2) ~= npts) 
    disp(['Error: Data epochs are not all the same length']);
    exit
  end
end

maxLag = max(lagList);

nptsEff = npts-maxLag;

LHSdata = [];
RHSdata = zeros(numEpochs*nptsEff,1+numLags*numChannels);
RHSdata(:,1) = ones(numEpochs*nptsEff,1);


for i=1:numEpochs

  newData = data{i}(:,maxLag+1:npts);

  LHSdata = [LHSdata;newData'];

  for j=1:numLags

    newData = data{i}(:,maxLag+1-lagList(j):npts-lagList(j));

    RHSdata((i-1)*nptsEff+1:i*nptsEff,(j-1)*numChannels+2:j*numChannels+1)=newData';

  end

end

A = inv(RHSdata'*RHSdata)*RHSdata'*LHSdata;

R = LHSdata-RHSdata*A;

end
