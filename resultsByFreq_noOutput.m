function [PC,HD,X,Xi,f,specFig_handl,PGCfig_handl] = resultsByFreq_noOutput(SW,dataList,index,cparams);

numComps = size(SW{index},1);

% Create the HD data

[M,N] = size(dataList);

if(M>N)
    dataList = dataList';
end

lagList = cparams.lagList;

[T,~,~,L] = pca(dataList);

numEpochs = size(T,2)/cparams.epochPts;

epochPts = cparams.epochPts;

for i=1:numEpochs
    PC{i} = T(1:numComps,(i-1)*epochPts+1:i*epochPts);
    HD{i} = SW{index}*dataList(:,(i-1)*epochPts+1:i*epochPts);
end

% Form an array for each

for i=1:numEpochs
    dataCube(:,:,i) = PC{i}';
    HDcube(:,:,i) = HD{i}';
end

cparams.trialave = 1;

dataCube = permute(dataCube,[1 3 2]);
HDcube = permute(HDcube,[1 3 2]);

tempfpass = cparams.fpass;

cparams.fpass = [0 55];

for i=1:numComps
    [SPC{i},f] = mtspectrumc(dataCube(:,:,i),cparams);
    [SHD{i},f] = mtspectrumc(HDcube(:,:,i),cparams);
end

fout=f;

ttlString = strcat('frequency',' ',num2str(index-1),'-',num2str(index+1));
specFig_handl=figure('position',[100 100 2000 2000],'name',ttlString,'NumberTitle','off','Visible','off');


for i=1:numComps
    ax1=subplot(numComps,4,4*i-3)
    plot(f,10*log10(SPC{i}));
    subplot(numComps,4,4*i-2)
    topoplot(abs(L(i,:)),'bipolar.loc')
    ax2=subplot(numComps,4,4*i-1)
    plot(f,10*log10(SHD{i}));
    subplot(numComps,4,4*i)
    topoplot(abs(SW{index}(i,:)),'bipolar.loc')
    
    axis(ax1,[cparams.fpass(1) cparams.fpass(2) -Inf Inf])
    axis(ax2,[cparams.fpass(1) cparams.fpass(2) -Inf Inf])
end
Xi = eye(numComps);
Xm = eye(numComps);

cparams.fpass = tempfpass;

[X,Xi,f] = PGC(HD,cparams);

PGCfig_handl = figure('Visible','off');
imagesc(X);
axis('image')
colorbar;

%saveas(f1,'standardFigure.jpg');

Xm = X;
end
