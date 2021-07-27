function success = plotResults(dataArray,HDtrans_file,cparams)

% Set up some strings so that everything is plotted correctly. 

base = HDtrans_file(1:strfind(HDtrans_file,'.dat')-1);
subject = base(1:6);
visit = base(8:9);
freqRange = base(11:end);
lowFreq = str2num(freqRange(1:strfind(freqRange,'_')-1));
highFreq = str2num(freqRange(strfind(freqRange,'_')+1:end));

cparams.fpass = [lowFreq highFreq];

% this function plots the results for a given transformation.

% First step is to plot the main figure

% just for compatability
SW{1} = load(HDtrans_file,'-ascii');

% I need epoch points, but I am going to assume it is in cparams.



[PC,HD,X,XofF,freq,spec_handl,PGC_handl] = resultsByFreq_noOutput(SW,dataArray,1,cparams);
% going to need epochPts, perhaps assume it is part fo cparams?
% I am going to want to rewrite RBF to accomodate this stuff. 

saveas(spec_handl,strcat(subject,'_',visit,'_',num2str(lowFreq),'_',num2str(highFreq),'_specAndHead.jpg'));
saveas(PGC_handl,strcat(subject,'_',visit,'_',num2str(lowFreq),'_',num2str(highFreq),'_PGC.jpg'));

% Plot a couple of random epochs, both in PC and HD.
numEpochs = size(dataArray,2)/cparams.epochPts;

% Choose some random segments:
numRandomSegments = 2;
randomSegments = randi([1 numEpochs],numRandomSegments,1);

f1 = figure('Visible','off');

for iter = 1:numRandomSegments
    
    subplot(2,numRandomSegments,iter)
    plotComps(PC{iter});
    subplot(2,numRandomSegments,iter+numRandomSegments);
    plotComps(HD{iter});
    
end

saveas(f1,strcat(subject,'_',visit,'_',num2str(lowFreq),'_',num2str(highFreq),'_samples.jpg'));

% From XofF obtain the by frequency value of the GC

GCints_handl = plot_GC_integrals(XofF,freq);

saveas(GCints_handl,strcat(subject,'_',visit,'_',num2str(lowFreq),'_',num2str(highFreq),'_GCintegrals.jpg'));

end

