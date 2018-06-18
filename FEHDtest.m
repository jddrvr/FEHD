% FEHDtest

% Runs the example shown in figures 2 and 3 of the paper.

% Load the data

load modelTimeSeries.mat

% Set the sampling rate for the data.
% This is the same for each of the examples.
cparams.Fs=1;

% ------------------------------------------------------
% REPRODUCE FIGURE 3

X = Xhrcl;

% Generate random linear combinations of the data

mixMatrix = rand(16,3);

Y = mixMatrix*X;

% Combine into principal components

% pca code written by Mike Repucci, and is also included in the supplementary
% material to the paper:
% Repucci, Schiff, Victor, General strategy for hierarchical decomposition
% of multivariate time series: Implications for temporal lobe seizures,
% Annals of Biomedical Engineering, 29:1135-1149: 2001
[T,~,~,L] = pca(Y);

% Keep the first 3 PCs

PCs{1} = T(1:3,:);

% Use 5 lags for AR model

lagList = [1:2];

% Call FEHD for frequency range [0.15 0.25]

cparams.fpass = [0.0 0.5];     

[HDfigure3,SW_figure3] = FEHD(PCs,lagList,cparams,6);

% eegplot.m written by Mike Repucci, and is also included in the supplementary
% material to the paper:
% Repucci, Schiff, Victor, General strategy for hierarchical decomposition
% of multivariate time series: Implications for temporal lobe seizures,
% Annals of Biomedical Engineering, 29:1135-1149: 2001
fig3 = eegplot(HDfigure3{1});
set(fig3,'name','Figure 3 hierarchical Generators');

% ----------------------------------------------------------------
% REPRODUCE FIGURE 4

X = XtwoToOne;

% Generate random linear combinations of the data

mixMatrix = rand(16,3);

Y = mixMatrix*X;

% Combine into principal components

[T,~,~,L] = pca(Y);

% Keep the first 3 PCs

PCs{1} = T(1:3,:);

% Use 2 lags for AR model

lagList = [1:2];

% Call FEHD for frequency range [0.0 0.5]

cparams.fpass = [0.0 0.5];     

[HDfigure4,SW_figure4] = FEHD(PCs,lagList,cparams,6);

fig4 = eegplot(HDfigure4{1});
set(fig4,'name','Figure 4 Hierarchical Generators');

[fig4PGC,~,f] = PGC(HDfigure4,cparams,lagList);

figure('name','Figure 4 Pairwise GC','NumberTitle','off')
imagesc(fig4PGC);colorbar;
%-------------------------------------
% REPRODUCE FIGURE 5

X = XtwoNets;

% Generate random linear combinations of the data

mixMatrix = rand(16,4);

Y = mixMatrix*X;

% Combine into principal components

[T,~,~,L] = pca(Y);

% Keep the first 3 PCs

PCs{1} = T(1:4,:);

% Use 2 lags for AR model

lagList = [1:2];

% Call FEHD for frequency range [0.0 0.5]

cparams.fpass = [0.0 0.5];     

[HDfigure5,SW_figure5] = FEHD(PCs,lagList,cparams,6);

fig5 = eegplot(HDfigure5{1});
set(fig5,'name','Figure 5 Hierarchical Generators');

[fig5PGC,~,f] = PGC(HDfigure5,cparams,lagList);

figure('name','Figure 5 Pairwise GC','NumberTitle','off');
imagesc(fig5PGC);colorbar;

%--------------------------------------
% REPRODUCE FIGURE 6 (CYCLIC GRAPH)

X = Xcyclic;

% Generate random linear combinations of the data

mixMatrix = rand(16,3);

Y = mixMatrix*X;

% Combine into principal components

[T,~,~,L] = pca(Y);

% Keep the first 3 PCs

PCs{1} = T(1:3,:);

% Use 2 lags

% Call FEHD for frequency range [0.0 0.5]

cparams.fpass = [0.0 0.5];

[HDfigure6,SW_figure6] = FEHD(PCs,lagList,cparams,6);

fig6 = eegplot(HDfigure6{1});
set(fig6,'name','Figure 6 Hierarchical Generators');

[fig6PGC,~,f] = PGC(HDfigure6,cparams,lagList);

figure('name','Figure 6 Pairwise GC','NumberTitle','off');
imagesc(fig6PGC);colorbar;

% -----------------------------------
% REPRODUCE FIGURE 8

X = XtwoWay;

% Generate random linear combinations of the data

mixMatrix = rand(16,3);

Y = mixMatrix*X;

% Combine into principal components

[T,~,~,L] = pca(Y);

% Keep the first 3 PCs

PCs{1} = T(1:3,:);

% Use 5 lags for AR model

lagList = [1:5];

% Call FEHD for frequency range [0.15 0.25]

cparams.fpass = [0.15 0.25];    

[HD_15_25,SW_15_25] = FEHD(PCs,lagList,cparams,6);

% Call FEHD for freuqnecy range [0.35 0.45];

cparams.fpass = [0.35 0.45];

[HD_35_45,SW_35_45] = FEHD(PCs,lagList,cparams,6);

f1 = eegplot(X);
set(f1,'name','Figure 7 Generators');

%f2 = eegplot(Y);
%set(f2,'name','Random linear combinations');

%f3 = eegplot(PCs{1});
%set(f3,'name','First 3 Principal components');

f4 = eegplot(HD_15_25{1});
set(f4,'name','Figure 7 FEHD, 0.15-0.25');

f5 = eegplot(HD_35_45{1});
set(f5,'name','Figure 7 FEHD, 0.35-0.45');

