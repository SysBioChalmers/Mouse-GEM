
%%%%%%%%%%%%%%%%%%%%%%%%
%% Fig2A
%%%%%%%%%%%%%%%%%%%%%%%%

% Load models
load('Human-GEM.mat');
load('Mouse-GEM.mat')
load('Rat-GEM.mat')
load('Zebrafish-GEM.mat')
load('Fruitfly-GEM.mat')
load('Worm-GEM.mat')


models{1}=ihuman;
models{2}=mouseGEM;
models{3}=ratGEM;
models{4}=zebrafishGEM;
models{5}=fruitflyGEM;
models{6}=wormGEM;

modelIds = {'Human1', 'Mouse1','Rat1','Zebrafish1','Fruitfly1','Worm1'};

% compare GEMs
compStruct = compareMultipleModels(models,false,false);

% check deviation in subsystems
subCoverage = (compStruct.subsystems.matrix - mean(compStruct.subsystems.matrix,2))./mean(compStruct.subsystems.matrix,2)*100;
subs_ind = any(abs(subCoverage) > 35,2);   % filter subsystems to plot
cmap = custom_cmap('redblue');
genHeatMap(subCoverage(subs_ind,:), 'colNames', modelIds, 'rowNames', compStruct.subsystems.ID(subs_ind),...
        'clusterDim','both','clusterDist','euclidean','colorMap', cmap,...
        'colorBounds', [-100, 100],'gridColor', 'k');
set(gca,'FontSize',14);
xtickangle(45)


%%%%%%%%%%%%%%%%%%%%%%%%%
sturucal comparison of the tINIT GEMs
%%%%%%%%%%%%%%%%%%%%%%%%%


% Load models
load('../data/tissueCelllineSpecificGEMs.mat');

% extract out models and ids
modelIds = arrayfun(@(i) tINIT_GEMs{i}.id,(1:numel(tINIT_GEMs))','UniformOutput',false);


%% compare models

% provide annotation for tissue- and cell type-specific GEMs
annotation = {
'GSE80465/Forebrain/Male/2';...     % APP23
'GSE80465/Forebrain/Male/6';...     % APP23
'GSE80465/Forebrain/Male/18';...    % APP23
'GSE80465/Forebrain/Male/24';...    % APP23
'GSE63943/HippocampalCA1/NA/10';... % APPPS1
'GSE100070/DentateGyrus/Male/4';... % APPPS1
'GSE110741/Hippocampus/Male/10';... % APPPS1
'GSE110741/Hippocampus/Male/4';...  % APPPS1
'GSE131659/Brain/Female/8';...      % APPswe/PSEN1dE9 (C57BL6)
'GSE131659/Brain/Male/8';...        % APPswe/PSEN1dE9 (C57BL6)
'GSE93678/Hippocampus/Female/13';...% APPswe/PSEN1dE9 (line 85)
'GSE104424/Hippocampus/Male/10';... % APPswe/PSEN1dE9 (line 85)
'GSE104424/Hippocampus/Male/6';...  % APPswe/PSEN1dE9 (line 85)
'GSE124266/Microglia/Female/11';... % TREM2KOJAX
'GSE124266/Microglia/Female/13';... % TREM2KOJAX
'GSE124266/Microglia/Male/11';...   % TREM2KOJAX
'GSE124266/Microglia/Male/13';...   % TREM2KOJAX
'GSE134031/Astrocyte/Female/2';...  % TREM2KOJAX
'GSE134031/Microglia/Female/2';...  % TREM2KOJAX
'GSE134031/Astrocyte/Female/16';... % TREM2KOJAX
'GSE134031/Microglia/Female/16';... % TREM2KOJAX
'GSE104381/Hippocampus/Male/8';...  % Trem2 KO (KOMP) x APPPS1
'GSE104381/Hippocampus/Male/4';...  % Trem2 KO (KOMP) x APPPS1
'GSE104381/Cortex/Male/8';...       % Trem2 KO (KOMP) x APPPS1
'GSE104381/Cortex/Male/4'...        % Trem2 KO (KOMP) x APPPS1
};

% compare models
compStruct = compareMultipleModels(tINIT_GEMs,false,false);


%% Fig S2. generate heatmap to compare model similarity (based on reaction content) with `structComp` field

% comapre rxn content using clustergram
cg = clustergram(compStruct.structComp, 'Symmetric', false,...
    'Colormap', flipud(bone), 'RowLabels', annotation,...
    'DisplayRatio', [0.1500 0.1500],...
    'ColumnLabels', annotation);
plot(cg);


%% Fig S3. generate heatmap comparing model subsystem coverage

subCoverage = (compStruct.subsystems.matrix - mean(compStruct.subsystems.matrix,2))./mean(compStruct.subsystems.matrix,2)*100;
subs_ind = any(abs(subCoverage) > 25,2);  % filter subsystems to plot

% compare model subsystem coverage
subSystemCoverage = clustergram(subCoverage(subs_ind,:), 'Colormap',...
    redbluecmap, 'DisplayRange', 100, 'rowLabels',...
    compStruct.subsystems.ID(subs_ind), 'columnLabels', annotation,...
    'ShowDendrogram', 'OFF');
plot(subSystemCoverage);


%% Fig 3C, generate 2D tSNE plot

modelType = {
'APP23';...     % 
'APP23';...     % 
'APP23';...    % 
'APP23';...    % 
'APPPS1';... % 
'APPPS1';... % 
'APPPS1';... % 
'APPPS1';...  % 
'APPswe/PSEN1dE9';... %  (C57BL6)
'APPswe/PSEN1dE9';... %  (C57BL6)
'APPswe/PSEN1dE9';... %  (line 85)
'APPswe/PSEN1dE9';... %  (line 85)
'APPswe/PSEN1dE9';... %  (line 85)
'TREM2KOJAX';... % 
'TREM2KOJAX';... % 
'TREM2KOJAX';... % 
'TREM2KOJAX';... % 
'TREM2KOJAX';... % 
'TREM2KOJAX';... % 
'TREM2KOJAX';... % 
'TREM2KOJAX';... % 
'Trem2KOAPPPS1';...  % Trem2 KO (KOMP) x APPPS1
'Trem2KOAPPPS1';...  % Trem2 KO (KOMP) x APPPS1
'Trem2KOAPPPS1';...  % Trem2 KO (KOMP) x APPPS1
'Trem2KOAPPPS1'...   % Trem2 KO (KOMP) x APPPS1
};

labels = {
'GSE80465 M 2';...     % APP23
'GSE80465 M 6';...     % APP23
'GSE80465 M18';...    % APP23
'GSE80465 M 24';...    % APP23
'GSE63943 10';... % APPPS1
'GSE100070 M 4';... % APPPS1
'GSE110741 M 10';... % APPPS1
'GSE110741 M 4';...  % APPPS1
'GSE131659 F 8';...      % APPswe/PSEN1dE9 (C57BL6)
'GSE131659 M 8';...        % APPswe/PSEN1dE9 (C57BL6)
'GSE93678 F 13';...% APPswe/PSEN1dE9 (line 85)
'GSE104424 M 10';... % APPswe/PSEN1dE9 (line 85)
'GSE104424 M 6';...  % APPswe/PSEN1dE9 (line 85)
'GSE124266 F 11';... % TREM2KOJAX
'GSE124266 F 13';... % TREM2KOJAX
'GSE124266 M 11';...   % TREM2KOJAX
'GSE124266 M 13';...   % TREM2KOJAX
'GSE134031 F 2';...  % TREM2KOJAX
'GSE134031 F 2';...  % TREM2KOJAX
'GSE134031 F 16';... % TREM2KOJAX
'GSE134031 F 16';... % TREM2KOJAX
'GSE104381 M 8';...  % Trem2 KO (KOMP) x APPPS1
'GSE104381 M 4';...  % Trem2 KO (KOMP) x APPPS1
'GSE104381 M 8';...       % Trem2 KO (KOMP) x APPPS1
'GSE104381 M 4'...        % Trem2 KO (KOMP) x APPPS1
};


% information to be presented: model type; tissue source; gender; month

% shape for tissue sources:
% Brain:         'o'  Circle
% cortex:        's'  Square
% hippocampus:   'd'  Diamond
% microglia:     'p'  Five-pointed star (pentagram)
% Astrocyte:     '^'  Upward-pointing triangle
% Dentate gyrus: 'v'  Downward-pointing triangle

% use colors for differenting models
cmap = [0.180 0.415 0.562;...  % 1. dark blue     - TREM2KOJAX
        0.891 0.312 0.100;...  % 2. orange red    - APPswe/PSEN1dE9
        0.975 0.674 0.239;...  % 3. orange yellow - Trem2KO (KOMP) x APPPS1
        0.440 0.679 0.157;...  % 4. grass green   - APPPS1
        1.800 0.000 1.000];    % 5. pink          - APP23

colorBasic = zeros(length(modelIds),1);
colorBasic(contains(modelType, 'TREM2KOJAX'))      = 1;
colorBasic(contains(modelType, 'APPswe/PSEN1dE9')) = 2;
colorBasic(contains(modelType, 'APPPS1'))          = 4;
colorBasic(contains(modelType, 'APP23'))           = 5;
colorBasic(contains(modelType, 'Trem2KOAPPPS1'))   = 3;

isBrain        = contains(annotation,'Brain','IgnoreCase',true);  
isCortex       = contains(annotation,'Cortex','IgnoreCase',true); 
isHippocampus  = contains(annotation,'Hippo','IgnoreCase',true); 
isMicroglia    = contains(annotation,'Microglia','IgnoreCase',true); 
isAstrocyte    = contains(annotation,'Astrocyte','IgnoreCase',true); 
isDentategyrus = contains(annotation,'Dentate','IgnoreCase',true); 


% plot 2D tSNE
tsnePlot = tsne(double(compStruct.reactions.matrix'),'Algorithm','exact','Distance','hamming','NumDimensions',2,'Perplexity',5);  % need to re-run tSNE for 2D

figure;
rng(20200814);

xlabel('tSNE 1');ylabel('tSNE 2');

hold on
scatter(tsnePlot(isBrain,1),tsnePlot(isBrain,2),140,cmap(colorBasic(isBrain),:),'filled','o'); % Brain
scatter(tsnePlot(isCortex,1),tsnePlot(isCortex,2),140,cmap(colorBasic(isCortex),:),'filled','s'); % Cortex
scatter(tsnePlot(isHippocampus,1),tsnePlot(isHippocampus,2),140,cmap(colorBasic(isHippocampus),:),'filled','d'); % hippocampus
scatter(tsnePlot(isMicroglia,1),tsnePlot(isMicroglia,2),260,cmap(colorBasic(isMicroglia),:),'filled','p'); % microglia
scatter(tsnePlot(isAstrocyte,1),tsnePlot(isAstrocyte,2),140,cmap(colorBasic(isAstrocyte),:),'filled','^'); % astrocyte
scatter(tsnePlot(isDentategyrus,1),tsnePlot(isDentategyrus,2),140,cmap(colorBasic(isDentategyrus),:),'filled','v'); % dentage gyrus
hold on

% use text label to present GEO id, gender and month
text(tsnePlot(:,1)-90, tsnePlot(:,2)+40, labels);

legTissue=legend('Brain','Cortex','Hippocampus','Microglia','Astrocyte','Dentate gyrus','Location','Southwest');
set(legTissue,'FontSize',15);


%% Fig S4, plot GEM performance base on the metabolic tasks that differ

taskFileName = '../data/metabolicTasks_Full.xlsx';
taskPerformance = compareMultipleModels(tINIT_GEMs,false,false,[],true,taskFileName);


isDiff = ~all(taskPerformance.funcComp.matrix == 0, 2) & ~all(taskPerformance.funcComp.matrix == 1, 2);
diffTasks = taskPerformance.funcComp.tasks(isDiff)

% visualize the matrix
spy(taskPerformance.funcComp.matrix(isDiff,:), 30);

% add labels
set(gca, 'XTick', 1:numel(tINIT_GEMs), 'XTickLabel', modelIds, 'XTickLabelRotation', 90, ...
    'YTick', 1:numel(diffTasks), 'YTickLabel', diffTasks, 'YAxisLocation', 'right');
xlabel(gca, '');
set(gca,'FontSize',12);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fig 5C proteomics data heatmap
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% study1: Mol Neurodegener. 2018 13(1):2
% Molecular and functional signatures in a novel Alzheimer’s disease mouse model assessed by quantitative proteomics

progression = {};
progression.conditions = {'APP/PS1-4month';'ADLP-4month';'APP/PS1-7month';'ADLP-7month';'APP/PS1-10month';'ADLP-10month'};
progression.genes = {'Gm2a';'Hexa';'Hexb';'Ctsb';'Ctsc';'Ctsd';'Ctss';'Ctsz'};

% import matrix from: 13024_2017_234_MOESM5_ESM.xlsx
progression.FC = [...
0	0	0	0	1.417925478	1.382074522;
0	0	1.251180359	1.260812087	1.538563536	1.527734807;
1.275993338	1.366642874	1.607705407	1.603235419	1.964457253	1.897694524;
0	0	0	0	1.277180841	1.177745402;
0	0	1.456410256	1.440615385	1.810240964	1.598934198;
0	0	1.394253493	1.41552826	1.430974478	1.413766435;
1.465475616	1.547802994	2.115828334	2.065989848	2.239342806	2.029973357;
1.298776098	1.296856251	1.605569948	1.628238342	1.974926254	2.016224189;
];


plotHeatmap = clustergram(progression.FC,...
                   'RowLabels',    progression.genes,...
                   'ColumnLabels', progression.conditions,...
                   'RowPdist', 'euclidean',...
                   'ColumnPdist', 'euclidean',...
                   'Displayrange', 2.5,...
                   'Standardize','Row',...
                   'Colormap', redbluecmap,...
                   'ImputeFun', @knnimpute);


plotHeatmap = clustergram(progression.FC,...
                   'RowLabels',    progression.genes,...
                   'ColumnLabels', progression.conditions,...
                   %'Cluster', 'Row',...
                   %'Cluster', 'Column',...
                   'RowPdist', 'euclidean',...
                   'ColumnPdist', 'euclidean',...
                   'Displayrange', 2.5,...
                   'Standardize','Row',...
                   'Colormap', redbluecmap,...
                   'ImputeFun', @knnimpute);


cgAxes =plot(plotHeatmap);

set(cgAxes, 'Clim', [80,90]) % 
set(gca,'FontSize',20);


%% study2: Exp Mol Med. 51(11):1-17
% Deep proteome profiling of the hippocampus in the 5XFAD mouse model reveals biological process alterations and a novel biomarker of Alzheimer’s disease

study2 = {};
study2.conditions = {'5XFAD-5month';'5XFAD-10month'};
study2.genes = {'Ctss';'Ctsb';'Ctsd';'Hexb';'Hexa';'Ctsz';'Gm2a';'Ctsf'};

% import matrix from: 12276_2019_326_MOESM4_ESM.xlsx
study2.FC = [...
5.212054712	4.84139003;
1.46113377	0.847124569;
1.845440234	3.058627178;
0	7.08353576;
2.462775624	4.808346022;
0	5.066205365;
0	3.276028084;
0	4.242094273;
];

plotHeatmap = clustergram(study2.FC,...
                   'Displayrange', 8,...
                   'Cluster', 'Column',...
                   'RowLabels',    study2.genes,...
                   'ColumnLabels', study2.conditions,...
                   'Colormap', redbluecmap,...
                   'ImputeFun', @knnimpute);
               
study2Heatmap =plot(plotHeatmap);
set(study2Heatmap,'ydir','reverse')
set(study2Heatmap,'FontSize',40)

