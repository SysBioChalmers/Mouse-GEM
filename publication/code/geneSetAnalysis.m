%
% code for gene set analysis
%

%% Load the KEGG_GSC.gmt gene set collection (GSC) file using the importGSC function

gscKEGG = importGSC('../data/MSigDB/c2.cp.kegg.v7.0.symbols.mouse.gmt');
gscReactome = importGSC('../data/MSigDB/c2.cp.reactome.v7.0.symbols.mouse.gmt');
gscHallmark = importGSC('../data/MSigDB/h.all.v7.0.symbols.mouse.gmt');
gscCombine = importGSC('../data/MSigDB/KeggReactiomeHallmark.mouse.gmt');


% load mouse model
load('Mouse-GEM.mat')


% Extract a metabolite-based GSC from the model using the extractMetaboliteGSC function.
gsc_met = extractMetaboliteGSC(mouseGEM, true);

% Extract a subsystem-based GSC from the model using the extractMetaboliteGSC function.
gsc_subSys = extractSubsystemGSC(mouseGEM);



%% gene-set analysis of selected models using Mouse1

% APPPS1, APPswe/PSEN1dE9, APP23, Trem2 KO (KOMP) x APPPS1, Trem2 KO (JAX)

% specify data annotations
data_accessions = {...
'GSE63943/HippocampalCA1/NA/10'...          % APPPS1
'GSE100070/DentateGyrus/Male/4'...          % APPPS1
'GSE110741/Hippocampus/Male/4'...           % APPPS1
'GSE110741/Hippocampus/Male/10'...          % APPPS1
'GSE131659/C57BL6/Brain/Female/8'...        % APPswe/PSEN1dE9
'GSE131659/C57BL6/Brain/Male/8'...          % APPswe/PSEN1dE9
'GSE104424/line85/Hippocampus/Male/10'...   % APPswe/PSEN1dE9
'GSE104424/line85/Hippocampus/Male/6'...    % APPswe/PSEN1dE9
'GSE93678/line85/Hippocampus/Female/13'...  % APPswe/PSEN1dE9
'GSE80465/Forebrain/Male/2'...              % APP23
'GSE80465/Forebrain/Male/6'...              % APP23
'GSE80465/Forebrain/Male/18'...             % APP23
'GSE80465/Forebrain/Male/24'...             % APP23
'GSE104381/Cortex/Male/8'...                % Trem2 KO (KOMP) x APPPS1
'GSE104381/Hippocampus/Male/4'...           % Trem2 KO (KOMP) x APPPS1
'GSE104381/Hippocampus/Male/8'...           % Trem2 KO (KOMP) x APPPS1
};


% specify list of DE data files to analyze
data_types = {...
'GSE63943_APPPS1_HippocampalCA1MixedNA10'...
'GSE100070_APPPS1_DentateGyrusMixedMale4'...
'GSE110741_APPPS1_HippocampusMixedMale4'...
'GSE110741_APPPS1_HippocampusMixedMale10'...
'GSE131659_APPPS1C57BL6_BrainMixedFemale8'...
'GSE131659_APPPS1C57BL6_BrainMixedMale8'...
'GSE104424_APPPS1line85_HippocampusMixedMale10'...
'GSE104424_APPPS1line85_HippocampusMixedMale6'...
'GSE93678_APPPS1line85_HippocampusMixedFemale13'...
'GSE80465_APP23_ForebrainMixedMale2'...
'GSE80465_APP23_ForebrainMixedMale6'...
'GSE80465_APP23_ForebrainMixedMale18'...
'GSE80465_APP23_ForebrainMixedMale24'...
'GSE104381_Trem2KOKOMPAPPPS1_CortexMixedMale8'...
'GSE104381_Trem2KOKOMPAPPPS1_HippocampusMixedMale4'...
'GSE104381_Trem2KOKOMPAPPPS1_HippocampusMixedMale8'...
};
fnames = strcat(data_types, '.txt');


% load and analyze each dataset
for i = 1:numel(fnames)
    DEdata = readtable(fnames{i});
    geneIDs   = DEdata.Var1;     % gene IDs
    log2FC    = DEdata.Var3;     % log2(fold-changes) (desease vs. control)
    pvals     = DEdata.Var6;     % fold-change significance (p-values)

    % reporter metabolite
    GSAresReMet{i} = geneSetAnalysis(geneIDs, pvals, log2FC, gsc_met, ...
                           'method', 'Reporter', 'nPerm', 50000, 'gsSizeLim', [20, 200]);
    GSAresReMet{i}.Properties.Description = data_accessions{i};  

    % subsystems
    GSAresSubsys{i} = geneSetAnalysis(geneIDs, pvals, log2FC, gsc_subSys, ...
                        'method', 'Reporter', 'nPerm', 50000, 'gsSizeLim', [20, 200]);
    GSAresSubsys{i}.Properties.Description = data_accessions{i};  

    % combined KEGG, Reactome, Hallmark
    GSAresCombine{i} = geneSetAnalysis(geneIDs, pvals, log2FC, gscCombine, 'method', 'Wilcoxon', 'nPerm', 50000, 'gsSizeLim', [20, Inf]);
    GSAresCombine{i}.Properties.Description = data_accessions{i};  

end
save('geneSetAnalysisResults.mat','GSAresCombine','GSAresReMet','GSAresSubsys');


%% Fig S5, using MiSGD

GSAheatmap(GSAresCombine, 'filterMethod', 'top each', 'filterThreshold', 5, 'colorMax', 5);
set(gca,'FontName','Arial','FontSize',14,'FontWeight','Bold');
xtickangle(45)


%% GEM-based analysis

% Fig S6, Subsystem
GSAheatmap(GSAresSubsys, 'filterMethod', 'pval', 'filterThreshold', 0.05, 'colorMax', 5);
set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold');
xtickangle(45)


% Fig 3D, reporter metabolites
GSAheatmap(GSAresReMet, 'filterMethod', 'top each', 'filterThreshold', 8, 'colorMax', 5);
%GSAheatmap(GSAresReMet, 'filterMethod', 'pval', 'filterThreshold', 0.0001, 'colorMax', 5);
set(gca,'FontName','Arial','FontSize',14,'FontWeight','Bold');
xtickangle(45)


