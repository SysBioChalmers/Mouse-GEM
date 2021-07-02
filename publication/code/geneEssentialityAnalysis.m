% script for evaluating gene essentiality predictions by mouse & rat GEMs

% SPECIFY MINIMUM THRESHOLD BIOMASS FLUX
BMthresh = 1;  % default = 1


%% Load gene essentiality data

% Mouse
f = fopen('../essentiality_data/OGEE_essentiality_data_mouse.txt');
temp = textscan(f, '%s%s');
mouseData = [temp{1}, temp{2}];
fclose(f);

% Fly
f = fopen('../essentiality_data/OGEE_essentiality_data_fly.txt');
temp = textscan(f, '%s%s');
flyData = [temp{1}, temp{2}];
fclose(f);

% Worm
f = fopen('../essentiality_data/OGEE_essentiality_data_worm.txt');
temp = textscan(f, '%s%s');
wormData = [temp{1}, temp{2}];
fclose(f);




%% MouseGEM

% load the model
tmp = load('../GEMs/Mouse-GEM.mat');
model = tmp.mouseGEM;

% add boundary metabolites
model = addBoundaryMets(model);

% load task structure
taskStruct = parseTaskList('../metabolic_tasks/metabolicTasks_Biomass_MA-GEMs.xlsx');
taskStruct.LBout(2) = BMthresh;  % set biomass production requirement

% check biomass task
checkTasks(model, [], true, false, false, taskStruct);

% run essentiality analysis
[metrics_mouseGEM, modelEssential_mouseGEM] = ...
    evalGeneEssentialityPred(model, mouseData, taskStruct);

% RESULTS:
%     sensitivity: 0.0628
%     specificity: 0.9950
%        accuracy: 0.5180
%              F1: 0.1176
%             MCC: 0.1582
%        p_enrich: 1.0156e-14


%% iMM1865

% load the model
tmp = load('../GEMs/iMM1865.mat');
model = tmp.iMM1865;

% load task structure
taskStruct = parseTaskList('../metabolic_tasks/metabolicTasks_Biomass_iMM1865.xlsx');
taskStruct.LBout(2) = BMthresh;  % set biomass production requirement

% rename peroxisome compartment from [x] to [p]
model.comps(strcmp(model.comps, 'x')) = {'p'};
model.mets = regexprep(model.mets, '_x$', '_p');

% add boundary metabolites
model = addBoundaryMets(model);

% add biomass metabolites
metsToAdd = {};
metsToAdd.mets = {'temp_biomass_c';'temp_biomass_e';'temp_biomass_x'};
metsToAdd.metNames = {'biomass';'biomass';'biomass'};
metsToAdd.compartments = {'c';'e';'x'};
metsToAdd.unconstrained = [false; false; true];
model = addMets(model,metsToAdd);
        
% add biomass transport rxns
rxnsToAdd = {};
rxnsToAdd.rxns = {'biomass_export'; 'biomass_exchange'};
rxnsToAdd.equations = {'biomass[c] => biomass[e]'; 'biomass[e] <=> biomass[x]'};
model = addRxns(model, rxnsToAdd, 3);
        
% set biomass rxn to produce biomass[c]
bm_met_ind = getIndexes(model, 'biomass[c]', 'metcomps');
bm_rxn_ind = find(model.c > 0);
model.S(bm_met_ind, bm_rxn_ind) = 1;

% get gene ID association information (entrez to symbol)
f = fopen('/Users/jonrob/Documents/PostDoc/HMA_Sandbox/mouseRatGEMs/data/MGI_ID_mapping.tsv');
gene_data = textscan(f, '%s %s %s %s', 'Delimiter', '\t', 'Headerlines', 2);
fclose(f);

% organize columns and remove duplicate rows
gene_data = [gene_data{3}, gene_data{2}];  % 3=entrez, 2=symbol
[~,gene_data_num] = ismember(gene_data, gene_data);
[~,uniq_ind] = unique(gene_data_num, 'rows');
gene_data = gene_data(uniq_ind,:);

% convert model grRules, genes, and rxnGeneMat to gene symbols
[model.grRules, model.genes, model.rxnGeneMat] = translateGrRules(model.grRules, gene_data);

% check biomass task
checkTasks(model, [], true, false, false, taskStruct);

% run essentiality analysis
[metrics_iMM1865, modelEssential_iMM1865] = ...
    evalGeneEssentialityPred(model, mouseData, taskStruct);

% RESULTS
%     sensitivity: 0.0207
%     specificity: 0.9982
%        accuracy: 0.5162
%              F1: 0.0405
%             MCC: 0.0899
%        p_enrich: 0.0026



%% MMR

% load the model
tmp = load('../GEMs/MMR.mat');
model = tmp.MMR;

% % adjust exchange reactions so that they are all facing the same direction
% % ([inside] <==> [outside]) i.e., negative flux implies consumption
% S = full(model.S);
% [~,exch_ind] = ismember(getExchangeRxns(model), model.rxns);
% backward_ind = exch_ind(any(S(:, exch_ind) > 0, 1));
% model.S(:, backward_ind) = -model.S(:, backward_ind);
% lb = model.lb(backward_ind);
% ub = model.ub(backward_ind);
% model.lb(backward_ind) = -ub;
% model.ub(backward_ind) = -lb;

% load task structure
taskStruct = parseTaskList('../metabolic_tasks/metabolicTasks_Biomass_MMR.xlsx');
taskStruct.LBout(2) = BMthresh;  % set biomass production requirement

% add boundary metabolites
model = addBoundaryMets(model);

% add biomass metabolites
metsToAdd = {};
metsToAdd.mets = {'temp_biomass_c';'temp_biomass_s';'temp_biomass_x'};
metsToAdd.metNames = {'biomass';'biomass';'biomass'};
metsToAdd.compartments = {'c';'s';'x'};
metsToAdd.unconstrained = [false; false; true];
model = addMets(model,metsToAdd);
        
% add biomass transport rxns
rxnsToAdd = {};
rxnsToAdd.rxns = {'biomass_export'; 'biomass_exchange'};
rxnsToAdd.equations = {'biomass[c] => biomass[s]'; 'biomass[s] <=> biomass[x]'};
model = addRxns(model, rxnsToAdd, 3);

% set biomass rxn to produce biomass[c]
bm_met_ind = getIndexes(model, 'biomass[c]', 'metcomps');
bm_rxn_ind = getIndexes(model, 'biomass_components', 'rxns');
model.S(bm_met_ind, bm_rxn_ind) = 1;

% check biomass task
checkTasks(model, [], true, false, false, taskStruct);

% run essentiality analysis
[metrics_MMR, modelEssential_MMR] = ...
    evalGeneEssentialityPred(model, mouseData, taskStruct);

% RESULTS
%     sensitivity: 0.0044
%     specificity: 0.9990
%        accuracy: 0.4789
%              F1: 0.0088
%             MCC: 0.0328
%        p_enrich: 0.1321



%% Fly-GEM

% load the model
tmp = load('../GEMs/Fruitfly-GEM.mat');
model = tmp.fruitflyGEM;

% add boundary metabolites
model = addBoundaryMets(model);

% load task structure
taskStruct = parseTaskList('../metabolic_tasks/metabolicTasks_Biomass_MA-GEMs.xlsx');
taskStruct.LBout(2) = BMthresh;  % set biomass production requirement

% check biomass task
checkTasks(model, [], true, false, false, taskStruct);


% need to convert the fly gene IDs from the essentiality data into symbols
f = fopen('../ID_conversion/gene_ID_conversion_fly.txt');
temp = textscan(f, '%s%s');
flyGeneID2sym = [temp{1}, temp{2}];
fclose(f);

[has_match, ind] = ismember(flyData(:,1), flyGeneID2sym(:,1));
flyData = flyData(has_match, :);
flyData(:,1) = flyGeneID2sym(ind(has_match), 2);
[~, uniq_ind] = unique(flyData(:,1));
flyData = flyData(uniq_ind, :);

% gene symbols in the GEM contain greek symbols, but the translation file
% contains latin spellings
greek2latin = {'α', 'alpha'
               'β', 'beta'
               'γ', 'gamma'
               'δ', 'delta'
               'ε', 'epsilon'
               'ζ', 'zeta'
               'η', 'eta'
               'θ', 'theta'
               'ι', 'iota'};

% safer to translate greek to latin because phrases like "eta" are prone to
% accidental matches.
for i = 1:size(greek2latin, 1)
    model.genes = regexprep(model.genes, greek2latin{i,1}, greek2latin{i,2});
    model.grRules = regexprep(model.grRules, greek2latin{i,1}, greek2latin{i,2});
end

% run essentiality analysis
[metrics_flyGEM, modelEssential_flyGEM] = ...
    evalGeneEssentialityPred(model, flyData, taskStruct);

% RESULTS
%     sensitivity: 0.0909
%     specificity: 0.9278
%        accuracy: 0.9086
%              F1: 0.0437
%             MCC: 0.0108
%         p_hyper: 0.3980


%% BMID000000141998 (fruitfly)

% load the model
tmp = load('../GEMs/BMID000000141998_fly.mat');
model = tmp.model;

% extract metabolite compartments from the met IDs
metCompsAbbrev = regexprep(model.mets, '^[^\[]+\[|\]$', '');
[~, model.metComps] = ismember(metCompsAbbrev, model.comps);

% adjust exchange reactions so that they are all facing the same direction
% ([inside] <==> [outside]) i.e., negative flux implies consumption
S = full(model.S);
[~, exch_ind] = ismember(getExchangeRxns(model), model.rxns);
backward_ind = exch_ind(any(S(:, exch_ind) > 0, 1));
model.S(:, backward_ind) = -model.S(:, backward_ind);
lb = model.lb(backward_ind);
ub = model.ub(backward_ind);
model.lb(backward_ind) = -ub;
model.ub(backward_ind) = -lb;

% add reversibility field
model.rev = (model.lb < 0) & (model.ub > 0);


% modify the gene IDs to remove invalid pieces
genes = model.genes;
for i = 1:numel(genes)
    gene_parts = strsplit(genes{i}, '_');
    genes(i) = gene_parts(2);
end

% convert flybase identifiers to symbols
[has_match, ind] = ismember(lower(genes), lower(flyGeneID2sym(:,1)));
genes(has_match) = flyGeneID2sym(ind(has_match), 2);

% genes contain another type of identifier as well (CG####)
f = fopen('../ID_conversion/gene_ID_conversion_flyCG.txt');
temp = textscan(f, '%s%s');
flyGeneCG2sym = [temp{1}, temp{2}];
fclose(f);
[has_match, ind] = ismember(lower(genes), lower(flyGeneCG2sym(:,1)));
genes(has_match) = flyGeneCG2sym(ind(has_match), 2);

% generate and convert grRules
model = generateGrRules(model);
model.grRules = translateGrRules(model.grRules, [model.genes, genes]);

% update gene and rxnGeneMat fields
[model.genes, model.rxnGeneMat] = getGenesFromGrRules(model.grRules);
model = rmfield(model, 'rules');  % this field will be regenerated later


% add boundary metabolites
model = addBoundaryMets(model);

% load task structure
taskStruct = parseTaskList('../metabolic_tasks/metabolicTasks_Biomass_BMID000000141998.xlsx');
taskStruct.LBequ = BMthresh;  % set biomass production requirement

% check biomass task
checkTasks(model, [], true, false, false, taskStruct);


% run essentiality analysis
[metrics_BMID000000141998, modelEssential_BMID000000141998] = ...
    evalGeneEssentialityPred(model, flyData, taskStruct);

% RESULTS
%     sensitivity: 0
%     specificity: 1
%        accuracy: 0.9780
%              F1: 0
%             MCC: NaN
%         p_hyper: 1
%
% Note: No genes were predicted as essential, thus the NaN MCC



%% FlySilico

% load the model
tmp = load('../GEMs/FlySilico_v1_corrected.mat');
model = tmp.model;

% extract metabolite compartments from the met IDs
metCompsAbbrev = regexprep(model.mets, '^[^\[]+\[|\]$', '');
[~, model.metComps] = ismember(metCompsAbbrev, model.comps);

% add reversibility field
model.rev = (model.lb < 0) & (model.ub > 0);

% convert genes to symbols
genes = model.genes;
f = fopen('../ID_conversion/gene_ID_conversion_flyCG.txt');
temp = textscan(f, '%s%s');
flyGeneCG2sym = [temp{1}, temp{2}];
fclose(f);
[has_match, ind] = ismember(lower(genes), lower(flyGeneCG2sym(:,1)));
genes(has_match) = flyGeneCG2sym(ind(has_match), 2);

% generate and convert grRules
model = generateGrRules(model);
model.grRules = translateGrRules(model.grRules, [model.genes, genes]);

% update gene and rxnGeneMat fields
[model.genes, model.rxnGeneMat] = getGenesFromGrRules(model.grRules);
model = rmfield(model, 'rules');  % this field will be regenerated later


% add boundary metabolites
model = addBoundaryMets(model);

% load task structure
taskStruct = parseTaskList('../metabolic_tasks/metabolicTasks_Biomass_FlySilico.xlsx');
taskStruct.LBequ = BMthresh;  % set biomass production requirement

% check biomass task
checkTasks(model, [], true, false, false, taskStruct);


% run essentiality analysis
[metrics_FlySilico, modelEssential_FlySilico] = ...
    evalGeneEssentialityPred(model, flyData, taskStruct);

% RESULTS
%     sensitivity: 0
%     specificity: 0.9796
%        accuracy: 0.9562
%              F1: 0
%             MCC: -0.0223
%         p_hyper: 1



%% Worm-GEM

% load the model
tmp = load('../GEMs/Worm-GEM.mat');
model = tmp.wormGEM;

% add boundary metabolites
model = addBoundaryMets(model);

% load task structure
taskStruct = parseTaskList('../metabolic_tasks/metabolicTasks_Biomass_MA-GEMs.xlsx');
taskStruct.LBout(2) = BMthresh;  % set biomass production requirement

% check biomass task
checkTasks(model, [], true, false, false, taskStruct);


% need to convert the worm gene IDs from the essentiality data into symbols
f = fopen('../ID_conversion/gene_ID_conversion_worm.txt');
temp = textscan(f, '%s%s%s');
wormGeneID2sym = [temp{1}, temp{2}, temp{3}];
fclose(f);

[has_match, ind] = ismember(wormData(:,1), wormGeneID2sym(:,1));
wormData = wormData(has_match, :);
wormData(:,1) = wormGeneID2sym(ind(has_match), 2);


% run essentiality analysis
[metrics_wormGEM, modelEssential_wormGEM] = ...
    evalGeneEssentialityPred(model, wormData, taskStruct);

% RESULTS
%     sensitivity: 0.2340
%     specificity: 0.9559
%        accuracy: 0.8544
%              F1: 0.3113
%             MCC: 0.2574
%        p_enrich: 4.9682e-12



%% iCEL1314 (worm)

% load the model
tmp = load('../GEMs/iCEL1314.mat');
model = tmp.model;

% extract metabolite compartments from the met IDs
metCompsAbbrev = regexprep(model.mets, '^[^\[]+\[|\]$', '');
model.comps = unique(metCompsAbbrev);
[~, model.metComps] = ismember(metCompsAbbrev, model.comps);

% add compartment information to the model
comp_abbrevs = {'e', 'extracellular'
                'c', 'cytosol'
                'm', 'mitochondria'};
[~,ind] = ismember(model.comps, comp_abbrevs(:,1));
model.compNames = comp_abbrevs(ind, 2);

% add reversibility field
model.rev = (model.lb < 0) & (model.ub > 0);

% add boundary metabolites
model = addBoundaryMets(model);


% load task structure
taskStruct = parseTaskList('../metabolic_tasks/metabolicTasks_Biomass_iCEL.xlsx');
taskStruct.LBequ = BMthresh;  % set biomass production requirement

% check biomass task
checkTasks(model, [], true, false, false, taskStruct);


% run essentiality analysis
[metrics_iCEL1314, modelEssential_iCEL1314] = ...
    evalGeneEssentialityPred(model, wormData, taskStruct);

% RESULTS
%     sensitivity: 0.2549
%     specificity: 0.9603
%        accuracy: 0.8584
%              F1: 0.3421
%             MCC: 0.2949
%         p_hyper: 5.2562e-11



%% iCEL1273 (worm)

% load the model
tmp = load('../GEMs/iCEL1273.mat');
model = tmp.model;

% extract metabolite compartments from the met IDs
metCompsNames = regexprep(model.mets, '^[^\[]+\[|\]$', '');
metCompsAbbrev = lower(cellstr(cellfun(@(x) x(1), metCompsNames)));
model.comps = unique(metCompsAbbrev);
[~, model.metComps] = ismember(metCompsAbbrev, model.comps);

% replace compartment names with abbreviations in met IDs
model.mets = regexprep(model.mets, '\[\w+\]$', '');
model.mets = strcat(model.mets, '[', metCompsAbbrev, ']');

% add compartment information to the model
model.compNames = unique(metCompsNames);

% add reversibility field
model.rev = (model.lb < 0) & (model.ub > 0);

% add boundary metabolites
model = addBoundaryMets(model);

% add grRules and rxnGeneMat fields
model = generateGrRules(model);
[~, model.rxnGeneMat] = getGenesFromGrRules(model.grRules);
model = rmfield(model, 'rules');  % this field will be regenerated later


% load task structure
taskStruct = parseTaskList('../metabolic_tasks/metabolicTasks_Biomass_iCEL.xlsx');
taskStruct.LBequ = BMthresh;  % set biomass production requirement

% check biomass task
checkTasks(model, [], true, false, false, taskStruct);


% run essentiality analysis
[metrics_iCEL1273, modelEssential_iCEL1273] = ...
    evalGeneEssentialityPred(model, wormData, taskStruct);

% RESULTS
%     sensitivity: 0.0753
%     specificity: 0.9859
%        accuracy: 0.8578
%              F1: 0.1296
%             MCC: 0.1429
%         p_hyper: 0.0022



