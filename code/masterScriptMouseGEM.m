%
% FILE NAME:    masterScriptMouseGEM.m
%
%
% PURPOSE: This script is for reconstruction of the Mouse-GEM, by using
%          the Human-GEM as template and taking into account mouse-specific
%          pathways/reactions.
%
%


%% Load Human-GEM as template
load('Human-GEM.mat');


% convert gene identifiers from Ensembl ids to gene symbols
[grRules,genes,rxnGeneMat] = translateGrRules(ihuman.grRules,'Name','ENSG');
ihuman.grRules    = grRules;
ihuman.genes      = genes;
ihuman.rxnGeneMat = rxnGeneMat;



%% Use MA reactions identifiers 

% load reaction annotaiton files
rxnAssoc = jsondecode(fileread('humanGEMRxnAssoc.JSON'));

%replace reaction identifiers with MA ids if available
ind = getNonEmptyList(rxnAssoc.rxnMAID);
ihuman.rxns(ind) = rxnAssoc.rxnMAID(ind);



%% Generate Mouse-GEM using Human-GEM as template

% get ortholog pairs from human to mouse
mouseOrthologPairs = extractAllianceGenomeOrthologs('human2MouseOrthologs.json');
mouseGEM = getModelFromOrthology(ihuman, mouseOrthologPairs);



%% Incorporate mouse-specific reactions

% get metabolic networks based on the KEGG annoation using RAVEN function
KEGG_human=getKEGGModelForOrganism('hsa');
KEGG_mouse=getKEGGModelForOrganism('mmu');

% remove reactions shared with human
MouseSpecificRxns=setdiff(KEGG_mouse.rxns,KEGG_human.rxns);

% remove reactions included in Human-GEM
MouseSpecificRxns=setdiff(MouseSpecificRxns,rxnAssoc.rxnKEGGID);

% get species-specific network for manual inspection and then
% organize species-specific pathways into two tsv files:
mouseSpecificNetwork=removeReactions(KEGG_mouse, setdiff(KEGG_mouse.rxns,MouseSpecificRxns), true, true, true);

% "mouseSpecificMets.tsv" contains new metabolites
metsToAdd = importTsvFile('mouseSpecificMets.tsv');

% "mouseSpecificRxns.tsv" contains new reactions
rxnsToAdd = importTsvFile('mouseSpecificRxns.tsv');
rxnsToAdd.subSystems = cellfun(@(s) {{s}}, rxnsToAdd.subSystems);

% integrate mouse-specific metabolic network
[mouseGEM, modelChanges] = addMetabolicNetwork(mouseGEM, rxnsToAdd, metsToAdd);


%% Gap-filling for biomass formation
[mouseGEM, gapfillNetwork]=gapfill4EssentialTasks(mouseGEM,ihuman);
% Added 0 reactions for gap-filling


%% Save the model into mat, yml, and xml

mouseGEM.id = 'Mouse-GEM';
save('../model/Mouse-GEM.mat', 'mouseGEM');
writeHumanYaml(mouseGEM, '../model/Mouse-GEM.yml');
exportModel(mouseGEM, '../model/Mouse-GEM.xml');

