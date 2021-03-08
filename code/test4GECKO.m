%
% code for generating a test model for GECKO automator
% 2020-3-8
%

% load Human-GEM v1.6 template
ihuman = importYaml('Human-GEM.yml');

% re-generate mouse GEM
mouseOrthologPairs = extractAllianceGenomeOrthologs('human2MouseOrthologs.json');

% add mouse-specific rxns and mets
metsToAdd = importTsvFile('mouseSpecificMets.tsv');
rxnsToAdd = importTsvFile('mouseSpecificRxns.tsv');

% try out new function updateAnimalGEM
[mouseGEM, speciesSpecNetwork, gapfillNetwork]=updateAnimalGEM(mouseOrthologPairs,rxnsToAdd,metsToAdd,'Mouse-GEM');

% fix grRules
[newRules, skipped] = simplifyGrRules(mouseGEM.grRules,true);
%grRule #5011 was too long (1397 characters), and therefore skipped.
%grRule #5093 was too long (2684 characters), and therefore skipped.
%grRule #5095 was too long (1628 characters), and therefore skipped.
%Failed to expand grRule #8181 due to excess length. Rule will remain in simplified form.

% update relevant fields
[genes,rxnGeneMat]  = getGenesFromGrRules(newRules);
mouseGEM.grRules    = newRules;
mouseGEM.genes      = genes;
mouseGEM.rxnGeneMat = rxnGeneMat;

% save model
save('../model/Mouse-GEM.mat','mouseGEM');

