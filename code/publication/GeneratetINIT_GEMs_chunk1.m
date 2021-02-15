%
% code for generation of tissue- and cell type-specific GEMs, chunk1
%


%% Add path of RAVEN and Human-GEM for loading functions
setRavenSolver('gurobi');


% Load processed TPM
d = load('../data/Models/TPMsPerCond_chunk1.mat');

f = readtable('../data/Models/TPMsPerCondGenes.txt', 'ReadVariableNames',false, 'ReadRowNames', false, 'Delimiter', '\t');
genes = strtrim(table2cell(f(:, 1)));

data = d.TPMs;
sampleIds = d.sampleNames;
data(isnan(data)) = 0; % for now, treat NaNs as 0.


%check that the data columns sum up to a number reasonably close to 10^6:
mean(sum(data,1))

arrayData = preData;
arrayData.levels = data;
arrayData.genes = genes;
arrayData.tissues = sampleIds;


%% Run tINIT algorithm

            
INIT_output = {};  % initialize INIT output structure

params.TimeLimit = 5000;%make sure tINIT gets enough time to run

for j = 1:length(arrayData.tissues)
    % run tINIT
    [init_model, metProduction, essentialRxnsForTasks, addedRxnsForTasks, deletedDeadEndRxns, deletedRxnsInINIT, taskReport] = ...
        getINITModel2(model,arrayData.tissues{j},[],[],arrayData,[],true,[],true,true,taskStruct, params);

    % organize results
    init_model.id = arrayData.tissues{j};
    INIT_output.model{j} = init_model;
    INIT_output.metProduction{j} = metProduction;
    INIT_output.essentialRxnsForTasks{j} = essentialRxnsForTasks;
    INIT_output.addedRxnsForTasks{j} = addedRxnsForTasks;
    INIT_output.deletedDeadEndRxns{j} = deletedDeadEndRxns;
    INIT_output.deletedRxnsInINIT{j} = deletedRxnsInINIT;
    INIT_output.taskReport{j} = taskReport;

    % save results each time in the loop if it fails
    save("models_chunk1.mat",'INIT_output');
end
   



