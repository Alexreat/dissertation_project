%[text] ## Initial setup Path COBRA Toolbox and Gurobi solver

addpath('/Users/valentinreateguirangel/Documents/MATLAB/libSBML-5')
savepath
%%
addpath(genpath('/Users/valentinreateguirangel/Documents/MATLAB/MyProject/COBRA_project/cobratoolbox'));
%%
%  COBRA initialization 
initCobraToolbox(false) %[output:670d9238]
%%
% Point COBRA at Gurobi
changeCobraSolver('gurobi', 'LP'); %[output:00c001d0]
%%
%[text] ## Load the Human-GEM model
%%
% Read the model 
model = readCbModel('/Users/valentinreateguirangel/Downloads/Human-GEM-1.19.0/model/Human-GEM.mat'); %[output:91542363]
%%
% Show a quick summary 
if isfield(model,'description')
    desc = model.description;
else
    desc = '--no description--';
end
fprintf('Loaded model:  %s\n', desc); %[output:6fb254ff]
fprintf('Reactions:     %d\nMetabolites:   %d\nGenes:         %d\n',... %[output:group:984bb587] %[output:2f41bf2b]
        numel(model.rxns), numel(model.mets), numel(model.genes)); %[output:group:984bb587] %[output:2f41bf2b]
%%
%[text] ## Load Expression vector for tumor
exprFile = "/Users/valentinreateguirangel/dissertation_project/notebooks/Zscore_vector_tumor.txt";
%%
% reading the table 
exprTbl = readtable(exprFile,          ...
                    'FileType','text', ...
                    'Delimiter','\t',  ...
                    'ReadVariableNames',false);
%%
% giving columns explicit names 
exprTbl.Properties.VariableNames = {'gene','expr'};
%%
% Table display 
disp(exprTbl(1:min(5,height(exprTbl)), :));        % first 5 rows %[output:8131d854]
fprintf('Rows in table: %d genes\n', height(exprTbl)); %[output:6b4e0370]
fprintf('Expression range: %.3g  –  %.3g\n', min(exprTbl.expr), max(exprTbl.expr)); %[output:4204a27f]
%%
%[text] ## Map expression onto model gem
% Build a fast lookup table:  key = gene ID,  value = expression level
exprMap = containers.Map(exprTbl.gene, exprTbl.expr);
%%
% Pre-allocate vector (same length & order as model.genes)
exprValues = nan(numel(model.genes),1);
%%
% Fill it where we have a match
for k = 1:numel(model.genes)
    g = model.genes{k};
    if isKey(exprMap, g)
        exprValues(k) = exprMap(g);
    end
end
%%
% Quick stats – how many genes got an expression value?
numMapped = sum(~isnan(exprValues));
fprintf('Mapped expression for %d / %d model genes (%.1f%%)\n', ... %[output:group:733f2b8b] %[output:2c341a1a]
        numMapped, numel(model.genes), 100*numMapped/numel(model.genes)); %[output:group:733f2b8b] %[output:2c341a1a]
%%
% Peek at the first few mapped entries
ix = find(~isnan(exprValues));
preview = table(model.genes(ix(1:min(5,end)))', exprValues(ix(1:min(5,end)))', ...
                'VariableNames', {'gene','expr'});
disp(preview); %[output:677d0e34]
%%
%[text] ## Map gene → reaction expression
%%
% Build the *expressionData* structure the function expects:
expressionData           = struct();
expressionData.gene      = exprTbl.gene;   % list of gene IDs
expressionData.value     = exprTbl.expr;   % their expression values
%  ➜  (If you ever have p-values, you'd add expressionData.sig)

% Map gene expression to reactions
%    • The function returns NaN for reactions with no mapped gene.
[expressionRxns, ~, ~] = mapExpressionToReactions(model, expressionData);

% Replace NaNs with –1 (iMAT expects –1 for "no data")
expressionRxns(isnan(expressionRxns)) = -100;

% Quick sanity check
mapped = sum(expressionRxns ~= -100);
fprintf('Mapped expression for %d / %d reactions (%.1f%%)\n', ... %[output:group:88944c52] %[output:17115b4b]
        mapped, numel(expressionRxns), 100*mapped/numel(expressionRxns)); %[output:group:88944c52] %[output:17115b4b]
%%
%[text] ## Build tumor specific model with iMat
%lowThr  = 2;     % ≤ lowThr  → force reaction *off*
%highThr = 4;     % ≥ highThr → force reaction *on*

loDummy =  -0.5;   % safely below all data
hiDummy =  0.5;   % safely above all data
%%
% Pack everything in an *options* structure ( API for COBRA )
options                 = struct();
options.solver          = 'iMAT';
options.expressionRxns  = expressionRxns;
options.threshold_lb    = loDummy;
options.threshold_ub    = hiDummy;
options.task_specific = {{'MAR13082'}, [1], [0.001]};
options.core = {'MAR10065', 'MAR10062', 'MAR10063', 'MAR10064', 'MAR13082'};
%%
%[text] ## Setting up Nutrients Medium
%%
% Look up nutrient name 
getHumanName = @(metID) model.metNames{find(strcmp(model.mets, metID))};
printRxnFormula(model, 'MAR09045') %[output:1d37b27d] %[output:936ceba3]
getHumanName('MAM03089e')  %[output:92c5011b]
%%
% Nutrients list 
model = changeRxnBounds(model, 'MAR09034', -10, 'l');   % glucose
model = changeRxnBounds(model, 'MAR09036', -10, 'l');   % glutamine
model = changeRxnBounds(model, 'MAR09035', -10, 'l');   % glutamate
model = changeRxnBounds(model, 'MAR09037', -10, 'l');   % aspartate

% Oxygen
model = changeRxnBounds(model, 'MAR09048', -1000, 'l'); % oxygen

% Inorganic essentials
model = changeRxnBounds(model, 'MAR09100', -1000, 'l'); % H2O
model = changeRxnBounds(model, 'MAR09101', -1000, 'l'); % NH4+
model = changeRxnBounds(model, 'MAR09102', -1000, 'l'); % phosphate
model = changeRxnBounds(model, 'MAR09103', -1000, 'l'); % sulfate
model = changeRxnBounds(model, 'MAR09104', -1000, 'l'); % bicarbonate

% Lipids & precursors
model = changeRxnBounds(model, 'MAR09091', -1, 'l');    % choline
model = changeRxnBounds(model, 'MAR09092', -1, 'l');    % ethanolamine
model = changeRxnBounds(model, 'MAR09093', -1, 'l');    % glycerol
model = changeRxnBounds(model, 'MAR09094', -1, 'l');    % palmitate
model = changeRxnBounds(model, 'MAR09095', -1, 'l');    % acetate

% Vitamins & cofactors
model = changeRxnBounds(model, 'MAR09120', -1, 'l');    % riboflavin
model = changeRxnBounds(model, 'MAR09121', -1, 'l');    % thiamine
model = changeRxnBounds(model, 'MAR09122', -1, 'l');    % folate
model = changeRxnBounds(model, 'MAR09123', -1, 'l');    % nicotinamide (NAD precursor)
model = changeRxnBounds(model, 'MAR09050', -1, 'l');    % histidine

% Optional: iron, magnesium, etc.
model = changeRxnBounds(model, 'MAR09110', -1, 'l');    % Fe2+
model = changeRxnBounds(model, 'MAR09111', -1, 'l');    % Mg2+
%%
% Runing iMat
tumorModel = createTissueSpecificModel(model, options); %[output:015b5e74]
%%
% Summary
fprintf('\nTumour-specific model created!\n'); %[output:5afff8f6]
fprintf(' Reactions retained: %d  (of %d)\n', numel(tumorModel.rxns), numel(model.rxns)); %[output:86daeb16]
fprintf(' Metabolites:        %d\n',          numel(tumorModel.mets)); %[output:7abc6e55]
fprintf(' Genes:              %d (of %d)\n',  numel(tumorModel.genes), numel(model.genes)); %[output:0d9b3e0d]
%%
save('tumorModel.mat', 'tumorModel');
%%
% Test biomass production in tumor model
findRxnIDs(tumorModel, 'MAR13082')  % human biomass %[output:03b53197]
findRxnIDs(tumorModel, 'MAR09034')  % glucose %[output:66276739]
findRxnIDs(tumorModel, 'MAR09048')  % oxygen %[output:8c3ca194]
%%
printRxnFormula(tumorModel, 'MAR13082') %[output:640a26a0] %[output:7fc4a0ae]
tumorModel.rxnNames{find(strcmp(tumorModel.rxns, 'MAR13082'))} %[output:3d4aa83b]
%%
printRxnFormula(tumorModel, 'MAR03907') %[output:9e8ecbef] %[output:48757a1d]
tumorModel.rxnNames{find(strcmp(tumorModel.rxns, 'MAR03907'))} %[output:088a0dc5]
%%
printRxnFormula(model, 'MAR13082')  % on full GEM model %[output:8e722da4] %[output:743ae127]
model.rxnNames{find(strcmp(model.rxns, 'MAR13082'))} %[output:9af3f170]
%%
printRxnFormula(tumorModel, 'MAR03905')  % ATP maintenance %[output:1a6e3778] %[output:6c0cc305]
printRxnFormula(tumorModel, 'MAR08636')  % Pyruvate kinase %[output:80dc1aa3] %[output:9f7f533d]
printRxnFormula(tumorModel, 'MAR08709')  % Hexokinase %[output:06db00a9] %[output:4dac4f5b]
printRxnFormula(tumorModel, 'MAR08680')  % Glycolysis intermediate %[output:0543f74a] %[output:0ab1090c]
%%
printRxnFormula(tumorModel, 'MAR08709')  % Hexokinase %[output:7504ea89] %[output:3dc036ef]
printRxnFormula(tumorModel, 'MAR03905')  % ATP maintenance %[output:7fc5282a] %[output:086d36fd]
printRxnFormula(tumorModel, 'MAR08636')  % Pyruvate kinase %[output:0a81e5c2] %[output:3e9d4990]
%%
rxns = findRxnsFromMets(model, 'MAM10012c');
printRxnFormula(model, rxns) %[output:3a6b358b] %[output:51a8e257]
%[text] 

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":32.1}
%---
%[output:670d9238]
%   data: {"dataType":"text","outputData":{"text":"\n\n      _____   _____   _____   _____     _____     |\n     \/  ___| \/  _  \\ |  _  \\ |  _  \\   \/ ___ \\    |   COnstraint-Based Reconstruction and Analysis\n     | |     | | | | | |_| | | |_| |  | |___| |   |   The COBRA Toolbox - 2025\n     | |     | | | | |  _  { |  _  \/  |  ___  |   |\n     | |___  | |_| | | |_| | | | \\ \\  | |   | |   |   Documentation:\n     \\_____| \\_____\/ |_____\/ |_|  \\_\\ |_|   |_|   |   <a href=\"http:\/\/opencobra.github.io\/cobratoolbox\">http:\/\/opencobra.github.io\/cobratoolbox<\/a>\n                                                  | \n\n > Checking if git is installed ...  Done (version: 2.39.5).\n > Checking if the repository is tracked using git ...  Done.\n > Checking if curl is installed ...  Done.\n > Checking if remote can be reached ...  Done.\n > Initializing and updating submodules (this may take a while)... Done.\n > Adding all the files of The COBRA Toolbox ...  Done.\n > Define CB map output... set to svg.\n > TranslateSBML is installed and working properly.\n > Configuring solver environment variables ...\n   - [*---] ILOG_CPLEX_PATH: --> set this path manually after installing the solver ( see <a href=\"https:\/\/opencobra.github.io\/cobratoolbox\/docs\/solvers.html\">instructions<\/a> )\n   - [-*--] GUROBI_PATH: \/Library\/gurobi1202\/macos_universal2\/matlab\n   - [*---] TOMLAB_PATH: --> set this path manually after installing the solver ( see <a href=\"https:\/\/opencobra.github.io\/cobratoolbox\/docs\/solvers.html\">instructions<\/a> )\n   - [*---] MOSEK_PATH: --> set this path manually after installing the solver ( see <a href=\"https:\/\/opencobra.github.io\/cobratoolbox\/docs\/solvers.html\">instructions<\/a> )\n   Done.\n > Checking available solvers and solver interfaces ...Could not find installation of mosek, so it cannot be tested\nCould not find installation of tomlab_snopt, so it cannot be tested\n Done.\n > Setting default solvers ...Could not find installation of mosek, so it cannot be tested\nCould not find installation of mosek, so it cannot be tested\n Done.\n > Saving the MATLAB path ... Done.\n   - The MATLAB path was saved in the default location.\n\n > Summary of available solvers and solver interfaces\n\n\t\t\tSupport \t   LP \t MILP \t   QP \t MIQP \t  NLP \t   EP \t  CLP\n\t------------------------------------------------------------------------------\n\tdqqMinos     \tactive        \t    1 \t    - \t    1 \t    - \t    - \t    - \t    -\n\tglpk         \tactive        \t    1 \t    1 \t    - \t    - \t    - \t    - \t    -\n\tgurobi       \tactive        \t    1 \t    1 \t    1 \t    1 \t    - \t    - \t    -\n\tlp_solve     \tlegacy        \t    1 \t    - \t    - \t    - \t    - \t    - \t    -\n\tmatlab       \tactive        \t    1 \t    - \t    - \t    - \t    1 \t    - \t    -\n\tmosek        \tactive        \t    0 \t    - \t    0 \t    - \t    - \t    0 \t    0\n\tpdco         \tactive        \t    1 \t    - \t    1 \t    - \t    - \t    1 \t    -\n\tqpng         \tpassive       \t    - \t    - \t    1 \t    - \t    - \t    - \t    -\n\tquadMinos    \tactive        \t    1 \t    - \t    - \t    - \t    - \t    - \t    -\n\ttomlab_snopt \tpassive       \t    - \t    - \t    - \t    - \t    0 \t    - \t    -\n\t------------------------------------------------------------------------------\n\tTotal        \t-             \t    7 \t    2 \t    4 \t    1 \t    1 \t    1 \t    0\n\n + Legend: - = not applicable, 0 = solver not compatible or not installed, 1 = solver installed.\n\n\n > You can solve LP problems using: 'gurobi' - 'pdco' \n > You can solve MILP problems using: 'gurobi' \n > You can solve QP problems using: 'gurobi' - 'pdco' \n > You can solve MIQP problems using: 'gurobi' \n > You can solve NLP problems using: \n > You can solve EP problems using: 'pdco' \n > You can solve CLP problems using: \n\n> Checking for available updates ... skipped\nremoving: \/Users\/valentinreateguirangel\/Documents\/MATLAB\/MyProject\/COBRA_project\/cobratoolbox\/src\/analysis\/thermo\/componentContribution\/new\nremoving: \/Users\/valentinreateguirangel\/Documents\/MATLAB\/MyProject\/COBRA_project\/cobratoolbox\/src\/analysis\/thermo\/groupContribution\/new\nremoving: \/Users\/valentinreateguirangel\/Documents\/MATLAB\/MyProject\/COBRA_project\/cobratoolbox\/src\/analysis\/thermo\/inchi\/new\nremoving: \/Users\/valentinreateguirangel\/Documents\/MATLAB\/MyProject\/COBRA_project\/cobratoolbox\/src\/analysis\/thermo\/molFiles\/new\nremoving: \/Users\/valentinreateguirangel\/Documents\/MATLAB\/MyProject\/COBRA_project\/cobratoolbox\/src\/analysis\/thermo\/protons\/new\nremoving: \/Users\/valentinreateguirangel\/Documents\/MATLAB\/MyProject\/COBRA_project\/cobratoolbox\/src\/analysis\/thermo\/trainingModel\/new\n","truncated":false}}
%---
%[output:00c001d0]
%   data: {"dataType":"text","outputData":{"text":"\n > changeCobraSolver: Gurobi interface added to MATLAB path.\n","truncated":false}}
%---
%[output:91542363]
%   data: {"dataType":"text","outputData":{"text":"Each model.subSystems{x} has been changed to a character array.\n","truncated":false}}
%---
%[output:6fb254ff]
%   data: {"dataType":"text","outputData":{"text":"Loaded model:  Human-GEM.mat\n","truncated":false}}
%---
%[output:2f41bf2b]
%   data: {"dataType":"text","outputData":{"text":"Reactions:     12971\nMetabolites:   8455\nGenes:         2887\n","truncated":false}}
%---
%[output:8131d854]
%   data: {"dataType":"text","outputData":{"text":"           <strong>gene<\/strong>            <strong>expr<\/strong>\n    <strong>___________________<\/strong>    <strong>____<\/strong>\n\n    {'ENSG00000000003'}     1  \n    {'ENSG00000000005'}     0  \n    {'ENSG00000000419'}     1  \n    {'ENSG00000000457'}     1  \n    {'ENSG00000000460'}     0  \n\n","truncated":false}}
%---
%[output:6b4e0370]
%   data: {"dataType":"text","outputData":{"text":"Rows in table: 60660 genes\n","truncated":false}}
%---
%[output:4204a27f]
%   data: {"dataType":"text","outputData":{"text":"Expression range: -1  –  1\n","truncated":false}}
%---
%[output:2c341a1a]
%   data: {"dataType":"text","outputData":{"text":"Mapped expression for 2884 \/ 2887 model genes (99.9%)\n","truncated":false}}
%---
%[output:677d0e34]
%   data: {"dataType":"text","outputData":{"text":"                                                         <strong>gene<\/strong>                                                                  <strong>expr<\/strong>         \n    <strong>_______________________________________________________________________________________________________________<\/strong>    <strong>_____________________<\/strong>\n\n    {'ENSG00000000419'}    {'ENSG00000001036'}    {'ENSG00000001084'}    {'ENSG00000001630'}    {'ENSG00000002549'}    1    1    1    0    1\n\n","truncated":false}}
%---
%[output:17115b4b]
%   data: {"dataType":"text","outputData":{"text":"Mapped expression for 8040 \/ 12971 reactions (62.0%)\n","truncated":false}}
%---
%[output:1d37b27d]
%   data: {"dataType":"text","outputData":{"text":"MAR09045\tMAM03089e \t<=>\t\n","truncated":false}}
%---
%[output:936ceba3]
%   data: {"dataType":"textualVariable","outputData":{"header":"1×1 cell array","name":"ans","value":"    {'MAM03089e  <=> '}\n"}}
%---
%[output:92c5011b]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"'tryptophan'"}}
%---
%[output:015b5e74]
%   data: {"dataType":"text","outputData":{"text":"Set parameter Username\nSet parameter LicenseID to value 2674995\nSet parameter TimeLimit to value 7200\nSet parameter IntFeasTol to value 1e-09\nSet parameter MIPGap to value 1e-12\nSet parameter LogFile to value \"MILPlog\"\nAcademic license - for non-commercial use only - expires 2026-06-04\nGurobi Optimizer version 12.0.2 build v12.0.2rc0 (mac64[arm] - Darwin 24.5.0 24F74)\n\nCPU model: Apple M1 Pro\nThread count: 10 physical cores, 10 logical processors, using up to 10 threads\n\nNon-default parameters:\nTimeLimit  7200\nIntFeasTol  1e-09\nMIPGap  1e-12\nLogToConsole  0\n\nOptimize a model with 21897 rows, 26413 columns and 82326 nonzeros\nModel fingerprint: 0x1edfd9d6\nVariable types: 12971 continuous, 13442 integer (13442 binary)\nCoefficient statistics:\n  Matrix range     [1e-04, 8e+04]\n  Objective range  [1e+00, 1e+00]\n  Bounds range     [1e+00, 1e+03]\n  RHS range        [1e+03, 1e+03]\nPresolve removed 13812 rows and 13868 columns\nPresolve time: 0.10s\nPresolved: 8085 rows, 12545 columns, 42163 nonzeros\nVariable types: 6633 continuous, 5912 integer (5911 binary)\nFound heuristic solution: objective 3.0000000\nFound heuristic solution: objective 4.0000000\nFound heuristic solution: objective 5.0000000\nFound heuristic solution: objective 6.0000000\n\nRoot relaxation: objective 7.235142e+03, 7399 iterations, 0.24 seconds (0.35 work units)\n\n    Nodes    |    Current Node    |     Objective Bounds      |     Work\n Expl Unexpl |  Obj  Depth IntInf | Incumbent    BestBd   Gap | It\/Node Time\n\n     0     0 7235.14237    0 1438    6.00000 7235.14237      -     -    1s\nH    0     0                    5737.0000000 7235.14237  26.1%     -    1s\nH    0     0                    5740.0000000 7235.14237  26.0%     -    1s\n     0     0 5883.10979    0  397 5740.00000 5883.10979  2.49%     -    1s\nH    0     0                    5741.0000000 5883.10979  2.48%     -    1s\nH    0     0                    5773.0000000 5883.10979  1.91%     -    1s\n     0     0 5849.04306    0  362 5773.00000 5849.04306  1.32%     -    1s\n     0     0 5849.04306    0  360 5773.00000 5849.04306  1.32%     -    1s\n     0     0 5849.02316    0  272 5773.00000 5849.02316  1.32%     -    1s\nH    0     0                    5775.0000000 5849.02290  1.28%     -    2s\nH    0     0                    5788.0000000 5849.02290  1.05%     -    2s\n     0     0 5849.02256    0  272 5788.00000 5849.02256  1.05%     -    2s\nH    0     0                    5789.0000000 5849.02256  1.04%     -    2s\n     0     0 5849.01062    0  258 5789.00000 5849.01062  1.04%     -    2s\nH    0     0                    5790.0000000 5849.01062  1.02%     -    2s\n     0     0 5849.01062    0  258 5790.00000 5849.01062  1.02%     -    2s\nH    0     0                    5793.0000000 5849.01062  0.97%     -    2s\nH    0     0                    5794.0000000 5849.01062  0.95%     -    2s\nH    0     0                    5795.0000000 5849.01062  0.93%     -    2s\nH    0     0                    5796.0000000 5849.01062  0.91%     -    2s\nH    0     0                    5798.0000000 5849.01062  0.88%     -    2s\nH    0     0                    5799.0000000 5849.01062  0.86%     -    2s\nH    0     0                    5803.0000000 5849.01062  0.79%     -    2s\nH    0     0                    5804.0000000 5849.01062  0.78%     -    2s\nH    0     0                    5805.0000000 5849.00768  0.76%     -    2s\n     0     0 5849.00768    0  251 5805.00000 5849.00768  0.76%     -    2s\nH    0     0                    5806.0000000 5849.00768  0.74%     -    2s\nH    0     0                    5807.0000000 5849.00768  0.72%     -    2s\nH    0     0                    5808.0000000 5849.00768  0.71%     -    2s\n     0     0 5849.00768    0   27 5808.00000 5849.00768  0.71%     -    3s\nH    0     0                    5832.0000000 5849.00762  0.29%     -    3s\nH    0     2                    5835.0000000 5849.00762  0.24%     -    3s\n     0     2 5849.00762    0   23 5835.00000 5849.00762  0.24%     -    3s\nH    1     4                    5837.0000000 5849.00741  0.21%  13.0    3s\nH    3     8                    5839.0000000 5849.00735  0.17%   5.7    3s\nH    5     8                    5841.0000000 5849.00735  0.14%   6.6    3s\nH   17    26                    5842.0000000 5849.00644  0.12%   7.2    3s\nH   53    58                    5843.0000000 5849.00614  0.10%   8.3    3s\nH  120   148                    5848.0000000 5849.00452  0.02%   5.3    3s\n*  156   151              26    5849.0000000 5849.00411  0.00%   4.7    3s\n\nCutting planes:\n  Learned: 1352\n  Gomory: 28\n  Cover: 6\n  Implied bound: 6\n  MIR: 53\n  Relax-and-lift: 10\n\nExplored 273 nodes (23026 simplex iterations) in 3.54 seconds (4.53 work units)\nThread count was 10 (of 10 available processors)\n\nSolution count 10: 5849 5848 5843 ... 5808\n\nOptimal solution found (tolerance 1.00e-12)\nBest objective 5.849000000000e+03, best bound 5.849000000000e+03, gap 0.0000%\n","truncated":false}}
%---
%[output:5afff8f6]
%   data: {"dataType":"text","outputData":{"text":"\nTumour-specific model created!\n","truncated":false}}
%---
%[output:86daeb16]
%   data: {"dataType":"text","outputData":{"text":" Reactions retained: 9902  (of 12971)\n","truncated":false}}
%---
%[output:7abc6e55]
%   data: {"dataType":"text","outputData":{"text":" Metabolites:        6508\n","truncated":false}}
%---
%[output:0d9b3e0d]
%   data: {"dataType":"text","outputData":{"text":" Genes:              2887 (of 2887)\n","truncated":false}}
%---
%[output:03b53197]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"9787"}}
%---
%[output:66276739]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"6086"}}
%---
%[output:8c3ca194]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"6100"}}
%---
%[output:640a26a0]
%   data: {"dataType":"text","outputData":{"text":"MAR13082\t45 MAM01371c + 0.0267 MAM01721n + 45 MAM02040c + 0.1124 MAM02847c + 0.4062 MAM03161c + 0.0012 MAM10012c + 5.3375 MAM10013c + 0.2212 MAM10014c + 0.4835 MAM10015c \t->\t45 MAM01285c + 45 MAM02039c + 45 MAM02751c + MAM03970c \n","truncated":false}}
%---
%[output:7fc4a0ae]
%   data: {"dataType":"textualVariable","outputData":{"header":"1×1 cell array","name":"ans","value":"    {'45 MAM01371c + 0.0267 MAM01721n + 45 MAM02040c + 0.1124 MAM02847c + 0.4062 MAM03161c + 0.0012 MAM10012c + 5.3375 MAM10013c + 0.2212 MAM10014c + 0.4835 MAM10015c  -> 45 MAM01285c + 45 MAM02039c + 45 MAM02751c + MAM03970c '}\n"}}
%---
%[output:3d4aa83b]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"'Generic human cell biomass reaction'"}}
%---
%[output:9e8ecbef]
%   data: {"dataType":"text","outputData":{"text":"MAR03907\tMAM01796c + MAM02554c \t->\tMAM01249c + MAM02039c + MAM02555c \n","truncated":false}}
%---
%[output:48757a1d]
%   data: {"dataType":"textualVariable","outputData":{"header":"1×1 cell array","name":"ans","value":"    {'MAM01796c + MAM02554c  -> MAM01249c + MAM02039c + MAM02555c '}\n"}}
%---
%[output:088a0dc5]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"'Ethanol:NADP+ oxidoreductase'"}}
%---
%[output:8e722da4]
%   data: {"dataType":"text","outputData":{"text":"MAR13082\t45 MAM01371c + 0.0267 MAM01721n + 45 MAM02040c + 0.1124 MAM02847c + 0.4062 MAM03161c + 0.0012 MAM10012c + 5.3375 MAM10013c + 0.2212 MAM10014c + 0.4835 MAM10015c \t->\t45 MAM01285c + 45 MAM02039c + 45 MAM02751c + MAM03970c \n","truncated":false}}
%---
%[output:743ae127]
%   data: {"dataType":"textualVariable","outputData":{"header":"1×1 cell array","name":"ans","value":"    {'45 MAM01371c + 0.0267 MAM01721n + 45 MAM02040c + 0.1124 MAM02847c + 0.4062 MAM03161c + 0.0012 MAM10012c + 5.3375 MAM10013c + 0.2212 MAM10014c + 0.4835 MAM10015c  -> 45 MAM01285c + 45 MAM02039c + 45 MAM02751c + MAM03970c '}\n"}}
%---
%[output:9af3f170]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"'Generic human cell biomass reaction'"}}
%---
%[output:1a6e3778]
%   data: {"dataType":"text","outputData":{"text":"MAR03905\tMAM01796c + MAM02552c \t->\tMAM01249c + MAM02039c + MAM02553c \n","truncated":false}}
%---
%[output:6c0cc305]
%   data: {"dataType":"textualVariable","outputData":{"header":"1×1 cell array","name":"ans","value":"    {'MAM01796c + MAM02552c  -> MAM01249c + MAM02039c + MAM02553c '}\n"}}
%---
%[output:80dc1aa3]
%   data: {"dataType":"text","outputData":{"text":"MAR08636\tMAM01632c \t<=>\tMAM01632e \n","truncated":false}}
%---
%[output:9f7f533d]
%   data: {"dataType":"textualVariable","outputData":{"header":"1×1 cell array","name":"ans","value":"    {'MAM01632c  <=> MAM01632e '}\n"}}
%---
%[output:06db00a9]
%   data: {"dataType":"text","outputData":{"text":"MAR08709\tMAM00290c + MAM02039c + MAM02555c \t<=>\tMAM00291c + MAM02554c \n","truncated":false}}
%---
%[output:4dac4f5b]
%   data: {"dataType":"textualVariable","outputData":{"header":"1×1 cell array","name":"ans","value":"    {'MAM00290c + MAM02039c + MAM02555c  <=> MAM00291c + MAM02554c '}\n"}}
%---
%[output:0543f74a]
%   data: {"dataType":"text","outputData":{"text":"MAR08680\tMAM01614c \t<=>\tMAM01614m \n","truncated":false}}
%---
%[output:0ab1090c]
%   data: {"dataType":"textualVariable","outputData":{"header":"1×1 cell array","name":"ans","value":"    {'MAM01614c  <=> MAM01614m '}\n"}}
%---
%[output:7504ea89]
%   data: {"dataType":"text","outputData":{"text":"MAR08709\tMAM00290c + MAM02039c + MAM02555c \t<=>\tMAM00291c + MAM02554c \n","truncated":false}}
%---
%[output:3dc036ef]
%   data: {"dataType":"textualVariable","outputData":{"header":"1×1 cell array","name":"ans","value":"    {'MAM00290c + MAM02039c + MAM02555c  <=> MAM00291c + MAM02554c '}\n"}}
%---
%[output:7fc5282a]
%   data: {"dataType":"text","outputData":{"text":"MAR03905\tMAM01796c + MAM02552c \t->\tMAM01249c + MAM02039c + MAM02553c \n","truncated":false}}
%---
%[output:086d36fd]
%   data: {"dataType":"textualVariable","outputData":{"header":"1×1 cell array","name":"ans","value":"    {'MAM01796c + MAM02552c  -> MAM01249c + MAM02039c + MAM02553c '}\n"}}
%---
%[output:0a81e5c2]
%   data: {"dataType":"text","outputData":{"text":"MAR08636\tMAM01632c \t<=>\tMAM01632e \n","truncated":false}}
%---
%[output:3e9d4990]
%   data: {"dataType":"textualVariable","outputData":{"header":"1×1 cell array","name":"ans","value":"    {'MAM01632c  <=> MAM01632e '}\n"}}
%---
%[output:3a6b358b]
%   data: {"dataType":"text","outputData":{"text":"MAR10065\t0.0345 MAM00291c + 0.0172 MAM00427c + 0.0172 MAM00427r + 0.0172 MAM00611c + 0.0172 MAM00611m + 0.0172 MAM00622c + 0.0172 MAM00622m + 0.0345 MAM00766c + 0.0172 MAM01008c + 0.0172 MAM01008r + 0.0172 MAM01022c + 0.0172 MAM01022r + 0.0172 MAM01024c + 0.0172 MAM01024r + 0.0345 MAM01026r + 0.0172 MAM01028c + 0.0172 MAM01028r + 0.0172 MAM01030c + 0.0172 MAM01030r + 0.0172 MAM01031c + 0.0172 MAM01031r + 0.0345 MAM01032c + 0.0345 MAM01060c + 0.0345 MAM01232c + 0.0345 MAM01401c + 0.0345 MAM01600m + 0.0345 MAM01803c + 0.0345 MAM01924c + 0.0345 MAM02049m + 0.0345 MAM02348c + 0.0345 MAM02394c + 0.0345 MAM02624c + 0.0345 MAM02814c + 0.0345 MAM02842c + 0.0345 MAM02978c + 0.0345 MAM02980c + 0.0345 MAM02984c + 0.0345 MAM03102m \t->\tMAM10012c \nMAR13082\t45 MAM01371c + 0.0267 MAM01721n + 45 MAM02040c + 0.1124 MAM02847c + 0.4062 MAM03161c + 0.0012 MAM10012c + 5.3375 MAM10013c + 0.2212 MAM10014c + 0.4835 MAM10015c \t->\t45 MAM01285c + 45 MAM02039c + 45 MAM02751c + MAM03970c \n","truncated":false}}
%---
%[output:51a8e257]
%   data: {"dataType":"matrix","outputData":{"columns":1,"header":"2×1 cell array","name":"ans","rows":2,"type":"cell","value":[["'0.0345 MAM00291c + 0.0172 MAM00427c + 0.0172 MAM00427r + 0.0172 MAM00611c + 0.0172 MAM00611m + 0.0172 MAM00622c + 0.0172 MAM00622m + 0.0345 MAM00766c + 0.0172 MAM01008c + 0.0172 MAM01008r + 0.0172 MAM01022c + 0.0172 MAM01022r + 0.0172 MAM01024c + 0.0172 MAM01024r + 0.0345 MAM01026r + 0.0172 MAM01028c + 0.0172 MAM01028r + 0.0172 MAM01030c + 0.0172 MAM01030r + 0.0172 MAM01031c + 0.0172 MAM01031r + 0.0345 MAM01032c + 0.0345 MAM01060c + 0.0345 MAM01232c + 0.0345 MAM01401c + 0.0345 MAM01600m + 0.0345 MAM01803c + 0.0345 MAM01924c + 0.0345 MAM02049m + 0.0345 MAM02348c + 0.0345 MAM02394c + 0.0345 MAM02624c + 0.0345 MAM02814c + 0.0345 MAM02842c + 0.0345 MAM02978c + 0.0345 MAM02980c + 0.0345 MAM02984c + 0.0345 MAM03102m  -> MAM10012c '"],["'45 MAM01371c + 0.0267 MAM01721n + 45 MAM02040c + 0.1124 MAM02847c + 0.4062 MAM03161c + 0.0012 MAM10012c + 5.3375 MAM10013c + 0.2212 MAM10014c + 0.4835 MAM10015c  -> 45 MAM01285c + 45 MAM02039c + 45 MAM02751c + MAM03970c '"]]}}
%---
