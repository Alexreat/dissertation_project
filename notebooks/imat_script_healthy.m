%[text] ## Initial setup Path COBRA Toolbox and Gurobi solver

addpath('/Users/valentinreateguirangel/Documents/MATLAB/libSBML-5')
savepath
%%
addpath(genpath('/Users/valentinreateguirangel/Documents/MATLAB/MyProject/COBRA_project/cobratoolbox'));
%%
%  COBRA initialization 
initCobraToolbox(false) %[output:0a1fd161]
%%
% Point COBRA at Gurobi
changeCobraSolver('gurobi', 'LP'); %[output:1bee066c]
%%

% 4️⃣  Sanity-check Gurobi by solving a 1-variable LP
toy         = struct();
toy.S       = sparse(1,1,1,1,1);    % 1 × 1 stoichiometric matrix (just "1")
toy.mets    = {'m1[c]'};            % single metabolite
toy.rxns    = {'v1'};               % single reaction
toy.lb      = 0;                    % 0 ≤ v ≤ 10
toy.ub      = 10;
toy.c       = 1;                    % maximise v
toy.b       = 0;                    % S·v = 0  (empty mass-balance)

sol = optimizeCbModel(toy);

if sol.stat == 1 %[output:group:9c44b45c]
    fprintf('✔️  Gurobi solved the toy LP (objective = %.2f)\n', sol.f); %[output:06e16cbf]
else
    error(' Gurobi did NOT solve the toy problem – check the solver setup.');
end %[output:group:9c44b45c]
%%
%[text] ## Load the Human-GEM model
% Read the model 
model = readCbModel('/Users/valentinreateguirangel/Downloads/Human-GEM-1.19.0/model/Human-GEM.mat'); %[output:4bf79fb9]

%%
% Show a quick summary 
if isfield(model,'description')
    desc = model.description;
else
    desc = '--no description--';
end
fprintf('Loaded model:  %s\n', desc); %[output:95f14e56]
fprintf('Reactions:     %d\nMetabolites:   %d\nGenes:         %d\n',... %[output:group:6ce8ebf3] %[output:0fc84f3d]
        numel(model.rxns), numel(model.mets), numel(model.genes)); %[output:group:6ce8ebf3] %[output:0fc84f3d]
%%
%[text] ## Load Expression vector for healthy
exprFile = "/Users/valentinreateguirangel/dissertation_project/notebooks/Zscore_vector_healthy.txt";
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
disp(exprTbl(1:min(5,height(exprTbl)), :));        % first 5 rows %[output:79d0477b]
fprintf('Rows in table: %d genes\n', height(exprTbl)); %[output:28d6e76f]
fprintf('Expression range: %.3g  –  %.3g\n', min(exprTbl.expr), max(exprTbl.expr)); %[output:367d7a01]
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
fprintf('Mapped expression for %d / %d model genes (%.1f%%)\n', ... %[output:group:1b4e8257] %[output:35e81730]
        numMapped, numel(model.genes), 100*numMapped/numel(model.genes)); %[output:group:1b4e8257] %[output:35e81730]
%%
% Peek at the first few mapped entries
ix = find(~isnan(exprValues));
preview = table(model.genes(ix(1:min(5,end)))', exprValues(ix(1:min(5,end)))', ...
                'VariableNames', {'gene','expr'});
disp(preview); %[output:811ec1d6]
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
fprintf('Mapped expression for %d / %d reactions (%.1f%%)\n', ... %[output:group:45c8b541] %[output:6cd4867e]
        mapped, numel(expressionRxns), 100*mapped/numel(expressionRxns)); %[output:group:45c8b541] %[output:6cd4867e]
%%
%[text] ## Build healthy specific model with iMat
lowThr  = -0.5;     % ≤ lowThr  → force reaction *off*
highThr = 0.5;     % ≥ highThr → force reaction *on*
%%
% Pack everything in an *options* structure ( API for COBRA )
options                 = struct();
options.solver          = 'iMAT';
options.expressionRxns  = expressionRxns;
options.threshold_lb    = lowThr;
options.threshold_ub    = highThr;
options.task_specific = {{'MAR13082'}, [1], [0.001]};
options.core = {'MAR10065', 'MAR10062', 'MAR10063', 'MAR10064', 'MAR13082'};
%%
%[text] ## Setting up Nutrients Medium
% Show all exchange reactions in the model
excRxns = model.rxns(findExcRxns(model));
disp(excRxns(1:20)); % Just show first 20 for now %[output:4c441c1e]
%%
% Look up nutrient name 
getHumanName = @(metID) model.metNames{find(strcmp(model.mets, metID))};
printRxnFormula(model, 'MAR09045') %[output:872cc467] %[output:763da8b0]
getHumanName('MAM03089e')  %[output:955a0640]
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
%[text] ## Runing iMat
healthyModel = createTissueSpecificModel(model, options); %[output:3c73afac]
%%
% Summary
fprintf('\nHealhty-specific model created!\n'); %[output:3274ab40]
fprintf(' Reactions retained: %d  (of %d)\n', numel(healthyModel.rxns), numel(model.rxns)); %[output:2f71ef64]
fprintf(' Metabolites:        %d\n',          numel(healthyModel.mets)); %[output:0ad53d11]
fprintf(' Genes:              %d (of %d)\n',  numel(healthyModel.genes), numel(model.genes)); %[output:24333cde]
%%
save('healthyModel.mat', 'healthyModel');
%%
% Test biomass production in tumor model
findRxnIDs(healthyModel, 'MAR13082') %[output:4139d6e9]
findRxnIDs(healthyModel, 'MAR09034')  % glucose %[output:0358030b]
findRxnIDs(healthyModel, 'MAR09048')  % oxygen %[output:098ab1c7]

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":30.2}
%---
%[output:0a1fd161]
%   data: {"dataType":"text","outputData":{"text":"\n\n      _____   _____   _____   _____     _____     |\n     \/  ___| \/  _  \\ |  _  \\ |  _  \\   \/ ___ \\    |   COnstraint-Based Reconstruction and Analysis\n     | |     | | | | | |_| | | |_| |  | |___| |   |   The COBRA Toolbox - 2025\n     | |     | | | | |  _  { |  _  \/  |  ___  |   |\n     | |___  | |_| | | |_| | | | \\ \\  | |   | |   |   Documentation:\n     \\_____| \\_____\/ |_____\/ |_|  \\_\\ |_|   |_|   |   <a href=\"http:\/\/opencobra.github.io\/cobratoolbox\">http:\/\/opencobra.github.io\/cobratoolbox<\/a>\n                                                  | \n\n > Checking if git is installed ...  Done (version: 2.39.5).\n > Checking if the repository is tracked using git ...  Done.\n > Checking if curl is installed ...  Done.\n > Checking if remote can be reached ...  Done.\n > Initializing and updating submodules (this may take a while)... Done.\n > Adding all the files of The COBRA Toolbox ...  Done.\n > Define CB map output... set to svg.\n > TranslateSBML is installed and working properly.\n > Configuring solver environment variables ...\n   - [*---] ILOG_CPLEX_PATH: --> set this path manually after installing the solver ( see <a href=\"https:\/\/opencobra.github.io\/cobratoolbox\/docs\/solvers.html\">instructions<\/a> )\n   - [-*--] GUROBI_PATH: \/Library\/gurobi1202\/macos_universal2\/matlab\n   - [*---] TOMLAB_PATH: --> set this path manually after installing the solver ( see <a href=\"https:\/\/opencobra.github.io\/cobratoolbox\/docs\/solvers.html\">instructions<\/a> )\n   - [*---] MOSEK_PATH: --> set this path manually after installing the solver ( see <a href=\"https:\/\/opencobra.github.io\/cobratoolbox\/docs\/solvers.html\">instructions<\/a> )\n   Done.\n > Checking available solvers and solver interfaces ...Could not find installation of mosek, so it cannot be tested\nCould not find installation of tomlab_snopt, so it cannot be tested\n Done.\n > Setting default solvers ...Could not find installation of mosek, so it cannot be tested\nCould not find installation of mosek, so it cannot be tested\n Done.\n > Saving the MATLAB path ... Done.\n   - The MATLAB path was saved in the default location.\n\n > Summary of available solvers and solver interfaces\n\n\t\t\tSupport \t   LP \t MILP \t   QP \t MIQP \t  NLP \t   EP \t  CLP\n\t------------------------------------------------------------------------------\n\tdqqMinos     \tactive        \t    1 \t    - \t    1 \t    - \t    - \t    - \t    -\n\tglpk         \tactive        \t    1 \t    1 \t    - \t    - \t    - \t    - \t    -\n\tgurobi       \tactive        \t    1 \t    1 \t    1 \t    1 \t    - \t    - \t    -\n\tlp_solve     \tlegacy        \t    1 \t    - \t    - \t    - \t    - \t    - \t    -\n\tmatlab       \tactive        \t    1 \t    - \t    - \t    - \t    1 \t    - \t    -\n\tmosek        \tactive        \t    0 \t    - \t    0 \t    - \t    - \t    0 \t    0\n\tpdco         \tactive        \t    1 \t    - \t    1 \t    - \t    - \t    1 \t    -\n\tqpng         \tpassive       \t    - \t    - \t    1 \t    - \t    - \t    - \t    -\n\tquadMinos    \tactive        \t    1 \t    - \t    - \t    - \t    - \t    - \t    -\n\ttomlab_snopt \tpassive       \t    - \t    - \t    - \t    - \t    0 \t    - \t    -\n\t------------------------------------------------------------------------------\n\tTotal        \t-             \t    7 \t    2 \t    4 \t    1 \t    1 \t    1 \t    0\n\n + Legend: - = not applicable, 0 = solver not compatible or not installed, 1 = solver installed.\n\n\n > You can solve LP problems using: 'gurobi' - 'pdco' \n > You can solve MILP problems using: 'gurobi' \n > You can solve QP problems using: 'gurobi' - 'pdco' \n > You can solve MIQP problems using: 'gurobi' \n > You can solve NLP problems using: \n > You can solve EP problems using: 'pdco' \n > You can solve CLP problems using: \n\n> Checking for available updates ... skipped\nremoving: \/Users\/valentinreateguirangel\/Documents\/MATLAB\/MyProject\/COBRA_project\/cobratoolbox\/src\/analysis\/thermo\/componentContribution\/new\nremoving: \/Users\/valentinreateguirangel\/Documents\/MATLAB\/MyProject\/COBRA_project\/cobratoolbox\/src\/analysis\/thermo\/groupContribution\/new\nremoving: \/Users\/valentinreateguirangel\/Documents\/MATLAB\/MyProject\/COBRA_project\/cobratoolbox\/src\/analysis\/thermo\/inchi\/new\nremoving: \/Users\/valentinreateguirangel\/Documents\/MATLAB\/MyProject\/COBRA_project\/cobratoolbox\/src\/analysis\/thermo\/molFiles\/new\nremoving: \/Users\/valentinreateguirangel\/Documents\/MATLAB\/MyProject\/COBRA_project\/cobratoolbox\/src\/analysis\/thermo\/protons\/new\nremoving: \/Users\/valentinreateguirangel\/Documents\/MATLAB\/MyProject\/COBRA_project\/cobratoolbox\/src\/analysis\/thermo\/trainingModel\/new\n","truncated":false}}
%---
%[output:1bee066c]
%   data: {"dataType":"text","outputData":{"text":"\n > changeCobraSolver: Gurobi interface added to MATLAB path.\n","truncated":false}}
%---
%[output:06e16cbf]
%   data: {"dataType":"text","outputData":{"text":"✔️  Gurobi solved the toy LP (objective = 0.00)\n","truncated":false}}
%---
%[output:4bf79fb9]
%   data: {"dataType":"text","outputData":{"text":"Each model.subSystems{x} has been changed to a character array.\n","truncated":false}}
%---
%[output:95f14e56]
%   data: {"dataType":"text","outputData":{"text":"Loaded model:  Human-GEM.mat\n","truncated":false}}
%---
%[output:0fc84f3d]
%   data: {"dataType":"text","outputData":{"text":"Reactions:     12971\nMetabolites:   8455\nGenes:         2887\n","truncated":false}}
%---
%[output:79d0477b]
%   data: {"dataType":"text","outputData":{"text":"           <strong>gene<\/strong>            <strong>expr<\/strong>\n    <strong>___________________<\/strong>    <strong>____<\/strong>\n\n    {'ENSG00000000003'}     1  \n    {'ENSG00000000005'}     0  \n    {'ENSG00000000419'}     1  \n    {'ENSG00000000457'}     1  \n    {'ENSG00000000460'}     0  \n\n","truncated":false}}
%---
%[output:28d6e76f]
%   data: {"dataType":"text","outputData":{"text":"Rows in table: 60660 genes\n","truncated":false}}
%---
%[output:367d7a01]
%   data: {"dataType":"text","outputData":{"text":"Expression range: 0  –  1\n","truncated":false}}
%---
%[output:35e81730]
%   data: {"dataType":"text","outputData":{"text":"Mapped expression for 2884 \/ 2887 model genes (99.9%)\n","truncated":false}}
%---
%[output:811ec1d6]
%   data: {"dataType":"text","outputData":{"text":"                                                         <strong>gene<\/strong>                                                                  <strong>expr<\/strong>         \n    <strong>_______________________________________________________________________________________________________________<\/strong>    <strong>_____________________<\/strong>\n\n    {'ENSG00000000419'}    {'ENSG00000001036'}    {'ENSG00000001084'}    {'ENSG00000001630'}    {'ENSG00000002549'}    1    1    1    0    1\n\n","truncated":false}}
%---
%[output:6cd4867e]
%   data: {"dataType":"text","outputData":{"text":"Mapped expression for 8040 \/ 12971 reactions (62.0%)\n","truncated":false}}
%---
%[output:4c441c1e]
%   data: {"dataType":"text","outputData":{"text":"    {'MAR00659'}\n    {'MAR07108'}\n    {'MAR07110'}\n    {'MAR07112'}\n    {'MAR07114'}\n    {'MAR07116'}\n    {'MAR07118'}\n    {'MAR07120'}\n    {'MAR07122'}\n    {'MAR07124'}\n    {'MAR07126'}\n    {'MAR09023'}\n    {'MAR09024'}\n    {'MAR09032'}\n    {'MAR09808'}\n    {'MAR09809'}\n    {'MAR09810'}\n    {'MAR09811'}\n    {'MAR09812'}\n    {'MAR09813'}\n\n","truncated":false}}
%---
%[output:872cc467]
%   data: {"dataType":"text","outputData":{"text":"MAR09045\tMAM03089e \t<=>\t\n","truncated":false}}
%---
%[output:763da8b0]
%   data: {"dataType":"textualVariable","outputData":{"header":"1×1 cell array","name":"ans","value":"    {'MAM03089e  <=> '}\n"}}
%---
%[output:955a0640]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"'tryptophan'"}}
%---
%[output:3c73afac]
%   data: {"dataType":"text","outputData":{"text":"Set parameter Username\nSet parameter LicenseID to value 2674995\nSet parameter TimeLimit to value 7200\nSet parameter IntFeasTol to value 1e-09\nSet parameter MIPGap to value 1e-12\nSet parameter LogFile to value \"MILPlog\"\nAcademic license - for non-commercial use only - expires 2026-06-04\nGurobi Optimizer version 12.0.2 build v12.0.2rc0 (mac64[arm] - Darwin 24.5.0 24F74)\n\nCPU model: Apple M1 Pro\nThread count: 10 physical cores, 10 logical processors, using up to 10 threads\n\nNon-default parameters:\nTimeLimit  7200\nIntFeasTol  1e-09\nMIPGap  1e-12\nLogToConsole  0\n\nOptimize a model with 21761 rows, 26277 columns and 82054 nonzeros\nModel fingerprint: 0xf7b6665f\nVariable types: 12971 continuous, 13306 integer (13306 binary)\nCoefficient statistics:\n  Matrix range     [1e-04, 8e+04]\n  Objective range  [1e+00, 1e+00]\n  Bounds range     [1e+00, 1e+03]\n  RHS range        [1e+03, 1e+03]\nPresolve removed 13835 rows and 13897 columns\nPresolve time: 0.10s\nPresolved: 7926 rows, 12380 columns, 41817 nonzeros\nVariable types: 6612 continuous, 5768 integer (5767 binary)\nFound heuristic solution: objective 3.0000000\nFound heuristic solution: objective 4.0000000\nFound heuristic solution: objective 5.0000000\nFound heuristic solution: objective 6.0000000\n\nRoot relaxation: objective 7.133322e+03, 6637 iterations, 0.17 seconds (0.25 work units)\n\n    Nodes    |    Current Node    |     Objective Bounds      |     Work\n Expl Unexpl |  Obj  Depth IntInf | Incumbent    BestBd   Gap | It\/Node Time\n\n     0     0 7133.32218    0 1410    6.00000 7133.32218      -     -    1s\nH    0     0                    5664.0000000 7133.32218  25.9%     -    1s\nH    0     0                    5665.0000000 7133.32218  25.9%     -    1s\nH    0     0                    5667.0000000 7133.32218  25.9%     -    1s\nH    0     0                    5717.0000000 7133.32218  24.8%     -    1s\n     0     0 5775.18393    0  378 5717.00000 5775.18393  1.02%     -    1s\n     0     0 5775.18393    0  298 5717.00000 5775.18393  1.02%     -    1s\n     0     0 5775.18393    0  298 5717.00000 5775.18393  1.02%     -    1s\nH    0     0                    5721.0000000 5775.18393  0.95%     -    2s\n     0     0 5775.18021    0  236 5721.00000 5775.18021  0.95%     -    2s\n     0     0 5775.17994    0  236 5721.00000 5775.17994  0.95%     -    2s\n     0     0 5775.17655    0  234 5721.00000 5775.17655  0.95%     -    2s\nH    0     0                    5775.0000000 5775.17655  0.00%     -    2s\n     0     0 5775.17655    0  234 5775.00000 5775.17655  0.00%     -    2s\n\nCutting planes:\n  Learned: 1328\n  Gomory: 61\n  Cover: 13\n  Implied bound: 7\n  MIR: 98\n  Relax-and-lift: 15\n\nExplored 1 nodes (14713 simplex iterations) in 2.23 seconds (2.94 work units)\nThread count was 10 (of 10 available processors)\n\nSolution count 10: 5775 5721 5717 ... 3\n\nOptimal solution found (tolerance 1.00e-12)\nBest objective 5.775000000000e+03, best bound 5.775000000000e+03, gap 0.0000%\n","truncated":false}}
%---
%[output:3274ab40]
%   data: {"dataType":"text","outputData":{"text":"\nHealhty-specific model created!\n","truncated":false}}
%---
%[output:2f71ef64]
%   data: {"dataType":"text","outputData":{"text":" Reactions retained: 9776  (of 12971)\n","truncated":false}}
%---
%[output:0ad53d11]
%   data: {"dataType":"text","outputData":{"text":" Metabolites:        6481\n","truncated":false}}
%---
%[output:24333cde]
%   data: {"dataType":"text","outputData":{"text":" Genes:              2887 (of 2887)\n","truncated":false}}
%---
%[output:4139d6e9]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"9660"}}
%---
%[output:0358030b]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"5980"}}
%---
%[output:098ab1c7]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"5994"}}
%---
