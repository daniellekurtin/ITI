%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup FOR LOCAL OR HPC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FOR RUNNING ON HPC
addpath(genpath('/rds/general/project/dan_data_tide/live/IndivTIStim/Behaviour'));
% For outputs
outputpath='/rds/general/project/dan_data_tide/live/IndivTIStim/Behaviour/Outputs';

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assess which neural metrics has the strongest
% relationship to EE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numParcels=600; 
load(strcat('TrialMetrics_',num2str(numParcels),'Parcels.mat'));
load(strcat('df_',num2str(numParcels),'ROIs.mat'));

[n_subj, n_sess]=size(df);
num_rois=size(df{1,1},2);
subjList={"S01","S02","S03","S04","S05","S06","S08"};
TMPerSubj={};

EvUpd=readtable('Ev_upd_df_PVL.csv');
EvUpd(:,1)=[];

%%
% First need to combine runs 1 and 2
for sesh=1:n_sess
for s=1:n_subj 

saveName=strcat(outputpath,'/random_forest_results_EVUpdate_',subjList{1,s},'V',num2str(sesh),'_',num2str(numParcels),'Parcels.mat');

fprintf('Running for subj: ')
fprintf(subjList{1,s})
fprintf(', V')
fprintf(num2str(sesh))
fprintf('\n *********************************************** \n')

% Get EvUpdate for this subj and sesh
EE=table2array(EvUpd(:,(s*sesh)));
if sum(isnan(EE))>0
toDel=isnan(EE);
EE(toDel,:)=[];
neuralMetrics(toDel,:)=[];
end
EE=zscore(EE);

% Get subj neural metrics
neuralMetrics=TrialMetrics{s,sesh};

% Zscore each col of NeuralMetrics except for the state timeseries
for ii=1:size(neuralMetrics,2)
neuralMetrics(:,ii)=zscore(neuralMetrics(:,ii));
end
[numTrials,numMetrics] = size(neuralMetrics);

%%%%%%%%%%%%%%%%%%%%%
% Start of random forest parameter setting
%%%%%%%%%%%%%%%%%%%%%
rng(42); % For reproducibility

% Parameters
numTrees = 100;           % Number of trees in the forest
numPermutations = 1000;   % Number of permutations for significance testing
testFraction = 0.3;       % Fraction of data used for testing
kFold = 5;                % Number of folds for cross-validation

% Initialize results storage
results = struct();
results.r2 = zeros(numMetrics, 1);
results.importance = cell(numMetrics, 1);
results.p_values = zeros(numMetrics, 1);
results.permR2 = zeros(numMetrics, numPermutations);

%% Random Forest Regression for each neural metric
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Computing Random Forest regression for each neural metric...\n');

for m = 1:numMetrics
    fprintf('Processing metric %d of %d\n', m, numMetrics);
    
    % Extract current neural metric
    currentMetric = neuralMetrics(:, m);
    
    % Initialize cross-validation
    cv = cvpartition(numTrials, 'KFold', kFold);
    
    % Storage for cross-validated predictions
    cv_predictions = zeros(numTrials, 1);
    
    % Cross-validated predictions
    for k = 1:kFold
        % Get training and testing indices
        trainIdx = cv.training(k);
        testIdx = cv.test(k);
        
        % Train random forest on training data
        Mdl = TreeBagger(numTrees, currentMetric(trainIdx), EE(trainIdx), ...
            'Method', 'regression', 'OOBPrediction', 'on', ...
            'PredictorSelection', 'curvature');
        
        % Predict on test data
        cv_predictions(testIdx) = predict(Mdl, currentMetric(testIdx));
    end
    
    % Calculate R² for cross-validated predictions
    SSE = sum((EE - cv_predictions).^2);
    SST = sum((EE - mean(EE)).^2);
    results.r2(m) = 1 - SSE/SST;
    
    % Train final model on all data to get feature importance
    finalMdl = TreeBagger(numTrees, currentMetric, EE, ...
        'Method', 'regression', 'OOBPrediction', 'on','OOBPredictorImportance', 'on', ...
        'PredictorSelection', 'curvature');
    
    % Get variable importance
    results.importance{m} = finalMdl.OOBPermutedPredictorDeltaError;
    
    % Permutation testing
    disp('Starting permutation testing')
    for p = 1:numPermutations
        % Permute the task performance data
        permEE = EE(randperm(numTrials));
        
        % Cross-validated predictions for permuted data
        perm_predictions = zeros(numTrials, 1);
        
        for k = 1:kFold
            % Get training and testing indices
            trainIdx = cv.training(k);
            testIdx = cv.test(k);
            
            % Train random forest on permuted training data
            permMdl = TreeBagger(numTrees, currentMetric(trainIdx), permEE(trainIdx), ...
                'Method', 'regression', 'PredictorSelection', 'curvature');
            
            % Predict on test data
            perm_predictions(testIdx) = predict(permMdl, currentMetric(testIdx));
        end
        
        % Calculate and store R² for permuted data
        SSE_perm = sum((permEE - perm_predictions).^2);
        SST_perm = sum((permEE - mean(permEE)).^2);
        results.permR2(m, p) = 1 - SSE_perm/SST_perm;
    end
    
    % Calculate p-value from permutation distribution
    results.p_values(m) = mean(results.permR2(m, :) >= results.r2(m));
end

%% Multiple comparison correction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% False Discovery Rate correction (Benjamini-Hochberg)
[~, sortedIndices] = sort(results.p_values);
numTests = length(results.p_values);
results.p_adjusted = zeros(numTests, 1);

for i = 1:numTests
    results.p_adjusted(sortedIndices(i)) = results.p_values(sortedIndices(i)) * numTests / i;
end

% Ensure no adjusted p-values exceed 1
results.p_adjusted = min(results.p_adjusted, 1);

% Save Results
save(saveName,'-v7.3');
disp('*******************************************************************')
disp('*****************    END OF THIS SUBJ!      *********************')
disp('*******************************************************************')


%% Compute target neural metric
[TMRankingSort_R2(:,1),TMRankingSort_R2(:,2)]=sort(abs(cell2mat(results.importance)),'descend');
TMPerSubj{sesh}(s,1)=TMRankingSort_R2(1,1);
TMPerSubj{sesh}(s,1)=TMRankingSort_R2(1,1);

end
end

%% Save everything
cd(outputpath)
save('BrainBehavAllSubj.mat')
