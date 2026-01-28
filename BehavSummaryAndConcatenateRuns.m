%% Load data
BehavDir='C:\Users\dk818\OneDrive - Imperial College London\ICRF\Modelling\Inputs\IGTBehavData';
addpath(BehavDir);
OutDir=('C:\Users\dk818\OneDrive - Imperial College London\ICRF\Modelling\Outputs\BehavCombinedRuns');
% For spm
addpath('C:\Users\dk818\OneDrive - Imperial College London\LilThingsAndMATLABPrograms\spm12')
addpath(genpath('C:\Users\dk818\OneDrive - Imperial College London\ICRF\Modelling'))
% Example subj ID = ITI_DK01_V1_Run1. The ITI is the study prefix; Initals
% and the "01" are just in case there's another DK, but V1 is "first visit"
% (i.e., non-stim sessions), and the Run1 and Run2 are the first and second
% IGT runs in the scanner

Subj={"S01","S02","S03","S04","S05","S06","S08"};
Visit={'V1','V2'};
% Visit={'V1'};
Run={'Run1','Run2'};

%% Start loop for each subject and visit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for subj=1:size(Subj,2)
for vis=1:size(Visit,2)

SubjRun1=dir(string(strcat(BehavDir,'\ITI_S0',num2str(subj),'_',Visit(1,vis),'_Run1*.csv')));
SubjRun2=dir(string(strcat(BehavDir,'\ITI_S0',num2str(subj),'_',Visit(1,vis),'_Run2*.csv')));
Run1=readtable(SubjRun1(1).name);
Run2=readtable(SubjRun2(1).name);
Runs=[Run1;Run2];

end

%% Compute task performance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Total amount won, mean RT, proportion of choosing each deck
disp('*****************************************************')
disp(strcat('For participant S0',num2str(subj+1),'_',Visit(1,vis)))
disp(strcat('Total points:',num2str(table2array(Runs(end,'total')))))

% Participant payout is points/200. Max can win is extra £10. 
disp(strcat('Participant payout: £',num2str((table2array(Runs(end,'total')))/200)))
disp(strcat('Mean reaction time: ',num2str(table2array(nanmean(Runs(:,'rt_select'))))))

% Make a table with subj as row, cols = V1 points, V1 payout, V2 points, V2
% payout, V1-V2 points
if vis==1
BehavResTbl(subj,1)=table2array(Runs(end,'total'));
BehavResTbl(subj,2)=table2array(Runs(end,'total'))/200;
else
BehavResTbl(subj,3)=table2array(Runs(end,'total'));
BehavResTbl(subj,4)=table2array(Runs(end,'total'))/200;
end

%% Concatenate data for PVL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ID, loss, win, deck_select, and n_moves
Runs=Runs(:,["ID","deck_select","loss","win","init_mark","n_moves"]);
ds=table2array(Runs(:,2));

% Let's help map data to V1, to make life easy
if vis==1
Runs(:,7)=Runs(:,2);

elseif vis==2
tmp=ds;
tmp(ds == 1) = 4;
tmp(ds == 2) = 1;
tmp(ds == 3) = 2;
tmp(ds == 4) = 3;
Runs(:,7)=array2table(tmp);

elseif vis==3
tmp=ds;
tmp(ds == 1) = 3;
tmp(ds == 2) = 4;
tmp(ds == 3) = 1;
tmp(ds == 4) = 2;
Runs(:,7)=array2table(tmp);

elseif vis==4
tmp=ds;
tmp(ds == 1) = 2;
tmp(ds == 2) = 3;
tmp(ds == 3) = 4;
tmp(ds == 4) = 1;
Runs(:,7)=array2table(tmp);

else
disp("not prepared for this!!")

end

Runs.Properties.VariableNames(7)="deck_select_map_to_V1";

%% Save outputs
FileName=strcat(OutDir,'\S0',num2str(subj),'_',Visit(1,vis),'.csv');
% writetable(Runs,string(FileName));

end

%% Visualise difference between visit 1 and 2
% Make another column for the difference in 
BehavResTbl(:,end+1)=BehavResTbl(:,3)-BehavResTbl(:,1);
BehavResTbl(:,end+1)=BehavResTbl(:,end)/200;

BehavResTbl=array2table(BehavResTbl);
BehavResTbl.Properties.VariableNames={'V1 Points','V1 Payout','V2 Points','V2 Payout','V1-V2 Points','V1-V2 Payount'};

%% Diff Total Points?
TP=table2array(BehavResTbl(:,[1,3]));

disp('Paired T-test for TP')
[h,p,ci,stats] = ttest(TP(:,1),TP(:,2))

%% Compute KL divergence between EvU
%{
General guidelines:
0-0.1: Very similar distributions
0.1-0.5: Moderate difference
0.5+: Substantial difference
1+: Very different distributions
%}
n_subj=size(EvUpd,2)/2;
for subj=1:n_subj
p=table2array(EvUpd(:,subj));
q=table2array(EvUpd(:,subj+n_subj));

% Your original KL divergence
P = p / sum(p);
Q = q / sum(q);

epsilon = 1e-10;
P = P + epsilon; Q = Q + epsilon;
P = P / sum(P); Q = Q / sum(Q);

% Jensen-Shannon Divergence (symmetric, bounded 0-1)
M = (P + Q) / 2;
observed_JS_div = 0.5 * sum(P .* log(P ./ M)) + 0.5 * sum(Q .* log(Q ./ M));

% Permutation test
n_permutations = 10000;
combined = [P, Q];
null_JS_div = zeros(n_permutations, 1);

for i = 1:n_permutations
    % Randomly shuffle assignment of values to columns
    shuffled = combined(randperm(length(P)), :);
    P_perm = shuffled(:, 1);
    Q_perm = shuffled(:, 2);
    
    % Normalize
    P_perm = P_perm / sum(P_perm);
    Q_perm = Q_perm / sum(Q_perm);
    
    M_perm = (P_perm + Q_perm) / 2;
    null_JS_div(i) = 0.5 * sum(P_perm .* log(P_perm ./ M_perm)) + 0.5 * sum(Q_perm .* log(Q_perm ./ M_perm));
end

% Calculate p-value
p_value(subj,1) = mean(null_JS_div >= observed_JS_div);



end

%% 
% --------------------------------------
% MAKE EVUPDATE PLOTS PER VISIT!!!
% --------------------------------------
% Make a distribution 
addpath('C:\Users\dk818\OneDrive - Imperial College London\LilThingsAndMATLABPrograms\ColorBrewer');
colors = brewermap(n_subj, 'YlGnBu');

% Option 2 - make a voilin plot with lines between the means
figure;
hold on;

% Position violins at x=1 and x=2
positions = [1, 2];

for subj = 1:n_subj
    % Your data: subjectData is 79x2 (79 trials, 2 visits)
    
    p=zscore(table2array(EvUpd(:,subj)));
    q=zscore(table2array(EvUpd(:,subj+n_subj)));
    subjectData=[p,q];

    % Create violins with transparency and specific color
    v = violinplot(subjectData, [], ...
        'ViolinColor', colors(subj,:), ...
        'ViolinAlpha', 0.2, ...
        'ShowMean', true, ...
        'ShowData', false, ...
        'Width', 0.15);
end

set(gca,'FontName','Arial','FontSize',18)
set(gca, 'XTick', [1 2], 'XTickLabel', {'Visit 1', 'Visit 2'});
ylabel('z-score EvU');


TP=table2array(BehavResTbl(:,[1,3]));
Payout=table2array(BehavResTbl(:,[2,4]));

%% Make TP violins 
figure()
hold on
v=violinplot(TP,[],'ViolinColor',colors(7,:),'ShowData',false);
% v(1).ViolinColor=colors(7,:);
% v(2).ViolinColor=colors(7,:);
% v(1).ShowData='false';
% v(2).ShowData='false';
for subj=1:n_subj
plot(positions, TP(subj,:), '-o', ...
    'Color', colors(subj,:), ...
    'LineWidth', 1.5, ...
    'MarkerFaceColor', colors(subj,:));
end
xticklabels({'Visit 1','Visit 2'})
set(gca,'FontName','Arial','FontSize',18)
ylabel('Total Points')

%% Make distribution
figure()
set(gcf,'Color','white')
hold on
for subj=1 %:n_subj
% Estimate PDF
[f, xi] = ksdensity(zscore(table2array(EvUpd(:,subj))));
% Create filled area plot
fill(xi, f, 'b', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
hold on;
plot(xi, f, 'LineWidth', 2, 'Color', 'k');
% xlabel('z-score EvU');
% ylabel('EvUpdate');
% title('Smooth Probability Density Function');
cmap = brewermap([], 'YlGnBu');
fill(xi, f, colors(subj,:), 'FaceAlpha', 0.6, 'EdgeColor', 'k', 'LineWidth', 1.5);
set(gca,'FontName','Arial','FontSize',18)
end


