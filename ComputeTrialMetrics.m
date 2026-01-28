%% Add paths
clear all
addpath(genpath('C:\Users\dk818\OneDrive - Imperial College London\ICRF\Modelling'));
% For outputs
outputpath='C:\Users\dk818\OneDrive - Imperial College London\ICRF\Modelling\Outputs';
% For miFC
addpath(genpath('C:\Users\dk818\OneDrive - Imperial College London\NCORE\miFCdFC\Functions'));
% For Brain Connectivity Toolbox
addpath('C:\Users\dk818\OneDrive - Imperial College London\LilThingsAndMATLABPrograms\BCT-main\BCT\2019_03_03_BCT')

%% Load data and set parameters
% Template networks and task timeseries
numParcels=200;
load(strcat('templateMatAndNetName_',num2str(numParcels),'Parcels.mat'));
% load(strcat('df_',num2str(numParcels),'ROIs.mat'));
load(strcat('df_',num2str(numParcels),'ROIs_8mmSmooth.mat'))

[n_subj, n_sess]=size(df);
[Tmax,num_rois]=size(df{1,1});
load TaskTimeseries.mat
num_metrics=size(templateMat,3);
TR=1.25;

%% Compute iFC metrics
for V=1:n_sess
for s=1:n_subj 
    disp(strcat('Running HBSD subj',num2str(s)))
    
    % Get the BOLD signals from this subject in this task
    BOLD = table2array(df{s,V})';
    Phase_BOLD=zeros(num_rois,Tmax); 

    % Get the BOLD phase using the Hilbert transform
    for seed=1:num_rois
        Phase_BOLD(seed,:)=angle(hilbert(BOLD(seed,:)));
    end
    for t=1:size(Phase_BOLD,2) %for each time point      
        iFC=zeros(num_rois); 
        for n=1:num_rois 
            for p=1:num_rois
                iFC(n,p)=cos(Phase_BOLD(n,t)-Phase_BOLD(p,t));
            end
        end

        % COMPUTE THE MEAN INTER-NETWORK CONNECTIVITY 
        for ii=1:num_metrics
        tmpMat=squeeze(templateMat(:,:,ii));
        BtwnNetworkConn{s,V}(ii,t)=mean(iFC(tmpMat==1));
        end
    end

    % Create the network metrics table
    % Store the mean between-network connectivity across vols per trial
    for ii=1:size(start_indices,1)-1
        start=start_indices(ii+1,1);       
        TrialMetrics{s,V}(ii,1:num_metrics)=mean(BtwnNetworkConn{s,V}(:,start:start+9),2)';
    end
end
end
%%
cd(outputpath)
% saveName=strcat('TrialMetrics_',num2str(numParcels),'Parcels.mat');
saveName=strcat('TrialMetrics_',num2str(numParcels),'Parcels_8mmSmooth.mat');
clearvars -except TrialMetrics BtwnNetworkConn saveName
save(saveName,'-v7.3')
