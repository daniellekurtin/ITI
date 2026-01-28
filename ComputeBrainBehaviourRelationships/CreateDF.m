% ==========================================
% Description:
%   This script performs the following processing steps on time series data:
%     1. Compute the set of valid template ROIs within the FOV for each subject;
%     3. Select a subset of ROIs (cortical and subcortical regions) from this intersection 
%        and extract the corresponding time series for each subject.
%
% Inputs:
%   - Atlas label and code, atlas image
%   - Time series, brain mask of BOLD for each subject
%
% Outputs:
%   - Selected intersection of valid ROIs across subjects
%   - Filtered time series for the selected ROIs in each subject
%
% Author: Huaxin Fan
% Date: July 1, 2025
% 
% Edited by Danielle Kurtin
% 25/08/2025
% ==========================================
clear all
close all
homeDir='C:\Users\dk818\OneDrive - Imperial College London\ICRF\Modelling';
addpath(genpath(homeDir))

subs = {'01','02','03','04','05','06','08'};
ts_subs = cell(numel(subs),1);
ROInm_existed_inFOV_subs = cell(numel(subs),1);

%% Start loop for each parcellation
for numParcels=[200,400,600]
% get the code and label of each brain region in the atlas
t_label = strcat('Schaefer2018_',num2str(numParcels),'Parcels_7Networks_order_LUT.txt');
fid = fopen(t_label);
data = textscan(fid, '%d %s %*d %*d %*d %*d', 'CommentStyle', '#');
fclose(fid);
ID = data{1};
Label = data{2}; 
ID = ID(2:end); 
Label = Label(2:end);

%% Start loop for each visit, subject, and run
for vis=1:2
for i = 1:numel(subs)
    subID = subs{i};
    Sub = ['sub-S', subID,'V',num2str(vis)];

for run=1:2
    % read timeseries
    ts_orig = readtable([homeDir,'\Inputs\Timeseries\TimeseriesFromHuaxinAtlasTransform\',Sub,'_ses-1_task-IGT_run-0',num2str(run),'_space-T1w_Schaefer2018_',num2str(numParcels),'Parcels_7Networks_SubT1Vol_TimeSeries.csv']);

    % read atlas image
    gunzip([homeDir,'\Inputs\Timeseries\TimeseriesFromHuaxinAtlasTransform\SubjMaskAndAtlas\',Sub,'_ses-1_task-IGT_run-0',num2str(run),'_space-T1w_Schaefer2018_',num2str(numParcels),'Parcels_7Networks_SubT1Vol_resamp3mm.nii.gz']);
    img = niftiread([homeDir,'\Inputs\Timeseries\TimeseriesFromHuaxinAtlasTransform\SubjMaskAndAtlas\',Sub,'_ses-1_task-IGT_run-0',num2str(run),'_space-T1w_Schaefer2018_',num2str(numParcels),'Parcels_7Networks_SubT1Vol_resamp3mm.nii']);
    delete([homeDir,'\Inputs\Timeseries\TimeseriesFromHuaxinAtlasTransform\SubjMaskAndAtlas\',Sub,'_ses-1_task-IGT_run-0',num2str(run),'_space-T1w_Schaefer2018_',num2str(numParcels),'Parcels_7Networks_SubT1Vol_resamp3mm.nii']); 

    unique_values = unique(img);   % number of unique values = regions with specific intensity 
    ROIidx_existed = unique_values(2:end); % for each ROI (row), get the intensity that labels it
    
    % Keep only the ROI we care about
    del=size(ROIidx_existed,1)-numParcels;
    ROIidx_existed(1:del,:)=[];

    % rename timeseries table as brain regions
    [tf, idx] = ismember(ROIidx_existed, ID);
    ROInm_existed = Label(idx);
 
    ts_existed = ts_orig(:,ROIidx_existed); % For subj 1, there are 241 regions with a unique intensity
    ts_existed.Properties.VariableNames = ROInm_existed;
        
    %% Make same format as other df
    SubCorName={'l_hip','l_amy','l_thl','l_nac','l_pal','l_put','l_cau','r_hip','r_amy','r_thl','r_nac','r_pal','r_put','r_cau'};
    for rr=1:size(SubCorName,2)
        ts_existed(:,end+1) = readtable([homeDir,'\Inputs\Timeseries\Subcortex\IGT_run-0',num2str(run),'\',Sub,'_',SubCorName{1,rr},'.csv']);
        ts_existed.Properties.VariableNames{1,end}=SubCorName{1,rr};
    end

    tmpdf{i,run}=ts_existed;
end

    df{i,vis}=[tmpdf{i,run};tmpdf{i,run}];
end
    
end

%% Save
cd(strcat(homeDir,'\Outputs'))
saveName=strcat('df_',num2str(numParcels),'ROIs.mat');
save(saveName,'df')



end

% %% Check if I actually have NaNs
% for subj=1:size(df,1)
% for vis=1:size(df,2)
% for r=1:size(df{1,1},2)
% 
%     if sum(isnan(table2array(df{subj,vis}(:,r))))>0
%         disp(strcat('NaN for S0',num2str(subj),'V',num2str(vis),'ROI',num2str(r)))
%     end
% 
% end
% end
% end
