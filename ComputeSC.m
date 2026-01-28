% ==========================================
% Description:
%   This script reads in normalised and non-normalised SC matrices and does
%   some quality control
%   
% NOTE: THIS SCRIPT DOES NOT INCLUDE SUB-S08 BECAUSE DWI PREPROC DID NOT
% FINISH ON THEM. ANNOYING, BUT WHATEVS.
%
% Inputs:
%   - SC matrices, normalised and not, from the HPC
%
% Outputs:
%   - SC 
%
% Author: Danielle Kurtin
% Date: September 24, 2025
% 
% ==========================================

%% Setup 
clear all
close all
homeDir='C:\Users\dk818\OneDrive - Imperial College London\ICRF\Modelling';
addpath(genpath(homeDir))
subs = {'01','02','03','04','05','06'};

%% Make a big array with normalised and non-normalised SC matrices
for sub=1:size(subs,2)

SC(:,:,sub)=table2array(readtable(strcat(homeDir,'/SC/sub-S',subs{1,sub},'V1/sub-S',subs{1,sub},'V1_ses-1_Schaefer2018_200Parcels_SC_matrix.csv')));

% Nor all SCNorm are 214x214; pad them out if needed
clear tmpSCnorm
tmpSCnorm=table2array(readtable(strcat(homeDir,'/SC/sub-S',subs{1,sub},'V1/sub-S',subs{1,sub},'V1_ses-1_Schaefer2018_200Parcels_SC_normalized.csv')));
[x,y]=size(tmpSCnorm);
if x<214
tmpSCnorm(end+1:214,:)=zeros(214-x,y);
end

if y<214
tmpSCnorm(:,end+1:214)=zeros(214,214-y);
end

SCnorm(:,:,sub)=tmpSCnorm;
end

SCmean=mean(SC,3);
SCnormmean=mean(SCnorm,3);

%% Assess outputs
figure()
imagesc(SCmean);
colorbar()

figure()
imagesc(SCnormmean);
colorbar()

% Compare to Jakub's 
load("C:\Users\dk818\OneDrive - Imperial College London\ICRF\ModellingOld\SC_schaefer200_17Networks_32fold_groupconnectome_2mm_symm.mat")
figure()
imagesc(SC)
colorbar()

figure()
imagesc(SC_D)
colorbar()

%% Save
% Takeaway - the SCmean and SC look comparable. I'm happy enough.
% Let's save it! 
cd(strcat(homeDir,'/Outputs'))
clearvars SC
SC=SCmean;
save('SC.mat','SC');