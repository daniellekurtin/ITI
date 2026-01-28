clear all
close all
addpath(genpath('/rds/general/project/dan_data_tide/live/IndivTIStim/PhenomModels'))
outputpath='/rds/general/project/dan_data_tide/live/IndivTIStim/PhenomModels/Outputs';

%% Set parameters 
% fMRI data params

load df_200ROIs.mat % dataframe with subj's time*ROI 
NSUB=size(df,1);
[Tmax,N]=size(df{1,1}); % num regions
indexN=1:N;
Isubdiag = find(tril(ones(N),-1));
TR=1.25;  % Repetition Time (seconds)
load BrainBehav_200Parcels.mat % Target metric value and code per subj
load TaskTimeseries.mat
load TrialMetrics_200Parcels.mat
load templateMatAndNetName_200Parcels.mat

% SC data params
load('SC.mat')
C = SC;
C=C(indexN,indexN);
maxC=1;
C = C/max(max(C))*maxC;

% Frequency of each node
load hopf_freq_sch200.mat
sigma=0.02;
a=-0.02*ones(N,2);
Tau=1;
epsFC=0.005; 

% NEW PARAMS FOR NL
omega=repmat(2*pi*f_diff(1:N)',1,2); 
omega(:,1)=-omega(:,1);
dsig = sqrt(dt)*sigma;

OptimStrat={"20FC80TM"};

for OS=1:size(OptimStrat,2)

% Before I used to run everything, but why would I need to do that? Can't I
% just load what I had before? 

load('Ceff_Optim'+OptimStrat{1,OS}+'_200Parcels.mat')

%% Perturbation
% sync pert --> move a towards positive
% noise pert --> move a towards negative
A=-0.02:0.001:0.02; %sync

for nsub=1:size(df,1)
for vis=1:size(df,2)

% TrialMetrics = TrialMetrics{nsub,vis}(numTrials,numMetrics) = miFC per trial for the target metric
% TMPerSubj{1,vis}(nsub,[pval, targetMetricCode, targetMetricmiFC]
neuralMetrics=TrialMetrics{nsub,vis}(:,TMPerSubj{1,vis}(nsub,2)); 
clear ROIList 
% Get list of regions in the target networks to loop through
% ROIList=IdentifyRegionsForPerturbation(TMPerSubj(nsub,2));
% Try to get list of ROIs using 
tmpMat=squeeze(templateMat(:,:,TMPerSubj{1,vis}(nsub,2)));
[rows,~,~]=find(tmpMat>0);
ROIList(:,1)=unique(rows);

% For each region . . . 
% ---------------------------------------
for nodestep=1:size(ROIList,1)
node=ROIList(nodestep,1);

% . . . and each value of A. . .
for astep=1:size(A,2)

a(node)=A(1,astep);
disp(strcat('Subj #',num2str(nsub),': running node',num2str(node),' and astep #',num2str(astep), 'which is ',num2str(a(node))))

% NL hopf
wC = squeeze(Ceff_sub{1,vis}(nsub,:,:));
sumC = repmat(sum(wC,2),1,2); % for sum Cij*xj
xs=zeros(Tmax,N);
z = 0.1*ones(N,2); % --> x = z(:,1), y = z(:,2)
nn=0;

% discard first 3000 time steps
for t=0:dt:1000
suma = wC*z - sumC.*z; % sum(Cij*xi) - sum(Cij)*xj
zz = z(:,end:-1:1); % flipped z, because (x.*x + y.*y)
z = z + dt*(a.*z + zz.*omega - z.*(z.*z+zz.*zz) + suma) + dsig*randn(N,2);
end

% actual modeling (x=BOLD signal (Interpretation), y some other oscillation)
for t=0:dt:((Tmax-1)*TR)
suma = wC*z - sumC.*z; % sum(Cij*xi) - sum(Cij)*xj
zz = z(:,end:-1:1); % flipped z, because (x.*x + y.*y)
z = z + dt*(a.*z + zz.*omega - z.*(z.*z+zz.*zz) + suma) + dsig*randn(N,2);
if abs(mod(t,TR))<0.01
nn=nn+1;
xs(nn,:)=z(:,1)';
end
end

% Compute FCsim, which is the same as iFC
BOLD=xs';
clear signal_filt Phase_BOLD;
for seed=1:N
Phase_BOLD(seed,:)= angle(hilbert(BOLD(seed,:)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NORMALLY, the error is found in the following steps: 
% FCsim=corrcoef(signal_filt');
% errorFC{nsub}(node,jj)=mean(mean((FCemp-FCsim).^2)); 
% jj=1+jj;

% HOWEVER, now we're gonna get dynamical. 
% We have to compute the iFC per timepoint
% then compute connectivity in our target network per trial
% then find error between simulated and empirical trial metric per trial

    % For each timepoint, compute the iFC
    for t=1:size(Phase_BOLD,2) %for each time point      
        iFC=zeros(N); 
        for n=1:N 
            for p=1:N
                iFC(n,p)=cos(Phase_BOLD(n,t)-Phase_BOLD(p,t));
            end
        end
        tmpMat=squeeze(templateMat(:,:,TMPerSubj{1,vis}(nsub,2)));

        % iFC*tmpMat gives a matrix, but it's weird format
        % iFC(tmpMat==1) is a N*1 array with the correct vals. . . 
        IFC(t,1)=nanmean(iFC(tmpMat==1));

    end % time

    % Downsample from TRs/vols to per trial 
    for ii=1:size(start_indices,1)-1
  
        start=start_indices(ii+1,1);    

        % Store the mean between-network connectivity across vols per trial
        % IFCperTrial{nsub,vis}(ii-1,nodestep,astep)=nanmean(squeeze(IFC(start:start+9,:,:)),'all');
        IFCperTrial{nsub,vis}(ii,nodestep,astep)=IFC(start:start+9,1);

       % FINALLY - Compute the error!
       errorFromEmpPerTrial{nsub,vis}(ii,nodestep,astep)=neuralMetrics(ii,1)-IFCperTrial{nsub,vis}(ii,nodestep,astep);
    end

end % a
end % node

cd(outputpath)
tmpsaveName=strcat('Perturb_Optim', OptimStrat{1,OS},'_Subj',num2str(nsub) ,'.mat');
save(tmpsaveName)
end % vis
end % subj


%% Find optimum
ComboNames=readtable('ComboNames.xlsx');
ComboNames=table2array(ComboNames);

for nsub=1:size(errorFromEmpPerTrial,2)
for vis=1:2

% Get list of regions in the target networks to loop through
tmpMat=squeeze(templateMat(:,:,TMPerSubj{1,vis}(nsub,2)));
[rows,~,~]=find(tmpMat>0);
ROIList(:,1)=unique(rows);

% Make a long form of TM for all node and a combos
clear miFC_data
h=0;
for nodestep=1:size(IFCperTrial{nsub,vis}(:,:,:),2)
for astep=1:size(a,2)
h=h+1;
tmp=squeeze(IFCperTrial{nsub,vis}(:,nodestep,astep));

% Use this to find maximum val
% NOTE - SubjTMList = ROI ID, ASTEP, TM
SubjTMList(h,1)=ROIList(nodestep,1);
SubjTMList(h,2)=a(1,astep);
SubjTMList(h,3)=nanmean(tmp);

% Use this to make surf plot
miFC_data(astep, nodestep)=nanmean(squeeze(IFCperTrial{nsub,vis}(:,nodestep,astep)));

end
end

% Create continuous indices for ROI IDs for plotting
len=1:1:size(ROIList,1);

% Create meshgrid for surface plot
[A, ROI_cont] = meshgrid(a, len);

% Find optimal val based on tstat
% If tstat is pos, then we want to maximise TM to maximise EEw. Neg Tstat =
% minimise TM for max EEw.
% NOTE - SUBJOPTIMAL = TM, ASTEP, ROI ID
if model_stats(nsub).tstat>0
[SubjOptimal(nsub,1),index]=max(SubjTMList(:,3));
else
[SubjOptimal(nsub,1),index]=min(SubjTMList(:,3));
end
SubjOptimal(nsub,2)=SubjTMList(index,2);
SubjOptimal(nsub,3)=SubjTMList(index,1);
OptimalROIName(nsub)=ComboNames(SubjOptimal(nsub,3),1);
end
end

cd(outputpath)
saveName=strcat('Perturb_Optim', OptimStrat{1,OS},'.mat');
save(saveName)
end