%{
The learning and fitting is done with the entire target metric matrix
No early stopping
%}

addpath(genpath('/rds/general/project/dan_data_tide/live/IndivTIStim/PhenomModels'))
outputpath='/rds/general/project/dan_data_tide/live/IndivTIStim/PhenomModels/Outputs';

% Conditions
%OptimStrategy={'20FC80TM','80FC20TM','50FC50TM'};
OptimStrategy={'80FC20TM'};

for OS=1:size(OptimStrategy,2)

clearvars -except OS OptimStrategy outputpath

disp('-------------------------------------------------------------------')
disp('-------------------------------------------------------------------')
disp('-------------------------------------------------------------------')
disp(strcat('NOW RUNNING OptimStrat',num2str(OS)))
disp('-------------------------------------------------------------------')
disp('-------------------------------------------------------------------')
disp('-------------------------------------------------------------------')

%% Set parameters 
numParcels=200;
load(strcat('templateMatAndNetName_',num2str(numParcels),'Parcels.mat'));
load(strcat('TrialMetrics_',num2str(numParcels),'Parcels_8mmSmooth.mat'));
load(strcat('df_',num2str(numParcels),'ROIs_8mmSmooth.mat'));
load(strcat('BrainBehav_',num2str(numParcels),'Parcels_8mmSmooth.mat'));

% fMRI data params
NSUB=size(df,1);
[Tmax,N]=size(df{1,1}); % num regions
indexN=1:N;
Isubdiag = find(tril(ones(N),-1));
TR=1.25;  % Repetition Time (seconds)

% SC data params
load('SC.mat')
C = SC;
C=C(indexN,indexN);
maxC=1;
C = C/max(max(C))*maxC;

% Frequency of each node
load(strcat('hopf_freq_sch',num2str(numParcels),'_8mmSmooth.mat'));
sigma=0.02;
a=-0.02*ones(N,2);
Tau=1;
epsFC=0.005; 
 
% NEW PARAMS FOR NL
omega=repmat(2*pi*f_diff(1:N)',1,2); 
omega(:,1)=-omega(:,1);
dt=0.1*TR/2; % ms
dsig = sqrt(dt)*sigma;

% for outputs
BOLD4All=zeros(NSUB,N,Tmax);
FCempAllSubj=zeros(NSUB,N,N);
FCsimAllSubj=zeros(NSUB,N,N); 
fittFC_sub_iter={};
errorFCtmAllSubj={};

%% Group FC emp
for sesh=1:size(df,2)
for nsub=1:NSUB
    ts=table2array(df{nsub,sesh})';
    FCemp=corrcoef(ts');
    FC(nsub,:,:)=FCemp;
end
end
FCemp=squeeze(mean(FC));

%% Fit Group/Global GEC
G=1;
Cnew=G*C;
olderror=100000;
disp('Fitting Group GEC for subj')

for iter=1:5000

    wC = Cnew;
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

    BOLD=xs';
    for seed=1:N
	ts33(seed,:)=BOLD(seed,:);
    end

    FCsim=corrcoef(ts33');
    errorFC(iter)=mean(mean((FCemp-FCsim).^2));

    % Stop if error is small and does not grow
    if mod(iter,100)<0.1
        errornow=mean(mean((FCemp-FCsim).^2));
        if  (olderror-errornow)/errornow<0.0001
            break;
        end
        if  olderror<errornow
            break;
        end
        olderror=errornow;
    end
    
    % Update C (i.e., GEC)
    for i=1:N 
        for j=1:N
            if (C(i,j)>0 || (abs(j-i)==100 && j<=N && i<=N))
            Cnew(i,j)=Cnew(i,j)+epsFC*(FCemp(i,j)-FCsim(i,j));
                if Cnew(i,j)<0
                    Cnew(i,j)=0;
                end
            end
        end
    end
    Cnew = Cnew/max(max(Cnew))*maxC;
end
Ceffgroup=Cnew;

%% Individual
BOLD4All={};
FCempAllSubj={};
FCsimAllSubj={};
Ceff_sub={};

for nsub=1:NSUB
for sesh=1:size(df,2)
    disp(strcat('Fitting individual GEC for S',num2str(nsub),'V',num2str(sesh)))

    % Compute FC
    ts=table2array(df{nsub,sesh})'; % get subj timeseries
    FCempSubj=corrcoef(ts'); % Compute FC	
    tmempSubj=TMPerSubj{1,sesh}(nsub,3);

    % Find subject's GEC
    f_diff_s=squeeze(f_diff_sub(nsub,sesh,:));
    omega=repmat(2*pi*squeeze(f_diff_sub(nsub,sesh,1:N)),1,2); omega(:,1)=-omega(:,1);
    Cnew=Ceffgroup; % if we dont start with groupo, would just use a subject´s C as starting point
    olderror=100000;

    for iter=1:5000
        wC = Cnew;
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

        % Filter simulated timeseries
        BOLD=xs';
        for seed=1:N
            ts33(seed,:)=BOLD(seed,:);
        end

        FCsim=corrcoef(ts33');
        
    	% TM
        tmpMat=squeeze(templateMat(:,:,TMPerSubj{1,sesh}(nsub,2)));
        tm=mean(FCsim(tmpMat==1));
 
        errorFCtm(iter)=nanmean((tmempSubj-tm).^2);  % we want to change this for each patient

        errorFC(iter)=mean(mean((FCempSubj-FCsim).^2));  % monitoring how FC evolves

        fittFC_sub_iter{nsub}(iter)=corr2(FCempSubj(Isubdiag),FCsim(Isubdiag));

        for i=1:N  %% learning
            for j=1:N
            if (C(i,j)>0 || (abs(j-i)==100 && j<=N && i<=N))
                % Cnew(i,j)=Cnew(i,j)+epsFC*(FCempSubj(i,j)-FCsim(i,j));
                if str2num(string(extract(OptimStrategy{1,OS},1))) ==2
                    Beta=0.20;
                elseif OS==8
                    Beta=0.80;
                else
                   Beta=0.50;
                end 
                
                WeightFC=Beta*(epsFC*(FCempSubj(i,j)-FCsim(i,j)));
                WeightTM=(1-Beta)*(epsFC*(tmempSubj-tm));
                Cnew(i,j)=Cnew(i,j)+WeightFC+WeightTM;

                if Cnew(i,j)<0
                    Cnew(i,j)=0;
                end
            end
            end
        end
        Cnew = Cnew/max(max(Cnew))*maxC;
    end
    BOLD4All{sesh}(nsub,:,1:size(BOLD,2))=BOLD;
    FCempAllSubj{sesh}(nsub,:,:)=FCempSubj;   % the FC per subj - can also compute group FC as the mean across subjs
    FCsimAllSubj{sesh}(nsub,:,:)=FCsim;
    Ceff_sub{sesh}(nsub,:,:)=Cnew;   % the fitted GEC
    Error(nsub,sesh)=errornow;  % save final error
    errorFCtmAllSubj{nsub, sesh}=errorFCtm; % save error across iterations - saves as a cell so that way we see the # iterations that occurred if early stopping is left  
    fittFC_sub(nsub,sesh)=corr2(FCempSubj(Isubdiag),FCsim(Isubdiag));  % save correlation between sim and emp FC

% Temporary save
cd(outputpath)
tmpsaveName=strcat('/rds/general/project/dan_data_tide/live/IndivTIStim/PhenomModels/Outputs/Ceff_Optim', OptimStrategy{1,OS},'_S',num2str(nsub) ,'V',num2str(sesh),'_',num2str(numParcels),'Parcels_8mmSmooth.mat');
save(tmpsaveName)

end
end

%% Save outputs
cd(outputpath)
saveName=strcat('/rds/general/project/dan_data_tide/live/IndivTIStim/PhenomModels/Outputs/Ceff_Optim', OptimStrategy{1,OS},'_',num2str(numParcels), 'Parcels_8mmSmooth.mat');
save(saveName,'-v7.3')

end
