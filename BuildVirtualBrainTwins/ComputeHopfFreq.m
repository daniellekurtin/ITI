clear all
close all
addpath(genpath('C:\Users\dk818\OneDrive - Imperial College London\ICRF\Modelling'))
outputpath='C:\Users\dk818\OneDrive - Imperial College London\ICRF\Modelling\Outputs';

% Load brain data
for numParcels=[200]
% load(strcat('df_',num2str(numParcels),'ROIs.mat'))
load(strcat('df_',num2str(numParcels),'ROIs_8mmSmooth.mat'))
[NSUB,NSESH]=size(df);
[Tmax,NPARCELLS]=size(df{1,1});
TR=1.25;  % Repetition Time (seconds)
countAbove=0;
countBelow=0;
Isubdiag = find(tril(ones(NPARCELLS),-1));
clearvars f_diff_sub fce f_diff
fce1=zeros(NSUB,NSESH,NPARCELLS,NPARCELLS);

for sub=1:NSUB
for sesh=1:NSESH
    
    disp(strcat('Running S',num2str(sub),'V',num2str(sesh)))

    ts=table2array(df{sub,sesh})';
    TT=Tmax;
    tss=zeros(NPARCELLS,Tmax);
    Ts = TT*TR;
    freq = (0:TT/2-1)/Ts;
    nfreqs=length(freq);
    
    for seed=1:NPARCELLS
        x(seed,:)=detrend(ts(seed,:)-mean(ts(seed,:)));
        pw = abs(fft(x(seed,:)));
        PowSpect2(:,seed) = pw(1:floor(TT/2)).^2/(TT/TR);
        PowSpect(:,seed,sub) = PowSpect2(:,seed); 
    end

    fce1(sub,sesh,:,:)=corrcoef(x','rows','pairwise');

    for seed=1:NPARCELLS
        Power_Areassub(:,seed)=gaussfilt(freq,PowSpect2(:,seed)',0.01);
    end

    [maxpowdata,index]=max(Power_Areassub);
    f_diff_sub2 = freq(index);
    f_diff_sub2(find(f_diff_sub2==0))=mean(f_diff_sub2(find(f_diff_sub2~=0)));

    % CUT OFF FREQS THAT ARE TOO LOW OR HIGH
    for seed=1:NPARCELLS
        if f_diff_sub2(seed) <0.02
            f_diff_sub2(seed)=0.02;
            countBelow=countBelow+1;
        elseif f_diff_sub2(seed)>0.07
            f_diff_sub2(seed)=0.07;
            countAbove=countAbove+1;
        end
    end
    f_diff_sub(sub,sesh,:)=f_diff_sub2;

end
end
fce=squeeze(mean(fce1,'all'));

Power_Areas=squeeze(mean(PowSpect,3));
for seed=1:NPARCELLS
    Power_Areas(:,seed)=gaussfilt(freq,Power_Areas(:,seed)',0.01);
end

% WAY TO CHECK - check peak for each subj for regions that have weird power
% freq
% plot(Power_Areassub(:,106))

[maxpowdata,index]=max(Power_Areas);
f_diff = freq(index);
f_diff(find(f_diff==0))=mean(f_diff(find(f_diff~=0)));

cd(outputpath)
save(strcat('hopf_freq_sch',num2str(numParcels),'_8mmSmooth.mat'),'f_diff_sub','fce','f_diff')

end

%% Lil QA plot
% x=1:1:size(f_diff_sub,3);
% figure()
% hold on
% for s=1:NSUB
% % plot(f_diff_sub(s,:))
% scatter(x,squeeze(f_diff_sub(s,1,:)),'filled')
% end
