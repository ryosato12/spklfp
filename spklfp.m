%%
clear;
close all;

rng('default'); 

%%
clc 
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaulttextInterpreter','latex','defaultAxesFontSize',16) 
format long

%% 
load('classifyUnits.mat')

% disp(Tall); 
% head(Tall);

%%
rats = {
% 'tk0056-171014-04', 'Buzsaki256.mat',[0 1200]+20;
% 'tk0056-171015-04', 'Buzsaki256.mat',[0 1200]+15;
% 'tk0062-171129-04', 'Buzsaki256.mat',[0 1200]+10;
% 'tk0062-171130-04', 'Buzsaki256.mat',[0 1200]+30;
% 'tk0064-180206-04','Buzsaki256.mat',[0 1200]+15;
% 'tk0064-180207-04','Buzsaki256.mat',[0 1200]+15;
% 'tk0067-180426-04', 'Buzsaki256.mat',[0 1200]+10;
% 'tk0067-180427-04','Buzsaki256.mat',[0 1200]+20;
% 'tk0068-180530-04','A8x32Edge.mat', [0 1200]+20;
% 'tk0068-180531-04','A8x32Edge.mat', [0 1200]+15;
% 'tk0069-180726-04','A8x32_5mm_35_300_160.mat',[0 1200]+10;
% 'tk0069-180727-04','A8x32_5mm_35_300_160.mat',[0 1200]+40;
'tk0070-180829-04','A8x32_5mm_35_300_160.mat',[0 1200]+12; % best performing
% 'tk0070-180830-04','A8x32_5mm_35_300_160.mat',[0 1200]+10;
% 'tk0072-181025-04','A8x32Edge.mat',[0 1200]+10;
% 'tk0072-181026-04','A8x32Edge.mat',[0 1200]+12;
% 'tk0074-181219-04','A8x32Edge.mat',[0 1200]+10;
% 'tk0074-181220-04','A8x32Edge.mat',[0 1200]+20;
% 'tk0075-190326-04','A8x32Edge.mat',[0 1200]+10;
% 'tk0075-190327-04','A8x32Edge.mat',[0 1200]+15;
% 'tk0076-190424-04','A8x32Edge.mat',[0 1200]+10;
% 'tk0076-190425-04','A8x32Edge.mat',[0 1200]+15
};

basePath = 'data';   
savePath = 'spk_LFP'; 

%%
for k = 1:size(rats,1)
    %%
    session = rats{k,1};  % linear track data
    mapfile = rats{k,2};  % channelmap ('Buzsaki256.mat', 'A8x32Edge.mat', or 'A8x32_5mm_35_300_160.mat')
    trange = rats{k,3};   % time range (sec)

    %%  
    idx = strfind(session,'-');
    rat = session(1:idx(1)-1);
    day = session(idx(1)+1:idx(2)-1);
    sess= session(idx(2)+1:end);

    %% 
    T = readtable('histology_test.xlsx','Sheet',rat,'Range','A1:H32');

    TT = T.Variables;

    lfpSUB = TT == 1;
    lfpSUBmol = TT == 2; 
    lfpCA1 = TT == 3;
    lfpCA3 = TT == 5; 

    % indices for ROI
    lfpSUB = reshape(fliplr(lfpSUB),[],1); 
    lfpSUBmol = reshape(fliplr(lfpSUBmol),[],1);
    lfpCA1 = reshape(fliplr(lfpCA1),[],1); 
    lfpCA3 = reshape(fliplr(lfpCA3),[],1);

    %% Local-field potential (1250-Hz)
    % lfp file path
    lfpPath = fullfile(basePath,rat,[day '-' sess]);

    % load lfp
    [lfp,t,~] = loadLfp(lfpPath,[],mapfile);
    lfp = lfp(1:256,:); 
    idx = trange(1)<=t & t<=trange(2);
    lfp = lfp(:,idx);
    t   = t(idx)-trange(1);

    %% Decimate by indlfp factor  
    indlfp = 8;
    Fs = 1250/indlfp;
    X = []; 
    for i = 1:size(lfp,1)
        X(i,:) = decimate(lfp(i,:),indlfp);
    end

    %% 
    X_wt = [];
    X_1 = X(:,1:end-1);
    for i = 1:size(X_1,1)
        X_wt(i,:,:) = cwt(X_1(i,:), ...
            'amor', ...
            Fs, ...
            VoicesPerOctave=2, ...
            FrequencyLimits=[1 Fs/2]); 
    end 

    %% Spike counts 
    % spike file path
    spkPath = fullfile(basePath,rat,[day '_sorted']);

    % load spike data 
    [ts,clu,cluOri,~,roi]=loadSpkHist(spkPath,[],session);

    % trim time range
    idx = trange(1)<ts & ts<trange(2);
    ts  = ts(idx)-trange(1);
    clu = clu(idx);

    %%
    tmp = ceil(ts*(1250/indlfp)); 
    nbins{1} = 1:size(X_1,2);
    nbins{2} = 1:max(clu);
    [spk_count,~] = hist3([tmp clu],nbins);

    % idx=roi==3; 
    % spk = spk(:,idx);

    %%
    % spk_rate = [];
    % windowSize = 50;
    % for i = 1:size(spk_count,2)
    %     spk_rate(:,i) = smoothdata(spk_count(:,i), 'movmean', windowSize);
    % end

    %% 
    spk_wt = [];
    for i = 1:size(spk_count,2)
        [spk_wt(i,:,:),f] = cwt(spk_count(:,i), ...
            'amor', ...
            Fs, ...
            VoicesPerOctave=2, ...
            FrequencyLimits=[1 Fs/2]); 
    end

    %% Linearize position 
    % whl file path
    whlPath = fullfile(basePath,rat,[day '-' sess]);
    
    % position scaling
    % 0.5217 (230pixels = 120cm) for open field, T-maze, Zigzag maze
    % 0.4130 (569pixels = 235cm) for linear track
    pixel2cm = 0.4130;
    
    % load pos (cm)
    % position is sampled at 39.0625 Hz (25.6 ms/sample)
    [t,x1,y1,x2,y2] = loadWhl(whlPath,pixel2cm);
    
    % flip y data
    y1 = -y1;
    y2 = -y2;

    % trim position data
    if trange(1)<t(1) || t(end)<trange(2)
        error 'Error: trange setting is out of range!'
    end
    idx = trange(1)<=t & t<trange(2);
    x1(~idx)=[];
    y1(~idx)=[];
    x2(~idx)=[];
    y2(~idx)=[];
    t(~idx)=[];
    
    % use the better-tracking LED
    if sum(x1==-1)>sum(x2==-1)
        posx = x2;
        posy = y2;
    else
        posx = x1;
        posy = y1;
    end
    post = t;
 
    pos = [x1 y1 x2 y2];

    %% Interpolate position 
    tmp = [];
    for i = 1:size(pos,2)
        tmp(:,i) = interp(pos(:,i),32/indlfp); 
    end
    [~,rpos,~,trial_num] = fix_position(tmp); 

    %% Splitting data into left moving even/odd trials 
    left_trials = rpos<=median(rpos);
    right_trials = rpos>=median(rpos);
    
    left_even_idx = find(mod(trial_num(:,left_trials),2) == 0);
    left_odd_idx = find(mod(trial_num(:,left_trials),2) == 1);
    
    right_even_idx = find(mod(trial_num(:,right_trials),2) == 0);
    right_odd_idx = find(mod(trial_num(:,right_trials),2) == 1);

    %% 
    lambda=10.^(-5:0); 
    k_fold=5; 
    bootstrap=0;

    CC_train = [];
    MSE_train = []; 
    wts = cell(0); 
    for i = 1:size(spk_wt,2)
        [cc,mse,wt] = ridgeCross( ...
            squeeze(spk_wt(1:size(spk_wt,1),i,left_even_idx)).', ...
            squeeze(X_wt(:,i,left_even_idx)).', ...
            k_fold,lambda,bootstrap);
        CC_train(i,:,:) = cc; 
        MSE_train(i,:,:) = mse; 
        wts{i} = wt; 
    end

    %% 
    CC_left_even = []; 
    for i = 1:size(CC_train,3)
        for j = 1:size(CC_train,1)
            [~,idx] = min(MSE_train(j,:,i));
            CC_left_even(i,j) = CC_train(j,idx,i); 
        end
    end

    %% 
    y_left_odd = [];
    y_right_odd = []; 
    for i = 1:size(spk_wt,2)
        for j = 1:size(spk_wt,1)
            [~,idx] = min(MSE_train(i,:,j));

            y_left_odd(j,i,:) = squeeze(X_wt(:,i, ...
            left_odd_idx)).'*squeeze(wts{i}(idx,:,j)).';

            y_right_odd(j,i,:) = squeeze(X_wt(:,i, ...
            right_odd_idx)).'*squeeze(wts{i}(idx,:,j)).';
        end 
    end 
    
    %%
    CC_left_odd = [];
    CC_right_odd = []; 
    for i = 1:size(y_left_odd,2)
        for j = 1:size(y_left_odd,1)
            CC_left_odd(i,j) = diag(corr(squeeze(y_left_odd(j,i,:)), ...
                squeeze(spk_wt(j,i,left_odd_idx))));
        end
    end

    for i = 1:size(y_right_odd,2)
        for j = 1:size(y_right_odd,1)
            CC_right_odd(i,j) = diag(corr(squeeze(y_right_odd(j,i,:)), ...
                squeeze(spk_wt(j,i,right_odd_idx))));
        end
    end

    %% 
    % find(Tall.sess == [rat,'-',day]); 

    units = table2array(Tall(513:563,'SUBp'));
    
    %% Plot results 
    % CA1 place cells appear sooner than SUB CA1 place cells... 

    figure(1); clf; 

    subplot(2,1,1)
    plot(f,abs(CC_left_odd(:,table2array(Tall(513:563,'SUBp')))),'Color','b')
    hold on
    plot(f,abs(CC_left_odd(:,table2array(Tall(513:563,'CA1p')))),'Color','r')
    hold off
    xlabel('Frequency (Hz)')
    ylabel('Correlation')
    title([rat, ': leftwards trials'])
    xlim([f(end) f(1)])

    subplot(2,1,2)
    plot(f,abs(CC_right_odd(:,table2array(Tall(513:563,'SUBp')))),'Color','b')
    hold on
    plot(f,abs(CC_right_odd(:,table2array(Tall(513:563,'CA1p')))),'Color','r')
    hold off
    xlabel('Frequency (Hz)')
    ylabel('Correlation')
    title([rat, ': rightwards trials'])
    xlim([f(end) f(1)])

end

% parallel computing?

