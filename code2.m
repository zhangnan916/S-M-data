% Load EEG data using EEGLAB
% Author: Based on EEG.history
clc;
clear all;

% Load channel information
load('D:\EEG_128channels_ERP_lanzhou_2015\chan_info_egi_128.mat');

% Load subject list from Excel
[~, ~, sub] = xlsread('D:\EEG_128channels_ERP_lanzhou_2015\subjects_information.xlsx', 'Sheet1', 'A2:A54');

% Loop over all subjects
for i = 1:53
    % Load raw EEG file
    EEG = pop_readegi(['D:\EEG_128channels_ERP_lanzhou_2015\' cell2mat(sub(i)) '.raw'], [], [], 'auto');

    % Set channel locations
    EEG.chanlocs = chanlocs;
    EEG.chaninfo = chaninfo;

    % Band-pass filtering
    EEG = pop_eegfiltnew(EEG, 0.3, [], 2750, 0, [], 1);   % High-pass at 0.3 Hz
    EEG = pop_eegfiltnew(EEG, [], 100, 34, 0, [], 1);     % Low-pass at 100 Hz

    % Notch filter at 50 Hz
    EEG.data = notchfilter(double(EEG.data), EEG.srate, 50);

    % Re-reference
    EEG = pop_reref(EEG, []);

    % Epoch extraction: from -0.1s to 0.5s for all cue types
    EEG = pop_epoch(EEG, {'fcue', 'hcue', 'scue'}, [-0.1 0.5], ...
                    'newname', 'EGI file epochs', 'epochinfo', 'yes');

    % Baseline correction (-100 to 0 ms)
    EEG = pop_rmbase(EEG, [-100 0]);

    % Run ICA
    EEG = pop_runica(EEG, 'extended', 1, 'interupt', 'on');

    % Artifact removal using ADJUST
    [aic] = ADJUST_EEG(EEG);
    EEG = pop_subcomp(EEG, aic, 0);

    EEG = eeg_checkset(EEG);

    % Re-extract epochs for each cue separately (0 to 0.5s)
    xx  = pop_epoch(EEG, {'fcue'}, [0 0.5], 'newname', 'EGI file epochs', 'epochinfo', 'yes');
    fcue = double(xx.data);

    xx1 = pop_epoch(EEG, {'hcue'}, [0 0.5], 'newname', 'EGI file epochs', 'epochinfo', 'yes');
    hcue = double(xx1.data);

    xx2 = pop_epoch(EEG, {'scue'}, [0 0.5], 'newname', 'EGI file epochs', 'epochinfo', 'yes');
    scue = double(xx2.data);

    % Save the processed data
    save(['D:\EEG_channels_2015\' cell2mat(sub(i)) '.mat'], 'fcue', 'hcue', 'scue');
end
