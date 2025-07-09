clear all; clc;

%% part1: compute group-level ERP

Subj = [1:10]; %% subject numbers
for i = 1:length(Subj)
    setname = strcat(num2str(i),'_LH.set'); %% filename of set file
    setpath = 'D:\day23\day23\Example_data\'; %% filepath of set files (need to be changed)
    EEG = pop_loadset('filename',setname,'filepath',setpath); %% load the data
    EEG = eeg_checkset( EEG );
    EEG_avg(i,:,:) = squeeze(mean(EEG.data,3)); %% single-subject ERPs; EEG_avg dimension: subj*channel*time
end
save('Group_level_ERP.mat','EEG_avg');  %% save the data of subjects

%% part2: plot group-level ERP

Cz=13; % select the channel to plot (display maximum response)
figure;
mean_data=squeeze(mean(EEG_avg(:,Cz,:),1)); %% select the data at Cz, average across subjects, mean_data: 1*3000
plot(EEG.times, mean_data,'k','linewidth', 1.5); %% plot the waveforms
set(gca,'YDir','reverse'); %% reverse the direction of Y axis
axis([-500 1000 -15 10]);  %% define the region to display
title('Group-level at Cz','fontsize',16); %% specify the figure name
xlabel('Latency (ms)','fontsize',16); %% name of X axis
ylabel('Amplitude (uV)','fontsize',16);  %% name of Y axis

%% part3: plot the scalp maps at dominant peaks

N2_peak = 207; P2_peak=374; %% dominant peaks on waveforms
N2_interval = find((EEG.times>=197)&(EEG.times<=217)); %% define the N2 intervals [peak-10 peak+10]
P2_interval = find((EEG.times>=364)&(EEG.times<=384)); %% define the P2 intervals [peak-10 peak+10]

N2_amplitude = squeeze(mean(EEG_avg(:,:,N2_interval),3));  %% N2 amplitude for each subject and each channel
P2_amplitude = squeeze(mean(EEG_avg(:,:,P2_interval),3));   %% P2 amplitude for each subject and each channel

Central = [13 14];
N2_measures = squeeze(mean(mean(EEG_avg(:,Central,N2_interval),2),3));  %% N2 amplitude for each subject  (for statistics)
P2_measures = squeeze(mean(mean(EEG_avg(:,Central,P2_interval),2),3));   %% P2 amplitude for each subject  (for statistics)

figure;
subplot(121); topoplot(mean(N2_amplitude),EEG.chanlocs,'maplimits',[-15 15]); title('N2 Amplitude','fontsize',16); %% N2 scalp map (group-level)
subplot(122); topoplot(mean(P2_amplitude),EEG.chanlocs,'maplimits',[-15 15]); title('P2 Amplitude','fontsize',16); %% P2 scalp map (group-level)

%% part4: series of scalp mps

time_interval = [0:100:500]; %% specify the time intervals to display (to be changed)
figure;
for i = 1:length(time_interval)
    latency_range = [time_interval(i) time_interval(i)+100]; %% lower and upper limits
    latency_idx = find((EEG.times>=latency_range(1))&(EEG.times<=latency_range(2))); %% interval of the specific regions
    Amplitude = squeeze(mean(mean(EEG_avg(:,:,latency_idx),1),3)); %% 1*channel (averaged across subjects and interval)
    subplot(2,3,i);
    topoplot(Amplitude,EEG.chanlocs,'maplimits',[-10 10]); %% topoplot(Amplitude,EEG.chanlocs);
    setname=strcat(num2str(latency_range(1)),'--',num2str(latency_range(2)),'ms'); %% specify the name of subplots
    title(setname,'fontsize',16); %% display the names of subplots
end
