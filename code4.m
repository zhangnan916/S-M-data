function [aic] = ADJUST_EEG(EEG)
% ADJUST_EEG - Automatically identify ICA components related to eye artifacts
% This script identifies blink, vertical/horizontal eye movement, and discontinuity-related ICs
% The returned output 'aic' is a list of component indices to be removed.

% Author: Based on the original ADJUST algorithm

%% Determine number of epochs
if length(size(EEG.data)) == 3
    num_epoch = size(EEG.data, 3);
else
    num_epoch = 0;
end

%% Check ICA activations
if isempty(EEG.icaact)
    disp('EEG.icaact not present. Recomputed from data.');
    if length(size(EEG.data)) == 3
        EEG.icaact = reshape(EEG.icaweights * EEG.icasphere * ...
                     reshape(EEG.data, size(EEG.data,1), []), ...
                     size(EEG.data,1), size(EEG.data,2), size(EEG.data,3));
    else
        EEG.icaact = EEG.icaweights * EEG.icasphere * EEG.data;
    end
end

%% Normalize IC topographies
topografie = EEG.icawinv';
for i = 1:size(EEG.icawinv, 2)
    ScalingFactor = norm(topografie(i, :));
    topografie(i, :) = topografie(i, :) / ScalingFactor;
    if num_epoch > 0
        EEG.icaact(i, :, :) = ScalingFactor * EEG.icaact(i, :, :);
    else
        EEG.icaact(i, :) = ScalingFactor * EEG.icaact(i, :);
    end
end

%% Channel position check
nopos_channels = [];
for el = 1:length(EEG.chanlocs)
    loc = EEG.chanlocs(el);
    if all([isempty(loc.X), isempty(loc.Y), isempty(loc.Z), ...
            isempty(loc.theta), isempty(loc.radius)])
        nopos_channels = [nopos_channels el];
    end
end
if ~isempty(nopos_channels)
    warning(['Channels ' num2str(nopos_channels) ...
        ' have incomplete location information. They will NOT be used to compute ADJUST spatial features']);
end
pos_channels = setdiff(1:length(EEG.chanlocs), nopos_channels);

%% Feature extraction
GDSF = compute_GD_feat(topografie, EEG.chanlocs(1,pos_channels), size(EEG.icawinv,2));
[SED, medie_left, medie_right] = computeSED_NOnorm(topografie, EEG.chanlocs(1,pos_channels), size(EEG.icawinv,2));
[SAD, var_front, var_back, mean_front, mean_back] = computeSAD(topografie, EEG.chanlocs(1,pos_channels), size(EEG.icawinv,2));
diff_var = var_front - var_back;

%% Epoch-based temporal features
K = zeros(num_epoch, size(EEG.icawinv,2));    % Kurtosis
Vmax = zeros(num_epoch, size(EEG.icawinv,2)); % Variance

for i = 1:size(EEG.icawinv,2)
    for j = 1:num_epoch
        Vmax(j,i) = var(EEG.icaact(i,:,j));
        K(j,i) = kurt(EEG.icaact(i,:,j));
    end
end

%% TK - Temporal Kurtosis
meanK = zeros(1, size(EEG.icawinv,2));
for i = 1:size(EEG.icawinv,2)
    if num_epoch > 100
        meanK(i) = trim_and_mean(K(:,i));
    else
        meanK(i) = mean(K(:,i));
    end
end

%% MEV - Maximum Epoch Variance
maxvar = zeros(1, size(EEG.icawinv,2));
meanvar = zeros(1, size(EEG.icawinv,2));
for i = 1:size(EEG.icawinv,2)
    if num_epoch > 100
        maxvar(i) = trim_and_max(Vmax(:,i)');
        meanvar(i) = trim_and_mean(Vmax(:,i)');
    else
        maxvar(i) = max(Vmax(:,i));
        meanvar(i) = mean(Vmax(:,i));
    end
end
nuovaV = maxvar ./ meanvar;

%% Threshold estimation using EM algorithm
[soglia_K, ~, ~]     = EM(meanK);
[soglia_SED, ~, ~]   = EM(SED);
[soglia_SAD, ~, ~]   = EM(SAD);
[soglia_GDSF, ~, ~]  = EM(GDSF);
[soglia_V, ~, ~]     = EM(nuovaV);

%% Artifact classification
horiz = intersect(intersect(find(SED >= soglia_SED), find(medie_left .* medie_right < 0)), ...
                  find(nuovaV >= soglia_V));
vert = intersect(intersect(find(SAD >= soglia_SAD), find(medie_left .* medie_right > 0)), ...
                 intersect(find(diff_var > 0), find(nuovaV >= soglia_V)));
blink = intersect(intersect(find(SAD >= soglia_SAD), find(medie_left .* medie_right > 0)), ...
                  intersect(find(meanK >= soglia_K), find(diff_var > 0)));
disc = intersect(find(GDSF >= soglia_GDSF), find(nuovaV >= soglia_V));

%% Combine all detected ICs
aic = unique([blink, horiz, vert]);  % Note: 'disc' not included in final rejection
return
