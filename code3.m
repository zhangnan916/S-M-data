function [filt] = notchfilter(dat, Fs, Fl, N)
% NOTCHFILTER - Line noise reduction filter for EEG/MEG data
%
% Syntax:
%   filt = notchfilter(dat, Fs, Fl)
%   filt = notchfilter(dat, Fs, Fl, N)
%
% Inputs:
%   dat - Data matrix (channels ¡Á time) [Nchans ¡Á Nsamples]
%   Fs  - Sampling frequency (Hz)
%   Fl  - Line noise frequency (Hz), e.g., 50
%         If a scalar (e.g., 50), filters out [48 52] Hz
%         If a two-element vector [low high], filters that band
%   N   - (Optional) Filter order, default is 4
%
% Output:
%   filt - Filtered data matrix (same size as input)
%
% References:
%   Original      (c) 2003, Pascal Fries
%   Modifications (c) 2003, Robert Oostenveld
%   This file is part of FieldTrip: http://www.ru.nl/neuroimaging/fieldtrip

if nargin < 4
    N = 4; % Default filter order
end

% Get data dimensions
Nchans   = size(dat, 1);
Nsamples = size(dat, 2);

% Nyquist frequency
Fn = Fs / 2;

% Define stopband frequency
if length(Fl) == 1
    % Default notch width of ¡À2 Hz around target frequency
    Fl = [Fl - 2, Fl + 2];
end

% Design Butterworth notch filter
[B, A] = butter(N, [min(Fl)/Fn, max(Fl)/Fn], 'stop');

% Apply zero-phase filtering (filtfilt for no phase distortion)
filt = filtfilt(B, A, dat')';  % Transpose for filtering along time, then back

end
