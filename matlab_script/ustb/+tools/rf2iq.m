function [iq, true_fs_new, f_cutoff] = rf2iq(rf, fs, fdemod, fs_new, f_cutoff, varargin)
% RF2IQ Convert RF data to IQ and downsample.
%   Syntax:
%       >> iq = rf2iq(rf, fs, fdemod, fs_new, f_cutoff)
%   Inputs:
%       rf        -- N-D array of real-valued rf traces in column orientation
%       fs        -- sampling frequency [Hz]
%       fdemod    -- demodulation frequency [Hz]
%       fs_new    -- new sampling freq after resampling
%       f_cutoff  -- cutoff frequency above which the signal will be filtered before
%                    resampling
%
%   Outputs
%       iq        -- N-D array of downsampled demodulated iq signal
%       f_cutoff  -- Cutoff filtering frequency for iq signal (= fs for non-filtered)
%

% 2015-03-26 ?yvind K.-V. Standal

%% input parsing
if nargin < 5,
    % use factor 2 oversampling compared to Nyquist
    f_cutoff = fs_new/2;
end
doresample = false; % need to interpolate (not just downsample by integer rate)
if nargin < 4,
    Ndown = 1; % no downsampling
else
    Ndown = fs/fs_new;
    if abs(Ndown - round(Ndown)) < 0.01, 
        Ndown = round(Ndown); % snap resampling to nearest integer
    else
        doresample = true;
        warning('Downsampling by a non-integer rate is highly inefficient.');
    end
end
true_fs_new=fs/Ndown;

%% preallocation and stuff
siz = size(rf);
N = prod(siz(3:end)); % collapse dimensions beyond 3
ind_new = 1:Ndown:siz(1);
L = length(ind_new); % length of downsampled array

% preallocate final array
iq = zeros([L,siz(2:end)], 'like', 1i);

% demodulation mixing vector
mix = exp(-2i*pi*fdemod/fs*(0:siz(1)-1)');

% otherwise not really need for low-pass filtering
dofilter = (Ndown > 1 && f_cutoff < fdemod);
if dofilter, % make low-pass filter
    % calculate window parameters
    [M, Wn, beta, typ] = kaiserord([0.9,1]*f_cutoff, [1,0], [1e-3,1e-3], fs);
    b = fir1(M, Wn, typ, kaiser(M+1,beta), 'noscale'); % filter design
end

%% main loop
for nn = 1:N, % loop over higher dimensions
    H = bsxfun(@times, hilbert(rf(:,:,nn)), mix);
    if dofilter,
        H(:,:) = filtfilt(b, 1, H); % lowpass filter
    end
    if doresample,
        iq(:,:,nn) = interp1(H, ind_new, 'linear', 0);
    elseif Ndown >= 2,
        iq(:,:,nn) = H(ind_new,:); % downsample
    else
        iq(:,:,nn) = H;
    end
end
