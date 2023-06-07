function channel_data = read_CPWC(h,N_planewaves,N_frames)

% Load the filename
load([h.data_folder,'/Sequence.mat']);
files_in_folder = dir(h.data_folder);

if nargin < 2
    N_planewaves = numel(Tx);   % number of plane waves
end


idx = 1;
for i = 1:length(files_in_folder)
    if isempty(strfind(files_in_folder(i).name,'_BDATA_RF')) == 0 && strfind(files_in_folder(i).name,'_BDATA_RF') && isempty(strfind(files_in_folder(i).name,'._'))
        frame_filename{idx} =  files_in_folder(i).name;
        split_string = strsplit(frame_filename{idx},'_layer0');
        frame_idx_list(idx) = str2num(split_string{1});
        idx = idx+1;
    end
end

[dummy,idx_sorted] = sort(frame_idx_list);
frame_filename_sorted = frame_filename(idx_sorted);

%Find number of frames
if nargin < 3 %Number of frames can be a parameter
    N_frames = length(frame_filename_sorted)
end

% Reading paramteres
channel_data = uff.channel_data();
channel_data.sampling_frequency = double(System.Parameters.sampleFreqMHz*10^6);;
channel_data.sound_speed = double(Parameter{1}.speedOfSoundMps);
channel_data.initial_time = 0;

% Create probe
channel_data.probe = uff.linear_array();
channel_data.probe.N = double(System.Transducer.elementCnt);
channel_data.probe.pitch = double(System.Transducer.elementPitchCm)/100;

% Save Pulse
channel_data.pulse = uff.pulse();
channel_data.pulse.center_frequency = double(Tw{1}.freqMHz*10^6);

% Read data
load([h.data_folder,'/',frame_filename_sorted{1}]);
var_names = who('AdcData_frame*'); % Names on variables containing RF data
data_initial = eval(var_names{1});

%% SEQUENCE GENERATION
angles = deg2rad(double(Roi{1}.steerAngleDegA));
for n=1:N_planewaves
    seq(n)=uff.wave();
    seq(n).probe=channel_data.probe;
    seq(n).source.azimuth=angles(n);
    seq(n).source.distance=Inf;
    seq(n).sound_speed=channel_data.sound_speed;
end
channel_data.sequence = seq;

no_samples = size(data_initial,2); % Number of samples in each sequence
time = [0:1/channel_data.sampling_frequency:(no_samples-1)/channel_data.sampling_frequency]'; %Set up initial time

%% Buffer for all data
all_data = single(zeros(no_samples,channel_data.N_elements,size(seq,2),N_frames));
for frame = 1:N_frames
    load([h.data_folder,'/',frame_filename_sorted{frame}]);
    data = eval(var_names{1});
    for transmission = 1:N_planewaves
        %Read RF data
        rfData = single([data(:,:,1+((transmission-1)*4))' data(:,:,2+((transmission-1)*4))' ...
            data(:,:,3+((transmission-1)*4))' data(:,:,4+((transmission-1)*4))']);
        
        %Compensate for steering regarding t0
        D = abs(channel_data.probe.geometry(1,1)-channel_data.probe.geometry(end,1));
        q = abs((D/2)*sin(channel_data.sequence(transmission).source.azimuth));
        
        channel_data.sequence(transmission).delay = -(q./channel_data.sound_speed-System.Transducer.delayOffsetUsec*10^-6);
        
        % build the dataset
        all_data(:,:,transmission,frame)=rfData;
    end
end

if channel_data.pulse.center_frequency > 10*10^6 %Only filter for L8-17 probe
    fc = channel_data.pulse.center_frequency;
    fs = channel_data.sampling_frequency;
    F = [fc-8*10^6 fc-7.5*10^6 fs/2-10^6 fs/2-0.5*10^6]
    [all_data] = tools.band_pass(all_data,channel_data.sampling_frequency,F);
else
    all_data = h.bandpass_filter_data(all_data,channel_data,0);
end
channel_data.data = single(all_data);
end