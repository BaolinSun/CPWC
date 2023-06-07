
clear all;
close all;

N_channels = 128;
N_pixels = 65536;
N_samples = 2500;
N_waves = 15;
N_frames = 2;

time = linspace(0,5e-05,N_samples);
data = random('normal',0,1,[N_samples N_channels N_waves N_frames]);
delay = random('unif',0,5e-05,[N_pixels N_channels N_waves]);

%% interp1
tic
for n_wave=1:N_waves
    for n_rx=1:N_channels
        temp = interp1(time,data(:,n_rx,n_wave,:),delay(:,n_rx,n_wave),'linear',0);
    end
end
[user,sys] = memory;
fprintf('interp1: %0.2f s, %0.2f GB\n',toc,(sys.PhysicalMemory.Total-sys.PhysicalMemory.Available)/2^30);

%% interpn
channel_idx_in = 1:N_channels;
frames_idx_in = 1:N_frames;
delay_in = repmat(delay, [1 1 1 N_frames]);
channel_idx_out = repmat(channel_idx_in, [N_pixels, 1, N_frames]);
frames_idx_out = repmat(permute(frames_idx_in,[1 3 2]), [N_pixels, N_channels, 1]);

tic
for n_wave=1:N_waves
    temp = interpn(time,channel_idx_in,frames_idx_in,...
               squeeze(data(:,:,n_wave,:)),...
               squeeze(delay_in(:,:,n_wave,:)),channel_idx_out,frames_idx_out,...
               'linear',0);
end
[user,sys] = memory;
fprintf('interpn: %0.2f s, %0.2f GB\n',toc,(sys.PhysicalMemory.Total-sys.PhysicalMemory.Available)/2^30);

%% griddedInterpolant 
F = griddedInterpolant({time,1:N_channels,1:N_waves,1:N_frames},data,'linear');
delay_in = repmat(delay, [1 1 1 N_frames]);
channel_in = repmat(1:N_channels, [N_pixels 1 N_waves N_frames]);
waves_in = repmat(permute(1:N_waves,[1 3 2]), [N_pixels, N_channels, 1, N_frames]);
frames_in = repmat(permute(1:N_frames,[1 3 4 2]), [N_pixels, N_channels, N_waves, 1]);

tic
temp = F(delay_in,channel_in,waves_in,frames_in);
[user,sys] = memory;
fprintf('griddedInterpolant: %0.2f s, %0.2f GB\n',toc,(sys.PhysicalMemory.Total-sys.PhysicalMemory.Available)/2^30);
