
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

%% interp1 GPU
time_gpu= gpuArray(time);
data_gpu= gpuArray(data);
delay_gpu= gpuArray(delay);

tic
for n_wave=1:N_waves
    for n_rx=1:N_channels
        temp = interp1(time_gpu,data_gpu(:,n_rx,n_wave,:),delay_gpu(:,n_rx,n_wave),'linear',0);
    end
end
[user,sys] = memory;
fprintf('interp1 GPU: %0.2f s, %0.2f GB\n',toc,(sys.PhysicalMemory.Total-sys.PhysicalMemory.Available)/2^30);
