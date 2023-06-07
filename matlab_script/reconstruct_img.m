function res = reconstruct_img(raw_data_dir, usimage_path)
    F_number = 1.75;

    Na = 75;

    CPW = zeros(6494, 128, 75, 1);  % impulse response channel data

    for n=1:Na
        %  Load the result

%        cmd=['load raw_data/raw_angle', num2str(n),'.mat']
%        eval(cmd);
        load([raw_data_dir,'/raw_angle_', num2str(n),'.mat'], 'raw_data')

        CPW(1:size(raw_data.data,1), :, n, 1) = raw_data.data;

        seq(n) = uff.wave();
        seq(n).probe = raw_data.sequence.probe;
        seq(n).source.azimuth = raw_data.sequence.source.azimuth;
        seq(n).source.distance = raw_data.sequence.source.distance;
        seq(n).sound_speed = raw_data.sequence.sound_speed;
        seq(n).delay = raw_data.sequence.delay;

    end

    channel_data = uff.channel_data();
    channel_data.sampling_frequency = raw_data.sampling_frequency;
    channel_data.sound_speed = raw_data.sound_speed;
    channel_data.initial_time = raw_data.initial_time;
    channel_data.pulse = raw_data.pulse;
    channel_data.probe = raw_data.probe;
    channel_data.sequence = seq;
    channel_data.data = CPW./max(CPW(:));

   sca=uff.linear_scan('x_axis',linspace(-19e-3,19e-3,256).', 'z_axis', linspace(5e-3,50e-3,256).');
    % sca=uff.linear_scan('x_axis',linspace(-50e-3,50e-3,512).', 'z_axis', linspace(5e-3,105e-3,512).');

    pipe=pipeline();
    pipe.channel_data=channel_data;
    pipe.scan=sca;

    pipe.receive_apodization.window=uff.window.tukey25;
    pipe.receive_apodization.f_number=F_number;

    b_data=pipe.go({midprocess.das() postprocess.coherent_compounding()});

    % if exist(save_dir,'dir')==0
    %     mkdir(save_dir);
    % end

    b_data.plot();
    b_data.save_as_png(usimage_path);

    res = 0;
end