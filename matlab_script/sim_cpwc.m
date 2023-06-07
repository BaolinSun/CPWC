function res = sim_cpwc(start_angle, end_angle, phantom_path, save_dir)

    if exist(save_dir,'dir')==0
        mkdir(save_dir)
    end    

    c0=1540;     % Speed of sound [m/s]
    fs=20.832e6;    % Sampling frequency [Hz]
    dt=1/fs;     % Sampling step [s] 
    
    %% field II initialisation
    % 
    % Next, we initialize the field II toolbox. Again, this only works if the 
    % Field II simulation program (<field-ii.dk>) is in MATLAB's path. We also
    % pass our set constants to it.
    path(path, 'C:\Users\Administrator\Documents\MATLAB\CPWC\matlab_script\Field_II_ver_3_30_windows');

%    path(path, '/home/hrzy/Desktop/matlab/Field_II_ver_3_30_linux')
    field_init(-1);
    set_field('c', c0);              % Speed of sound [m/s]
    set_field('fs', fs);             % Sampling frequency [Hz]
    set_field('use_rectangles', 1);  % use rectangular elements
    % set_field('no_ascii_output', 1);


    probe = uff.linear_array();
    f0                      = 5.208e+06;      % Transducer center frequency [Hz] 5.208MHz 5.1333e+06MHz
    lambda                  = c0/f0;           % Wavelength [m]
    probe.element_height    = 5e-3;            % Height of element [m]
    probe.pitch             = 0.300e-3;        % probe.pitch [m]
    kerf                    = 0.03e-03;        % gap between elements [m]
    probe.element_width     = probe.pitch-kerf;% Width of element [m]
    lens_el                 = 20e-3;           % position of the elevation focus
    probe.N                 = 128;             % Number of elements
    pulse_duration          = 2.5;             % pulse duration [cycles]



    pulse = uff.pulse();
    pulse.fractional_bandwidth = 0.67;        % probe bandwidth [1]
    pulse.center_frequency = f0;
    t0 = (-1/pulse.fractional_bandwidth/f0): dt : (1/pulse.fractional_bandwidth/f0);
    impulse_response = gauspuls(t0, f0, pulse.fractional_bandwidth);
    impulse_response = impulse_response-mean(impulse_response); % To get rid of DC

    te = (-pulse_duration/2/f0): dt : (pulse_duration/2/f0);
    excitation = 30 * square(2*pi*f0*te+pi/2);
    one_way_ir = conv(impulse_response,excitation);
    two_way_ir = conv(one_way_ir,impulse_response);
    lag = length(two_way_ir)/2+1;


    % figure;
    % plot((0:(length(two_way_ir)-1))*dt -lag*dt,two_way_ir); hold on; grid on; axis tight
    % plot((0:(length(two_way_ir)-1))*dt -lag*dt,abs(hilbert(two_way_ir)),'r')
    % plot([0 0],[min(two_way_ir) max(two_way_ir)],'g');
    % legend('2-ways pulse','Envelope','Estimated lag');
    % title('2-ways impulse response Field II');
    
    %% Aperture Objects
    % Next, we define the the mesh geometry with the help of Field II's
    % *xdc_linear_array* function.

    noSubAz=round(probe.element_width/(lambda/8));        % number of subelements in the azimuth direction
    noSubEl=round(probe.element_height/(lambda/8));       % number of subelements in the elevation direction
    Th = xdc_focused_array (probe.N, probe.element_width, probe.element_height, kerf, lens_el, noSubAz, noSubEl, [0 0 Inf]); 
    Rh = xdc_focused_array (probe.N, probe.element_width, probe.element_height, kerf, lens_el, noSubAz, noSubEl, [0 0 Inf]);
    % Th = xdc_linear_array (probe.N, probe.element_width, probe.element_height, kerf, noSubAz, noSubEl, [0 0 Inf]); 
    % Rh = xdc_linear_array (probe.N, probe.element_width, probe.element_height, kerf, noSubAz, noSubEl, [0 0 Inf]); 

    % We also set the excitation, impulse response and baffle as below:
    xdc_excitation (Th, excitation);
    xdc_impulse (Th, impulse_response);
    xdc_baffle(Th, 0);
    xdc_center_focus(Th,[0 0 0]);
    xdc_impulse (Rh, impulse_response);
    xdc_baffle(Rh, 0);
    xdc_center_focus(Rh,[0 0 0]);
    
    %% Define plane wave sequence
    % Define the start_angle and number of angles
    F_number = 1.75;
    alpha_max = atan(1/2/F_number);
    Na=75;                                      % number of plane waves
    F=1;                                        % number of frames
    alpha=linspace(-alpha_max,alpha_max,Na);    % vector of angles [rad]


    %path_phantom = '/home/xuepeng/ultrasound/plane_wave_beamforming/PICMUS_SIM/database/simulation/resolution_distorsion/resolution_distorsion_simu_phantom.hdf5';
    %% Define phantom
    % Define some points in a phantom for the simulation
    phantom = us_phantom();
    phantom.read_file(phantom_path);

    % point_position_bak = phantom.sca;

%    step = 1;
%    point_position = phantom.sca(1:step:length(phantom.sca), :);
%    point_amplitudes = phantom.amp(1:step:length(phantom.amp));

    point_position = phantom.sca;
    point_amplitudes = phantom.amp;

%    load(phamtom_path,'sca')
%    load(phamtom_path,'amp')

%    point_position = sca;
%    point_amplitudes = amp;

    %% output data
    cropat=round(2*50e-3/c0/dt);    % maximum time sample, samples after this will be dumped
    CPW=zeros(cropat,probe.N,Na,F);  % impulse response channel data

    
    %% Compute CPW signals
    raw_data = uff.channel_data();
    for n = start_angle : end_angle
        cmd = ['Calculating angle ',num2str(n),' of ',num2str(Na)];
        disp(cmd);
        % transmit aperture
        xdc_apodization(Th,0,ones(1,probe.N));
        xdc_times_focus(Th,0,probe.geometry(:,1)'.*sin(alpha(n))/c0);
        
        % receive aperture
        xdc_apodization(Rh, 0, ones(1,probe.N));
        xdc_focus_times(Rh, 0, zeros(1,probe.N));

        % do calculation
        [v,t]=calc_scat_multi(Th, Rh, point_position, point_amplitudes);
        
        % build the dataset
        CPW(1:size(v,1), :, n, 1)=v;
        
        % Save transmit sequence
        seq=uff.wave();
        seq.probe=probe;
        seq.source.azimuth=alpha(n);
        seq.source.distance=Inf;
        seq.sound_speed=c0;
        seq.delay = -lag*dt+t;

        raw_data.sampling_frequency = fs;
        raw_data.sound_speed = c0;
        raw_data.initial_time = 0;
        raw_data.pulse = pulse;
        raw_data.probe = probe;
        raw_data.sequence = seq;
        raw_data.data = v;
        
        name = 'raw_data';

        save([save_dir, '/raw_angle_', num2str(n),'.mat'], name);

    end

    disp(['Calculate angle from ', num2str(start_angle),' to ', num2str(end_angle), ' finished...']);
    field_end();

% Return nothing
res = 99;

end