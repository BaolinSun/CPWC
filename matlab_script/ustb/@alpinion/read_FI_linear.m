 function channel_data = read_FI_linear(h,N_frames,start_frame)
            
            if nargin < 3
                start_frame = 1;
            end
            
            % Load the filename
            load([h.data_folder,'/Sequence.mat']);
            files_in_folder = dir(h.data_folder);
            
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
            
            %%
            
            % Reading parameters
            channel_data = uff.channel_data();
            channel_data.sampling_frequency = double(System.Parameters.sampleFreqMHz*10^6);
            channel_data.sound_speed = double(Parameter{1}.speedOfSoundMps);
            channel_data.initial_time = 0;
            
            % Create probe
            channel_data.probe = uff.linear_array();
            channel_data.probe.N = double(System.Transducer.elementCnt);
            channel_data.probe.pitch = double(System.Transducer.elementPitchCm)*10^-2;
            channel_data.probe.element_width = double(System.Transducer.elementWidthCm)*10^-2;
            
            % Save Pulse
            channel_data.pulse = uff.pulse();
            channel_data.pulse.center_frequency = double(Tw{1}.freqMHz*10^6);
            
            %Need this to figure out how many scan lines
            load([h.data_folder,'/',frame_filename_sorted{1}]);
            var_names = who('AdcData_scanline*'); % Names on variables containing RF data
            
            %Find number of frames
            if nargin < 2 %Number of frames can be a parameter
                N_frames = length(frame_filename_sorted)
            end
            
            %% SEQUENCE GENERATION
            N=length(var_names);        % number of focused beams
            seq=uff.wave();
            for n=1:N
                seq(n)=uff.wave();
                seq(n).probe=channel_data.probe;
                seq(n).source.xyz=[double(Tx{n}.originPosCm(1))*10^-2 0 double(Tx{n}.focalPosCm)*10^-2];
                seq(n).sound_speed=channel_data.sound_speed;
            end
            channel_data.sequence = seq;
            
            
            %% Read data
            %Calculate the first offset so we know how many "extra" samples we need
            hardcoded_offset = (sqrt(seq(1).source.x^2+seq(1).source.z^2) - seq(1).source.z)/channel_data.sound_speed;
            no_samples = double(AlignedSampleNum+ceil(hardcoded_offset*channel_data.sampling_frequency)); % Number of samples in each sequence
            t_out = [0:1/channel_data.sampling_frequency:(no_samples-1)/channel_data.sampling_frequency]'; %Set up initial time
            N_receive_channels = length(Rx{1}.aperture);
            channel_data.N_active_elements = N_receive_channels;
            N_sc = double(ScanlineNum);
            RxMux = calculate_rx_mux(h,N_receive_channels,channel_data,N_sc);
            % Buffer for all data
            all_data = single(zeros(no_samples,channel_data.N_elements,size(seq,2),1));
            RF_Ch_data = zeros(no_samples, channel_data.probe.N_elements);
            %
            tools.workbar()
            frame_idx = 1;
            
            f_start = channel_data.pulse.center_frequency-2*channel_data.pulse.center_frequency/3;
            f_stop = channel_data.pulse.center_frequency+2*channel_data.pulse.center_frequency/3;
            f_transition = channel_data.pulse.center_frequency/5;
            %f_transition = 1.5e6;
            %interpolation_factor = 5;
            F = [f_start f_start+f_transition f_stop f_stop+f_transition];
            %%
            for frame = start_frame:N_frames+start_frame-1
                load([h.data_folder,'/',frame_filename_sorted{frame}]);
                var_names = who('AdcData_scanline*'); % Names on variables containing RF data
                for scan_line = 1:N_sc
                    tools.workbar((scan_line+(frame_idx-1)*N_sc)/(N_frames*N_sc),sprintf('Reading %d frame(s) of FI data from Alpinion.',N_frames),'Reading FI data from Alpinion.')
                    %Read RF data
                    cmd = ['data_raw = ',var_names{scan_line},';'];
                    
                    eval(cmd);
                    
                    data_reordered = h.reorder_chdata(scan_line, channel_data.probe.N_elements, N_sc, N_receive_channels, double(data_raw));
                    data_reordered = data_reordered';
                    
                    MuxTbl = RxMux(scan_line,:);
                    idx = find((MuxTbl > 0) & (MuxTbl <= channel_data.probe.N_elements));
                    idx1 = MuxTbl(idx);
                    
                    RF_Ch_data = RF_Ch_data*0; % Set to zero
                    RF_Ch_data(1:size(data_reordered,1),idx1) = data_reordered(:,idx);

                    clear data_reordered;

                    %Calculate t0 compensation
                    hardcoded_offset = (sqrt(channel_data.sequence(scan_line).source.x^2+channel_data.sequence(scan_line).source.z^2) - channel_data.sequence(scan_line).source.z);
                    channel_data.sequence(scan_line).delay = -(-hardcoded_offset./channel_data.sound_speed + System.Transducer.delayOffsetUsec*10^-6);

                    %build the dataset
                    all_data(:,:,scan_line,frame_idx) = RF_Ch_data;
        
                    clear(var_names{scan_line});
                end
                frame_idx = frame_idx + 1;
            end
            tools.workbar(1)
            if channel_data.pulse.center_frequency > 10*10^6 
                fc = channel_data.pulse.center_frequency;
                fs = channel_data.sampling_frequency;
                F = [fc-8*10^6 fc-7.5*10^6 fs/2-10^6 fs/2-0.5*10^6]
                [all_data] = tools.band_pass(all_data,channel_data.sampling_frequency,F);
            else
                all_data = h.bandpass_filter_data(all_data,channel_data,0);
            end
            channel_data.data = single(all_data);
        end