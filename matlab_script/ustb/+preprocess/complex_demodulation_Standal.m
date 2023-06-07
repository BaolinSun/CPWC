classdef complex_demodulation_Standal < preprocess
    %COMPLEX_DEMODULATION   Matlab implementation of complex demodulation
    %
    %   authors: Øyvind K.-V. Standal
    %
    %   $Last updated: 2017/09/10$
    
    %% constructor
    methods (Access = public)
        function h=complex_demodulation_Standal()
            h.name='Complex demodulation MATLAB';
            h.reference='www.ustb.no';
            h.implemented_by={'Øyvind K.-V. Standal'};
            h.version='v1.0.5';
        end
    end
    
    properties (Access = public)
        plot_on                     % plot intermediate graphs
        modulation_frequency        % modulation frequency [Hz]
        downsample_frequency        % sampling frequency after downsampling [Hz]
    end
    
    methods
        function output=go(h)
            
            % check if we can skip calculation
            if h.check_hash()
                output= h.output;
                return;
            end
            
            % modulation frequency
            if isempty(h.modulation_frequency)||(h.modulation_frequency<eps)
                warning('The modulation frequency is not specified. The estimated central frequency will be used.');
                
                % power spectrum
                [fx pw] = tools.power_spectrum(h.input.data,h.sampling_frequency);
                assert(sum(pw)>0,'Dataset is zero');
                
                % computing central frequency and bandwidth
                disp('Estimating power spectrum');
                fpw=filter(ones(1,26)./26,1,pw);fpw=[fpw(13:end); zeros(12,1)];
                [dc ic]=max(fpw.*(fx>0).'); fc=fx(ic);
                bw_up=min(fx((fx>fc)&(fpw<dc/2).')); % -6dB upper limit
                bw_do=max(fx((fx<fc)&(fpw<dc/2).')); % -6dB down limit
                fc=(bw_up+bw_do)/2;                  % center frequency
                bw=2*(bw_up-fc);
                
                % set modulation ferquency
                h.modulation_frequency=fc;
            end
            
            doresample = false; % need to interpolate (not just downsample by integer rate)
            if isempty(h.downsample_frequency),
                Ndown = 1; % no downsampling
            else
                Ndown = round(h.sampling_frequency/h.downsample_frequency);
                %if abs(Ndown - round(Ndown)) < 0.01,
                %    Ndown = round(Ndown);
                %else
                %    doresample = true;
                %    warning('Downsampling by a non-integer rate is highly inefficient.');
                %end
            end
            
            siz = size(h.input.data);
            N = prod(siz(3:end)); % collapse dimensions beyond 3
            ind_new = 1:Ndown:siz(1);
            L = length(ind_new); % length of downsampled array
            
            % preallocate final array
            iq = zeros([L,siz(2:end)], 'like', 1i);
            
            % demodulation mixing vector
            mix = exp(-2i*pi*h.modulation_frequency/h.sampling_frequency*(0:siz(1)-1)');
            
            % otherwise not really need for low-pass filtering
            dofilter = (Ndown > 1);
            if dofilter, % make low-pass filter
                % calculate window parameters
                [M, Wn, beta, typ] = kaiserord([0.9,1]*h.modulation_frequency, [1,0], [1e-2,1e-2], h.sampling_frequency);
                b = fir1(M, Wn, typ, kaiser(M+1,beta), 'noscale'); % filter design
                assert(size(h.input.data,1)>3*(length(b)-1),'Too many low-pass filter coefficients! Try relaxing the filter design, or increasing the samples in the dataset.');
            end
            
            % preallocate temp array for Hilbert
            % H = zeros(siz(1:2), 'like', 1i);
            for nn = 1:N, % loop over higher dimensions
                H = bsxfun(@times, hilbert(h.input.data(:,:,nn)), mix);
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
            
            % downsampling
            dt=Ndown/h.sampling_frequency;
            t=(h.input.time(1):dt:h.input.time(end));
            
            % copy results
            h.output=uff.channel_data(h.input);
            h.output.modulation_frequency=h.modulation_frequency;
            h.output.initial_time = t(1);
            h.output.sampling_frequency = 1/(t(2)-t(1));
            h.output.data=iq;
            
            % pass reference
            output = h.output;
            
            % update hash
            h.save_hash();            
        end
    end
    
        %% set methods
    methods
%         function h=set.input(h,in_data)
%             assert(isa(in_data,'uff.channel_data'), 'The input channel_data is not a CHANNEL_DATA class. Check HELP CHANNEL_DATA.');
%             assert(in_data.modulation_frequency==0,sprintf('The input channel_data is already demodulated with %0.2 MHz',in_data.modulation_frequency/1e6));
%             assert(~isempty(in_data.data),'The input channel_data is empty');
%             assert(any(in_data.data(:)>0),'The input channel_data is zero');
%             
%             timediff=in_data.time(2)-in_data.time(1);
%             assert(timediff>0,'The time interval of the time vector in the channel_data is zero');
%             
%             h.sampling_frequency=1/timediff;
%             h.channel_data=in_data;
%         end
        function h=set.plot_on(h,in_plot_on)
            assert(isa(in_plot_on,'logical'), 'The input plot_on is not a LOGICAL class (true/false). Check HELP LOGICAL.');
            h.plot_on=in_plot_on;
        end
        function h=set.modulation_frequency(h,in_modulation_frequency)
            assert(numel(in_modulation_frequency)==1, 'The modulation_frequency must be a escalar');
            h.modulation_frequency=in_modulation_frequency;
        end
        function h=set.downsample_frequency(h,in_downsample_frequency)
            assert(numel(in_downsample_frequency)==1, 'The downsample_frequency must be a escalar');
            h.downsample_frequency=in_downsample_frequency;
        end
    end

end
