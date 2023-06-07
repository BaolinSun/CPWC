classdef modified_autocorrelation_displacement_estimation < postprocess
    % MODIFIED AUTOCORRELATION DISPLACEMENT ESTIMATION   
    % 
    % Process to estimate displacement. This was originally introduced
    % estimate blood velocity, we have however modified it to estimate
    % displacement in stead. This implementation also estimates the center
    % frequency as suggested by 
    %
    % Excerpt from my master thesis http://www.olemarius.net/Thesis/ :
    %
    % In 1985 Barber et al. published A New Time Domain Technique for 
    % Velocity Measurements Using Doppler Ultrasound. This technique was 
    % able to estimate the Doppler frequency using the autocorrelation function. 
    % The estimation is done in the time domain, allowing a very fast implementation. 
    % This technique was the basis of the first real-time blood flow imaging system 
    % demonstrated by Kasai et al. (1985). The thoroughly mathematical analysis of 
    % the technique was done by Angelsen and Kristoffersen (1983).
    %
    % The center frequency is attenuated proportional to the depth in ultrasound imaging.
    % In 1995 Loupas et al. suggested to include an estimate of the center frequency, 
    % and called this the 2D autocorrelator. We will call it the modified autocorrelation method. 
    %
    % Credits to Thomas B?rstad for providing the original implementation
    % in his masters thesis:
    % B?rstad, T. K. (2010). Comparison of three ultrasound velocity estimators
    % for strain imaging of the brain. NTNU.
    %
    %   Authors: Ole Marius Hoel Rindal <olemarius@olemarius.net>
    %   $Last updated: 2017/08/15$
    
    %% constructor
    methods (Access = public)
        function h=modified_autocorrelation_displacement_estimation()
            h.name='Autocorrelation Displacement Estimation';
            h.reference=['Loupas, T., Powers, J. T., & Gill, R. W. (1995). An Axial Velocity Estimator for Ultrasound Blood Flow Imaging, Based on a Full Evaluation of the Doppler Equation by Means of a Two-Dimensional Autocorrelation Approach. IEEE Transactions on Ultrasonics, Ferroelectrics and Frequency Control, 42(4).'...
                         'Barber, W. D., Eberhard, J. W., & Karr, S. G. (1985). A new time domain technique for velocity measurements using Doppler ultrasound. IEEE Transaction on Biomedical Engineering, 32(3)'...
                         'Kasai, C., Namekawa, K., Koyano, A., & Omoto, R. (1985). Real-Time Two-Dimensional Blood Flow Imaging Using an Autocorrelation Technique. IEEE Transactions on Sonics and Ultrasonics, 32(3).'...
                         'Angelsen, B.A.J., & Kristoffersen, K. (1983). Discrete time estimation of the mean Doppler frequency in ultrasonic blood velocity measurements. IEEE Transaction on Biomedical Engineering, (4).']; 
            h.implemented_by='Ole Marius Hoel Rindal <olemarius@olemarius.net>, Thomas B?rstad <thomas.borstad@gmail.com>';
            h.version='v1.0.0';
        end
    end
    
    %% Additional properties
    properties
        z_gate = 4
        x_gate = 2
        packet_size = 6
        estimated_center_frequency
        channel_data
    end
    
    methods
        function output=go(h)
            % check if we can skip calculation
            if h.check_hash()
                output = h.output; 
                return;
            end 
            
            [N_pixels Nrx Ntx N_frames]=size(h.input.data);
            
            assert(N_frames>h.packet_size,'The number of frames needs to be higher than the packet size');
            assert(Nrx==1,'The pulsed doppler speckle traking can only be used between frames');
            assert(Ntx==1,'The pulsed doppler speckle traking can only be used between frames');
            assert(mod(h.z_gate,2)==0,'Please use an even number for the z_gate');
            assert(mod(h.x_gate,2)==0,'Please use an even number for the x_gate');
            
            % declare output structure
            output=uff.beamformed_data(h.input); % ToDo: instead we should copy everything but the data
            
            % save scan
            output.scan = h.input.scan;
            
            % calculate sampling frequency in image
            h.input.calculate_sampling_frequency(h.channel_data.sound_speed);
            
            % get images in matrix format
            images = h.input.get_image('none-complex');
            
            % create a buffer for the output
            displacement_data = zeros(size(h.input.data,1),size(h.input.data,2),...
                        size(h.input.data,3),size(h.input.data,4)-h.packet_size+1);
            temp_fc_hat = zeros(size(images));
            
            tools.workbar()
            for i = h.packet_size:h.input.N_frames
                 tools.workbar(i/h.input.N_frames,'Calculating modified autocorr displacement estimation','Modified autocorr displacement estimation');        
                
                % buffer
                temp_disp = zeros(size(images(:,:,1,1,1)));
                
                %Calculate displacement
                [temp_disp(h.z_gate/2:end-h.z_gate/2-1,h.x_gate/2:end-h.x_gate/2),dummy,temp_fc_hat(h.z_gate/2:end-h.z_gate/2-1,h.x_gate/2:end-h.x_gate/2,i-h.packet_size+1)] = modified_pulsed_doppler_displacement_estimation(h,images(:,:,i-h.packet_size+1:i));
                displacement_data(:,1,1,i-h.packet_size+1) = temp_disp(:);
            end
            tools.workbar(1)
            h.estimated_center_frequency = temp_fc_hat;
            output.data = displacement_data;
        end
    end
    
    methods (Access = private)
        function [d_2,f_hat,fc_hat,C] = modified_pulsed_doppler_displacement_estimation(h,X)
            
            c   = h.channel_data.sound_speed;           %Speed of sound
            fs  = h.input.sampling_frequency; %Sampling frequency
            U = h.z_gate;
            V = h.x_gate;
            O = h.packet_size;
            
            X_conj = conj(X); %Complex conjugate of X
            
            %Autocorrelation lag
            R_0_1 = sum(X(1:end-1,:,1:end-1).*X_conj(1:end-1,:,2:end),3);
            R_0_1 = conv2(R_0_1, ones(U,1), 'valid');
            R_0_1 = conv2(R_0_1, ones(1,V), 'valid');
            
            %Autocorrelation lag
            R_1_0 = sum(X(1:end-1,:,:).*X_conj(2:end,:,:),3);
            R_1_0 = conv2(R_1_0, ones(U,1), 'valid');
            R_1_0 = conv2(R_1_0, ones(1,V), 'valid');
            
            %Estimated doppler frequency
            f_hat = (angle((R_0_1))/(2*pi));
            
            %Estimated centeral frequency
            fc_hat = abs((angle(R_1_0))/(2*pi*1/fs));
            
            %Correlation coefficient estimation qualityindicator.
            C = sum(X(1:end-1,:,:).*X_conj(1:end-1,:,:),3);
            C = conv2(C, ones(U,1), 'valid');
            C = conv2(C, ones(1,V), 'valid');
            C = (O/(O-1))*abs(R_0_1)./C;
            
            %Modified autocorrelation method
            d_2 = c*f_hat./(2*fc_hat);
        end
        
        
    end
end