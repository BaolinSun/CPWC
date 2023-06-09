function [filtered_p,h,w] = high_pass(p,Fs,F)
    % function filtered_p = high_pass(p,Fs,[upper_freq_on upper_freq_off])
    
    % filter specification
    A=[0 1];                % band type: 0='stop', 1='pass'
    dev=[1e-3 1e-3];        % ripple/attenuation spec
    [M,Wn,beta,typ]= kaiserord(F,A,dev,Fs);  % window parameters
    b=fir1(M,Wn,typ,kaiser(M+1,beta),'noscale'); % filter design

    [h,w] = freqz(b);
    
    % filtering
    filt_delay=round((length(b)-1)/2);
    filtered_p=filter(b,1,[p; zeros(filt_delay,size(p,2),size(p,3),size(p,4))],[],1);

    % correcting the delay
    filtered_p=filtered_p((filt_delay+1):end,:,:);
    
end

