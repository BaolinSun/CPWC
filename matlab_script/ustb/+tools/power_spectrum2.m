function [f, F] = power_spectrum2(data,fs)
    % function [fx, F] = power_spectrum(data,fs)
        
    % max size of array before we loop instead to save memory. Don't know what is a
    % sensible size here.
    maxsize = 100e6;
    
    N = size(data,1);
    
    % for very large N-D arrays consider looping 
    w = whos('data');
    if w.bytes < maxsize,
        F = mean(abs(fft(data(:,:), N)), 2);
        F = F/max(F);
    else 
        F = zeros(N,1);
        [~,~,J] = size(data); % collapses dimensions beyond 3
        for jj = 1:J,
            F = F + mean(abs(fft(data(:,:,jj))), 2);
        end
        F = F/max(F);
    end
    f = (0:N-1)'/N * fs;
end

