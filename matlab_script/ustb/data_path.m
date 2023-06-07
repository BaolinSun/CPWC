function path_out = data_path()
%DATA_PATH Returns the path (root) of data for USTB.
% Returns the value of environment variable 'DATA_PATH' if it exists,
% otherwise it returns path to `data` folder in USTB
%
%   Usage: out = ustb_data_path()
%
    path_out = getenv('DATA_PATH'); 
    
    if isempty(path_out),
        path_out = [ustb_path(), filesep, 'data'];
    end
    
    if ~exist(path_out, 'dir'),
        % we must duplicate filesep to get the message through warning
        ind_sep=find(path_out==filesep);
        for n=1:numel(ind_sep)
            path_out = [path_out(1:ind_sep(n)) path_out(ind_sep(n):end)];
            if(n<numel(ind_sep)) ind_sep(n+1:end)=ind_sep(n+1:end)+1; end
        end
        
        warning('USTB:DataPathNotExist', ['"', path_out, '"', ' does not exist!']);
    end
end
