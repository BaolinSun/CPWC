function a=build_mex(filename)

try
    c_flags = '';
    d_flags = '';
    ld_flags = '';
    
    sys_str = lower(computer('arch'));
    
    if ~isempty(strfind(sys_str,'win64'))
        obj_ext = 'obj';
        d_flags = [d_flags '-D_WIN_ '];
    elseif ~isempty(strfind(sys_str,'glnxa64'))
                
        [main_ver, minor_ver, rev_ver] = get_gcc_version();
        
        if main_ver < 4
            error('Must have gcc version 4 or higher');
        end
        
        d_flags=['-D_UNIX_ '];
        c_flags='-I/usr/include/tbb';
        ld_flags='-L/usr/lib64 -ltbb';
     
        obj_ext = 'o';

    elseif ~isempty(strfind(sys_str,'maci64'))
        
        [main_ver, minor_ver, rev_ver] = get_gcc_version();
        
        if main_ver < 4
            error('Must have gcc version 4 or higher');
        end
        
        d_flags=['-D_UNIX_ '];
        c_flags='-I/usr/local/Cellar/tbb/2017_U5/include/tbb/';
        ld_flags='-L/usr/local/Cellar/tbb/2017_U5/lib -ltbb';
     
        obj_ext = 'o';
        
    else
        error('Only MAC, WINDOWS and LINUX is supported.');
    end
    
    mex_str = ['mex ' d_flags ' ' c_flags ' ' ld_flags ' source' filesep filename '.cpp'];
    
    disp(mex_str)
    eval(mex_str);
catch err
    %cd(currdir);
    rethrow(err);
end
delete('*.o');
    

function [main_ver, minor_ver, rev_ver] = get_gcc_version()

main_ver = 0;
minor_ver = 0; 
rev_ver = 0;
[ret str] = system('gcc --version');
if ret == 127
    error('gcc not present on system');
end

remain = str;
while true 
    [tok remain] = strtok(remain,' ');
    if isempty(tok)
        break;
    end
    ver_nums = sscanf(tok,'%d.%d.%d');
    if isempty(ver_nums) || length(ver_nums) ~= 3
        continue
    end
    
    main_ver = ver_nums(1);
    minor_ver = ver_nums(2);
    rev_ver = ver_nums(3);
    break;
end
