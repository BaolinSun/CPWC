function out = ustb_path()
%USTB_PATH Returns the path of the UltraSound ToolBox
%   Usage: out = ustb_path()
    
    [out,name,ext] = fileparts(which('ustb_path'));
end