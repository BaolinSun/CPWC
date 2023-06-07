function button = dialog_timeout(delay, varargin)
% questdlg function with timeout property
%
% Based on timeoutDlg by MathWorks Support Team
% https://uk.mathworks.com/matlabcentral/answers/96229-how-can-i-have-a-dialog-box-or-user-prompt-with-a-time-out-period
%
% button = questdlgtimeout(delay,'qstring')
% button = questdlgtimeout(delay,qstring,title)
% button = questdlgtimeout(delay,qstring,title,default)
% button = questdlgtimeout(delay,qstring,title,str1,str2,default)
% button = questdlgtimeout(delay,qstring,title,str1,str2,str3,default)
% button = questdlgtimeout(delay,qstring,title, ..., options)
%
% INPUT ARGUMENTS
% delay       Duration in second during withich the dialog is maintained
%
% var1,var2,...
%             Accepts input arguments for builtin questdlg. See
%             documentation of questdlg
%
% OUTPUT ARGUMENTS
% button       The dialog has three default buttons, Yes, No, and Cancel.
%              If the user presses one of these three buttons, button is
%              set to the name of the button pressed. If the user presses
%              the close button on the dialog without making a choice,
%              button returns as an empty character vector (''). If the
%              user presses the Return key, button returns with a value of
%              'Yes'.
%
%              If you provide default or options, button will be the
%              default value.
%
%
% See also
% questdlg, timer, scr20170308_154424_questdlgtimeout
%
% Written by Kouichi C. Nakamura Ph.D.
% MRC Brain Network Dynamics Unit
% University of Oxford
% 08-Mar-2017 16:06:58

f1 = findall(0, 'Type', 'figures');
t = timer('TimerFcn', {@closeit f1}, 'StartDelay', delay);
start(t);

dlg = @questdlg;

% Call the dialog
button = dlg(varargin{:});

if isempty(button)
  if length(varargin) >= 3
      if  isstruct(varargin{end})
          button = varargin{end}.Default;
      elseif ischar(varargin{end})

            button = varargin{end};

        else
            error('unexpected syntax')
        end
    else % no default provided
        % leave button empty
    end
end

% Delete the timer
if strcmp(t.Running, 'on')
  stop(t);
end
delete(t);

function closeit(src, event, f1)
fprintf(1,'Time out! ');
f2 = findall(0, 'Type', 'figure');
fnew = setdiff(f2, f1);
if ishandle(fnew)
  close(fnew);
end
