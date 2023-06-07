function next_msg=text_progress_bar(c,previous_msg,barsize)
    if nargin<1||isempty(c) c=0; end
    if nargin<2||isempty(previous_msg) previous_msg=''; end
    if nargin<3||isempty(barsize) barsize=20; end
    
    % the message
    n_completed=round(c/100*barsize);
    msg= ['[' repmat('=',[1 n_completed]) repmat('.',[1 barsize-n_completed]) ']' sprintf(' %0.0f',c)];
    
    % the display
    fprintf([previous_msg, msg '%%']);
    next_msg = repmat(sprintf('\b'), 1, length(msg)+1);
end
