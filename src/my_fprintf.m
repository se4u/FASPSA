function my_fprintf(channel, string, varargin)
global disable_print;
if isempty(disable_print) || disable_print ~= channel
    fprintf(channel, string, varargin{:}); 
end