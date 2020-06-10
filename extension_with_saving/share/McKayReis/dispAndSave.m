function dispAndSave( str, fid )
%DISPANDSAVE prints a string to screen and writes it in file.
%   In both cases (screen and file) a new line is entered at the end of
%   string.
%
%   str (input) is the string
%   fid (input) is the file id created by fopen
%   
%   usage example:
%   fid = fopen([filename],'w');
%   dispAndSave('hello world', fid);
%   fclose(fid);

disp(str);
fprintf(fid,'%s\n',str);



end

