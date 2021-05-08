function coordData = GetData(line)
% obtain the float data from a string of character
%inputs: 
%             line  -  a string of character
%outputs:
%        coordData  -  float data 

%---------------------------------------------------------
line(1:2) = [];
coordData = str2num(line);