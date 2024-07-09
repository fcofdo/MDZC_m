function mValues = IMDC_Read_Matrix(strFileToRead, strFormat, strDelimiter);
%

if ~ischar(strFileToRead)
    error('Input argument must be a character array!')
end

% ********************************************************
% Replace substrings '\' and '\\' with '/' for convinience
% ********************************************************

strFileToRead = strrep(strFileToRead,'\','/');
strFileToRead = strrep(strFileToRead,'//','/');

% *********
% Open file
% *********

fid = fopen(strFileToRead, 'r');
if (fid == -1)
   % The file did not open
   error(['Problems opening file ', strFileToRead])
end

% *********************************
% Read first valid line in the File
% *********************************

NumLines = 0;
while (1)
   line = fgets(fid);
   if (line == -1)
      % End of file.
      fclose(fid);
      error(['No data in file ' strFileToRead])
   end
   firstChar = sscanf(line, '%s', 1);
   if (firstChar ~= '*')
       % This is the first line with values
       % I need to check how many columns I need to read
       vParts = strread(line, '%s', 'delimiter', ' ');
       nCols = size(vParts,1);
       break;
  else
      NumLines = NumLines + 1;
   end
end

fclose(fid);

% Check that the following line is in fact a line with data in columns and
% that the number of columns coincides with the defined number of columns in strColsFormat 

if nargin==2
    strDelimiter = ' ';
end

mTmp = strread(line,'%s','delimiter',strDelimiter);

nCols = size(mTmp,1);
mCols = strread(strFormat,'%s','delimiter',' ');
nExpectedCols = size(mCols,1);
if (nCols ~= nExpectedCols)
    % Expected and read number of columns do not coincide
    error('The expected number of columns does not coincide the read values in the file')
end

% Prepare command to read the file
strTmp = '[';
for I=1:nCols
    strTmp = [strTmp ' mValues(:,' num2str(I) ')'];
end
strTmp = [strTmp ']'];

% Read the file
%[mValues(:,1) mValues(:,2) mValues(:,3) mValues(:,4)] = textread(strFileToRead, strFormat, 'headerlines', NumLines);
eval([strTmp '= textread(strFileToRead, strFormat, ''headerlines'', NumLines, ''delimiter'', strDelimiter);']);

% *********************************
% Check if it really read something
% *********************************

if isempty(mValues)
    error(['Problems reading data in file ' strFileToRead])
end

return;
