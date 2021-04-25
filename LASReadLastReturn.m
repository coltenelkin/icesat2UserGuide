function outfile = LASReadLastReturn(infilename, outfilename)
% LASREAD reads in files in LAS 1.1 format and outputs comma delimited text files
% note: In order to make the program work for LAS 1.2 format go to 
% line 40 and change it to "if VersionMajor ~= 1"
%
% INPUT
% infilename:   input file name in LAS 1.1 format 
%               (for example, 'myinfile.las') 
% outfilename:  output file name in text format 
%               (for example, 'myoutfile.txt')
% nFields:      default value of 1 outputs X, Y and Z coordinates of the 
%               point - [X Y Z]. 
%               A value of 2 gives Intensity as an additional attribute - [X Y Z I].
%               A value of 3 gives the Return number and the Number of returns 
%               in addition to the above - [X Y Z I R N].                
%           
% OUTPUT
% outfile:      the output matrix
% 
% EXAMPLE
% A = LASRead ('infile.las', 'outfile.txt', 3)
%
% Cici Alexander
% September 2008 (updated 26.09.2008)
% 
% Modifications 4 October 2020 by Colten Elkin:
% changed bottom few lines of script to only save last return.
% The object is to save only terrain (not canopy) elevations
% in order to create a digital terrain model using the downsampling.mlx
% live script. Contact coltenelkin@u.boisestate.edu for info. 

% Open the file
fid =fopen(infilename);

% Check whether the file is valid
if fid == -1
    error('Error opening file')
end

% Check whether the LAS format is 1.1
fseek(fid, 24, 'bof');
VersionMajor = fread(fid,1,'uchar');
VersionMinor = fread(fid,1,'uchar');
if VersionMajor ~= 1 %|| VersionMinor ~= 1
    error('LAS format is not 1.1')
end

% Read in the offset to point data
fseek(fid, 96, 'bof');
OffsetToPointData = fread(fid,1,'uint32');

% Read in the scale factors and offsets required to calculate the coordinates
fseek(fid, 131, 'bof');
XScaleFactor = fread(fid,1,'double');
YScaleFactor = fread(fid,1,'double');
ZScaleFactor = fread(fid,1,'double');
XOffset = fread(fid,1,'double');
YOffset = fread(fid,1,'double');
ZOffset = fread(fid,1,'double');

% The number of bytes from the beginning of the file to the first point record
% data field is used to access the attributes of the point data
%
c = OffsetToPointData;

% Read in the X coordinates of the points
%
% Reads in the X coordinates of the points making use of the
% XScaleFactor and XOffset values in the header.
fseek(fid, c, 'bof');
X1=fread(fid,inf,'int32',24);
X=X1*XScaleFactor+XOffset;
X = X(:); % turn vectors into column vectors for matrix generation below

% Read in the Y coordinates of the points
fseek(fid, c+4, 'bof');
Y1=fread(fid,inf,'int32',24);
Y=Y1*YScaleFactor+YOffset;
Y = Y(:); % column vector

% Read in the Z coordinates of the points
fseek(fid, c+8, 'bof');
Z1=fread(fid,inf,'int32',24);
Z=Z1*ZScaleFactor+ZOffset;
Z = Z(:); % column vector


% Read in the Intensity values of the points
fseek(fid, c+12, 'bof');
Int=fread(fid,inf,'uint16',26);
Int = Int(:); % column vector


% Read in the Return Number of the points. The first return will have a
% return number of one, the second, two, etc.
fseek(fid, c+14, 'bof');
Rnum=fread(fid,inf,'bit3',221);
Rnum = Rnum(:); % column vector
Rnum = round(Rnum); % turn into integers

% Read in the Number of Returns for a given pulse.
fseek(fid, c+14, 'bof');
fread(fid,1,'bit3');
Num=fread(fid,inf,'bit3',221);
Num = Num(:); % column vector
Num = round(Num); % turn into integers

counter = 1:length(X); % create indexing column for cutdown purposes
counter = counter(:); % column vector

matrix  = [counter, X, Y, Z, Int, Rnum, Num]; % create matrix from values

% use logical indexing to eliminate rows where Rnum does not equal num
% (i.e. where the return is not the last return)

cutdown = matrix(matrix(:,6) == matrix(:,7)); % rows with wanted data

matrix = matrix(cutdown, [2 3 4]); % shrink matrix down to just X, Y, and Z of last return

X = matrix(:,1); % create new X column vector with only wanted data
Y = matrix(:,2); % new Y column vector
Z = matrix(:,3); % new Z column vector

% Write out the file with X, Y and Z coordinates, representing only last
% returns

outfileheader = ['X' 'Y' 'Z'];
outfile = [X Y Z];



dlmwrite(outfilename,outfileheader);
dlmwrite(outfilename,outfile, '-append','precision','%.2f','newline','pc');