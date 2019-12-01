function [lan_data] = lanread(lan_filename,thepath)
% LANREAD Read Landsat data file type .lan
%
% Ex: im = lanread('montana.lan')
% size(im) = [m,n,d];
%
% Alt: im = lanread('montana.lan',thepath) where thepath is
%      path to the *.lan-files.  The image must either be located in the
%      current directory, in the one of the directories
%      .../data/ or .../data/protected/lan under the fms150-directory,
%      or be present on the Matlab path.
%
% Available Landsat images are:
%    littlecoriver.lan
%    mississippi.lan
%    montana.lan
%    paris.lan
%    rio.lan
%    tokyo.lan
%
% (From landsatdemo in the image analysis toolbox.)

% $Id: lanread.m 3325 2007-04-06 15:52:34Z finn $

if (nargin<2), thepath = []; end
if isempty(thepath)
  [p,n,e]=fileparts(which('fms150path'));
  thepath= {fullfile('.',filesep),...
            fullfile(p,filesep,'data',filesep),...
            fullfile(p,filesep,'data',filesep,'protected',filesep),...
            fullfile(p,filesep,'data',filesep,'protected',filesep, ...
                     'lan',filesep),...
            ''};
elseif ischar(thepath)
  thepath = {thepath};
end

fid = -1;
if (fid<0)
  for path_idx=1:length(thepath)
    filename = sprintf('%s%s',thepath{path_idx},lan_filename);
    fid = fopen(filename,'r');
    if (fid>=0), break; end
  end
  if (fid<0) % If not found anywhere
    error(sprintf('Could not open file: %s',filename));
  end
end

% find out how big the image is based on file size,
% assuming square image, 7 bands
nbands = 7;
fseek(fid,0,'eof');
file_bytes = ftell(fid);
nlines = floor(sqrt(file_bytes/nbands));
nsamples = nlines;

% skip header
nbytes_header = 128;
fseek(fid,nbytes_header,'bof');

% prepend * to read data into an array that has the same class as the data
A = fread(fid,[nsamples nlines*nbands],'*uint8'); 

fclose(fid);

% put data into a 3D array
A_3dim = reshape(A,nsamples,nbands,nlines);
lan_data = permute(A_3dim,[3 1 2]);
