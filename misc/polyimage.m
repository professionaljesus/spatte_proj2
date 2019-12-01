% POLYIMAGE Computes an aliased indicator image for a polygon.
%
%  indic = polyimage(marks,sz)
%
%  marks = (p-by-2)-matrix with landmark coordinates.
%  sz = image size in pixels.
%
%  Polyimage can behave slightly strange for some self intersecting polygons.
%
% SEE ALSO: indicshape, indicshape_old

% Implemented in a MATLAB mex file using:
%  polyimage.c
%  garbage.h
%  log.h
%  preprocess.h
%  rasterize.h
%  scanline.h
%  types.h
%
% $Id: polyimage.m 4586 2012-10-08 16:18:33Z johanl $
% Copyright 2006, Gustav Lindström (matlab interface Johan Lindström)

%If platform dependent mex file does not exist, revert to old
%function behaviour using roipoly.
indic = double(roipoly(ones(sz),marks(:,1),marks(:,2)));

% END
