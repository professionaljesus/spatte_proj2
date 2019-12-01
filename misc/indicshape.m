function indic=indicshape(marks,sz)
% INDICSHAPE Compute an indicator image for a deformed template.
%
%  indic=indicshape(marks,sz)
%
%  marks = (p-by-2)-matrix with landmark coordinates.
%  sz = image size in pixels.

% $Id: indicshape.m 3020 2006-10-11 15:43:08Z johanl $

indic = double(polyimage(marks,sz)>.5);

%indic = double(roipoly(ones(sz),marks(:,1),marks(:,2)));
