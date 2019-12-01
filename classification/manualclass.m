function cl=manualclass(x,K,method)
% MANUALCLASS Mark image pixels as belonging to different classes.
%
% cl=manualclass(x,K)
% cl=manualclass(x,K,method)
%
%  x: Image.
%  K: The number of classes
%  method: Which method to use:
%          'pixels': mark indiviual pixels (default).
%          'rectangles': mark rectangles.
%  cl: cl(i,j) = k for pixels marked as class k
%      cl(i,j) = 0 elsewhere
%
% Example:
%  x = lanread('rio.lan');
%  cl0 = manualclass(landsatimage(x),3);

% Copyright (c) 2002 Finn Lindgren
% $Revision: 1562 $  $Date: 2002-04-18 02:48:05 +0200 (tor, 18 apr 2002) $

% Parse input parameters
if nargin<3
  method = [];
end
if isempty(method)
  method = 'pixels';
end

m = size(x,1);
n = size(x,2);
imagesc(x);

switch (method)
 case 'pixels',
  cl = zeros(m,n);
  for k=1:K
    fprintf(['Click on pixels belonging to class %i. ',...
	     'Hit RETURN when you are done.'],k)
    [x_new,y_new]=ginput;
    cl_idx = max(1,round([min(m,y_new),min(n,x_new)]));
    fprintf(['%s%i pixels marked for class %i.',...
	     '                                            ','\n'],...
	    char(13),size(cl_idx,1),k)
    cl = cl+(sparse(cl_idx(:,1),cl_idx(:,2),...
		    ones(size(cl_idx,1),1),m,n)>0)*k;
  end
 case 'rectangles',
  cl = zeros(m,n);
  for k=1:K
    fprintf(['Click the two corners of a rectangle with ',...
	     'pixels belonging to class %i.'],k)
    [x_new,y_new]=ginput(2);
    cl_idx = sort(max(1,round([min(m,y_new),min(n,x_new)])),1);
    fprintf(['%s%i pixels marked for class %i.',...
	     '                                            ','\n'],...
	    char(13),prod((diff(cl_idx,1)+1)),k)
    cl(cl_idx(1,1):cl_idx(2,1),cl_idx(1,2):cl_idx(2,2)) = k;
  end
 otherwise,
  error(['Unknown method "', method, '".']);
end