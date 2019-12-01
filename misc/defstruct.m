function y=defstruct(x,def)
% DEFSTRUCT Fill struct with default values if needed.
%
%  y=defstruct(x,def)

% $Revision: 2952 $  $Date: 2006-09-22 18:45:42 +0200 (fre, 22 sep 2006) $
% Copyright (c) Finn Lindgren 2002

if isempty(x)
  y = def;
else
  y = x;
  fields = fieldnames(def);
  for i=1:length(fields)
    if ~isfield(y,fields{i});
      y = setfield(y,fields{i},getfield(def,fields{i}));
    end
  end
  % To achieve the same field ordering as in def,
  %   and discard unknown things:
  % y = def;
  % fields = fieldnames(def);
  % for i=1:length(fields)
  %   if isfield(x,fields{i});
  %     y = setfield(y,fields{i},getfield(x,fields{i}));
  %   end
  % end
end
