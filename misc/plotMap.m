function plotMap(map, data, edge)
% plotMap Plots and colours a map based on a cell array of polygons.
%
% plotMap(map)
% plotMap(map, data)
% plotMap(map, data, edge)
%
% map - a cell array containing cells with polygons as n-by-2 matrices.
% data - a vector of length(data)==length(map), used to determine colours
%        for each polygon. If empty only the polygon outlines are plotted.
% edge - colour of edge bordering each polygon. Should be:
%    colourtriplet: edge=[1 1 0]
%        character: edge='k'
%         no edges: edge='none'

if nargin<2, data=nan; end
if nargin<3, edge=[0.5 0.5 0.5]; end

%check sizes
if length(data)==1, 
  data=repmat(data, [length(map) 1]); 
elseif length(data)~=length(map)
  error('data should be of same length as map (%u)', length(map))
end

%check if edges defined
if all(isnan(data(:))) && strcmp(edge,'none')
  edge='k';
end
  
isHold = ishold();
if ~isHold, plot(nan,nan); end
%figure
hold on
for i=1:length(map)
  if isnan(data(i)) && ~strcmp(edge,'none')
    for j=1:length(map{i})
      line(map{i}{j}(:,1), map{i}{j}(:,2), 'color', edge)
    end
  else
    for j=1:length(map{i})
      patch(map{i}{j}(:,1), map{i}{j}(:,2), ...
        data(i), 'edgecolor', edge)
    end
  end
end
if ~isHold, hold off; end


