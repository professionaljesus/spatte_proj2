function fmsn20path
% FMSN20PATH Set path to fmsn20-subdirectories.
%
% Is called automatically by fmsn20 on the efd network.
% Call manually on other systems.

% $Id: fmsn20path.m 5090 2017-11-05 19:35:05Z johanl $

%find path to this file
[p,~,~]=fileparts(which(mfilename));
%generate a sub path
P = genpath(p);
%drop any .svn folders
I = strfind(P,[filesep '.svn']);
if ~isempty(I)
  J = strfind(P,';');
  for i=length(J)-1:-1:1
    if J(i)<I(end) && I(end)<J(i+1)
      I(end) = [];
      P(J(i)+1:J(i+1)) = [];
    end
  end
end
%add the path
addpath(P);

%set doublebuffer for the figures
set(0,'defaultfiguredoublebuffer','on')
