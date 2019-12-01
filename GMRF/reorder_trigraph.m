function [FV1,S1,ind0,ind1]=reorder_trigraph(FV0,S0,order)
% REORDER_TRIGRAPH Reorder nodes in a triangular graph for nice Cholesky
%
% [FV1,S1,ind0,ind1]=reorder_trigraph(FV0,S0)
%
%  ind0 = amd("typical precision matrix")
%  ind1 = the inverse ordering
%
%  S1 == S0(ind0,:)
%  S0 == S1(ind1,:)
%  FV1 == ind1(FV0)
%  FV0 == ind0(FV1)

% $Id: reorder_trigraph.m 4836 2014-12-10 11:09:32Z johanl $

if (nargin<3), order = []; end
if isempty(order), order = 2; end

n = size(S0,1);

VV = fv2vv(FV0,n).*1;
if (order==1)
  VV = speye(n)+VV;
else
  VV = speye(n)+VV*VV;
end
%test if amd is a built in function (it is as of 2006b)
if ( exist('amd','builtin')==0 )
  %nope, use old slow symamd
  ind0 = symamd(VV);
else
  ind0 = amd(VV);
end
ind0 = ind0(:);
ind1 = (1:n)';
ind1(ind0) = ind1;

S1 = S0(ind0,:);
FV1 = ind1(FV0);

%% Sanity check:
% [S1(ind1,:)-S0]
% [S1(2,:)-S0(ind0(2),:),S0(2,:)-S1(ind1(2),:)]


function VV=fv2vv(FV,n)
% FV2VV Calculate sparse VV graph matrix from FV graph.

if (nargin<2)
  n = max(FV(:));
end

VV=sparse([FV(:,1);FV(:,1);FV(:,2)],...
          [FV(:,2);FV(:,3);FV(:,3)],...
          ones(size(FV,1)*3,1),n,n);

VV=((VV+VV')>0); % Handle open triangulations by symmetrisation.
