function R=rotmat(phi)
% ROTMAT Compute a 2D rotation matrix
%
% R=rotmat(phi)

R = [cos(phi), -sin(phi); sin(phi), cos(phi)];
