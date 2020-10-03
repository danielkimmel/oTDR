function mm = minmax(M)
%
% mm = minmax(M)
%
% This replaces the MATLAB native MINMAX() function, which was moved to the
% DeepLearning toolbox. 
%
% Accepts the R x Q matrix M and returns an R x 2 matrix MM containing the
% minimum and maximum value of each row in M.
%
% Daniel Kimmel, 2019 Aug 02

mm = [min(M,[],2),max(M,[],2)];