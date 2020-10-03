function [out] = commonCondVar(RA,CC)
% [out] = commonCondVar(RA,CC)
% 
% Compute variance associated with the common condition response CC
% explained by each regression axis RA.
%
% ACCEPTS
% RA -- N x M matrix of coefficients to describe M axes in N-dimensional
% space. (neurons x axes). RA may be a cell embedding the NxM matrix. 
% CC -- T x N matrix of data values for each of N dimensions at T time
% points (times x neurons). CC may be a cell embedding the TxN matrix. 
%
% RETURNS
% out -- structure with following fields
%       .V_neu -- variance across time associated with each neuron (N x 1)
%       .V_tot -- total variance summed across neurons (1 x 1)
%       .V_RA -- variance across time explained by each axes (M x 1)
%       .VE_RA -- .V_RA normalized by total variance .V_tot
%
% copyright Daniel Kimmel, 2017 January 15.

if iscell(RA)
    RA = RA{1};
end
if iscell(CC)
    CC = CC{1};
end

% compute variance ACROSS TIME for each neuron
out.V_neu = var(CC,[],1)'; % store as neurons x 1

% total variance ACROSS NEURONS
out.V_tot = sum(out.V_neu);

% project CC response onto each RA (no mean subtraction):
proj = projData(CC,RA,false);

% compute variance ACROSS TIME (not conditions) for each RA:
out.V_RA = squeeze(var(proj,[],1)); % store as RAs x 1

% normalize as percent variance
out.VE_RA = out.V_RA ./ out.V_tot * 100;
