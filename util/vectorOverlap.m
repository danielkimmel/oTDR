function [out] = vectorOverlap(sRA,varargin)
%
% [dRA_lowD_prop,dRA_lowD_prop_bySRA,dRA_lowD_mag,dRA_lowD_mag_bySRA] = ...
%    vectorOverlap(sRA,varargin);
%
% Computes the magnitude and proportion of the high-D vectors dRA captured
% by the low-D space spanned by the vectors sRA.
% 
% ACCEPTS
% sRA -- N x J matrix specifying J N-dimensional vectors that define a
%        J-dimensional subspace into which dRA_orig is projected.
%
% THE FOLLOWING vars are passed as part of varargin, but are NOT optional
% dRA_orig -- N x K matrix specifying K N-dimensional vectors. Vectors
%        should full magnitude (i.e., NOT normalized). 
%
% RETURNS
% dRA_lowD_prop -- ...
%
% Daniel Kimmel, 2017 June 10

% collect varargin
dRA_orig = varargin{1};

% convert sRA from cell to vector if necessary (passed in as cell when
% vectorOverlap is called by sampleRandSubspaces)
if iscell(sRA)
    sRA = sRA{1};
end

% compute magnitude of dRA_orig and normalize dRA_orig vectors 
[dRA_mag,dRA] = normVects(dRA_orig);

% project high-D dRA into low-D space. This provides a proportion of dRA
% caputured by the sRAs because dRA defines unit vectors, thus projection
% of 1 indicates 100% capture.
dRA_lowD_prop_bySRA = abs(dRA' * sRA);

% % find the proportion of the dRA captured by each of the 3 sRAs. We
% % normalize this magnitude by the magnitude of the high-D dRA, even
% % though this is redudant since the high-D sRA has magnitude 1.
% dRA_lowD_prop_bySRA = bsxfun(@rdivide,abs(dRA_lowD),dRA_mag_unit);

% compute proportion of dRA captured by lowD space. Since the sRAs are
% orthogonal, we find the proportion of dRA captured by taking the
% magnitude of the projected vector (i.e., sqrt of sum of squared
% components, summed across sRAs). No need to normalize for the proportion,
% since dRA defined unit vectors. 
dRA_lowD_prop = sqrt(sum(dRA_lowD_prop_bySRA .^ 2,2));

% % check
% if any(bsxfun(@minus,sqrt(sum(dRA.^2,2)),sqrt(sum(dRA_lowD.^2,2))) < 0)
%     error('projection onto dRA cannot be more than projection magnitude into 3D space')
% end

% Backsolve for the actual magnitude of the projected dRA by
% scaling the high-D magnitude by the proportion captured
dRA_lowD_mag = dRA_mag .* dRA_lowD_prop;
% Do the same for the magnitude of each component
dRA_lowD_mag_bySRA = bsxfun(@times, dRA_lowD_prop_bySRA, dRA_mag);

% compute difference between magnitudes of original vector and projected
% vector:
dRA_lowD_magDiff = dRA_mag - dRA_lowD_mag;
dRA_lowD_magDiff_bySRA = bsxfun(@minus,dRA_mag,dRA_lowD_mag_bySRA);

% double-check no negative values
if any(dRA_lowD_magDiff(:) < 0) || any(dRA_lowD_magDiff_bySRA(:) < 0)
    error('Difference in vector magnitudes < 0, which is not expected.')
end

% store in output var:
out.dRA_lowD_prop = dRA_lowD_prop;
out.dRA_lowD_prop_bySRA = dRA_lowD_prop_bySRA;
out.dRA_lowD_mag = dRA_lowD_mag;
out.dRA_lowD_mag_bySRA = dRA_lowD_mag_bySRA;
out.dRA_lowD_magDiff  = dRA_lowD_magDiff;
out.dRA_lowD_magDiff_bySRA = dRA_lowD_magDiff_bySRA;
