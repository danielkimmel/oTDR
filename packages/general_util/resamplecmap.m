% Resample a colormap by interpolation or decimation
function cmap = resamplecmap(cmap, clevels, xi)

t = cmap;
if nargin < 3
    xi = linspace(1,clevels,size(t,1)); 
end
xi([1 end]) = [1 clevels]; % These need to be exact for the interpolation to 
% work and we don't want machine precision messing it up
cmap = [interp1(xi, t(:,1), 1:clevels);...
        interp1(xi, t(:,2), 1:clevels);...
        interp1(xi, t(:,3), 1:clevels)]';
end