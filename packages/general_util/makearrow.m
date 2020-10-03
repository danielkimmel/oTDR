function arrow = makearrow(pars)
% makearrow makes the vertices of a 2d arrow
%
% pars has fields:
% dist: the distance between the base of the arrow and the origin
% dir: the direction of the arrow
% length: the lenght of the arrow
% width: the width of the arrow
% lengthratio: the fraction of length taken up by the base
% widthratio: the fraction of width taken up by the base
% 
% arrow = makearrow(pars)

% The total length
L = pars.length;

% The total width
W = pars.width;

% The length of the asta
x1 = L * pars.lengthratio;

% The halfwidth of the asta
y1 = W * pars.widthratio / 2;

% The lenght of the punta
x2 = L - x1;

% The width of the ali della punta
y2 = W/2 - y1;


% The unit vector in direction of the arrow
u = [cos(pars.dir*pi/180); sin(pars.dir*pi/180)];

% The unit vector indicating the orthogonal direction to the right
ur = [u(2); -u(1)];


% Make all the points in the arrow
arrow = zeros(2,7);
arrow(:,1) = u*pars.dist + ur*y1;
arrow(:,2) = arrow(:,1) + u*x1;
arrow(:,3) = arrow(:,2) + ur*y2;
arrow(:,4) = u*pars.dist + u*L;
arrow(:,5) = arrow(:,3) - ur*W;
arrow(:,6) = arrow(:,2) - ur*y1*2;
arrow(:,7) = arrow(:,1) - ur*y1*2;


return

addpath('u:/code/toolboxes/mtools');

pars.dist = -1;
pars.dir = 270;
pars.length = 1;
pars.width = 0.7;
pars.lengthratio = 0.4;
pars.widthratio = 0.5;

arrow = makearrow(pars);
figure; plot(arrow(1,:),arrow(2,:));
set(gca,'dataaspectratio',[1 1 1]);
axis off
% cd('C:\Users\Valerio\talks\2004-VSS\figures');
% figname = 'arrow';
% print(gcf,'-depsc2',figname);




