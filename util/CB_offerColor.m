function [out] = CB_offerColor(varargin)
% [out] = CB_offerColor(cond)
%
% Define color scheme for offer size in CB experiments. 
% Returns OUT, which is a 5x3 RGB array for the 5 conditions.  Optionally
% pass vectors COND to specify which conditions (1:5) to include in OUT.
%
% Daniel Kimmel, April 7, 2013

% based on uoregon.edu below, but converted to print friendly according to
% illustrator and added grey for middle point

%% collect optional input vars

if length(varargin) > 0
    condN = varargin{1};
else
    condN = [];
end

%% define

% DK's old colors from thesis, circa 5/20/13
% out = [...
%     80 23 81
%     181 78 157
%     153 153 153
%     38 179 75
%     32 81 40]/255;

% Gamal's new colors circa 1/11/16
out = [...
             0    0.2745    0.4980
    0.4000    0.6000    1.0000
    0.3216    0.8588    1.0000
    1.0000    0.6235    0.0353
    1.0000    0.2118    0.1373];

if ~isempty(condN)
    out = out(condN,:);
end

return

out = jet(m);
out(1,:) = [66 103 177]/255;
out(2,:) = [49 208 208]/255; %[102 204 204]/255;
out(3,:) = [102 204 51]/255; %[0 222 0]/255; %[0 190 0]/255;
out(4,:) = [255 128 0]/255; %[255 153 51]/255; % [238 145 5]/255; %[255 170 0]/255;
out(5,:) = [255 103 103]/255; % [255 51 51]/255; % [255 65 65]/255;

% color brewer
out(1,:) = [123 50 148]/255;
out(2,:) = [194 165 207]/255;
out(3,:) = [100 100 100]/255; 
out(4,:) = [166 219 160]/255; 
out(5,:) = [0 136 55]/255; 

% hsv
out(1,:) = [0.0312 0 1];
out(2,:) = [0.969 0 1]; %[102 204 204]/255;
out(3,:) = [0.49 0.49 0.49];% gray % [0 150 71]/255; %gray green  % [0 0.812 1];% cyan %[0 222 0]/255; %[0 190 0]/255;
out(4,:) = [0.0312 1 0]; %[255 153 51]/255; % [238 145 5]/255; %[255 170 0]/255;
out(5,:) = [1 0 0.188]; % red-orange % [1 0 0.0938]; % red % [255 51 51]/255; % [255 65 65]/255;

% hsv - colorsafe
out(1,:) = [58 83 164]/255;
out(2,:) = [194 99 167]/255; %[102 204 204]/255;
out(3,:) = [101 103 43]/255; % army green % [0.4 0.4 0.4];% gray % [0 150 71]/255; %gray green  % [0 0.812 1];% cyan %[0 222 0]/255; %[0 190 0]/255;
out(4,:) = [106 189 69]/255; %[255 153 51]/255; % [238 145 5]/255; %[255 170 0]/255;
out(5,:) = [237 29 55]/255; % red-orange % [1 0 0.0938]; % red % [255 51 51]/255; % [255 65 65]/255;



% color brewer - colorsafe - light purple
out(1,:) = [118 42 131]/255;
out(2,:) = [175 141 195]/255;
out(3,:) = [231 212 232]/255;
out(4,:) = [127 191 123]/255;
out(5,:) = [27 120 55]/255;

% color brewer - colorsafe - light yellow
out(1,:) = [69 117 180]/255;
out(2,:) = [145 191 219]/255;
out(3,:) = [211 209 122]/255; % [254 224 144]/255; % color brewer version
out(4,:) = [252 141 89]/255;
out(5,:) = [215 48 39]/255;

% http://geography.uoregon.edu/datagraphics/color_scales.htm
% Green-to-Magenta, 16 Steps
temp = [...
0	0.316	0
0	0.526	0
0	0.737	0
0	0.947	0
0.316	1	0.316
0.526	1	0.526
0.737	1	0.737
1	1	1
1	0.947	1
1	0.737	1
1	0.526	1
1	0.316	1
0.947	0	0.947
0.737	0	0.737
0.526	0	0.526
0.316	0	0.316];

clear out
out = temp(fliplr([1 3 12 13 16]),:);
out(3,:) = [0.5 0.5 0.5];



    
    
    
