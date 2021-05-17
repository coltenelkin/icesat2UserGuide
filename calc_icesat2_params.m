function [params, range, SDev] = calc_icesat2_params(icesat2, tif, R2, offsets)
% Function COREGISTER_ICESAT2 coregisters icesat-2 data with a corresponding digital
% terrain model 
% INPUTS: icesat2 = a csv file with icesat 2 elevations created using the
%                       h5 to csv jupyter notebook
%             tif = the matrix created using geotiffread() of the terrain
%                       parameter (slope, aspect, etc)
%              R2 = the cell map refernce created as the second output in
%                       geotiffread()
%         offsets = a [2 1] vector that serves as the spatial offsets in
%                       the x and y directions (meters)
% OUTPUTS: params = the terrain parameter from the tif averaged over the
%                   bounding box of the icesat2 footprint
%           range = the range (max - min) of the parameter in the footprint
%            SDev = the standard deviation of the parameter in the
%            footprint

% Created 30 November 2020 by Colten Elkin (coltenelkin@u.boisestate.edu)
% last modified 30 November 2020

T = readtable(icesat2);

easts = T.Easting(:); % pull out the easting values
norths = T.Northing(:); % pull out the northings
footwidth = 11; % approx. width of icesat2 shot footprint in meters

theta = zeros(size(norths)); % initialize empty matrices
xs = {};
ys = {};
xpoly = nan([1,5]);
ypoly = nan([1,5]);


for r = 1:length(theta)-1
    if strcmp(icesat2(end-44:end-40), 'ATL08') == 1 % ATL08 commands
        theta(r) = abs(atan((norths(r+1) - norths(r))/(easts(r+1) - easts(r)))); % trig to get angle theta along-track
        
        % get the x and y vectors to form the polygon
        xpoly(1) = easts(r) + (footwidth/2) - footwidth/2*cos((pi/2) - theta(r)); % calculate the 4 corners in the x direction
        xpoly(2) = easts(r) + (footwidth/2) + footwidth/2*cos((pi/2) - theta(r));
        xpoly(3) = easts(r) - (footwidth/2) + footwidth/2*cos((pi/2) - theta(r));
        xpoly(4) = easts(r) - (footwidth/2) - footwidth/2*cos((pi/2) - theta(r));
        xpoly = xpoly+offsets(1); % adjust by the easting offset
        xs{r} = [xpoly(1), xpoly(2), xpoly(3), xpoly(4), xpoly(1)]; % save the corners as a vector in the x-es cell array
        
        ypoly(1) = norths(r) - 50 - footwidth/2*cos((pi/2) - theta(r)); % calculate the 4 corners in the y direction
        ypoly(2) = norths(r) - 50 + footwidth/2*cos((pi/2) - theta(r));
        ypoly(3) = norths(r) + 50 + footwidth/2*cos((pi/2) - theta(r));
        ypoly(4) = norths(r) + 50 - footwidth/2*cos((pi/2) - theta(r));
        ypoly = ypoly+offsets(2); % adjust by the nothing offset
        ys{r} = [ypoly(1), ypoly(2), ypoly(3), ypoly(4), ypoly(1)]; % save the corners as a vector in the y-s cell array
        
    elseif strcmp(icesat2(end-44:end-40), 'ATL06') == 1 % ATL06 commands
        theta(r) = abs(atan((norths(r+1) - norths(r))/(easts(r+1) - easts(r)))); % trig to get angle theta along-track
        
        % get the x and y vectors to form the polygon
        xpoly(1) = easts(r) + (footwidth/2) - footwidth/2*cos((pi/2) - theta(r)); % calculate the 4 corners in the x direction
        xpoly(2) = easts(r) + (footwidth/2) + footwidth/2*cos((pi/2) - theta(r));
        xpoly(3) = easts(r) - (footwidth/2) + footwidth/2*cos((pi/2) - theta(r));
        xpoly(4) = easts(r) - (footwidth/2) - footwidth/2*cos((pi/2) - theta(r));
        xpoly = xpoly+A(1); % adjust by the easting offset
        xs{r} = [xpoly(1), xpoly(2), xpoly(3), xpoly(4), xpoly(1)]; % save the corners as a vector in the x-es cell array
        
        ypoly(1) = norths(r) - 20 - footwidth/2*cos((pi/2) - theta(r)); % calculate the 4 corners in the y direction
        ypoly(2) = norths(r) - 20 + footwidth/2*cos((pi/2) - theta(r));
        ypoly(3) = norths(r) + 20 + footwidth/2*cos((pi/2) - theta(r));
        ypoly(4) = norths(r) + 20 - footwidth/2*cos((pi/2) - theta(r));
        ypoly = ypoly+A(2); % adjust by the nothing offset
        ys{r} = [ypoly(1), ypoly(2), ypoly(3), ypoly(4), ypoly(1)]; % save the corners as a vector in the y-s cell array
        
    end
end


x = R2.XWorldLimits(1)+0.5*R2.CellExtentInWorldX:R2.CellExtentInWorldX:R2.XWorldLimits(end)-0.5*R2.CellExtentInWorldX; % get a vector of x coords
y = R2.YWorldLimits(1)+0.5*R2.CellExtentInWorldY:R2.CellExtentInWorldY:R2.YWorldLimits(end)-0.5*R2.CellExtentInWorldY; % get a vector of y coords


[xgrid, ygrid] = meshgrid(x, y); % create grids of each of the x and y coords
parameter_report = zeros([1, length(xs)]);

for t = 1:length(xs)
    xv = xs{t}; % bounding box x vector
    yv = ys{t}; % bounding box y vector
    
    
    % first trimming
    in = inpolygon(xgrid, ygrid, xv, yv); % get logical array of in values
%     pointsinx = xgrid(in); % save x locations
%     pointsiny = ygrid(in); % save y locations
    in2 = flip(in); % create a flipped in-grid (need row, column instead of column, row)
    paramsin = tif(in2); % saveparameters
    parameter_report(t) = nanmean(paramsin);
end
params = parameter_report(:);
