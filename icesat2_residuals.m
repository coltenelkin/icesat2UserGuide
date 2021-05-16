function [rmsez, residuals, differences] = icesat2_residuals(icesat2, elevations, R2, Abest, offset, ChangeTif, TifYear)
% Function ICESAT2_RESIDUALS calculates residuals of a given ICESat-2 track
% INPUTS: icesat2 = a csv file with icesat 2 elevations created using the
%                       h5 to csv jupyter notebook
%      elevations = the matrix created using geotiffread()
%              R2 = the cell map refernce created as the second output in
%                       geotiffread()
%           Abest = a [2 1] vector that serves as the best spatial offsets in
%                       the easting and northing directions (meters)
%          offset = the offset in elevation values
%                       pulled from the function value using the fminsearch approach
%       ChangeTif = (OPTIONAL) a tif with the same coordinate information
%                   as the elevations tif that represents the annual rate 
%                   of change of a glacier's surface in meters. The sign
%                   convention is positive numbers indicate an accumulating
%                   surface, negative numbers indicate a melting surface,
%                   and areas outside of the glacier outline should be NaN
%                   values
%         TifYear = (OPTIONAL, REQUIRED IF USING A CHANGE TIF) the year of
%                   the input elevations tif used. This script multiplies the number
%                   of years since the input elevations by the average anual rate of
%                   change at each elevation point to get an estimated glacier
%                   surface map in the year of the ICESat-2 transects.
%                   Format is yyyy.
% OUTPUTS:  rmsez = the root mean squared difference in elevation values
%                       before offset
%       residuals = the calculated residuals (vector) between the icesat-2
%                       elevations and their corresponding (offset) DTM 
%                       elevations
%     differences = the vector of point-by-point elevation differences
%                       before elevation offset

% Created 30 October 2020 by Colten Elkin (coltenelkin@u.boisestate.edu)
% last modified 30 January 2021

% most recent update: added two lines to deal with various no data values
% used in DTM storage (-9999; 4.02e38, etc.);

elevations(elevations < -10) = nan; % throw out trash data
elevations(elevations > 10000) = nan; % more trash takeout
if nargin == 7 % check to see if there's a glacier change tif
    yyyy = icesat2(end-38:end-35); % pull the year from the filename
    yyyy = str2num(yyyy); % turn the year string into a number
    numyears = yyyy - TifYear; % get number of years since the input elevation tif
    EleChange = numyears * ChangeTif; % scale the change tif by the number of years since the elevations tif
    elevations = elevations + EleChange;
end

A = Abest; % save as the name that's in the loop below

T = readtable(icesat2); % read in the first icesat-2 atl08 csv file

zmod = T.Elevation(1:end-1); % save the 'model' elevations (icesat-2 elevations)

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
        xpoly = xpoly+A(1); % adjust by the easting offset
        xs{r} = [xpoly(1), xpoly(2), xpoly(3), xpoly(4), xpoly(1)]; % save the corners as a vector in the x-es cell array
        
        ypoly(1) = norths(r) - 50 - footwidth/2*cos((pi/2) - theta(r)); % calculate the 4 corners in the y direction
        ypoly(2) = norths(r) - 50 + footwidth/2*cos((pi/2) - theta(r));
        ypoly(3) = norths(r) + 50 + footwidth/2*cos((pi/2) - theta(r));
        ypoly(4) = norths(r) + 50 - footwidth/2*cos((pi/2) - theta(r));
        ypoly = ypoly+A(2); % adjust by the nothing offset
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
elevation_report = zeros([1, length(xs)]);

for t = 1:length(xs)
    xv = xs{t}; % bounding box x vector
    yv = ys{t}; % bounding box y vector
    
    
    % first trimming
    in = inpolygon(xgrid, ygrid, xv, yv); % get logical array of in values
    pointsinx = xgrid(in); % save x locations
    pointsiny = ygrid(in); % save y locations
    in2 = flip(in); % create a flipped in-grid (need row, column instead of column, row)
    elevationsin = elevations(in2); % save elevations
    elevation_report(t) = nanmean(elevationsin);
end
ztruth = elevation_report(:); % create column vector

differences = zmod-ztruth; % calculate differences
differences(differences > 80) = nan; % toss bad points
differences(differences < -80) = nan; % toss bad points

rmsez = sqrt(nanmean((differences).^2)); % calculate rmsez

residuals = (zmod+offset)-ztruth; % calculate residuals