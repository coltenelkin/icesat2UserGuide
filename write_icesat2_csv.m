function write_icesat2_csv(inputdir, outputdir, shapefile)
% function WRITE_ICESAT2_CSV reads in an h5 file and outputs csv files of
% individual beams and their pertinent data
% INPUTS: inputdir = directory pointing to .h5 files (ends with '/') (string)
%        outputdir = directory where you want to save csv files (also ends
%                      with '/') (string)
%        shapefile = directory and name of shapefile of the region of
%                      interest that serves for clipping the icesat-2 data 
%                      down (string)
%                    Note: the shapefile should be in a decimal-degree
%                       coordinate system

% created 21 December 2020 by Colten Elkin (coltenelkin@u.boisestate.edu)
% requires matlab function deg2utm.m (available here:
% https://www.mathworks.com/matlabcentral/fileexchange/10915-deg2utm )
% also requires that path to deg2utm is added in the add path line below:
addpath('/Users/glaciologygroup/Desktop/elkin/ms_code/matlab_scripts'); % add path for calling deg2utm later

% last edited: 19 April 2021
% most recent update: added the geoid elevation to get zRefined

if(~isfolder(outputdir)) % create the output directory if it doesnt already exist
   mkdir(outputdir) 
end

watershed = shaperead(shapefile);
h5files = dir([inputdir,'*.h5']); % pull out the h5 files
beams = {'gt1r', 'gt1l', 'gt2r', 'gt2l', 'gt3r', 'gt3l'}; % list of icesat2 beams for inner loop
for filecounter = 1:length(h5files) % loop through icesat2 files
    % check to see whether it's ATL08 or ATL06
    if strcmp(h5files(filecounter).name(1:5), 'ATL08') == 1 % ATL08 commands
        h5test = [h5files(filecounter).folder,'/',h5files(filecounter).name]; % get string pointing to n-th h5 file
        
        for beamcount = 1:length(beams)
            beam = beams{beamcount}; % set beam
            % read in data
            terrain = h5read(h5test, ['/',beam,'/land_segments/terrain/h_te_mean']); % read terrain means
            lat = h5read(h5test,  ['/',beam,'/land_segments/latitude']); % read lats
            lon = h5read(h5test,  ['/',beam,'/land_segments/longitude']); % read lons
            bright = h5read(h5test,  ['/',beam,'/land_segments/brightness_flag']); % read in brightness of shot
            std = h5read(h5test,  ['/',beam,'/land_segments/terrain/h_te_std']); % standard deviation
            can = h5read(h5test,  ['/',beam,'/land_segments/canopy/h_canopy']); % canopy elevation
            
            % crop data to area of interest
            % first, crude cropping
            lonlims = watershed.BoundingBox(:, 1); % save upper and lower longitudes of the watershed
            latlims = watershed.BoundingBox(:, 2); % save upper and lower latitudes of the watershed
            
            
            % note: this trimming process is ~4x faster than just using inpolygon
            Ix = find(lon > min(lonlims) & lon < max(lonlims)); % find longitudes between limits
            lon = lon(Ix); % cut down based on longitudes
            lat = lat(Ix);
            terrain = terrain(Ix);
            bright = bright(Ix);
            std = std(Ix);
            can = can(Ix);
            
            Iy = find(lat > min(latlims) & lat < max(latlims)); % find lats between limits
            lat = lat(Iy); % cut down based on latitudes
            lon = lon(Iy);
            terrain = terrain(Iy);
            bright = bright(Iy);
            std = std(Iy);
            can = can(Iy);
            
            % now do a final trip to the actual watershed bounds
            wats = [watershed.X', watershed.Y']; % save coordinate tuples from the waterhsed shapefile
            Ifinal = inpolygon(lon, lat, wats(:,1), wats(:,2));
            
            lat = lat(Ifinal); % save data as vectors after final clipping
            
            if ~isempty(lat) % continue if data is inside the region of interest
                lon = lon(Ifinal); % conitnue saving data
                terrain = terrain(Ifinal);
                bright = bright(Ifinal);
                std = std(Ifinal);
                can = can(Ifinal);
                can(can > 1000) = nan; % change canopy elevation no data value to nan
                
                
                % use deg2utm script to write easting and northing values
                [easts, norths] = deg2utm(lat, lon);
                
                % create a final structure with all of the data
                s = struct(); % set structure for fields
                s.Latitude = lat; % set lats
                s.Longitude = lon; % set lons
                s.Elevation = terrain; % fill additional fields
                s.Elevation(s.Elevation >= 10e20) = NaN; % set Nans (icesat2 default nan value is 4.028e38)
                s.Canopy = can;
                s.std = std;
                s.Easting = easts;
                s.Northing = norths;
                s.Brightness_Flag = bright;
                
                testgeoid = geoidheight(s.Latitude, s.Longitude); % get geoid offsets for lat/lons.
                % Note^ the orange error is just fine -- in our case negative lats are
                % west, and the geoid height function is equipped to deal with that
    
                zrefined = s.Elevation - testgeoid; % apply geoid shift
                s.zRefined = zrefined; % save in the structure
                s.Geoid = testgeoid;
                
                
                table = struct2table(s); % convert to a table
                h5filename = h5test(end-38:end-3); % save h5 filename
                outputname = [outputdir, h5filename, '_', beam, '.csv']; % save full filename
                writetable(table, outputname)
            end
        end
        % if not ATL08, ATL06?
    elseif strcmp(h5files(filecounter).name(1:5), 'ATL06') == 1 % enter ATL06 commands
        h5test = [h5files(filecounter).folder,'/',h5files(filecounter).name]; % get string pointing to n-th h5 file
        for beamcount = 1:length(beams)
            beam = beams{beamcount}; % set beam
            % read in data
            terrain = h5read(h5test, ['/',beam,'/land_ice_segments/h_li']); % read terrain means
            lat = h5read(h5test,  ['/',beam,'/land_ice_segments/latitude']); % read lats
            lon = h5read(h5test,  ['/',beam,'/land_ice_segments/longitude']); % read lons
            bright = h5read(h5test,  ['/',beam,'/land_ice_segments/sigma_geo_h']); % read in vertical geolocation error
            std = h5read(h5test,  ['/',beam,'/land_ice_segments/h_li_sigma']); % standard deviation
            
            % crop data to area of interest
            % first, crude cropping
            lonlims = watershed.BoundingBox(:, 1); % save upper and lower longitudes of the watershed
            latlims = watershed.BoundingBox(:, 2); % save upper and lower latitudes of the watershed
            
            
            % note: this trimming process is ~4x faster than just using inpolygon
            Ix = find(lon > min(lonlims) & lon < max(lonlims)); % find longitudes between limits
            lon = lon(Ix); % cut down based on longitudes
            lat = lat(Ix);
            terrain = terrain(Ix);
            bright = bright(Ix);
            std = std(Ix);
            
            Iy = find(lat > min(latlims) & lat < max(latlims)); % find lats between limits
            lat = lat(Iy); % cut down based on latitudes
            lon = lon(Iy);
            terrain = terrain(Iy);
            bright = bright(Iy);
            std = std(Iy);
            
            % now do a final trip to the actual watershed bounds
            wats = [watershed.X', watershed.Y']; % save coordinate tuples from the waterhsed shapefile
            Ifinal = inpolygon(lon, lat, wats(:,1), wats(:,2));
            
            lat = lat(Ifinal); % save data as vectors after final clipping
            
            if ~isempty(lat) % continue if data is inside the region of interest
                lon = lon(Ifinal); % conitnue saving data
                terrain = terrain(Ifinal);
                bright = bright(Ifinal);
                std = std(Ifinal);
                
                
                % use deg2utm script to write easting and northing values
                [easts, norths] = deg2utm(lat, lon);
                
                % create a final structure with all of the data
                s = struct(); % set structure for fields
                s.Latitude = lat; % set lats
                s.Longitude = lon; % set lons
                s.Elevation = terrain; % fill additional fields
                s.Elevation(s.Elevation >= 10e20) = NaN; % set Nans (icesat2 default nan value is 4.028e38)
                s.std = std;
                s.Easting = easts;
                s.Northing = norths;
                s.Vert_Geo_error = bright;

                testgeoid = geoidheight(s.Latitude, s.Longitude); % get geoid offsets for lat/lons.
                % Note^ the orange error is just fine -- in our case negative lats are
                % west, and the geoid height function is equipped to deal with that
    
                zrefined = s.Elevation - testgeoid; % apply geoid shift
                s.zRefined = zrefined; % save in the structure
                s.Geoid = testgeoid;
                
                
                table = struct2table(s); % convert to a table
                h5filename = h5test(end-38:end-3); % save h5 filename
                outputname = [outputdir, h5filename, '_', beam, '.csv']; % save full filename
                writetable(table, outputname)
            end
        end
    end
end
