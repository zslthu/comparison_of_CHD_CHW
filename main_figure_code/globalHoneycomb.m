function  [value,hc] = globalHoneycomb(cdata, data, class_num, data_min, data_max, grid, pi,col,ctitle,flag_text,font_size)

fullpath = mfilename('fullpath');
[path,~]=fileparts(fullpath);

S = shaperead([path '/map/AR6_NEW.shp']);

% data = data_diff_space; % "example" is the name of input array
% class_num = 10;  % Total number of categories
%
% data_min = -0.15;
% data_max = 0.15;
%
% grid = 0.5;

%% set coordinate
[lon, lat] = meshgrid(-180+grid/2:grid:180-grid/2,90-grid/2:-grid:-90+grid/2);
% [lon, lat] = meshgrid(-179.5:1:179.5,89.5:-1:-89.5);
[m,n] = size(data);
mask = nan(m,n);
[Sm,~] = size(S);

%% make mask
for i = 1:Sm
    [in,~] = inpolygon(lon,lat,S(i).X,S(i).Y);
    mask(in) = i;
end
% Sm

%% calculate
value = nan(Sm,1);
for i = 1:Sm
    value(i) = nanmean(data(mask==i));  % Average the selected area
    % value(i) = mean(data(mask==i),'omitnan');  % Average the selected area
    value(i) = median(data(mask==i),'omitnan');  % Average the selected area
end
% value
% cdata = colormap('jet');  % Colormap name
% cdata = flip(colormap(brewermap([],'RdBu')));
[m,~] = size(cdata);
ind =  ceil(linspace(1, m-1, class_num));
cdata_sel = cdata(ind,:);

% data_min = min(value);
% data_max = max(value);

boundary = linspace(data_min, data_max, class_num+1);
class_ind = discretize(value, boundary);

class_ind(value<=data_min) = 1;
class_ind(value>=data_max) = class_num;

color_list = cdata_sel(class_ind,:);


%% prepare for plot
area_list = strings(Sm,1);
for i = 1:Sm
    area_list(i) = S(i).Acronym;
end
[~,area_ind] = sort(area_list);
value = value(area_ind);
color_list = color_list(area_ind,:);
% area_list'

%% plot
if pi
    % 如果没有提供optionalParam1的值，则默认为10
    if nargin < 11
        font_size = 5;
    end
    honeycomb_map(roundn(value,-4) ,color_list, flag_text,font_size)

    %%%- set colorbar
    cmin = col(1);
    cmax = col(2);
    stepi = col(3);

    hc = colorbar('horiz');
    %%%- position
    pos_hc = get(hc,'Position');
    pos_hc(4) = pos_hc(4)*0.4;
    pos_hc(2) = pos_hc(2)-pos_hc(4)*7;
    set(hc,'Position',pos_hc)
    %%%- ticklabel
    caxis([cmin,cmax])
    set(hc,'ticks',cmin:stepi:cmax)
    title(hc,ctitle,'FontSize',12)
end
end

