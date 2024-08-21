clc;clear;

outputDir_obs = '../results/frequency/obs';
outputDir = '../results/frequency/cmip6_his/';
nameModel = {'ACCESS-CM2' 'ACCESS-ESM1-5' 'BCC-CSM2-MR' 'CanESM5' 'CESM2' 'CMCC-ESM2' 'EC-Earth3' ...
    'FGOALS-g3' 'GFDL-CM4' 'IPSL-CM6A-LR' 'MIROC6' 'MPI-ESM1-2-LR' 'MRI-ESM2-0' 'NorESM2-LM'};
im_n = length(nameModel);
nameV = {'shdi' 'shwi'};

%% Plot settings
Nh=3; Nw=3; % subplot
gap = [0.04 0.01];
marg_h = [0.08 0.08];
marg_w = [0.08 0.08];
[ha, pos] = tight_subplot(Nh, Nw, gap, marg_h, marg_w);
set(gcf,'Position',[1400 200 800 600])

ABlabel0 ={'a' 'b'};
ABlabel1 = {'c' 'd' 'e'};
ABlabel2 = {'f' 'g' 'h'};

LineS = {'-.' '--'};
cdata_line = [0.8 0 0; 0 0 0.8];
cdata = colormap(flip(brewermap([],'RdBu')));
Ylabel = {'CHD' 'CHW'};

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% temproal
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  obs
data_temp = [];
im = 1;
for iv = 1:2
    FileDir = [outputDir_obs,'/Frequency_',char(nameV(iv)),'_all_cru.mat'];
    load(FileDir)

    data_temp(:,im) = [nan(1,64-size(data_frequency,2)),mean(data_frequency,1,'omitnan')];
    data_temp_change = [mean(data_frequency(:,end-20+1:end)-data_frequency(:,1:20),2,'omitnan')];

    data_yearly = data_temp(:,im);

    diff_mean = mean(data_yearly(end-20+1:end))-mean(data_yearly(1:20));
    diff_std = std(data_temp_change,'omitnan');

    %% plot temporal
    loc_ha = iv;
    axes(ha(loc_ha))

    if iv ==1
        pos_c = cell2mat(pos(loc_ha));
        pos_c(1) = pos_c(1)+0.02;
        pos_c(3) = pos_c(3)+pos_c(3)*0.35;
        pos_c(4) = pos_c(4)*0.9;
    else
        pos_c = cell2mat(pos(loc_ha));
        pos_c(1) = pos_c(1)+pos_c(3)*0.35+0.07;
        pos_c(3) = pos_c(3)+pos_c(3)*0.35;
        pos_c(4) = pos_c(4)*0.9;
    end
    set(ha(loc_ha),'position',pos_c)

    h(im,iv) = plot(1951:2014,data_yearly,'Color',cdata_line(iv,:),'LineStyle',char(LineS(im)),'LineWidth',1);hold on

    nameCRU = {'\Delta{\itf}_{CHD}^{ cru}' '\Delta{\itf}_{CHW}^{ cru}'};
    text(0.02,0.9,[nameCRU{iv} '= ' num2str(diff_mean,'%.2f') '\pm' num2str(diff_std,'%.2f')],'units','normalized')
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  cmip6-his
data_temp = [];
data_temp_change = [];
data_space = [];
for iv = 1:2
    for im = 1:im_n
        FileDir = [outputDir,'/Frequency_',char(nameV(iv)),'_all_',char(nameModel(im)), '.mat'];
        load(FileDir)

        data_temp(:,im) = [nan(1,64-size(data_frequency,2)),mean(data_frequency,1,'omitnan')];
        data_temp_temp(:,im) = [mean(data_frequency(:,end-20+1:end)-data_frequency(:,1:20),2,'omitnan')];
    end
    data_yearly = smoothcurve(mean(data_temp,2,'omitnan'),1);
    data_yearly_range = range(data_temp,2);

    data_temp_change = data_temp_temp(:);

    diff_mean = mean(data_yearly(end-20+1:end))-mean(data_yearly(1:20));
    diff_std = std(data_temp_change,'omitnan');

    %% plot temporal
    loc_ha = iv;
    axes(ha(loc_ha))

    h(2,iv) =plot(1951:2014,data_yearly,'Color',cdata_line(iv,:),'LineStyle','-','LineWidth',1);hold on
    %%- range
    x = [1951:2014]';
    y25 = data_yearly - data_yearly_range/2;
    y75 = data_yearly + data_yearly_range/2;
    hl = fill([x',fliplr(x')],[y25',fliplr(y75')],cdata_line(iv,:),'FaceAlpha',0.1 ,'LineStyle', 'none');
    %%- range

    %%%%%%%%%%%% -format
    ylim([0 0.5]); xlim([1950 2015])
    ylabel(['{\itf}_{' char(Ylabel(iv)) '}'])
    legend([h(1,iv),h(2,iv)],{'cru'  'cmip6'},'box','off','Location','northeast')

    nameCMIP6 = {'\Delta{\itf}_{CHD}^{ cmip6}' '\Delta{\itf}_{CHW}^{ cmip6}'};
    text(0.02,0.75,[nameCMIP6{iv} '= ' num2str(diff_mean,'%.2f') '\pm' num2str(diff_std,'%.2f')],'units','normalized')

    text(-0.05,1.1,char(ABlabel0(iv)),'Units','normalized','FontWeight','bold','FontSize',14)

end
axes(ha(3));axis off


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% spatial
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  obs
data_temp =[];
data_space = [];
im = 1;
for iv = 1:2 % shdi % shwi

    FileDir = [outputDir_obs,'/Frequency_',char(nameV(iv)),'_all_cru.mat'];
    load(FileDir)

    data_temp = [mean(data_frequency(:,end-20+1:end)-data_frequency(:,1:20),2,'omitnan')];
    data_space(:,iv) = mean(data_temp,2,'omitnan');
    data_space_ma = flipud(reshape(data_space(:,iv),360,720));

    %% plot spatial
    loc_ha = 3+iv;
    axes(ha(loc_ha))

    data =  data_space_ma;
    data_min = -0.3;
    data_max = 0.3;
    data_step = 0.1;
    class_num = 10;
    grid = 0.5;
    col = [data_min data_max data_step];
    ctitle = ['\Delta{\itf}_{' Ylabel{iv} '}^{ cru}'];
    flag_text = 1;
    globalHoneycomb(cdata, data, class_num, data_min, data_max, grid, pi,col,ctitle,flag_text);
    text(-0.05,0.9,char(ABlabel1(iv)),'Units','normalized','FontWeight','bold','FontSize',14)
end
data_diff_space = flipud(reshape(data_space(:,1)-data_space(:,2),360,720));

%% plot spatial - difference
loc_ha = 6;
axes(ha(loc_ha))

data =  data_diff_space;
data_min = -0.15;
data_max = 0.15;
data_step = 0.05;
class_num = 10;
grid = 0.5;
col = [data_min data_max data_step];
ctitle = ['\Delta{\itf}_{' Ylabel{1} '}^{ cru} - \Delta{\itf}_{' Ylabel{2} '}^{ cru}'];
globalHoneycomb(cdata, data, class_num, data_min, data_max, grid, pi,col,ctitle,flag_text);
text(-0.05,0.9,ABlabel1(3),'Units','normalized','FontWeight','bold','FontSize',14)


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  cmip6-his
data_temp = [];
data_space = [];
for iv = 1:2 % shdi % shwi
    for im = 1:im_n
        FileDir = [outputDir,'/Frequency_',char(nameV(iv)),'_all_',char(nameModel(im)), '.mat'];
        load(FileDir)

        data_temp(:,im) = [mean(data_frequency(:,end-20+1:end)-data_frequency(:,1:20),2,'omitnan')];
    end
    data_space(:,iv) = mean(data_temp,2,'omitnan');
    data_space_mm = flipud(reshape(data_space(:,iv),180,360));

    %% plot spatial
    loc_ha = 6+iv;
    axes(ha(loc_ha))

    data =  data_space_mm;
    data_min = -0.3;
    data_max = 0.3;
    data_step = 0.1;
    class_num = 10;
    grid = 1;
    col = [data_min data_max data_step];
    ctitle = ['\Delta{\itf}_{' Ylabel{iv} '}^{ cmip6}'];
    globalHoneycomb(cdata, data, class_num, data_min, data_max, grid, pi,col,ctitle,flag_text);
    text(-0.05,0.9,char(ABlabel2(iv)),'Units','normalized','FontWeight','bold','FontSize',14)

end
data_diff_space = flipud(reshape(data_space(:,1)-data_space(:,2),180,360));

%% plot spatial - difference
loc_ha = 9;
axes(ha(loc_ha))

data =  data_diff_space;
data_min = -0.15;
data_max = 0.15;
data_step = 0.05;
class_num = 10;
grid = 1;
col = [data_min data_max data_step];
ctitle = ['\Delta{\itf}_{' Ylabel{1} '}^{ cmip6} - \Delta{\itf}_{' Ylabel{2} '}^{ cmip6}'];
globalHoneycomb(cdata, data, class_num, data_min, data_max, grid, pi,col,ctitle,flag_text);
text(-0.05,0.9,ABlabel2(3),'Units','normalized','FontWeight','bold','FontSize',14)


% % %%-save
% savefig_name = ['/Volumes/SYSU/Research/CompoundDrought/code/final_v1/figure_response_v2/new_Part1_1_obs_cmip6_v3'];
% exportgraphics(gcf,[savefig_name,'.pdf'],'ContentType','vector')
% exportgraphics(gcf,[savefig_name,'.jpg'],'Resolution',600)


