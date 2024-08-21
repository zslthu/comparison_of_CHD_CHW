clc;clear;

nameModel = {'ACCESS-CM2' 'ACCESS-ESM1-5' 'BCC-CSM2-MR' 'CanESM5' 'CESM2' 'CMCC-ESM2' 'EC-Earth3' ...
    'FGOALS-g3' 'GFDL-CM4' 'IPSL-CM6A-LR' 'MIROC6' 'MPI-ESM1-2-LR' 'MRI-ESM2-0' 'NorESM2-LM'}; 
im_n = length(nameModel);

%% read the region map
grid_1 = './map/AR6_NEW_1.tif';
grid_1 = double(imread(grid_1));
grid_1(grid_1==255) = nan;
grid_code = unique(grid_1(~isnan(grid_1)));

grid_05 = './map/AR6_NEW_05.tif';
grid_05 = double(imread(grid_05));
grid_05(grid_05==255) = nan;

nameRegion_shp = {'PAC' 'GIC' 'NWN' 'NEN' 'WNA' 'CNA' 'ENA' 'NCA' 'SCA' ...
    'CAR' 'NWS' 'NSA' 'NES' 'SAM' 'SWS' 'SES' 'SSA' 'NEU' ...
    'WCE' 'EEU' 'MED' 'SAH' 'WAF' 'CAF' 'NEAF' 'SEAF' 'WSAF' ...
    'ESAF' 'RAR' 'WSB' 'ESB' 'RFE' 'WCA' 'ECA' 'TIB' 'EAS' ...
    'ARP' 'SAS' 'SEA' 'NAU' 'CAU' 'EAU' 'SAU' 'NZ'};
nameRegion_label = {'NWN','NEN','GIC','WNA','CNA','ENA','NCA','SCA',...
    'CAR','NWS','NSA','SAM','NES','SWS','SES','SSA','NEU','RAR','WCE','EEU',...
    'WSB','ESB','RFE','MED','WCA','ECA','TIB','EAS','SAH','ARP','SAS','SEA',...
    'WAF','CAF','NEAF','WSAF','SEAF','ESAF','NAU','CAU','EAU','SAU','NZ','PAC'};

positions = zeros(size(nameRegion_label));
for i = 1:length(nameRegion_label)
    positions(i) = find(ismember(nameRegion_shp, nameRegion_label{i}));
end

%% figure set
fig = 1;
if fig
    figure;
    Nh=2; Nw=3; 
    gap = [0.1 0.03];
    marg_h = [0.2 0.05];
    marg_w = [0.35 0.06];
    [ha, pos] = tight_subplot(Nh, Nw, gap, marg_h, marg_w);
    set(gcf,'Position',[1200 300 820 600])
end

AB1 = {'a' 'b'};
AB2 = {'c' 'd'};
Ylabel = {'{\itS}_{CHD}' '{\itS}_{CHW}'};

load('../results/frequency/skill_score.mat')
for iv = 1:2
    for ic = 1:length(grid_code)
        skill_grid_1(grid_1 == grid_code(ic)) = mean(skill_score(ic,:,iv));
    end    

    %% figure map
    if fig ==1
        cdata = colormap(brewermap([],'YlGnBu'));

        loci = (iv-1)*3+1;
        axes(ha(loci))

        pos_change = cell2mat(pos(loci));
        pos_change(1) = pos_change(1) - pos_change(3)*1.2-0.02;
        pos_change(3) = pos_change(3)*2.2;
        set(gca,'position',pos_change)

        data =  skill_grid_1;
        data_min = 0.3;
        data_max = 1;
        class_num = 10;
        grid = 1;
        col = [data_min data_max 0.1];
        ctitle = Ylabel{iv};
        flag_text = 1;
        globalHoneycomb(cdata, data, class_num, data_min, data_max, grid, 1, col,ctitle, flag_text, 8);
        text(0,1.05,AB1(iv),'units','normalize','fontsize',16,'fontweight','bold')

        if iv ==1
            text(1,1.05,'c','units','normalize','fontsize',16,'fontweight','bold')
            text(1.6,1.05,'d','units','normalize','fontsize',16,'fontweight','bold')
        end

        if iv==2
            x = [0.43 0.47];
            y = [0.555 0.555];
            annotation('textarrow',x,y,'String','better performance')
            x = [0.43 0.47];
            y = [0.13 0.13];
            annotation('textarrow',x,y,'String','better performance')
        end
    end
end

%% heatmap for each model
for iv = 1:2
    skill_temp = skill_score(:,:,iv);
    skill_temp = skill_temp(positions,:); 

    if fig
        loci = iv+1;
        axes(ha(loci))

        pos_change = cell2mat(pos(loci));
        pos_change(2) = pos_change(2) - pos_change(4)*1.1;
        pos_change(4) = pos_change(4)*2.1;
        set(gca,'position',pos_change)

        h = heatmap(skill_temp,'ColorbarVisible','off');
        colormap(cdata)
        clim([0.3, 1])

        %-format
        ax = gca;
        ax.XData = nameModel;
        if iv ==1
            ax.YData = nameRegion_label;
        else
            ax.YDisplayLabels = nan(size(ax.YDisplayData));
        end
        title(Ylabel{iv})
        axes(ha(5)); axis off;
        axes(ha(6)); axis off;
    end
    cons_mode(:,iv) = abs(1-mean(skill_temp./mean(skill_temp),2));
end


fig_2 = 1;
if fig_2
    figure;
    %% figure set
    Nh=1; Nw=2; % subplot
    gap = [0.1 0.01];
    marg_h = [0.25 0.05];
    marg_w = [0.08 0.08];
    [ha, pos] = tight_subplot(Nh, Nw, gap, marg_h, marg_w);
    set(gcf,'Position',[1200 300 800 250])

    Ylabel2 = {'{\itCNS}_{CHD}' '{\itCNS}_{CHW}'};
    ABlab = {'a' 'b'};
    for iv = 1:2
        for ic = 1:length(grid_code)
            skill_grid_1(grid_1 == grid_code(ic)) = cons_mode(ic,iv);
        end
        median(skill_grid_1(:),'omitnan')

        cdata = colormap(brewermap([],'YlGnBu'));

        axes(ha(iv))
        data =  skill_grid_1;
        data_min = 0;
        data_max = .6;
        class_num = 10;
        grid = 1;
        col = [data_min data_max 0.1];
        ctitle = Ylabel2{iv};
        flag_text = 1;
        globalHoneycomb(cdata, data, class_num, data_min, data_max, grid, 1, col,ctitle, flag_text,8);
        if iv==2
            x = [0.16 0.13];
            y = [0.11 0.11];
            annotation('textarrow',x,y,'String','better consistency')
            x = [0.59 0.56];
            y = [0.11 0.11];
            annotation('textarrow',x,y,'String','better consistency')
        end
        text(-0.02,0.98,ABlab{iv},'units','normalized','fontsize',14,'fontweight','bold')
    end
end