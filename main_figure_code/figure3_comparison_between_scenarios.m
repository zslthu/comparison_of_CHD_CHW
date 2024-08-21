clc;clear;

nameScenario = {'his_nat' 'LS3MIP_pdlc'};
nameScenario_dir = {'hisnat' 'ls'};
nameModel_his_nat = {'CESM2' 'IPSL-CM6A-LR' 'MIROC6' 'ACCESS-CM2' 'ACCESS-ESM1-5' 'BCC-CSM2-MR' 'CanESM5' 'CNRM-CM6-1'  ...
    'FGOALS-g3' 'GFDL-CM4' 'MRI-ESM2-0' 'NorESM2-LM'};
nameModel_LS3MIP_pdlc = {'CESM2' 'IPSL-CM6A-LR' 'MIROC6' 'CMCC-ESM2' 'EC-Earth3' 'MPI-ESM1-2-LR'}; 
nameV = {'shdi' 'shwi'};

outputDir = '../results/frequency/cmip6_scenarios/';

%% plot set
figure;
Nh=2; Nw=3; % subplot
gap = [0.04 0.0];
marg_h = [0.08 0.2];
marg_w = [0.08 0.08];
[ha, pos] = tight_subplot(Nh, Nw, gap, marg_h, marg_w);
set(gcf,'Position',[1200 300 700 450])

cdata_line = [0.8 0 0; 0 0 0.8];
cdata = colormap(flip(brewermap([],'RdBu')));

legendS(1,:) = {'{\itf}_{CHD-hist-nat}' '{\itf}_{CHW-hist-nat}'};
legendS(2,:)= {'{\itf}_{CHD-LFMIP-pdLC}' '{\itf}_{CHW-LFMIP-pdLC}'};
nameT = {'Impact of anthropogenic climate change' 'Impact of land-atmosphere feedback'};
ABlabel1 = {'b' 'c' 'd'};

%% plot time series
for is = 1 % scenario
    eval(['nameModel = nameModel_' char(nameScenario(is))]);
    im_n = length(nameModel);
    for iv = 1:2 % shdi % shwi
        data_yearly = [];
        data_yearly_his = [];

        PR_mean_model = [];
        PR_std_model = [];
        for im = 1:im_n
            Dir = [outputDir '/' char(nameScenario_dir(is))];
            FileDir = [Dir,'/Frequency_',char(nameV(iv)),'_all_',char(nameModel(im)), '.mat'];
            data = load(FileDir);

            Dir = [outputDir '/his/'];
            FileDir = [Dir,'/Frequency_',char(nameV(iv)),'_all_',char(nameModel(im)), '.mat'];
            data_his = load(FileDir);

            data_yearly(:,im)= mean(data.data_frequency,2,'omitnan');
            data_yearly_his(:,im)= mean(data_his.data_frequency,2,'omitnan');

            PR_mean_model(im,iv) = median(1-1./(data_yearly_his(:,im)./data_yearly(:,im)),'omitnan');
            PR_std_model(im,iv) = std(1-1./(data_yearly_his(:,im)./data_yearly(:,im)),'omitnan');

            %% figure each model
            xloc = [0.1 -0.1];
            AB = {'a' 'b'};

            axes(ha(1))
            pos_change = cell2mat(pos(1));
            pos_change(1) = pos_change(1)+0.04;
            pos_change(3) = pos_change(3)*2.7;
            pos_change(4) = pos_change(4)*0.9;
            pos_change(2) = pos_change(2)+pos_change(4)*0.3;
            set(gca,'position',pos_change)

            h(im,iv) = plot(im-xloc(iv),PR_mean_model(im,iv),'*','Color',cdata_line(iv,:),'LineWidth',1.5);hold on
            ymin = PR_mean_model(im,iv)-PR_std_model(im,iv);
            ymax = PR_mean_model(im,iv)+PR_std_model(im,iv);
            plot(ones(2,1)*(im-xloc(iv)),[ymin ymax],'k-');hold on;

            %-format
            ylabel('{\itFAR}')
            ylim([-0.4 1.2])
            xlim([0 im_n+1])
            set(gca,'ytick',[-0.4:0.4:1.2])
            set(gca,'xticklabel',nameModel,'xtick',1:im_n,'xticklabelRotation',45)
            title(nameT(is))
            text(0,1.08,'a','unit','normalize','fontweight','bold','fontsize',14)
            if im==1
                plot([0 13],[0 0],'k--')
            end
        end

        data_yearly_mm = mean(data_yearly,2,'omitnan');   
        data_yearly_his_mm = mean(data_yearly_his,2,'omitnan');   

        PR_temp(:,iv) = data_yearly_his_mm./data_yearly_mm;
        PR_temp(:,iv) = 1-1./(data_yearly_his_mm./data_yearly_mm);

        PR_map  = flipud(reshape(PR_temp(:,iv),180,360));
        PR_mean(is,iv) = median(PR_temp(:,iv),'omitnan');

        %% figure ensemble map
        if is==1
            Ylabel = {'{\itFAR}_{CHD}^{acc}','{\itFAR}_{CHW}^{acc}'};
        else
            Ylabel = {'{\itFAR}_{CHD}^{laf}','{\itFAR}_{CHW}^{laf}'};
        end

        loc_ha = 3+iv;
        axes(ha(loc_ha))

        data =  PR_map;
        data_min = -1.2;
        data_max = 1.2;
        class_num = 10;
        grid = 1;
        col = [data_min data_max 0.3];
        ctitle = char(Ylabel(iv));
        globalHoneycomb(cdata, data, class_num, data_min, data_max, grid, 1, col,ctitle,1);
        text(-0.05,0.9,char(ABlabel1(iv)),'Units','normalized','FontWeight','bold','FontSize',14)
    end
    if is==1
        legend([h(1,1),h(1,2)],{'{\itFAR}_{CHD}^{acc}','{\itFAR}_{CHW}^{acc}'},'box','off','location','southwest','Orientation','horizontal')
    else
        legend([h(1,1),h(1,2)],{'{\itFAR}_{CHD}^{laf}','{\itFAR}_{CHW}^{laf}'},'box','off','location','northwest','Orientation','horizontal')
    end

    PR_diff_map  = flipud(reshape(PR_temp(:,1)-PR_temp(:,2),180,360));

    %% figure ensemble map - difference
    if is==1
        Ylabel = {'{\itFAR}_{CHD}^{acc}','{\itFAR}_{CHW}^{acc}'};
    else
        Ylabel = {'{\itFAR}_{CHD}^{laf}','{\itFAR}_{CHW}^{laf}'};
    end

    loc_ha = 6;
    axes(ha(loc_ha))

    data =  PR_diff_map;
    data_min = -0.5;
    data_max = 0.5;
    class_num = 10;
    grid = 1;
    col = [data_min data_max 0.1];
    ctitle = [Ylabel{1} ' - ' Ylabel{2}];
    globalHoneycomb(cdata, data, class_num, data_min, data_max, grid, 1, col,ctitle,1);
    text(-0.05,0.9,char(ABlabel1(3)),'Units','normalized','FontWeight','bold','FontSize',14)

end
axes(ha(2));axis off;
axes(ha(3));axis off;



