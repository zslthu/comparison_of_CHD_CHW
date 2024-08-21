clc;clear;

nameScenario = {'hist-nat' 'LFMIP-pdlc'}; 
nameScenario2 = {'hist_nat' 'LS3MIP_pdlc'}; 
nameModel_hist_nat = {'CESM2' 'IPSL-CM6A-LR' 'MIROC6' 'ACCESS-CM2' 'ACCESS-ESM1-5' 'BCC-CSM2-MR' 'CanESM5' 'CNRM-CM6-1'  ...
    'FGOALS-g3' 'GFDL-CM4' 'MRI-ESM2-0' 'NorESM2-LM'};
nameModel_LS3MIP_pdlc = {'CESM2' 'IPSL-CM6A-LR' 'MIROC6' 'CMCC-ESM2' 'EC-Earth3' 'MPI-ESM1-2-LR'};

outputDir = '../results/attribution';
nameVar = {'tas' 'pr' 'relation'};
nameCE_F = {'SP-TI_R' '-SP-TI_R' }; % 'shdi' 'shwi'
nameCE_F2 = {'SPI-STI' '-SPI-STI' }; % 'shdi' 'shwi'

%% plot set
fig_box =1;
if fig_box ==1
    Nh=1; Nw=4;
    gap = [0.04 0.04];
    marg_h = [0.15 0.1]; 
    marg_w = [0.05 0.05];
    [ha, pos] = tight_subplot(Nh, Nw, gap, marg_h, marg_w);
    set(gcf,'Position',[1200 200 750 350])

    nameX = {'{\itta}' '{\itpr}' '{\itcorr}'};
    pointSym = {'*' '+' 'o' 'x' 's' 'd' 'p' 'h' '^' 'v' '<' '>' 'x' 's'};
end

%% plot set
fig_spatial =0;
if fig_spatial ==1
    Nh=3; Nw=4; 
    gap = [0.06 0.0];
    marg_h = [0.15 0.1]; 
    marg_w = [0.05 0.05];
    [ha, pos] = tight_subplot(Nh, Nw, gap, marg_h, marg_w);
    set(gcf,'Position',[1200 200 800 600])

    nameV2 = {'ta' 'pr' 'corr'};
    ABlabel = {'a' 'd' 'g' 'j' 'b' 'e' 'h' 'k' 'c' 'f' 'i' 'l'};
end
nameV = {'CHD' 'CHW'};
nameT = {'anthropogenic climate change' 'land-atmosphere feedback'};

%% calculation
for is = 1:2 % scenario
    eval(['nameModel = nameModel_' char(nameScenario2(is))]);
    im_n = length(nameModel);

    data_im = [];
    diff_mean_im = [];
    for ic = 1:2
        for im = 1:im_n
            %- original
            FileDir = [outputDir '/Joint_probability_' nameScenario{is} '/original/' ...
                nameCE_F{ic} '/' nameCE_F2{ic} '*' nameModel{im} '*.npy'];

            temp_dir = dir(FileDir).folder;
            temp_name = dir(FileDir).name;
            fileinfo = [temp_dir '/' temp_name];

            data_base = readNPY(fileinfo);   
            data_base(data_base==0) = nan;
            for iv = 1:3 
                %- change
                FileDir = [outputDir '/Joint_probability_' nameScenario{is} '/change_' nameVar{iv},'/' ...
                    nameCE_F{ic} '/' nameCE_F2{ic} '*' nameModel{im} '*.npy'];

                temp_dir = dir(FileDir).folder;
                temp_name = dir(FileDir).name;
                fileinfo = [temp_dir '/' temp_name];

                data_change = readNPY(fileinfo);   
                data_change(data_change==0) = nan;

                data_diff = data_change - data_base;

                data_im(:,iv,im) = data_diff;
                data_mean_im(iv,im)  = mean(data_diff,'omitnan');
            end
        end
        data_iv = mean(data_im,3,'omitnan');

        %% boxplot
        if fig_box
            %%%%%%%%%%%%%%%%% figure
            loc_ha = (is-1)*2+ic;
            axes(ha(loc_ha))

            boxplot(data_iv,'Symbol','','Widths',0.3,'Colors','k');hold on;
            h = findobj(gca, 'Tag', 'Median');
            set(h, 'Color', 'red', 'LineWidth', 2);

            plot([0 4],[0 0],'k--');hold on;
            for im = 1:im_n
                hl{im} = plot([1 2 3],data_mean_im(:,im),'.','Marker',pointSym{im},'LineWidth',1); hold on
            end
            ylim([-0.015 0.035])

            %- format
            xticklabels(nameX)
            yticks([-0.01:0.01:0.03])
            ylabel(['Contribution to ' nameV{ic} ' (\Delta{\itP}^{' nameV{ic} '})'])
            if ic == 1
                hLegend = legend([hl{:}],nameModel,'box','off','location','northeast','fontsize',8);
                text(0.24,1.05, ['Impact of ' nameT{is}],'units','normalized','fontsize',12,'fontweight','bold')
            end
            set(gca,'TickLabelInterpreter','tex')
            ABlab = {'a' 'b' 'c' 'd'};
            text(-0.1,0.98,ABlab{loc_ha},'units','normalized','fontsize',14,'fontweight','bold')
        end

        %% spatial
        if fig_spatial
            cdata = colormap(flipud(brewermap([],'RdBu')));

            for iv = 1:3
                data_diff_mm = flipud(reshape(data_iv(:,iv),180,360));

                loc_ha = (is-1)*2+ic;
                loc_ha = (iv-1)*4+loc_ha;
                axes(ha(loc_ha))

                if is==2
                    pos_change = cell2mat(pos(loc_ha));
                    pos_change(1) = pos_change(1)+pos_change(3)*0.08;
                    set(gca,'Position',pos_change)
                end

                pi = 1;
                data =  data_diff_mm;

                data_min = -0.05;
                data_max = 0.05;
                step_i   = 0.01;
                if iv~=1
                    data_min = -0.02;
                    data_max = 0.02;
                    step_i   = 0.01/2.5;
                end
                class_num =10;
                grid = 1;
                col = [data_min data_max step_i*2];
                ctitle = ['\Delta{\itP}_{' nameV2{iv} '}^{' nameV{ic} '}'];
                [~,hc] = globalHoneycomb(cdata, data, class_num, data_min, data_max, grid, pi,col,ctitle,0);

                text(-0.05,0.9,char(ABlabel(loc_ha)),'Units','normalized','FontWeight','bold','FontSize',14)
                if iv==1 && ic==1
                    text(0.24,1.05, ['Impact of ' nameT{is}],'units','normalized','fontsize',13,'fontweight','bold')
                end
            end
        end
    end
end