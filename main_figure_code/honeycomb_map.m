function honeycomb_map(value,color_list,flag_text,font_size)


area=["NWN","NEN","GIC","WNA","CNA","ENA","NCA","SCA","CAR","NWS","NSA","SAM","NES","SWS","SES","SSA","NEU","RAR","WCE","EEU","WSB","ESB","RFE","MED","WCA","ECA","TIB","EAS","SAH","ARP","SAS","SEA","WAF","CAF","NEAF","WSAF","SEAF","ESAF","MDG","NAU","CAU","EAU","SAU","NZ" "PAC"];
x=[2,4,7,1,3,5,2,3,5,6,8,7,9,6,8,7,16,22,15,17,19,21,23,16,18,20,22,24,15,17,21,25,14,16,18,15,17,16,19.5,26,25,27,26,28.5,28.5];
y=[25,25,25,22,22,22,19,16,16,13,13,10,10,7,7,4,25,25,22,22,22,22,22,19,19,19,19,19,16,16,16,16,13,13,13,10,10,7,9,12,9,9,6,5.5,13];
x(17:end)=x(17:end)-2;  % Atlantic Regional Interval

idele = find(area=='MDG');
area(idele) = [];
x(idele) = [];
y(idele) = [];

[area,area_ind] = sort(area);
x = x(area_ind);
y = y(area_ind);


n=size(x,2);
A=nan(n,6);
B=nan(n,6);
for i=1:n
    A(i,:)=[x(i),x(i)+1,x(i)+1,x(i),x(i)-1,x(i)-1].*(sqrt(3)/2);
    B(i,:)=[y(i),y(i)-1,y(i)-3,y(i)-4,y(i)-3,y(i)-1].*(1/2);
end

hold on; 
% value
for i=1:n
    % txt=[strcat(area(i)),num2str(value(i))];
    txt=[strcat(area(i))];
    hh = patch(A(i,:), B(i,:), color_list(i,:),'edgecolor','w');
    if flag_text
        text(A(i,1), B(i,1)-1,txt,'Color','black','FontSize',font_size,'HorizontalAlignment','center')
    end
end
% -(sqrt(3)/4)
hold off; 
axis equal;
axis off;

%             %%%- set colorbar
%        colorbar('horiz');
%         %%%- position
%         pos_hc = get(hc,'Position');
%         pos_hc(4) = pos_hc(4)*0.4;
%         pos_hc(2) = pos_hc(2)-pos_hc(4)*7;
%         set(hc,'Position',pos_hc)
%         %%%- ticklabel
%         stepi = 0.1;
%         caxis([data_min,data_max])
%         set(hc,'ticks',[data_min:stepi:data_max])
%         title(hc,'\Delta{\itf}_{CHD} - \Delta{\itf}_{CHW}','FontSize',12)
end
