% %for figure placement depending on monitor configuration
% work = true;
% if work == true
%     anchor_v = -500;
%     anchor_h = -1900;
% else
%     anchor_v = 50;
%     anchor_h = 0;
% end
% 
% ND_name = strcat(h5read("..\candidates_data.h5","/candidates/label"),num2str(h5read("..\candidates_data.h5","/candidates/ND_No.")));
% NDs_List = unique(ND_name);
% NDs_List = NDs_List(NDs_List ~= "S3_A11");

%% NPL S2B Map
NPL_array = h5read("..\exemplary_map.h5","/exemplary/NPL_map");
imfilt = NPL_array > 7;
f = figure(4); clf
imagesc(NPL_array);
stats = regionprops("table",imfilt,"Centroid", ...
        "MajorAxisLength","MinorAxisLength","BoundingBox","ConvexImage", ...
        "Area");
stats = [stats(1:7,:);stats(7,:);stats(9:end,:)]; %insert duplicate ND detected to replicate that in NPLs map (caused by differences in object detection algorithm)
stats.Centroid(7,:) = 0.99.*stats.Centroid(7,:);
diameters = mean([stats.MajorAxisLength stats.MinorAxisLength],2);
radii = diameters/2;
txt_anchor_L = zeros(height(stats),2);
txt_anchor_L(:,:) = stats.Centroid + [2 -2];
hold on
viscircles(stats.Centroid([1:4,6:end],:),2.*radii([1:4,6:end]),'color','yellow','LineWidth',2)
txt_fsize = 24;
for i = 1:height(stats)
    if ismember(i,9:11)
    text(txt_anchor_L(i,1),txt_anchor_L(i,2),[' \leftarrow ',num2str(i)],'Color','yellow','FontSize',txt_fsize,'Rotation',15,'FontWeight','bold')
    elseif i == 5
        viscircles(stats.Centroid(i,:),2.*radii(i),'color',[0 0.75 0.2],'LineWidth',3)
        text(txt_anchor_L(i,1),txt_anchor_L(i,2),['  \leftarrow ',num2str(i),' (B)'],'Color',[0 0.75 0.2],'FontSize',txt_fsize,'Rotation',-15,'FontWeight','bold')
    elseif i == 14
    text(txt_anchor_L(i,1),txt_anchor_L(i,2),[' \leftarrow ',num2str(i)],'Color','yellow','FontSize',txt_fsize,'Rotation',-20,'FontWeight','bold')
    elseif ismember(i,[8,16])
    text(txt_anchor_L(i,1),txt_anchor_L(i,2),[' \leftarrow ',num2str(i)],'Color','yellow','FontSize',txt_fsize,'Rotation',-45,'FontWeight','bold')
    else
    text(txt_anchor_L(i,1),txt_anchor_L(i,2),[' \leftarrow ',num2str(i)],'Color','yellow','FontSize',txt_fsize,'Rotation',-15,'FontWeight','bold')
    end
end
hold off
colormap hot
xlabel("")
ylabel("")
cb_npl = colorbar("southoutside","Limits",[0.8 1e2],"TickDirection","out","FontSize",28,...
    "TickLabelInterpreter","tex","Ticks",[1 10 100],"ticklabels",{"10^{0}","10^{1}","10^{2}"},...
    "Box","on","LineWidth",2);
clim([0.8 1e2])
cb_npl.Color = [1 1 0];
ylabel(cb_npl,"Count rate, I (kcps)","FontWeight","bold","FontSize",28)
axis square
set(gca,"XTickLabel",[],"YTickLabel",[],"FontWeight","bold","ColorScale","log","FontSize",28,...
    "XTick",[],"YTick",[],"XColor",[1 1 0],"YColor",[1 1 0],"LineWidth",2)
f.Position = [anchor_h anchor_v 550 660];

%% figure(1)
title_offset = [0.5, -0.3, 0];
f = figure(1);clf
UoM_array = h5read("..\exemplary_map.h5","/exemplary/UoM_map");
t = tiledlayout(1,3,"TileSpacing","compact");
nexttile;
imagesc(h5readatt("..\exemplary_map.h5","/exemplary/UoM_map/","x (µm)"),...
    h5readatt("..\exemplary_map.h5","/exemplary/UoM_map/","y (µm)"),UoM_array)
colormap bone
xlabel("")
ylabel("")
cb = colorbar("westoutside","FontWeight","bold","FontSize",20,"Box","on","LineWidth",2,"TickDirection","out",...
    "TickLabelInterpreter","tex","Ticks",[0.1 1 10 100],"ticklabels",["10^{-1}","10^{0}","10^{1}","10^{2}"]);
ylabel(cb,"Count rate, I (kcps)","FontSize",20,"FontWeight","bold")
xlim([-22 50]) %crop to 65x65 µm
ylim([-50 22])
axis square
set(gca,"XTickLabel",[],"YTickLabel",[],"FontWeight","bold","ColorScale","log","Box","on","LineWidth",2)
clim([1e-1 1e2])
title("(a)", 'Units', 'normalized', 'Position', title_offset)

%B
uom_xData = h5read("..\candidates_data.h5","/candidates/Irradiation_Power,_P_(µW)/18");
uom_yData = h5read("..\candidates_data.h5","/candidates/ND_{pdep}_(cps)/18");
npl_xData = h5read("..\candidates_data.h5","/candidates/Irradiation_Power,_P_(µW)_[NPL]/18");
npl_yData = h5read("..\candidates_data.h5","/candidates/ND_{pdep}_(cps)_[NPL]/18");
[Psat_uom,psat_full_uom,psat_ND_uom,psat_laser_uom] = pdep_components(uom_xData,uom_yData);
[Psat_npl,psat_full_npl,psat_ND_npl,psat_laser_npl] = pdep_components(npl_xData,npl_yData);

ax2a = axes(t);
ax2a.Layout.Tile=2;
fplot(ax2a,psat_ND_npl,'Color','r','LineStyle','-.','LineWidth',2) %NPL ND emission

ax2b = axes(t);
ax2b.Layout.Tile=2;
fplot(ax2b,psat_full_npl,'Color','r','LineWidth',4) %NPL Full Psat
hold on
plot(ax2b,NaN,'Color','k','LineWidth',4,'Visible','on') %Invisible NaN for legend
hold off

ax2c = axes(t);
ax2c.Layout.Tile=2;
xline(ax2c,Psat_npl(2),'Color','r','LineStyle','-','LineWidth',2) %NPL saturation

ax2d = axes(t);
ax2d.Layout.Tile=2;
fplot(ax2d,psat_laser_npl,'Color','r','LineStyle','--','LineWidth',2) %NPL Laser

ax2e = axes(t);
ax2e.Layout.Tile=2;
fplot(ax2e,psat_laser_uom,'Color','k','LineStyle','--','LineWidth',2) %UOM Laser

ax2f = axes(t);
ax2f.Layout.Tile=2;
fplot(ax2f,psat_ND_uom,'Color','k','LineStyle','-.','LineWidth',2) %UOM ND emission

ax2g = axes(t);
ax2g.Layout.Tile=2;
fplot(ax2g,psat_full_uom,'Color','k','LineWidth',4) %UOM Full Psat

ax2h = axes(t);
ax2h.Layout.Tile=2;
plot(ax2h,uom_xData,uom_yData,'ko','MarkerSize',5,'LineWidth',2,'MarkerFaceColor',[1 1 1]) %UOM datapoints

ax2i = axes(t);
ax2i.Layout.Tile=2;
plot(ax2i,npl_xData,npl_yData,'ro','MarkerSize',5,'LineWidth',2,'MarkerFaceColor',[1 1 1]) %NPL datapoints

ax2j = axes(t);
ax2j.Layout.Tile=2;
xline(ax2j,Psat_uom(2),'Color','k','LineStyle','--','LineWidth',2) %UOM saturation


uomnum = sort((0:200:1000)*1e-3,2,"ascend");
uomstr = string(uomnum);

nplnum = sort((0:1000:8000)*1e-3,2,"ascend");
nplstr = string(nplnum);

totax = [ax2a,ax2b,ax2c,ax2d,ax2e,ax2f,ax2g,ax2h,ax2i,ax2j];
for i = 1:numel(totax)
    set(totax(i),'fontsize',20,'YTick',0:2e4:12e4,'YTickLabel',"");
    if i == 1 %set bottom layer axes to white background
        totax(i).Color = [1 1 1];
    else
        totax(i).Color = 'none';
    end
    totax(i).Box = 'off';
    totax(i).PlotBoxAspectRatio = [1 1 1];
    set(totax(i),'fontsize',20,'YTick',0:2e4:12e4,'YTickLabel',"",'YLim',[0 12e4],'XTickLabel',"",'linewidth',2);
    if ismember(i,[1:4,numel(totax)-1]) %NPL Axes
        totax(i).XAxisLocation = "top";
        totax(i).XLim = ([0 12*Psat_npl(2)]);
        totax(i).YAxisLocation = "right";
        totax(i).XTick = nplnum;
    else %UOM axes
        totax(i).XAxisLocation = "bottom";
        totax(i).XLim = ([0 12*Psat_uom(2)]);
        totax(i).YAxisLocation = "left";
        totax(i).XTick = uomnum;
    end
end

set(ax2d,'fontsize',18,'YTick',0:2e4:12e4,'YTickLabel',"",'XTick',1e3.*nplnum,'XTickLabel',"",'YColor',[1 1 1],'XColor',[1 1 1]);
set(ax2j,'fontsize',18,'YTick',0:2e4:12e4,'YTickLabel',string((0:20:120)'),'XTick',1e3.*uomnum,'XTickLabel',uomstr);
set(ax2i,'fontsize',18,'YTick',0:2e4:12e4,'YTickLabel',"",'YColor','r','XTick',1e3.*nplnum,'XTickLabel',nplstr,'XColor','r');
legend(ax2b,["NPL","UoM"],"Location","northwest")
xlabel(ax2j,"Power, P (mW)")
xlabel(ax2i,"Power, P (mW)")
ylabel(ax2j,"Count rate, I (kcps)")
title("(b)", 'Units', 'normalized', 'Position', title_offset)

ax3 = nexttile;
benyam_offset_manual = 3.5;
npl_t = h5read("..\candidates_data.h5","/candidates/Delay_Time,_\tau_(s)_[NPL]/18");
npl_g2 = h5read("..\candidates_data.h5","/candidates/g^{(2)}(t)_[NPL]/18");
uom_t = h5read("..\candidates_data.h5","/candidates/Delay_Time,_\tau_(s)/18");
uom_g2 = h5read("..\candidates_data.h5","/candidates/g^{(2)}(t)/18");

f1 = g2fit(npl_t,npl_g2);
Y1 = f1(npl_t);
f2 = g2fit(uom_t,uom_g2);
Y2 = f1(uom_t);
plot(1e9.*npl_t-benyam_offset_manual,Y1,"r-","LineWidth",4)
hold on
plot(1e9.*(uom_t-f2.e),Y2.',"k-","LineWidth",2)
plot(1e9.*npl_t-benyam_offset_manual,npl_g2,"ro","MarkerSize",3.5,'MarkerFaceColor',[1 1 1],'LineWidth',1.2)
plot(1e9.*(uom_t-f2.e),uom_g2,"ko","MarkerSize",3.5,'MarkerFaceColor',[1 1 1],'LineWidth',1.2)
hold off
figform
axis square
xlim([-100 100])
ylim([0 1.5])
yticks([0, 0.5, 1, 1.5]);yticklabels(string([0, 0.5, 1, 1.5]));
xlabel("Delay time, t (ns)"); ylabel('g^{(2)}(t)');
set(gca,'FontSize',20)
title("(c)", 'Units', 'normalized', 'Position', title_offset)
    %inset
    inax = axes('Position',[.8225 .22 .080 .27],'TickDir','in','YTick',[0,0.5,1],'YTickLabel',{'0','0.5','1'});
    plot(1e9.*npl_t-benyam_offset_manual,Y1.',"r-","LineWidth",3)
    hold on
    plot(1e9.*(uom_t-f2.e),Y2,'k-','LineWidth',2)
    hold off
    box on
    figform
    axis square
    xlim([-1000 1000])
    ylim([0 1.4])

f.Position = [anchor_h anchor_v 1600 600];
set(findobj(gcf,'type','axes'),"FontSize",18,"FontWeight","bold")
set(inax,'FontSize',14,'TickDir','in','FontWeight','bold')

%% figure(2) candidates parameters
%These are selected due to 3600s acquisition time for TTTR mode
selected_NDs = [...
    "S2_B",5,8;...
    "S2_B",18,17;...
    "S2_Ba",22,15;...
    "S2_Cd",2,6;...
    "S3_Ad",3,21;...
    "S3_D",2,22;...
    ];

NDs_table = table('Size',[height(h5read("..\candidates_data.h5","/candidates/manuscript_code")),23],...
    'VariableTypes',["string","string","double"...
    "double","double","double","double","double","double",...
    "double","double","double","double","double","double",...
    "double","double","double","double",...
    "double","double","double","double"]);

NDs_table.Properties.VariableNames = ["manuscript code","NPL code","Acquisition Power, P (µW)","g^{(2)}(0)","g^{(2)}(0) [NPL]","g^{(2)}(0) [nErr]","g^{(2)}(0) [nErr] [NPL]",...
    "I_{80} (cps)","I_{80} (cps) [NPL]","I_{80} (cps) [nErr]","I_{80} (cps) [nErr] [NPL]","I_{80} (cps) [pErr]","I_{80} (cps) [pErr] [NPL]",...
    "P_{sat} (µW)","P_{sat} (µW) [NPL]","P_{sat} (µW) [nErr]","P_{sat} (µW) [nErr] [NPL]","P_{sat} (µW) [pErr]","P_{sat} (µW) [pErr] [NPL]",...
    "k_{\infty}","k_{\infty} [NPL]","k_{\infty} [nErr]","k_{\infty} [nErr] [NPL]"];

%populate table
NDs_table.("manuscript code") = h5read("..\candidates_data.h5","/candidates/manuscript_code");
NDs_table.("NPL code") = h5read("..\candidates_data.h5","/candidates/NPL_code");
NDs_table.("Acquisition Power, P (µW)") = h5read("..\candidates_data.h5","/candidates/Acquisition_Power,_P_(µW)");

NDs_table.("g^{(2)}(0)") = h5read("..\candidates_data.h5","/candidates/g^{(2)}(0)");
NDs_table.("g^{(2)}(0) [NPL]") = h5read("..\candidates_data.h5","/candidates/g^{(2)}(0)_[NPL]");
NDs_table.("g^{(2)}(0) [nErr]") = h5read("..\candidates_data.h5","/candidates/g^{(2)}(0)_[nErr]");
NDs_table.("g^{(2)}(0) [nErr] [NPL]") = h5read("..\candidates_data.h5","/candidates/g^{(2)}(0)_[nErr]_[NPL]");

NDs_table.("I_{80} (cps)") = h5read("..\candidates_data.h5","/candidates/I_{80}_(cps)");
NDs_table.("I_{80} (cps) [NPL]") = h5read("..\candidates_data.h5","/candidates/I_{80}_(cps)_[NPL]");
NDs_table.("I_{80} (cps) [nErr]") = h5read("..\candidates_data.h5","/candidates/I_{80}_(cps)_[nErr]");
NDs_table.("I_{80} (cps) [nErr] [NPL]") = h5read("..\candidates_data.h5","/candidates/I_{80}_(cps)_[pErr]_[NPL]");
NDs_table.("I_{80} (cps) [pErr]") = h5read("..\candidates_data.h5","/candidates/I_{80}_(cps)_[nErr]");
NDs_table.("I_{80} (cps) [pErr] [NPL]") = h5read("..\candidates_data.h5","/candidates/I_{80}_(cps)_[pErr]_[NPL]");

NDs_table.("P_{sat} (µW)") = h5read("..\candidates_data.h5","/candidates/P_{sat}_(µW)");
NDs_table.("P_{sat} (µW) [NPL]") = h5read("..\candidates_data.h5","/candidates/P_{sat}_(µW)_[NPL]");
NDs_table.("P_{sat} (µW) [nErr]") = h5read("..\candidates_data.h5","/candidates/P_{sat}_(µW)_[nErr]");
NDs_table.("P_{sat} (µW) [nErr] [NPL]") = h5read("..\candidates_data.h5","/candidates/P_{sat}_(µW)_[nErr]_[NPL]");
NDs_table.("P_{sat} (µW) [pErr]") = h5read("..\candidates_data.h5","/candidates/P_{sat}_(µW)_[pErr]");
NDs_table.("P_{sat} (µW) [pErr] [NPL]") = h5read("..\candidates_data.h5","/candidates/P_{sat}_(µW)_[pErr]_[NPL]");

NDs_table.("k_{\infty}") = h5read("..\candidates_data.h5","/candidates/k_{\infty}");
NDs_table.("k_{\infty} [NPL]") = h5read("..\candidates_data.h5","/candidates/k_{\infty}_[NPL]");
NDs_table.("k_{\infty} [nErr]") = h5read("..\candidates_data.h5","/candidates/k_{\infty}_[nErr]");
NDs_table.("k_{\infty} [nErr] [NPL]") = h5read("..\candidates_data.h5","/candidates/k_{\infty}_[nErr]_[NPL]");

correlation_table = NDs_table(str2double(selected_NDs(:,3)),:);

correlation_table = sortrows(correlation_table,"manuscript code","ascend");
NDs_List = unique(ND_name);
f = figure(2);clf
colororder({'r','k'}) %set so Y-axis left (red) and right (black)
tiledlayout("vertical","TileSpacing","none","Padding","compact");

nexttile; %g(2)(0)
t = bar(1:height(NDs_List),[correlation_table.("g^{(2)}(0) [NPL]"),correlation_table.("g^{(2)}(0)")],'EdgeColor','none','FaceAlpha',0.6,'BarWidth',1,'GroupWidth',0.7);
hold on
errorbar(t(1).XEndPoints,...
    correlation_table.("g^{(2)}(0) [NPL]"),correlation_table.("g^{(2)}(0) [nErr] [NPL]"),...
    'ro',"LineWidth",2,"MarkerSize",3,"MarkerFaceColor",[1 1 1]);
errorbar(t(2).XEndPoints,...
    correlation_table.("g^{(2)}(0)"),correlation_table.("g^{(2)}(0) [nErr]"),...
    'ko',"LineWidth",2,"MarkerSize",3,"MarkerFaceColor",[1 1 1]);
yline(0.5,'b--','LineWidth',2)
legend(["NPL","UoM"],'Location','northwest',FontSize=18)
ylabel("g^{(2)}(0)")
hold off
figform
xticks(1:height(NDs_List)); xticklabels([])
axis padded
yticks(0:0.5:1);yticklabels(string(0:0.5:1));ylim([0 1.1])

nexttile; %Isat
t = bar(1:height(NDs_List),[correlation_table.("I_{80} (cps) [NPL]"),correlation_table.("I_{80} (cps)")],'EdgeColor','none','FaceAlpha',0.6,'BarWidth',1,'GroupWidth',0.7);
hold on
errorbar(t(1).XEndPoints,...
    correlation_table.("I_{80} (cps) [NPL]"),correlation_table.("I_{80} (cps) [nErr] [NPL]"),...
    'ro',"LineWidth",2,"MarkerSize",3,"MarkerFaceColor",[1 1 1]);
errorbar(t(2).XEndPoints,...
    correlation_table.("I_{80} (cps)"),correlation_table.("I_{80} (cps) [nErr]"),...
    'ko',"LineWidth",2,"MarkerSize",3,"MarkerFaceColor",[1 1 1]);
ylabel("I_{80} (kcps)")
hold off
figform
xticks(1:height(NDs_List)); xticklabels([])
axis padded
yticks(0:25e3:100e3);yticklabels(string(0:25:100));ylim([0 99e3])

nexttile; %Psat
t = bar(1:height(NDs_List),[correlation_table.("P_{sat} (µW) [NPL]"),correlation_table.("P_{sat} (µW)")],'EdgeColor','none','FaceAlpha',0.6,'BarWidth',1,'GroupWidth',0.7);
hold on
errorbar(t(1).XEndPoints,...
    correlation_table.("P_{sat} (µW) [NPL]"),correlation_table.("P_{sat} (µW) [nErr] [NPL]"),...
    'ro',"LineWidth",2,"MarkerSize",3,"MarkerFaceColor",[1 1 1]);
errorbar(t(2).XEndPoints,...
    correlation_table.("P_{sat} (µW)"),correlation_table.("P_{sat} (µW) [nErr]"),...
    'ko',"LineWidth",2,"MarkerSize",3,"MarkerFaceColor",[1 1 1]);
ylabel("P_{sat} (µW)")
hold off
figform
xticks(1:height(NDs_List)); xticklabels([])
axis padded
yticks(0:200:600);yticklabels(string(0:200:600));ylim([0 1.2*max(correlation_table.("P_{sat} (µW) [NPL]"))])

nexttile; %k_inf
t = bar(1:height(NDs_List),[correlation_table.("k_{\infty} [NPL]"),correlation_table.("k_{\infty}")],'EdgeColor','none','FaceAlpha',0.6,'BarWidth',1,'GroupWidth',0.7);
hold on
errorbar(t(1).XEndPoints,...
    correlation_table.("k_{\infty} [NPL]"),correlation_table.("k_{\infty} [nErr] [NPL]"),...
    'ro',"LineWidth",2,"MarkerSize",3,"MarkerFaceColor",[1 1 1]);
errorbar(t(2).XEndPoints,...
    correlation_table.("k_{\infty}"),correlation_table.("k_{\infty} [nErr]"),...
    'ko',"LineWidth",2,"MarkerSize",3,"MarkerFaceColor",[1 1 1]);
ylabel("k_{\infty} (kcps)")
hold off
figform
xticks(1:height(NDs_List)); xticklabels(correlation_table.("manuscript code"))
axis padded
yticks(0:50e3:200e3);yticklabels(string(0:50:200));ylim([0 199e3])

set(findobj(gcf,'type','axes'),"FontSize",18,'XLim',[0.4 6.6])
f.Position = [anchor_h anchor_v 600 850];

%% figure(3) Ratios
tet_offset = [0.5, -0.3, 0];
selected_NDs = [...
    "S2_B",5,8;...
    %"S2_B",18,17;...
    "S2_Ba",22,15;...
    "S2_Cd",2,6;...
    "S3_Ad",3,21;...
    %"S3_D",2,22;... %this is the one mislabelled at NPL (mislabelled as S3_B ND#2)
    ];

correlation_table = NDs_table(ismember(NDs_table.("NPL code"),strcat(selected_NDs(:,1),selected_NDs(:,2))),:);
mask = (correlation_table.("g^{(2)}(0)") < 0.5 & correlation_table.("Acquisition Power, P (µW)") < 800);
sub_table = correlation_table(mask,:);

unq_NDs = unique(sub_table.("manuscript code"));
Iratio = zeros(numel(unq_NDs),3);
Pratio = zeros(numel(unq_NDs),3);
sgl_iratio = [];
sgl_pratio = [];
j = 1; %counter
for i = 1:numel(unq_NDs)
    septable = sub_table(sub_table.("manuscript code") == unq_NDs(i),:);
    Iratio(i,1) = mean(double(septable.("I_{80} (cps)"))./double(septable.("I_{80} (cps) [NPL]")));
    Iratio(i,2) = mean(double(septable.("I_{80} (cps) [nErr]"))./double(septable.("I_{80} (cps) [nErr] [NPL]")));
    Iratio(i,3) = mean(double(septable.("I_{80} (cps) [pErr]"))./double(septable.("I_{80} (cps) [pErr] [NPL]")));

    Pratio(i,1) = mean(septable.("P_{sat} (µW)")./septable.("P_{sat} (µW) [NPL]"));
    Pratio(i,2) = mean(double(septable.("P_{sat} (µW) [nErr]"))./double(septable.("P_{sat} (µW) [nErr] [NPL]")));
    Pratio(i,3) = mean(double(septable.("P_{sat} (µW) [pErr]"))./double(septable.("P_{sat} (µW) [pErr] [NPL]")));

    sgl_iratio = [sgl_iratio,Iratio(i,1)];
    sgl_pratio = [sgl_pratio,Pratio(i,1)];
end

f = figure(3);clf
tiledlayout("vertical","TileSpacing","none","Padding","compact")

nexttile;
disp(strcat("Mean I(sat) ratio is ",num2str(mean(sgl_iratio)),"±",num2str(std(sgl_iratio))));
yline(mean(Iratio(:,1)),'k--',"LineWidth",2)
hold on
yp = [mean(sgl_iratio)-std(sgl_iratio) mean(sgl_iratio)-std(sgl_iratio) mean(sgl_iratio)+std(sgl_iratio) mean(sgl_iratio)+std(sgl_iratio)];
xp = [-10 10 10 -10];
patch(xp,yp,[0 0 0], 'edgecolor', 'none', 'FaceAlpha', 0.1)
errorbar(1:numel(unq_NDs),Iratio(:,1),...
    std(Iratio(:,1)),std(Iratio(:,1)),...
    'ro',"LineWidth",2,...
    "MarkerSize",8,"MarkerFaceColor",[1 1 1])
text(0.62,0.43,strcat("I_{80} ratio = ",num2str(round(mean(sgl_iratio),2)),"(",num2str(round(100*std(sgl_iratio),0)),")"),...
    "FontSize",12,"FontWeight","bold")
hold off
xticks(1:numel(unq_NDs)); xticklabels([])
yticks(0:0.2:0.8);yticklabels(0:0.2:0.8)
ylabel("I_{80}^{UoM}/I_{80}^{NPL}")
figform
axis padded
xlim([0.5 numel(unq_NDs)+0.5])

nexttile;
yline(mean(Pratio(:,1)),'k--',"LineWidth",2)
hold on
yp = [mean(sgl_pratio)-std(sgl_pratio) mean(sgl_pratio)-std(sgl_pratio) mean(sgl_pratio)+std(sgl_pratio) mean(sgl_pratio)+std(sgl_pratio)];
xp = [-10 10 10 -10];
patch(xp,yp,[0 0 0], 'edgecolor', 'none', 'FaceAlpha', 0.1)
errorbar(1:numel(unq_NDs),Pratio(:,1),...
    Pratio(:,2),Pratio(:,3),...
    'ro',"LineWidth",2,...
    "MarkerSize",8,"MarkerFaceColor",[1 1 1])
text(0.62,0.4,strcat("P_{sat} ratio = ",num2str(round(mean(sgl_pratio),2)),"(",num2str(10*round(std(sgl_pratio),1)),")"),...
    "FontSize",12,"FontWeight","bold")
hold off
xticks(1:numel(unq_NDs)); xticklabels(unq_NDs)
yticks(0:0.2:0.8);yticklabels(0:0.2:0.8)
disp(strcat("Mean P(sat) ratio is ",num2str(mean(sgl_pratio)),"±",num2str(std(sgl_pratio))));
ylabel("P_{sat}^{UoM}/P_{sat}^{NPL}")
figform
axis padded
xlim([0.5 numel(unq_NDs)+0.5])

disp("=====================================")
set(findobj(gcf,'type','axes'),"FontSize",18)
f.Position = [anchor_h anchor_v 550 400];

%% FUNCTIONS
function [Psat,psat_full,psat_ND,psat_laser] = pdep_components(X,Y)
    [fitresult, ~] = psat_fit(X,Y);
    coeffs = [coeffvalues(fitresult); confint(fitresult,0.68)]';
    coeffs = [coeffs(:,1)-coeffs(:,2),coeffs(:,1),coeffs(:,3)-coeffs(:,1)]; % sorting in to lower,val,upper
    k_inf = coeffs(1,:);
    Psat = coeffs(2,:);

    psat_full = @(P) k_inf(2)*P/(P+Psat(2)) + coeffs(3,2)*P;
    psat_ND = @(P) k_inf(2)*P/(P+Psat(2));
    psat_laser = @(P) coeffs(3,2)*P;
end

function rptplot(tbl,col,mstyle)
unq_NDs = unique(tbl.("NPL code"));
for i = 1:numel(unique(tbl.("NPL code")))
    septable = tbl(tbl.("NPL code") == unq_NDs(i),:);
    hold on
    errorbar(1:max(septable.("Repeat No.")),double(septable.(col))/mean(double(septable.(col))),...
        double(septable.(strcat(col," [nErr]")))/mean(double(septable.(col))),...
        double(septable.(strcat(col," [pErr]")))/mean(double(septable.(col))),...
        mstyle(i),'MarkerFaceColor',[1 1 1],'MarkerSize',4,'LineWidth',1.5,'LineStyle','none');
end
xticks(1:max(tbl.("Repeat No.")));xticklabels([]) 
hold off
ylabel("g^{(2)}(0)")
figform
axis padded
yline(1,'k--')
%ylim([0.5 1.5])
end