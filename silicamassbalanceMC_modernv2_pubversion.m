%%
%Silicon isotope mass balance Monte Carlo model
%Elizabeth Trower, June 2021
%this code was designed with Matlab 2018b

%This model is based on the bones of the code from Trower & Fischer 2019;
%it was redesigned to explore Si isotope mass balance in the modern silica
%cycle by Lizzy Trower and Shaily Rahman.

%This is the publication version that pulls data from a static .xlsx file.

%%
%This section loads d30Si datasets from the spreadsheet, converts the data
%from strings to double precision values, and generates inverse cumulative
%distribution functions. It also sets up distributions for the sizes of
%each sink/source.

clear

n = 10000;

%bulk silicate Earth
BSE = -0.29;
BSE_1sd = 0.04;
pdfBSE = makedist('Normal','mu',BSE,'sigma',BSE_1sd);
z00 =rand(1,n);
icdfBSE = icdf(pdfBSE,z00);

%SINKS
%diatoms
diatomdata = xlsread('d30Si_data_compilation_05072021.xlsx',3,'A:A');
diatomseddata = xlsread('d30Si_data_compilation_05072021.xlsx',4,'A:A');
alldiatomdata = cat(1,diatomdata,diatomseddata);
pdfd30Sidiatom = fitdist(diatomdata,'Kernel');
pdfd30Sidiatomsed = fitdist(diatomseddata,'Kernel');
pdfd30Sialldiatom = fitdist(alldiatomdata,'Kernel');
z0 = rand(1,n);
icdfd30Sidiatom = icdf(pdfd30Sidiatom,z0);
icdfd30Sidiatomall = icdf(pdfd30Sidiatom,z0);

Fdiatomavg = 9.2;
Fdiatom1sd = 1.6;
pdfFdiatom = truncate(makedist('Normal','mu',Fdiatomavg,'sigma',Fdiatom1sd),...
    0,inf); %these are truncated distributions so there are no negative sink sizes
z0b = rand(1,n);
icdfFdiatom = icdf(pdfFdiatom,z0b);
%this currently attributes the entire F_B sink to diatoms

%radiolarians
raddata = xlsread('d30Si_data_compilation_05072021.xlsx',6,'A:A');
pdfd30Sirad = fitdist(raddata,'Kernel');
z1 = rand(1,n);
icdfd30Sirad = icdf(pdfd30Sirad,z1);

%radiolarian data is currently not included in the marine silica sinks due
%to lack of constraints on radiolarian flux (we therefore assume that the
%radiolarian contribution is negligible compared to the diatom flux).

%sponges
spongedata = xlsread('d30Si_data_compilation_05072021.xlsx',7,'A:A');
pdfd30Sisponge = fitdist(spongedata,'Kernel');
z2 = rand(1,n);
icdfd30Sisponge = icdf(pdfd30Sisponge,z2);

Fspongeavg = 1.7;
Fsponge1sd = 1.6;
pdfFsponge = truncate(makedist('Normal','mu',Fspongeavg,'sigma',Fsponge1sd),...
    0,inf);
z2b = rand(1,n);
icdfFsponge = icdf(pdfFsponge,z2b);

%reverse weathering
%there is no existing data for direct measurements of the d30Si values of
%reverse weathering products

FRWavg = 4.7;
FRW1sd = 2.3;
pdfFRW = truncate(makedist('Normal','mu',FRWavg,'sigma',FRW1sd),...
    0,inf);
z3b = rand(1,n);
icdfFRW = icdf(pdfFRW,z3b);

%SOURCES
%aeolian dust
%this section uses extra column data to choose either loess or soil data
%for d30Si of aeolian dust. it does not make a substantial difference to
%the overall balance.
dustdata = xlsread('d30Si_data_compilation_05072021.xlsx',9,'A:A');
[num,txt,raw] = xlsread('d30Si_data_compilation_05072021.xlsx',9,'F:F');
dustdatatype = categorical(txt(2:end));
loessdata = dustdata(dustdatatype == 'loess');
soildata = dustdata(dustdatatype == 'soil');
pdfd30Sidust = fitdist(loessdata,'Kernel');
%pdfd30Sidust = fitdist(soildata,'Kernel');
z4 = rand(1,n);
icdfd30Sidust = icdf(pdfd30Sidust,z4);

Fdustavg = 0.5;
Fdust1sd = 0.5;
pdfFdust = truncate(makedist('Normal','mu',Fdustavg,'sigma',Fdust1sd),...
    0,inf);
z4b = rand(1,n);
icdfFdust = icdf(pdfFdust,z4b);

%hydrothermal fluids
hydrothermaldata = xlsread('d30Si_data_compilation_05072021.xlsx',10,'A:A');
pdfd30Sihydrothermal = fitdist(hydrothermaldata,'Kernel');
z5 = rand(1,n);
icdfd30Sihydrothermal = icdf(pdfd30Sihydrothermal,z5);

Fhydrothermalavg = 1.7;
Fhydrothermal1sd = 0.9;
pdfFhydrothermal = truncate(makedist('Normal','mu',Fhydrothermalavg,'sigma',...
    Fhydrothermal1sd),0,inf);
z5b = rand(1,n);
icdfFhydrothermal = icdf(pdfFhydrothermal,z5b);

%rivers
riverdata = xlsread('d30Si_data_compilation_05072021.xlsx',11,'A:A');
[num2,txt2,raw2] = xlsread('d30Si_data_compilation_05072021.xlsx',11,'F:F');
riverdatatype = categorical(txt2(2:end));
[num3,txt3,raw3] = xlsread('d30Si_data_compilation_05072021.xlsx',11,'G:G');
riverCOSCATtype = categorical(txt3(2:end));
rivermouthdata = riverdata(riverdatatype == 'mouth');
pdfd30Sirivermouth = fitdist(rivermouthdata,'Kernel');
pdfd30Siriver = fitdist(riverdata,'Kernel');
z6 = rand(1,n);
icdfd30Siriver = icdf(pdfd30Siriver,z6);
icdfd30Sirivermouth = icdf(pdfd30Sirivermouth,z6);

Friveravg = 6.3;
Friver1sd = 0.4;
pdfFriver = truncate(makedist('Normal','mu',Friveravg,'sigma',Friver1sd),...
    0,inf);
z6b = rand(1,n);
icdfFriver = icdf(pdfFriver,z6b);

%rivers_aSi
riverasidata = xlsread('d30Si_data_compilation_05072021.xlsx',13,'A:A');
pdfd30Siriverasi = fitdist(riverasidata,'Kernel');
z7 = rand(1,n);
icdfd30Siriverasi = icdf(pdfd30Siriverasi,z7);

Friverasiavg = 1.9;
Friverasi1sd = 1;
pdfFriverasi = truncate(makedist('Normal','mu',Friverasiavg,'sigma',Friverasi1sd),...
    0,inf);
z7b = rand(1,n);
icdfFriverasi = icdf(pdfFriverasi,z7b);

%groundwater
GWdata = xlsread('d30Si_data_compilation_05072021.xlsx',12,'A:A');
GWsalinity = xlsread('d30Si_data_compilation_05072021.xlsx',12,'F:F');
GWfresh = GWdata(GWsalinity<4);
GWmarine = GWdata(GWsalinity>25);
pdfd30SiGWfresh = fitdist(GWfresh,'Kernel');
z8 = rand(1,n);
icdfd30SiGWfresh = icdf(pdfd30SiGWfresh,z8);
pdfd30SiGWmarine = fitdist(GWmarine,'Kernel');
z8a = rand(1,n);
icdfd30SiGWmarine = icdf(pdfd30SiGWmarine,z8a);

FGWmarineavg = 2.9;
FGWmarine1sd = 1.1;
pdfFGWmarine = truncate(makedist('Normal','mu',FGWmarineavg,'sigma',FGWmarine1sd),...
    0,inf);
z8b = rand(1,n);
icdfFGWmarine = icdf(pdfFGWmarine,z8b);

FGWfreshavg = 0.7;
FGWfresh1sd = 0.3;
pdfFGWfresh = truncate(makedist('Normal','mu',FGWfreshavg,'sigma',...
    FGWfresh1sd),0,inf);
z8c = rand(1,n);
icdfFGWfresh = icdf(pdfFGWfresh,z8c);


%low temperature silicate weathering
%currently parameterized with same d30Si dataset as river aSi
z9 = rand(1,n);
icdfd30SilowTweath = icdf(pdfd30Siriverasi,z9);


FlowTweathavg = 1.9;
FlowTweath1sd = 0.7;
pdfFlowTweath = truncate(makedist('Normal','mu',FlowTweathavg,'sigma',...
    FlowTweath1sd),0,inf);
z9b = rand(1,n);
icdfFlowTweath = icdf(pdfFlowTweath,z9b);

%ice sheet meltwaters
ISMWdata = xlsread('d30Si_data_compilation_05072021.xlsx',15,'A:A');
pdfd30SiISMW = fitdist(ISMWdata,'Kernel');
z10 = rand(1,n);
icdfd30SiISMW = icdf(pdfd30SiISMW,z10);

FISMWavg = 0.3;
FISMW1sd = 0.3;
pdfFISMW = truncate(makedist('Normal','mu',FISMWavg,'sigma',FISMW1sd),...
    0,inf);
z10b = rand(1,n);
icdfFISMW = icdf(pdfFISMW,z10b);

%%
%Generate plots of d30Si distributions of sinks + sources
f1 = figure;
edges = -6:0.5:3;
xvals_pdf = -6:0.1:3;
s11 = subplot(1,3,1);
yyaxis left
histogram(diatomdata,edges)
hold on
histogram(alldiatomdata,edges)
yyaxis right
yvals_diatom = pdf(pdfd30Sidiatom,xvals_pdf);
yvals_alldiatom = pdf(pdfd30Sialldiatom,xvals_pdf);
yvals_diatomsed = pdf(pdfd30Sidiatomsed,xvals_pdf);
hold on
plot(xvals_pdf,yvals_diatom,'LineWidth',1.5)
plot(xvals_pdf,yvals_alldiatom,'LineWidth',1.5)
plot(xvals_pdf,yvals_diatomsed,'LineWidth',1.5)
xline(-0.35,'g','LineWidth',1.5)
xlabel('\delta^3^0Si')
title('diatoms')
xlim([-6 3])

s12 = subplot(1,3,2);
yyaxis left
histogram(raddata,edges)
yyaxis right
yvals_rad = pdf(pdfd30Sirad,xvals_pdf);
hold on
plot(xvals_pdf,yvals_rad,'LineWidth',1.5)
xline(-0.35,'g','LineWidth',1.5)
xlabel('\delta^3^0Si')
title('radiolarians')
xlim([-6 3])

s13 = subplot(1,3,3);
yyaxis left
histogram(spongedata,edges)
yyaxis right
yvals_sponge = pdf(pdfd30Sisponge,xvals_pdf);
hold on
plot(xvals_pdf,yvals_sponge,'LineWidth',1.5)
xline(-0.35,'g','LineWidth',1.5)
xlabel('\delta^3^0Si')
title('sponges')
xlim([-6 3])

f2 = figure;
s21 = subplot(3,3,1);
yyaxis left
histogram(loessdata,edges)
yyaxis right
yvals_dust = pdf(pdfd30Sidust,xvals_pdf);
hold on
plot(xvals_pdf,yvals_dust,'LineWidth',1.5)
xline(-0.35,'g','LineWidth',1.5)
xlabel('\delta^3^0Si')
title('dust')
xlim([-6 3])

s22 = subplot(3,3,2);
yyaxis left
histogram(hydrothermaldata,edges)
yyaxis right
yvals_hydrothermal = pdf(pdfd30Sihydrothermal,xvals_pdf);
hold on
plot(xvals_pdf,yvals_hydrothermal,'LineWidth',1.5)
xline(-0.35,'g','LineWidth',1.5)
xlabel('\delta^3^0Si')
title('hydrothermal')
xlim([-6 3])

s23 = subplot(3,3,3);
yyaxis left
histogram(riverdata,edges)
yyaxis right
yvals_river = pdf(pdfd30Siriver,xvals_pdf);
hold on
plot(xvals_pdf,yvals_river,'LineWidth',1.5)
xline(-0.35,'g','LineWidth',1.5)
xlabel('\delta^3^0Si')
title('river DSi')
xlim([-6 3])

s24 = subplot(3,3,4);
yyaxis left
histogram(riverasidata,edges)
yyaxis right
yvals_riverasi = pdf(pdfd30Siriverasi,xvals_pdf);
hold on
plot(xvals_pdf,yvals_riverasi,'LineWidth',1.5)
xline(-0.35,'g','LineWidth',1.5)
xlabel('\delta^3^0Si')
title('river aSi')
xlim([-6 3])

s25 = subplot(3,3,5);
yyaxis left
histogram(GWfresh,edges)
yyaxis right
yvals_GWfresh = pdf(pdfd30SiGWfresh,xvals_pdf);
hold on
plot(xvals_pdf,yvals_GWfresh,'LineWidth',1.5)
xline(-0.35,'g','LineWidth',1.5)
xlabel('\delta^3^0Si')
title('GW_f_r_e_s_h')
xlim([-6 3])

s26 = subplot(3,3,6);
yyaxis left
histogram(GWmarine,edges)
yyaxis right
yvals_GWmarine = pdf(pdfd30SiGWmarine,xvals_pdf);
hold on
plot(xvals_pdf,yvals_GWmarine,'LineWidth',1.5)
xline(-0.35,'g','LineWidth',1.5)
xlabel('\delta^3^0Si')
title('GW_m_a_r_i_n_e')
xlim([-6 3])

s27 = subplot(3,3,7);
yyaxis left
histogram(riverasidata,edges)
yyaxis right
yvals_riverasi = pdf(pdfd30Siriverasi,xvals_pdf);
hold on
plot(xvals_pdf,yvals_riverasi,'LineWidth',1.5)
xline(-0.35,'g','LineWidth',1.5)
xlabel('\delta^3^0Si')
title('low T weath')
xlim([-6 3])

s28 = subplot(3,3,8);
yyaxis left
histogram(ISMWdata,edges)
yyaxis right
yvals_ISMW = pdf(pdfd30SiISMW,xvals_pdf);
hold on
plot(xvals_pdf,yvals_ISMW,'LineWidth',1.5)
xline(-0.35,'g','LineWidth',1.5)
xlabel('\delta^3^0Si')
title('ISMW')
xlim([-6 3])

%%
%This section generates plots to visualize the river d30Si dataset,
%including a comparison of the full (unweighted) dataset, river mouth only
%dataset, and a flux-weighted estimate.

%Generate initial figure with all river data + river mouth only data
figure
xvals_pdf = -6:0.1:3;
yvals_river = pdf(pdfd30Siriver,xvals_pdf);
plot(xvals_pdf,yvals_river,'LineWidth',1.5)
hold on
yvals_rivermouth = pdf(pdfd30Sirivermouth,xvals_pdf);
plot(xvals_pdf,yvals_rivermouth,'LineWidth',1.5)
xlabel('\delta^3^0Si')
title('river DSi')
xlim([-1 3])

%Generate a flux-weighted estimate for river DSi d30Si weighting by COSCAT
%(see Durr et al., 2011).
COSCATs = categories(riverCOSCATtype);
nCOSCATs = length(COSCATs);
icdfd30Si_COSCATs = zeros(nCOSCATs,n);

%Sort river data by COSCAT and generate a range of d30Si estimates
%consistent with each COSCAT.
for m = 1:nCOSCATs
    cat = COSCATs(m);
    catriverdata = riverdata(riverCOSCATtype == cat);
    pdfd30Sicatriver = fitdist(catriverdata,'Kernel');
    zcat = rand(1,n);
    icdfd30Si_COSCATs(m,:) = icdf(pdfd30Sicatriver,zcat);
end

%Grab DSi flux for each COSCAT from spreadsheet (data from Durr et al.,
%2011).
DSi_flux_COSCATs = xlsread('COSCAT DSi flux data.xlsx',1,'B:B');
DSi_flux_sum = sum(DSi_flux_COSCATs);

%Generate flux-weighted river DSi d30Si estimates.
icdfd30Si_river_fluxweight = sum((DSi_flux_COSCATs./DSi_flux_sum)'*...
    icdfd30Si_COSCATs,1);

%Calculate pdf of flux-weighted river DSi d30Si estimates and add to
%figure.
pdfd30Si_river_fluxweight = fitdist(icdfd30Si_river_fluxweight','Kernel');
yvals_river_fluxweight = pdf(pdfd30Si_river_fluxweight,xvals_pdf);
plot(xvals_pdf,yvals_river_fluxweight,'LineWidth',1.5)

%Add subglacial meltwaters to figure.
plot(xvals_pdf,yvals_ISMW,'LineWidth',1.5)
legend ('all river data','river mouth data','flux-weighted river data',...
    'ISMW','Location','northwest')

%%
%Sink isotope mass balance
%This section runs the mass balance calculation for the silica sinks
%and output a distribution of values describing the offset from BSE.

F_sinks = icdfFdiatom + icdfFsponge;
d30Si_norm_sinks = (icdfFdiatom.*icdfd30Sidiatom + icdfFsponge.*...
    icdfd30Sisponge)./F_sinks;
figure
subplot(3,2,1)
histogram(d30Si_norm_sinks,-4:.1:4)
xline(-0.35,'g','LineWidth',1.5)
title('normalized \delta^3^0Si, sinks, excluding RW')
xlabel('\delta^3^0Si')
ylabel('count')
xlim([-4 4])

subplot(3,2,2)
histogram(icdfFdiatom,0:.5:30)
hold on
histogram(icdfFsponge,0:.5:30)
histogram(icdfFRW,0:.5:30)
histogram(icdfFRW + F_sinks,0:.5:30)
xlabel('flux (Tmol)')
ylabel('count')
title('sink total fluxes')
legend('diatom', 'sponge', 'RW', 'total BSi + RW')
xlim([0 30])

F_sources = icdfFdust + icdfFhydrothermal + icdfFriver + icdfFGWmarine + ...
    icdfFISMW + icdfFriverasi + icdfFlowTweath + icdfFGWfresh;
d30Si_norm_sources = (icdfFdust.*icdfd30Sidust + icdfFhydrothermal.*...
    icdfd30Sihydrothermal + icdfFriver.*icdfd30Siriver + icdfFGWmarine.*...
    icdfd30SiGWmarine + icdfFISMW.*icdfd30SiISMW + icdfFriverasi.*...
    icdfd30Siriverasi + icdfFlowTweath.*icdfd30SilowTweath +...
    icdfFGWfresh.*icdfd30SiGWfresh)./F_sources;

subplot(3,2,3)
histogram(d30Si_norm_sources,-4:.1:4)
xline(-0.35,'g','LineWidth',1.5)
title('normalized \delta^3^0Si, sources')
xlabel('\delta^3^0Si')
ylabel('count')
xlim([-4 4])

subplot(3,2,4)
% histogram(icdfFriver,0:.5:30)
% hold on
histogram(F_sources,0:.5:30)
xlabel('flux (Tmol)')
ylabel('count')
title('source total fluxes')
xlim([0 30])
%legend('river DSi','total')

d30SiRW_est = (d30Si_norm_sources.*(F_sinks + icdfFRW) - ...
    (icdfFdiatom.*icdfd30Sidiatom +icdfFsponge.*icdfd30Sisponge))...
    ./icdfFRW;
subplot(3,2,5)
yyaxis left
hRW = histogram(d30SiRW_est,-4:.1:4);
hRW.FaceColor = [0 0.5 0.5];
hold on
histogram(d30Si_norm_sinks,-4:.1:4)
title('\delta^3^0Si_R_W estimates')
hold on
yyaxis right
xvals = -2:.01:3.5;
yvals = pdf(pdfd30Sidiatom,xvals);
plot(xvals,yvals)
xline(-0.35,'g','LineWidth',1.5)
legend('model est','BSi sinks','diatom data','BSE','Location','northwest')
xlim([-4 4])
xlabel('\delta^3^0Si')

subplot(3,2,6)
histogram(F_sources-(F_sinks + icdfFRW),-15:.5:15)
xlabel('flux (Tmol)')
ylabel('count')
title('flux imbalance, sources - sinks')
xline(0,'k','LineWidth',1.5)
xlim([-15 15])

%%
%This section will run a series of sensitivity tests by changing key parts
%of the mass balance and generate a plot comparing the model output, i.e.
%predictions of RW d30Si.

%First plot the base case for visual comparison
figure
subplot(3,2,1)
h11 = histogram(d30SiRW_est,-4:.1:4);
ylabel('count')
h11.FaceColor = [0 0.5 0.5];
xlabel('\delta^3^0Si_R_W')
title('base model')
ylim([0 400])
xlim([-4 4])

%This will calculate d30Si_RW by first estimate F_RW assuming mass balance
FRW_mbest = F_sources - F_sinks;
FRW_mbest(FRW_mbest<0) = NaN;
d30SiRW_mbest = (d30Si_norm_sources.*(F_sinks + FRW_mbest) - ...
    (icdfFdiatom.*icdfd30Sidiatom +icdfFsponge.*icdfd30Sisponge))...
    ./FRW_mbest;

subplot(3,2,2)
h21 = histogram(d30SiRW_mbest,-4:.1:4);
hold on
h22 = histogram(d30SiRW_est,-4:.1:4);
ylabel('count')
h22.FaceColor = [0 0.5 0.5];
xlabel('\delta^3^0Si_R_W')
title ('require flux mass balance')
ylim([0 400])
xlim([-4 4])

%Alternative mass balance w/ smaller diatom sink (7 Tmol)
%currently forces flux mass balance, but can swap out for commented-out
%text to not use this assumption.

Fdiatomavg_alt = 7;
Fdiatom1sd_alt = 3.4;
pdfFdiatom_alt = truncate(makedist('Normal','mu',Fdiatomavg_alt,'sigma',Fdiatom1sd_alt),...
    0,inf); %these are truncated distributions so there are no negative sink sizes
z0b_alt = rand(1,n);
icdfFdiatom_alt = icdf(pdfFdiatom_alt,z0b_alt);

F_sinks_alt = icdfFdiatom_alt + icdfFsponge;
%the following lines force flux mass balance
FRW_mbest_alt = F_sources - F_sinks_alt;
FRW_mbest_alt(FRW_mbest_alt<0) = NaN;
d30SiRW_est_alt = (d30Si_norm_sources.*(F_sinks_alt + FRW_mbest_alt) - ...
    (icdfFdiatom_alt.*icdfd30Sidiatom +icdfFsponge.*icdfd30Sisponge))...
    ./FRW_mbest_alt;
%the previous lines could be substituted for the following equation to not
%force flux mass balance
% d30SiRW_est_alt = (d30Si_norm_sources.*(F_sinks_alt + icdfFRW) - ...
%     (icdfFdiatom_alt.*icdfd30Sidiatom +icdfFsponge.*icdfd30Sisponge))...
%     ./icdfFRW;

subplot(3,2,3)
h31 = histogram(d30SiRW_est_alt,-4:.1:4);
hold on
h32 = histogram(d30SiRW_est,-4:.1:4);
ylabel('count')
h32.FaceColor = [0 0.5 0.5];
title('F_d_i_a_t_o_m = 7 Tmol')
xlabel('\delta^3^0Si_R_W')
ylim([0 400])
xlim([-4 4])

%Alternative mass balance w/ larger diatom sink (14 Tmol)
%currently forces flux mass balance, but can swap out for commented-out
%text to not use this assumption.

Fdiatomavg_alt2 = 14;
Fdiatom1sd_alt2 = 3.4;
pdfFdiatom_alt2 = truncate(makedist('Normal','mu',Fdiatomavg_alt2,'sigma',Fdiatom1sd_alt2),...
    0,inf); %these are truncated distributions so there are no negative sink sizes
z0b_alt2 = rand(1,n);
icdfFdiatom_alt2 = icdf(pdfFdiatom_alt2,z0b_alt2);

F_sinks_alt2 = icdfFdiatom_alt2 + icdfFsponge;
%the following lines force flux mass balance
FRW_mbest_alt2 = F_sources - F_sinks_alt2;
FRW_mbest_alt2(FRW_mbest_alt2<0) = NaN;
d30SiRW_est_alt2 = (d30Si_norm_sources.*(F_sinks_alt2 + FRW_mbest_alt2) - ...
    (icdfFdiatom_alt2.*icdfd30Sidiatom +icdfFsponge.*icdfd30Sisponge))...
    ./FRW_mbest_alt2;
%the previous lines could be substituted for the following equation to not
%force flux mass balance
% d30SiRW_est_alt2 = (d30Si_norm_sources.*(F_sinks_alt2 + icdfFRW) - ...
%     (icdfFdiatom_alt2.*icdfd30Sidiatom +icdfFsponge.*icdfd30Sisponge))...
%     ./icdfFRW;

subplot(3,2,4)
h41 = histogram(d30SiRW_est_alt2,-4:.1:4);
hold on
h42 = histogram(d30SiRW_est,-4:.1:4);
ylabel('count')
h42.FaceColor = [0 0.5 0.5];
title('F_d_i_a_t_o_m = 14 Tmol')
xlabel('\delta^3^0Si_R_W')
ylim([0 400])
xlim([-4 4])

%Use only river mouth data for river DSi d30Si distribution.

d30Si_norm_sources_alt = (icdfFdust.*icdfd30Sidust + icdfFhydrothermal.*...
    icdfd30Sihydrothermal + icdfFriver.*icdfd30Sirivermouth + icdfFGWmarine.*...
    icdfd30SiGWmarine + icdfFISMW.*icdfd30SiISMW + icdfFriverasi.*...
    icdfd30Siriverasi + icdfFlowTweath.*icdfd30SilowTweath +...
    icdfFGWfresh.*icdfd30SiGWfresh)./F_sources;
d30SiRW_est_rivermouth = (d30Si_norm_sources_alt.*(F_sinks + icdfFRW) - ...
    (icdfFdiatom.*icdfd30Sidiatom +icdfFsponge.*icdfd30Sisponge))...
    ./icdfFRW;

subplot(3,2,5)
h51 = histogram(d30SiRW_est_rivermouth,-4:.1:4);
hold on
h52 = histogram(d30SiRW_est,-4:.1:4);
ylabel('count')
h52.FaceColor = [0 0.5 0.5];
title('river mouth DSi only')
xlabel('\delta^3^0Si_R_W')
ylim([0 400])
xlim([-4 4])

%Use flux-weighted river data for river DSi d30Si distribution.

%calculate d30SiRW_est_river_fluxweight
d30Si_norm_sources_alt2 = (icdfFdust.*icdfd30Sidust + icdfFhydrothermal.*...
    icdfd30Sihydrothermal + icdfFriver.*icdfd30Si_river_fluxweight + icdfFGWmarine.*...
    icdfd30SiGWmarine + icdfFISMW.*icdfd30SiISMW + icdfFriverasi.*...
    icdfd30Siriverasi + icdfFlowTweath.*icdfd30SilowTweath +...
    icdfFGWfresh.*icdfd30SiGWfresh)./F_sources;
d30SiRW_est_river_fluxweight = (d30Si_norm_sources_alt2.*(F_sinks + icdfFRW) - ...
    (icdfFdiatom.*icdfd30Sidiatom +icdfFsponge.*icdfd30Sisponge))...
    ./icdfFRW;

subplot(3,2,6)
h61 = histogram(d30SiRW_est_river_fluxweight,-4:.1:4);
hold on
h62 = histogram(d30SiRW_est,-4:.1:4);
ylabel('count')
h62.FaceColor = [0 0.5 0.5];
title('flux-weighted river DSi')
xlabel('\delta^3^0Si_R_W')
ylim([0 400])
xlim([-4 4])

%%
%This section uses the estimates of the flux-weighted marine silica sinks
%(diatom + sponges + reverse weathering) and d30Si datasets of continental
%silica sinks (phytoliths, pedogenic clays, freshwater diatoms) to invert
%for the size of the continental silica sink in order for the combination
%of marine + continental silica sinks to be in mass balance with BSE.

%Calculate size and flux-weighted d30Si distribution for marine sinks.
F_sinks_marine = icdfFdiatom + icdfFsponge + icdfFRW;
d30Si_norm_sinks_marine = (icdfFdiatom.*icdfd30Sidiatom + icdfFsponge.*...
    icdfd30Sisponge + icdfFRW.*d30SiRW_est)./F_sinks_marine;

%Load pedogenic clay and generate pdf
pedclaydata = xlsread('d30Si_data_compilation_05072021.xlsx',16,'A:A');
pdfpedclay = fitdist(pedclaydata,'Kernel');
zclay = rand(1,n);
icdfClay = icdf(pdfpedclay,zclay);

%Generate alternate estimate of true pedogenic clay d30Si based on a mean
%d30Si value of -2 permil
d30Sipedclay_rev = -2;
d30Sipedclay_rev_1sd = .5;
d30Sipdfpedclay_rev = makedist('Normal','mu',d30Sipedclay_rev,'sigma',...
    d30Sipedclay_rev_1sd);
zpedclay_rev =rand(1,n);
icdfd30Sipedclay_rev = icdf(d30Sipdfpedclay_rev,zpedclay_rev);

%Plot histogram of pedogenic clay data, which is interpreted here as a
%mixture of true pedogenic clays + clay-sized primary aluminosilicates
%based on sampling techniques, and alternative "unmixed" estimate of
%pedogenic clay d30Si centered on -2 permil.
figure
yyaxis left
histogram(pedclaydata,edges)
ylabel('count')
hold on
yyaxis right
yvals_clay = pdf(d30Sipdfpedclay_rev,-6:.1:3);
plot(-6:.1:3,yvals_clay)
ylabel('probability density')
xline(-0.35,'g','LineWidth',1.5);
xlabel('\delta^3^0Si')
xlim([-6 3])
legend('clay data','proposed unmixed fit')

%Grab phytolith data from spreadsheet and calculate pdf.
phytdata = xlsread('d30Si_data_compilation_05072021.xlsx',17,'A:A');
pdfd30Siphyt_data = fitdist(phytdata,'Kernel');
zphyt2 = rand(1,n);
icdfd30Siphyt_data = icdf(pdfd30Siphyt_data,zphyt2);

%Plot phyolith data and pdf.
figure
yyaxis left
histogram(phytdata,edges)
ylabel('count')
hold on
yyaxis right
yvals_phytdata = pdf(pdfd30Siphyt_data,-6:.1:3);
plot(-6:.1:3,yvals_phytdata)
ylabel('probability density')
xline(-0.35,'g','LineWidth',1.5);
xlabel('\delta^3^0Si')
xlim([-6 3])
legend('phytolith data','kernel pdf estimate')

%Grad freshwater diatom data from spreadsheet and calculate pdf
freshdiatomdata = xlsread('d30Si_data_compilation_05072021.xlsx',18,'A:A');
pdfd30Sifreshdiatom_data = fitdist(freshdiatomdata,'Kernel');
zfreshdiatom = rand(1,n);
icdfd30Sifreshdiatom_data = icdf(pdfd30Sifreshdiatom_data,zfreshdiatom);

%Plot freshwater diatom data and pdf.
figure
yyaxis left
histogram(freshdiatomdata,edges)
ylabel('count')
hold on
yyaxis right
yvals_freshdiatom = pdf(pdfd30Sifreshdiatom_data,-1:.1:3.5);
plot(-1:.1:3.5,yvals_freshdiatom)
ylabel('probability density')
xline(-0.35,'g','LineWidth',1.5);
xlabel('\delta^3^0Si')
xlim([-6 3])
legend('freshwater diatom data','kernel pdf estimate')

%Run continental sink mixure modeling with three sinks - this section
%basically explores all possible mixtures of the three continental sinks
%to generate mixture-weighted d30Si values and then inverts for the total
%continental sink size required for each mixture to be in mass balance with
%BSE + marine silica sinks.

%Generate vectors for each of the mixtures.
%xx = clay
%yy = phytoliths
%zz = freshwater diatoms

N = 101;
x = linspace(0, 1, N);
y = x;

[xx, yy] = meshgrid(x, y);
zz = 1 - (xx + yy);
zz(zz<0) = NaN;
sums = xx + yy + zz;

xx_vector = xx(:);
yy_vector = yy(:);
zz_vector = zz(:);

xx_vector(isnan(zz_vector) == 1) = [];
yy_vector(isnan(zz_vector) == 1) = [];
zz_vector(isnan(zz_vector) == 1) = [];

%Preallocate space for for-loop.
F_cont_mix3 = zeros(n,length(zz_vector));
F_cont_mix3_medians = zeros(length(zz_vector),1);
F_cont_mix3_NaNs = zeros(length(zz_vector),1);

%Run mass balance calculation for each mixture.
for counter21 = 1:length(zz_vector)
    d30Sicontmix3 = xx_vector(counter21)*icdfd30Sipedclay_rev + ...
        yy_vector(counter21)*icdfd30Siphyt_data + ...
        zz_vector(counter21)*icdfd30Sifreshdiatom_data;
    F_cont_mix3(:,counter21) = F_sinks_marine.*(d30Si_norm_sinks_marine-icdfBSE)./...
        (icdfBSE-d30Sicontmix3);
    F_cont_mix3(F_cont_mix3(:,counter21)<0,counter21) = NaN;
    F_cont_mix3_medians(counter21) = median(F_cont_mix3(isfinite(F_cont_mix3(:,counter21)),counter21));
    F_cont_mix3_NaNs(counter21) = sum(isnan(F_cont_mix3(:,counter21)))/n*100;
end

%Generate figure to visualize continental silica sink estimates using
%ternary plots.

figure
%Plot estimates of total continental silica sink size required for mass
%balance for each mixture.
subplot(3,2,1)
ternpcolor(xx_vector,yy_vector,F_cont_mix3_medians)
ternlabel('clay','phytolith','freshwater diatom')
caxis([0 80])
title('F_c_o_n_t (Tmol/yr)')

%Plot contours of estimates of total continental silica sink size.
subplot(3,2,2)
terncontour(xx_vector,yy_vector,zz_vector,F_cont_mix3_medians,5:5:60)

%Plot percentage of NaNs for each mixture, which illustrates that some
%mixtures are significantly less plausible than others.
subplot(3,2,3)
ternpcolor(xx_vector,yy_vector,F_cont_mix3_NaNs)
caxis([0 100])
title('% of NaNs')

%Plot contours of percentage of NaNs for each mixture.
subplot(3,2,4)
terncontour(xx_vector,yy_vector,zz_vector,F_cont_mix3_NaNs,[5 10 15 20 25])

%This section now uses the estimates of the size of the pedogenic clay
%silica sink (calculated as the corresponding fraction of the continental
%silica sink size estimate for each mixture) to determine a continental
%silica storage ratio, which relates to the average stoichiometry of
%silicate weathering reactions.

%Preallocate space for for-loop.
F_DSi_weathflux_mix3 = zeros(n,length(zz_vector));
cont_storage_ratio_mix3 = zeros(n,length(zz_vector));
cont_storage_ratio_mix3_medians = zeros(length(zz_vector),1);

%Run mass balance calculations.
for counter22 = 1:length(zz_vector)
   F_DSi_weathflux_mix3(:,counter22) = icdfFriver' + icdfFGWfresh' + icdfFriverasi' +...
       icdfFGWmarine' + icdfFISMW' + F_cont_mix3(:,counter22).*yy_vector(counter22) + ...
       F_cont_mix3(:,counter22).*zz_vector(counter22);
   cont_storage_ratio_mix3(:,counter22) = F_cont_mix3(:,counter22).*xx_vector(counter22)./...
       F_DSi_weathflux_mix3(:,counter22);
   cont_storage_ratio_mix3_medians(counter22) = median(cont_storage_ratio_mix3...
       (isfinite(cont_storage_ratio_mix3(:,counter22)),counter22));
end

%Plot continental storage ratios.
subplot(3,2,5)
ternpcolor(xx_vector,yy_vector,cont_storage_ratio_mix3_medians)
caxis([0 0.8])
title('clay:DSi ratio')

%Plot contours of continental storage ratios.
subplot(3,2,6)
terncontour(xx_vector,yy_vector,zz_vector,cont_storage_ratio_mix3_medians,...
    0.05:0.05:0.75)

%%
%This section uses seawater data to estimate d30Si diatom to assess whether 
%the diatom dataset is representative.

%redefine xvals_pdf in case this section is run out of sequence
xvals_pdf = -6:0.1:3;

%Grab seawater DSi d30Si data from spreadsheet and identify subsets of
%seawater dataset.
seawaterdata = xlsread('d30Si_data_compilation_05072021.xlsx',2,'A:A');
seawaterdepth = xlsread('d30Si_data_compilation_05072021.xlsx',2,'H:H');
seawater_shallow = seawaterdata(seawaterdepth<500);
[num4,txt4,raw4] = xlsread('d30Si_data_compilation_05072021.xlsx',2,'E:E');
seawaterlocation = categorical(txt4(2:end));
seawaterdataSO = seawaterdata(seawaterlocation == 'Southern Ocean');
seawaterdataSOshallow = seawaterdata(seawaterlocation == 'Southern Ocean'...
    & seawaterdepth<500);
pdfd30Siseawater = fitdist(seawaterdata,'Kernel');
pdfd30Siseawater_shallow = fitdist(seawater_shallow,'Kernel');

%Estimate diatom d30Si data assuming -1.1 permil offset from seawater DSi
%d30Si (i.e., no Rayleigh fractionation).
pdfd30Sidiatom_SWest = fitdist(seawaterdata-1.1,'Kernel');
pdfd30Sidiatom_SWest_shallow = fitdist(seawater_shallow-1.1,'Kernel');
pdfd30Sidiatom_SWestSO = fitdist(seawaterdataSO-1.1,'Kernel');
pdfd30Sidiatom_SWestSOshallow = fitdist(seawaterdataSOshallow-1.1,'Kernel');

%Generate datasets to be used in different versions of mass balance
%calculation.
zsea = rand(1,n);
icdfd30Siseawater = icdf(pdfd30Siseawater,zsea);
icdfd30Sidiatom_SWest = icdf(pdfd30Sidiatom_SWest,zsea);
icdfd30Sidiatom_SWest_shallow = icdf(pdfd30Sidiatom_SWest_shallow,zsea);
icdfd30Sidiatom_SWestSO = icdf(pdfd30Sidiatom_SWestSO,zsea);
icdfd30Sidiatom_SWestSOshallow = icdf(pdfd30Sidiatom_SWestSOshallow,zsea);

%Plot pdfs of observed diatom data vs. estimates based on seawater DSi
%d30Si.
figure
plotcolors = parula(5);
plot(xvals_pdf,yvals_diatom,'LineWidth',1.5)
hold on
plot(xvals_pdf,yvals_diatomsed,'LineWidth',1.5)
plot(xvals_pdf,yvals_alldiatom,'LineWidth',1.5)
yvals_diatom_SWest = pdf(pdfd30Sidiatom_SWest,xvals_pdf);
plot(xvals_pdf,yvals_diatom_SWest,'Color',plotcolors(1,:))
yvals_diatom_SWest_shallow = pdf(pdfd30Sidiatom_SWest_shallow,xvals_pdf);
plot(xvals_pdf,yvals_diatom_SWest_shallow,'Color',plotcolors(2,:))
yvals_diatom_SWestSO = pdf(pdfd30Sidiatom_SWestSO,xvals_pdf);
plot(xvals_pdf,yvals_diatom_SWestSO,'Color',plotcolors(3,:))
yvals_diatom_SWestSOshallow = pdf(pdfd30Sidiatom_SWestSOshallow,xvals_pdf);
plot(xvals_pdf,yvals_diatom_SWestSOshallow,'Color',plotcolors(4,:))
legend('water column diatom data (n = 240)','sediment diatom data (n = 24)',...
    'all diatom data (n = 264)',...
    'SW estimate, all data (n = 1261)','SW estimate, <500m (n = 579)',...
    'SW estimate, Southern Ocean (n = 240)',...
    'SW estimate, Southern Ocean < 500m (n = 105)',...
    'location','northwest')
title('diatom data vs. SW estimates')
xlim([-2 3])
