% this file plots the figures
clear all;
currdir = pwd;
restoredefaultpath;
% replicate Figure 1 in KM (2012)

fixed_rf=0;
label_R = 'cmy';
prior_type = 'minnesota';
ndraws  =1e4;
prior_only=0;
save_every=10;
spec= ['rgibbs_',label_R,'prior_only_',num2str(prior_only),'prior_',prior_type,'fixed_rf_',num2str(fixed_rf),'_ndraws',num2str(ndraws),'number_of_shocks8','save_every',num2str(save_every)];
%spec= ['rgibbs_',label_R,'prior_only_',num2str(prior_only),'prior_',prior_type,'fixed_rf_',num2str(fixed_rf),'_ndraws',num2str(ndraws),'number_of_shocks8','_thinning',num2str(save_every)];

%L =load('matfiles/rgibbs_cmyprior_only_0prior_minnesotafixed_rf_0_ndraws100.mat')

L =load(['matfiles/',filesep,spec,'.mat']);




nhorizons = 21;%size(L.L,1);2
nvar      = size(L.L,2);
nshocks   = size(L.L,3);



L_q50=nan(nhorizons,nvar,nshocks); % store conjF quantile 50th
L_q16=nan(nhorizons,nvar,nshocks); % store conjF quantile 16th
L_q84=nan(nhorizons,nvar,nshocks); % store conjF quantile 84th


nburn=0;

for ii=1:nhorizons
    for jj=1:nvar
        for kk=1:nshocks
            


            L_q50(ii,jj,kk) = quantile(L.L(ii,jj,kk,nburn+1:end),0.5);
            L_q16(ii,jj,kk) = quantile(L.L(ii,jj,kk,nburn+1:end),0.16);
            L_q84(ii,jj,kk) = quantile(L.L(ii,jj,kk,nburn+1:end),0.84);

       

        end
    end
end




% figure settings
ftsizeaxis   = 9;
ftsizexlabel = 9;
ftsizetitle  = 9;
ftsizelegend = 9;
ftlinewidth  = 1.0;
medianwidth  = 1.0;


plot_uv_scale = 1;
horizon       = nhorizons-1;

if ~exist('pngfiles', 'dir')
    mkdir('pngfiles');
    disp('Created directory: pngfiles');
else
    disp('Directory pngfiles already exists');
end


close all
%% figure 1: demand

hFig = figure('name','Demand');
set(hFig, 'Position', [20 20 700/2 900/2])

demand_variables_labels = {'Real GDP','PCE Price Index','Fed Funds Rate','Nonresidential Investment','Unemployment','Real Wage'};
demand_variables_index = [1,11,25,4,19,15];

posterior_bands = 'orange';
posterior_medians  = 'red';


for i=1:length(demand_variables_index)
subplot(3,2,i)
variable   = demand_variables_index(i);
shock      = 1;
auv16=L_q16(:,variable,shock)'*plot_uv_scale;
auv50=L_q50(:,variable,shock)'*plot_uv_scale;
buv84=L_q84(:,variable,shock)'*plot_uv_scale;


[~,~]=jbfill([0:horizon],auv16,buv84,rgb(posterior_bands),rgb(posterior_bands),0,0.5);
hold on
hold on
hmed=plot([0:horizon],auv50,'linestyle','-','linewidth',2,'color',rgb(posterior_medians));
hold on
hzero=plot([0:horizon],auv50*0,'k-.','linewidth',0.5);



xlabel('Quarters','FontSize',ftsizexlabel)
set(gca,'FontSize',ftsizeaxis)
set(gca,'LineWidth',ftlinewidth)
title(demand_variables_labels(i))
box off
legend off
grid on
% Force exactly 5 equally spaced ticks including zero
ax = gca;
yLimits = ylim;



% Determine the range and create exactly 5 ticks with zero included
yRange = yLimits(2) - yLimits(1);
step = yRange / 4;  % 4 intervals = 5 ticks

% Find which tick should be zero
zeroPosition = -yLimits(1) / yRange;  % Position of zero as fraction of range
zeroIndex = round(zeroPosition * 4) + 1;  % Which of the 5 ticks should be zero

% Create 5 ticks with zero at the appropriate position
tickValues = zeros(1, 5);
tickValues(zeroIndex) = 0;

% Fill in the other ticks
for j = 1:5
    if j ~= zeroIndex
        tickValues(j) = (j - zeroIndex) * step;
    end
end

% Round to 2 decimal places
roundedTicks = round(tickValues, 2);

% Set the ticks
yticks(roundedTicks);
ylim([min(yLimits(1),roundedTicks(1)),max(yLimits(end),roundedTicks(end))])

% Set x-ticks at 0, 5, 10, 15, 20
xticks([0 5 10 15 20]);
xlim([0,20])
end

set(gcf, 'PaperPositionMode', 'auto');
print(['epsfiles/',num2str(spec),'_shock_demand_new_format.eps'],'-depsc');
print(['pngfiles/',num2str(spec),'_shock_demand_new_format.png'],'-dpng');



%% figure 2: investment

hFig = figure('name','Invesment');
set(hFig, 'Position', [20 20 700/2 900/2])

demand_variables_labels = {'Real GDP','PCE Price Index','Fed Funds Rate','Nonresidential Investment','Unemployment','Real Wage'}
demand_variables_index = [1,11,25,4,19,15];

for i=1:length(demand_variables_index)
subplot(3,2,i)
variable   = demand_variables_index(i);
shock      = 2;
auv16=L_q16(:,variable,shock)'*plot_uv_scale;
auv50=L_q50(:,variable,shock)'*plot_uv_scale;
buv84=L_q84(:,variable,shock)'*plot_uv_scale;


[~,~]=jbfill([0:horizon],auv16,buv84,rgb(posterior_bands),rgb(posterior_bands),0,0.5);
hold on
hold on
hmed=plot([0:horizon],auv50,'linestyle','-','linewidth',2,'color',rgb(posterior_medians));
hold on
hzero=plot([0:horizon],auv50*0,'k-.','linewidth',0.5);




hold on
xlabel('Quarters','FontSize',ftsizexlabel)
set(gca,'FontSize',ftsizeaxis)
set(gca,'LineWidth',ftlinewidth)
title(demand_variables_labels(i))
box off
legend off
grid on
% Force exactly 5 equally spaced ticks including zero
ax = gca;
yLimits = ylim;

% Determine the range and create exactly 5 ticks with zero included
yRange = yLimits(2) - yLimits(1);
step = yRange / 4;  % 4 intervals = 5 ticks

% Find which tick should be zero
zeroPosition = -yLimits(1) / yRange;  % Position of zero as fraction of range
zeroIndex = round(zeroPosition * 4) + 1;  % Which of the 5 ticks should be zero

% Create 5 ticks with zero at the appropriate position
tickValues = zeros(1, 5);
tickValues(zeroIndex) = 0;

% Fill in the other ticks
for j = 1:5
    if j ~= zeroIndex
        tickValues(j) = (j - zeroIndex) * step;
    end
end

% Round to 2 decimal places
roundedTicks = round(tickValues, 2);

% Set the ticks
yticks(roundedTicks);
ylim([min(yLimits(1),roundedTicks(1)),max(yLimits(end),roundedTicks(end))])

% Set x-ticks at 0, 5, 10, 15, 20
xticks([0 5 10 15 20]);
xlim([0,20])
end

set(gcf, 'PaperPositionMode', 'auto');
print(['epsfiles/',num2str(spec),'_shock_investment_new_format.eps'],'-depsc');
print(['pngfiles/',num2str(spec),'_shock_investment_new_format.png'],'-dpng');






%% figure 3: financial

hFig = figure('name','Financial');
set(hFig, 'Position', [20 20 700/2 900/2])

demand_variables_labels = {'Real GDP','PCE Price Index','Fed Funds Rate','Nonresidential Investment','Unemployment','Real Wage'}
demand_variables_index = [1,11,25,4,19,15];

for i=1:length(demand_variables_index)
subplot(3,2,i)
variable   = demand_variables_index(i);
shock      = 3;
auv16=L_q16(:,variable,shock)'*plot_uv_scale;
auv50=L_q50(:,variable,shock)'*plot_uv_scale;
buv84=L_q84(:,variable,shock)'*plot_uv_scale;


[~,~]=jbfill([0:horizon],auv16,buv84,rgb(posterior_bands),rgb(posterior_bands),0,0.5);

hold on
hold on
hmed=plot([0:horizon],auv50,'linestyle','-','linewidth',2,'color',rgb(posterior_medians));
hold on
hzero=plot([0:horizon],auv50*0,'k-.','linewidth',0.5);
hold on


xlabel('Quarters','FontSize',ftsizexlabel)
set(gca,'FontSize',ftsizeaxis)
set(gca,'LineWidth',ftlinewidth)
title(demand_variables_labels(i))
box off
legend off
grid on
% Force exactly 5 equally spaced ticks including zero
ax = gca;
yLimits = ylim;

% Determine the range and create exactly 5 ticks with zero included
yRange = yLimits(2) - yLimits(1);
step = yRange / 4;  % 4 intervals = 5 ticks

% Find which tick should be zero
zeroPosition = -yLimits(1) / yRange;  % Position of zero as fraction of range
zeroIndex = round(zeroPosition * 4) + 1;  % Which of the 5 ticks should be zero

% Create 5 ticks with zero at the appropriate position
tickValues = zeros(1, 5);
tickValues(zeroIndex) = 0;

% Fill in the other ticks
for j = 1:5
    if j ~= zeroIndex
        tickValues(j) = (j - zeroIndex) * step;
    end
end

% Round to 2 decimal places
roundedTicks = round(tickValues, 2);

% Set the ticks
yticks(roundedTicks);
ylim([min(yLimits(1),roundedTicks(1)),max(yLimits(end),roundedTicks(end))])

% Set x-ticks at 0, 5, 10, 15, 20
xticks([0 5 10 15 20]);
xlim([0,20])
end

set(gcf, 'PaperPositionMode', 'auto');
print(['epsfiles/',num2str(spec),'_shock_financial_new_format.eps'],'-depsc');
print(['pngfiles/',num2str(spec),'_shock_financial_new_format.png'],'-dpng');








%% figure 4: monetary

hFig = figure('name','Monetary');
set(hFig, 'Position', [20 20 700/2 900/2])

demand_variables_labels = {'Real GDP','PCE Price Index','Fed Funds Rate','Nonresidential Investment','Unemployment','Real Wage'}
demand_variables_index = [1,11,25,4,19,15];

for i=1:length(demand_variables_index)
subplot(3,2,i)
variable   = demand_variables_index(i);
shock      = 4;
auv16=L_q16(:,variable,shock)'*plot_uv_scale;
auv50=L_q50(:,variable,shock)'*plot_uv_scale;
buv84=L_q84(:,variable,shock)'*plot_uv_scale;

[~,~]=jbfill([0:horizon],auv16,buv84,rgb(posterior_bands),rgb(posterior_bands),0,0.5);

hold on
hold on
hmed=plot([0:horizon],auv50,'linestyle','-','linewidth',2,'color',rgb(posterior_medians));
hold on
hzero=plot([0:horizon],auv50*0,'k-.','linewidth',0.5);
hold on

xlabel('Quarters','FontSize',ftsizexlabel)
set(gca,'FontSize',ftsizeaxis)
set(gca,'LineWidth',ftlinewidth)
title(demand_variables_labels(i))
box off
legend off
grid on
% Force exactly 5 equally spaced ticks including zero
ax = gca;
yLimits = ylim;

% Determine the range and create exactly 5 ticks with zero included
yRange = yLimits(2) - yLimits(1);
step = yRange / 4;  % 4 intervals = 5 ticks

% Find which tick should be zero
zeroPosition = -yLimits(1) / yRange;  % Position of zero as fraction of range
zeroIndex = round(zeroPosition * 4) + 1;  % Which of the 5 ticks should be zero

% Create 5 ticks with zero at the appropriate position
tickValues = zeros(1, 5);
tickValues(zeroIndex) = 0;

% Fill in the other ticks
for j = 1:5
    if j ~= zeroIndex
        tickValues(j) = (j - zeroIndex) * step;
    end
end

% Round to 2 decimal places
roundedTicks = round(tickValues, 2);

% Set the ticks
yticks(roundedTicks);
ylim([min(yLimits(1),roundedTicks(1)),max(yLimits(end),roundedTicks(end))])

% Set x-ticks at 0, 5, 10, 15, 20
xticks([0 5 10 15 20]);
end

set(gcf, 'PaperPositionMode', 'auto');
print(['epsfiles/',num2str(spec),'_shock_monetary_new_format.eps'],'-depsc');
print(['pngfiles/',num2str(spec),'_shock_monetary_new_format.png'],'-dpng');






%% figure 5: government Spending

hFig = figure('name','Government Spending');
set(hFig, 'Position', [20 20 700/2 900/2])

demand_variables_labels = {'Real GDP','PCE Price Index','Fed Funds Rate','Nonresidential Investment','Unemployment','Real Wage'}
demand_variables_index = [1,11,25,4,19,15];

for i=1:length(demand_variables_index)
subplot(3,2,i)
variable   = demand_variables_index(i);
shock      = 5;
auv16=L_q16(:,variable,shock)'*plot_uv_scale;
auv50=L_q50(:,variable,shock)'*plot_uv_scale;
buv84=L_q84(:,variable,shock)'*plot_uv_scale;

[~,~]=jbfill([0:horizon],auv16,buv84,rgb(posterior_bands),rgb(posterior_bands),0,0.5);

hold on
hold on
hmed=plot([0:horizon],auv50,'linestyle','-','linewidth',2,'color',rgb(posterior_medians));
hold on
hzero=plot([0:horizon],auv50*0,'k-.','linewidth',0.5);
hold on
xlabel('Quarters','FontSize',ftsizexlabel)
set(gca,'FontSize',ftsizeaxis)
set(gca,'LineWidth',ftlinewidth)
title(demand_variables_labels(i))
box off
legend off
grid on
hold on

% Force exactly 5 equally spaced ticks including zero
ax = gca;
yLimits = ylim;

% Determine the range and create exactly 5 ticks with zero included
yRange = yLimits(2) - yLimits(1);
step = yRange / 4;  % 4 intervals = 5 ticks

% Find which tick should be zero
zeroPosition = -yLimits(1) / yRange;  % Position of zero as fraction of range
zeroIndex = round(zeroPosition * 4) + 1;  % Which of the 5 ticks should be zero

% Create 5 ticks with zero at the appropriate position
tickValues = zeros(1, 5);
tickValues(zeroIndex) = 0;

% Fill in the other ticks
for j = 1:5
    if j ~= zeroIndex
        tickValues(j) = (j - zeroIndex) * step;
    end
end

% Round to 2 decimal places
roundedTicks = round(tickValues, 2);

% Set the ticks
yticks(roundedTicks);
ylim([min(yLimits(1),roundedTicks(1)),max(yLimits(end),roundedTicks(end))])

% Set x-ticks at 0, 5, 10, 15, 20
xticks([0 5 10 15 20]);
xlim([0,20])
end

set(gcf, 'PaperPositionMode', 'auto');
print(['epsfiles/',num2str(spec),'_shock_Gspending_new_format.eps'],'-depsc');
print(['pngfiles/',num2str(spec),'_shock_Gspending_new_format.png'],'-dpng');






%% figure 6: technology

hFig = figure('name','technology');
set(hFig, 'Position', [20 20 700/2 900/2])

demand_variables_labels = {'Real GDP','PCE Price Index','Fed Funds Rate','Nonresidential Investment','Unemployment','Real Wage'}
demand_variables_index = [1,11,25,4,19,15];

for i=1:length(demand_variables_index)
subplot(3,2,i)
variable   = demand_variables_index(i);
shock      = 6;
auv16=L_q16(:,variable,shock)'*plot_uv_scale;
auv50=L_q50(:,variable,shock)'*plot_uv_scale;
buv84=L_q84(:,variable,shock)'*plot_uv_scale;


[~,~]=jbfill([0:horizon],auv16,buv84,rgb(posterior_bands),rgb(posterior_bands),0,0.5);

hold on
hold on
hmed=plot([0:horizon],auv50,'linestyle','-','linewidth',2,'color',rgb(posterior_medians));
hold on
hzero=plot([0:horizon],auv50*0,'k-.','linewidth',0.5);
hold on

xlabel('Quarters','FontSize',ftsizexlabel)
set(gca,'FontSize',ftsizeaxis)
set(gca,'LineWidth',ftlinewidth)
title(demand_variables_labels(i))
box off
legend off
grid on
% Force exactly 5 equally spaced ticks including zero
ax = gca;
yLimits = ylim;

% Determine the range and create exactly 5 ticks with zero included
yRange = yLimits(2) - yLimits(1);
step = yRange / 4;  % 4 intervals = 5 ticks

% Find which tick should be zero
zeroPosition = -yLimits(1) / yRange;  % Position of zero as fraction of range
zeroIndex = round(zeroPosition * 4) + 1;  % Which of the 5 ticks should be zero

% Create 5 ticks with zero at the appropriate position
tickValues = zeros(1, 5);
tickValues(zeroIndex) = 0;

% Fill in the other ticks
for j = 1:5
    if j ~= zeroIndex
        tickValues(j) = (j - zeroIndex) * step;
    end
end

% Round to 2 decimal places
roundedTicks = round(tickValues, 2);

% Set the ticks
yticks(roundedTicks);
ylim([min(yLimits(1),roundedTicks(1)),max(yLimits(end),roundedTicks(end))])

% Set x-ticks at 0, 5, 10, 15, 20
xticks([0 5 10 15 20]);
xlim([0,20])
end

set(gcf, 'PaperPositionMode', 'auto');
print(['epsfiles/',num2str(spec),'_shock_technology_new_format.eps'],'-depsc');
print(['pngfiles/',num2str(spec),'_shock_technology_new_format.png'],'-dpng');







%% figure 6: technology

hFig = figure('name','technology');
set(hFig, 'Position', [20 20 700/2 900/2])

demand_variables_labels = {'Real GDP','PCE Price Index','Fed Funds Rate','Nonresidential Investment','Unemployment','Real Wage'}
demand_variables_index = [1,11,25,4,19,15];

for i=1:length(demand_variables_index)
subplot(3,2,i)
variable   = demand_variables_index(i);
shock      = 6;
auv16=L_q16(:,variable,shock)'*plot_uv_scale;
auv50=L_q50(:,variable,shock)'*plot_uv_scale;
buv84=L_q84(:,variable,shock)'*plot_uv_scale;


[~,~]=jbfill([0:horizon],auv16,buv84,rgb(posterior_bands),rgb(posterior_bands),0,0.5);

hold on
hold on
hmed=plot([0:horizon],auv50,'linestyle','-','linewidth',2,'color',rgb(posterior_medians));
hold on
hzero=plot([0:horizon],auv50*0,'k-.','linewidth',0.5);
hold on

xlabel('Quarters','FontSize',ftsizexlabel)
set(gca,'FontSize',ftsizeaxis)
set(gca,'LineWidth',ftlinewidth)
title(demand_variables_labels(i))
box off
legend off
grid on
% Force exactly 5 equally spaced ticks including zero
ax = gca;
yLimits = ylim;

% Determine the range and create exactly 5 ticks with zero included
yRange = yLimits(2) - yLimits(1);
step = yRange / 4;  % 4 intervals = 5 ticks

% Find which tick should be zero
zeroPosition = -yLimits(1) / yRange;  % Position of zero as fraction of range
zeroIndex = round(zeroPosition * 4) + 1;  % Which of the 5 ticks should be zero

% Create 5 ticks with zero at the appropriate position
tickValues = zeros(1, 5);
tickValues(zeroIndex) = 0;

% Fill in the other ticks
for j = 1:5
    if j ~= zeroIndex
        tickValues(j) = (j - zeroIndex) * step;
    end
end

% Round to 2 decimal places
roundedTicks = round(tickValues, 2);

% Set the ticks
yticks(roundedTicks);
ylim([min(yLimits(1),roundedTicks(1)),max(yLimits(end),roundedTicks(end))])

% Set x-ticks at 0, 5, 10, 15, 20
xticks([0 5 10 15 20]);
xlim([0,20])
end

set(gcf, 'PaperPositionMode', 'auto');
print(['epsfiles/',num2str(spec),'_shock_technology_new_format.eps'],'-depsc');
print(['pngfiles/',num2str(spec),'_shock_technology_new_format.png'],'-dpng');






%% figure 7: laborsupply

hFig = figure('name','laborsupply');
set(hFig, 'Position', [20 20 700/2 900/2])

demand_variables_labels = {'Real GDP','PCE Price Index','Fed Funds Rate','Nonresidential Investment','Unemployment','Real Wage'}
demand_variables_index = [1,11,25,4,19,15];

for i=1:length(demand_variables_index)
subplot(3,2,i)
variable   = demand_variables_index(i);
shock      = 7;
auv16=L_q16(:,variable,shock)'*plot_uv_scale;
auv50=L_q50(:,variable,shock)'*plot_uv_scale;
buv84=L_q84(:,variable,shock)'*plot_uv_scale;


[~,~]=jbfill([0:horizon],auv16,buv84,rgb(posterior_bands),rgb(posterior_bands),0,0.5);

hold on
hold on
hmed=plot([0:horizon],auv50,'linestyle','-','linewidth',2,'color',rgb(posterior_medians));
hold on
hzero=plot([0:horizon],auv50*0,'k-.','linewidth',0.5);
hold on

xlabel('Quarters','FontSize',ftsizexlabel)
set(gca,'FontSize',ftsizeaxis)
set(gca,'LineWidth',ftlinewidth)
title(demand_variables_labels(i))
box off
legend off
grid on
% Force exactly 5 equally spaced ticks including zero
ax = gca;
yLimits = ylim;

% Determine the range and create exactly 5 ticks with zero included
yRange = yLimits(2) - yLimits(1);
step = yRange / 4;  % 4 intervals = 5 ticks

% Find which tick should be zero
zeroPosition = -yLimits(1) / yRange;  % Position of zero as fraction of range
zeroIndex = round(zeroPosition * 4) + 1;  % Which of the 5 ticks should be zero

% Create 5 ticks with zero at the appropriate position
tickValues = zeros(1, 5);
tickValues(zeroIndex) = 0;

% Fill in the other ticks
for j = 1:5
    if j ~= zeroIndex
        tickValues(j) = (j - zeroIndex) * step;
    end
end

% Round to 2 decimal places
roundedTicks = round(tickValues, 2);

% Set the ticks
yticks(roundedTicks);
ylim([min(yLimits(1),roundedTicks(1)),max(yLimits(end),roundedTicks(end))])

% Set x-ticks at 0, 5, 10, 15, 20
xticks([0 5 10 15 20]);
xlim([0,20])
end

set(gcf, 'PaperPositionMode', 'auto');
print(['epsfiles/',num2str(spec),'_shock_laborsupply_new_format.eps'],'-depsc');
print(['pngfiles/',num2str(spec),'_shock_laborsupply_new_format.png'],'-dpng');








%% figure 7: wagebargaining

hFig = figure('name','wagebargaining');
set(hFig, 'Position', [20 20 700/2 900/2])
demand_variables_labels = {'Real GDP','PCE Price Index','Fed Funds Rate','Nonresidential Investment','Unemployment','Real Wage'}
demand_variables_index = [1,11,25,4,19,15];

for i=1:length(demand_variables_index)
subplot(3,2,i)
variable   = demand_variables_index(i);
shock      = 8;
auv16=L_q16(:,variable,shock)'*plot_uv_scale;
auv50=L_q50(:,variable,shock)'*plot_uv_scale;
buv84=L_q84(:,variable,shock)'*plot_uv_scale;


[~,~]=jbfill([0:horizon],auv16,buv84,rgb(posterior_bands),rgb(posterior_bands),0,0.5);

hold on
hold on
hmed=plot([0:horizon],auv50,'linestyle','-','linewidth',2,'color',rgb(posterior_medians));
hold on
hzero=plot([0:horizon],auv50*0,'k-.','linewidth',0.5);
hold on

xlabel('Quarters','FontSize',ftsizexlabel)
set(gca,'FontSize',ftsizeaxis)
set(gca,'LineWidth',ftlinewidth)
title(demand_variables_labels(i))
box off
legend off
grid on
% Force exactly 5 equally spaced ticks including zero
ax = gca;
yLimits = ylim;

% Determine the range and create exactly 5 ticks with zero included
yRange = yLimits(2) - yLimits(1);
step = yRange / 4;  % 4 intervals = 5 ticks

% Find which tick should be zero
zeroPosition = -yLimits(1) / yRange;  % Position of zero as fraction of range
zeroIndex = round(zeroPosition * 4) + 1;  % Which of the 5 ticks should be zero

% Create 5 ticks with zero at the appropriate position
tickValues = zeros(1, 5);
tickValues(zeroIndex) = 0;

% Fill in the other ticks
for j = 1:5
    if j ~= zeroIndex
        tickValues(j) = (j - zeroIndex) * step;
    end
end

% Round to 2 decimal places
roundedTicks = round(tickValues, 2);

% Set the ticks
yticks(roundedTicks);
ylim([min(yLimits(1),roundedTicks(1)),max(yLimits(end),roundedTicks(end))])

% Set x-ticks at 0, 5, 10, 15, 20
xticks([0 5 10 15 20]);
xlim([0,20])
end

set(gcf, 'PaperPositionMode', 'auto');
print(['epsfiles/',num2str(spec),'_shock_wagebargaining_new_format.eps'],'-depsc');
print(['pngfiles/',num2str(spec),'_shock_wagebargaining_new_format.png'],'-dpng');

% 
% close all;
% 
% for i=1:35
% var = i;
% 
% figure(i)
% shock=5;
% 
% % Get data
% data1 = squeeze(L.L(1,var,shock,:));
% data2 = squeeze(Lcmy.L(1,var,shock,:));
% 
% % Define common bins
% minVal = min([data1; data2]);
% maxVal = max([data1; data2]);
% binEdges = linspace(minVal, maxVal, 30); % 30 bins, adjust as needed
% 
% % Create normalized histograms
% histogram(data1, binEdges, 'Normalization', 'probability');
% hold on;
% histogram(data2, binEdges, 'Normalization', 'probability');
% pause
% end
