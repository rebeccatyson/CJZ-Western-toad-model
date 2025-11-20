% Created by Noah Marshall - Summer 2020
% Modified by Rebecca Tyson

%{
Plots (histograms of # of simulations in which the population was endangered) overlaid by (histograms of # of simulations in which the population was NOT endangered)

Endangerement was defined by COSEWIC standards for a local population. The population is considered enangered if < 250 breeding adults are present. 
As the model is female focussed, the threshold for endangerment was set at 125

Values are estimated through repeated sampling from a triangular or uniform dist. with values described in 'Supplementary Material' of Toad Draft
%}

% In this version of the code, the output is NaN if (1-zetag*mz)<0.

clear; clc;
%% Use the model to find the average sustainable level of highway mortality
% Idea of this is to run simulations over the parameter space to determine
% the average minimal sustainable survivorship level (maximal sustainable
% mortality) in the long run. Thus the model should be ran for a long time
% and some populatin level should be chosen that represents a unsustainable
% population level.
%% Initialize parameters sampled using LHS
K = 3000000;
runs=1000000;
threatened_thresh = 500;
endangered_thresh = 125;

mortz_dist = makedist('Uniform','lower',0,'upper',1); % basic road-crossing mortality
a_dist = makedist('Uniform','lower',0.053,'upper',0.24); % alpha
r_dist = makedist('Uniform','lower',3000,'upper',7000); % r
d_dist = makedist('Uniform','lower',0.0014,'upper',0.0021); % delta
dx_dist = makedist('Uniform','lower',0.00211,'upper',0.0333); % mu_x
dy_dist = makedist('Uniform','lower',0.0041,'upper',0.0211); % mu_y
dz_dist = makedist('Uniform','lower',0.0016,'upper',0.0041); % mu_z
zetag_dist = makedist('Uniform', 'lower', 1.0, 'upper', 3.0); % zeta_g
zetax_dist = makedist('Uniform', 'lower', 1.0, 'upper', 2.0); % zeta_x

mortzs = random(mortz_dist,runs,1);
as = random(a_dist,runs,1);
ds = random(d_dist,runs,1);
rs = random(r_dist,runs,1);
dxs = random(dx_dist,runs,1);
dys = random(dy_dist,runs,1);
dzs = random(dz_dist,runs,1);
zetags = random(zetag_dist,runs,1);
zetaxs = random(zetax_dist, runs, 1);

% open output file for parameter values
fileID = fopen('dataThreatEndangerment_vsMortality_RTzeta_NEW.txt', 'w');
fprintf(fileID, '%+8s, %+8s, %+8s, %+8s, %+8s, %+8s, %+8s, %+8s, %+8s, %+8s\n', ...
    'z1', 'mz', 'alpha', 'r', 'delta', 'zetax', 'zetag', 'mux', 'muy', 'muz')

%% Run Simulations
for i=1:runs     
    params = [mortzs(i), as(i), ds(i), rs(i), zetaxs(i), dxs(i), dys(i), dzs(i)]; 

    if ((1-zetags(i)*mortzs(i)) < 0)
        z1 = NaN;
    else
        z1 = getZ1star_mort_zeta_NEW(mortzs(i), as(i), rs(i), zetags(i), zetaxs(i), ...
            ds(i), dxs(i), dys(i), dzs(i), K);
        if z1 <= 0
            z1 = 0;
        end
    end

    popMax(i,:) = z1/(params(1)^2 * params(2));
    popMin(i,:) = z1/as(i);

    % print parameter sets leading to healthy populations
    if ((dxs(i) >= 0.027) && (z1 > threatened_thresh) && (mortzs(i) >= 0.1))
        % disp(['mz = ', num2str(mortzs(i)), ', alpha = ', num2str(as(i)), ...
        %       ', r = ', num2str(rs(i)), ', delta = ', num2str(ds(i)), ...
        %       ', zetax = ', num2str(zetaxs(i)), ', zetag = ', num2str(zetags(i)), ...
        %       ', mux = ', num2str(dxs(i)), ', muy = ', num2str(dys(i)),  ...
        %       ', muz = ', num2str(dzs(i)), ', z1 = ', num2str(z1)]);   
        fprintf(fileID, '%8.4f, %8.4f, %8.4f, %8.4f, %8.4f, %8.4f, %8.4f, %8.4f, %8.4f, %8.4f\n', ...
            z1, mortzs(i), as(i), rs(i), ds(i), zetaxs(i), zetags(i), dxs(i), dys(i), dzs(i));
    end
end

% close the data file
fclose(fileID)

%% Find simulations such that the final adult population is less than the required threshold
Epoor_mortz = mortzs(popMin <= endangered_thresh);
Tpoor_mortz = mortzs((popMin < threatened_thresh) & (popMin > endangered_thresh));
%Tpoor_mortz = mortzs(popMin < threatened_thresh);
good_mortz = mortzs(popMin >= threatened_thresh);
mortz = 0:0.02:1;
figure;
gh = histogram(good_mortz,0:0.02:1,'Facecolor','green'); % check hist or histcounts, to see if I can get the bars
hold on;
eh = histogram(Epoor_mortz,0:0.02:1,'Facecolor','red'); 
th = histogram(Tpoor_mortz,0:0.02:1,'Facecolor','yellow'); 
grid
ecrossing = mortz(find(eh.Values < th.Values,1));
tcrossing = mortz(find(th.Values < gh.Values,1));
% xline(ecrossing,'k--');
% xline(tcrossing,'k--');
yt = get(gca, 'YTick'); 
top = max([eh.Values,th.Values,gh.Values]);
axis([0 1 0 top])
% ytix = linspace(min(yt), max(yt), 10);  % Noah's code
ytix = linspace(0, top, 11);  % Rebecca's code
set(gca, 'YTick',ytix, 'YTickLabel',round(ytix/top,2))
% set(gca, 'YTick',ytix, 'YTickLabel',round(ytix/25000,2))  % Rebecca
set(gca,'FontSize',14)
set(findall(gcf,'type','text'),'FontSize',14)
crossings_at = sprintf('%s = %d and %d',"crossings at: ", ecrossing, ...
    tcrossing);
%legend('Minimum breeding Pop. < 125','Minimum breeding Pop. >125',cross_at)
% str = sprintf('%s %d)',"s_z (n = ",runs);
% xlabel(str)
xlabel('mortality ($m_z$)', 'interpreter', 'latex')
ylabel('Probability of Endangerment', 'interpreter', 'latex')
%title('Predicted Probability of Endangerment')

% experimenting with colours
% try these:
%    instead of "red", use "#A2142F" (dark red)
%    instead of "green", use "#4DBEE" (light blue)
%    instead of "yellow", use "#EDB129" (gold)

% Making a stacked plot
[NEpoor,edges_Epoor] = histcounts(Epoor_mortz, 0:0.02:1);
[NTpoor,edges_Tpoor] = histcounts(Tpoor_mortz, 0:0.02:1);
[Ngood,edges_good] = histcounts(good_mortz, 0:0.02:1);

% stack1 = figure(2);
% bEndg = bar(0.01:0.02:0.99, [NEpoor; NTpoor; Ngood], 'stacked', 'BarWidth', 1);
% bEndg(1).FaceColor = "#ffd700"; 
% bEndg(2).FaceColor = "#fa8775";
% bEndg(3).FaceColor = "#cd34b5";
% 
% stack2 = figure(3);
% bEndg = bar(0.01:0.02:0.99, [NEpoor; NTpoor; Ngood], 'stacked', 'BarWidth', 1);
% bEndg(1).FaceColor = "#ffd700"; 
% bEndg(2).FaceColor = "#ea5f94";
% bEndg(3).FaceColor = "#0000ff";
% 
% stack3 = figure(4);
% bEndg = bar(0.01:0.02:0.99, [NEpoor; NTpoor; Ngood], 'stacked', 'BarWidth', 1);
% bEndg(2).FaceColor = "#ffd700"; 
% bEndg(1).FaceColor = "#ea5f94";
% bEndg(3).FaceColor = "#0000ff";
% 
% stack4 = figure(5);
% bEndg = bar(0.01:0.02:0.99, [NEpoor; NTpoor; Ngood], 'stacked', 'BarWidth', 1);
% bEndg(2).FaceColor = "#ffd700"; 
% bEndg(1).FaceColor = "#fa8775";
% bEndg(3).FaceColor = "#cd34b5";

% % non-normalised counts, colourblind-safe palette
% stack3 = figure;
% bEndg = bar(0.01:0.02:0.99, [NEpoor; NTpoor; Ngood], 'stacked', 'BarWidth', 1);
% bEndg(2).FaceColor = "#ffd700"; 
% bEndg(1).FaceColor = "#ea5f94";
% bEndg(3).FaceColor = "#9d02d7";
% legend('Endangered', 'Threatened', 'Healthy');

% normalised counts, colourblind-safe palette
stack4 = figure;
Ntot = NEpoor + NTpoor + Ngood;
bEndg = bar(0.01:0.02:0.99, [NEpoor./Ntot; NTpoor./Ntot; Ngood./Ntot], ...
    'stacked', 'BarWidth', 1);
bEndg(2).FaceColor = "#ffd700"; 
bEndg(1).FaceColor = "#ea5f94";
bEndg(3).FaceColor = "#9d02d7";
legend('Endangered', 'Threatened', 'Healthy', 'interpreter', 'latex', ...
    'FontSize', 16, 'Location', 'southeast');
xlabel('mortality ($m_z$)', 'interpreter', 'latex', 'FontSize', 18)
ylabel('Probability of Status', 'interpreter', 'latex', 'FontSize', 18)
axis([0 1 0 1]);
