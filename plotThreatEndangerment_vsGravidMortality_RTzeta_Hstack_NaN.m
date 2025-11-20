% Created by Noah Marshall - Summer 2020
% Modified by Rebecca Tyson

%{
Plots (histograms of # of simulations in which the population was endangered) overlaid by (histograms of # of simulations in which the population was NOT endangered)

Endangerement was defined by COSEWIC standards for a local population. The population is considered enangered if < 250 breeding adults are present. 
As the model is female focussed, the threshold for endangerment was set at 125

Values are estimated through repeated sampling from a triangular or uniform dist. with values described in 'Supplementary Material' of Toad Draft
%}

% For this version of the code, if (1-zetag*mz)<0, the output is NaN.

clear; clc;
%% Use the model to find the average sustainable level of highway mortality
% Idea of this is to run simulations over the parameter space to determine
% the average minimal sustainable survivorship level (maximal sustainable
% mortality) in the long run. Thus the model should be ran for a long time
% and some populatin level should be chosen that represents a unsustainable
% population level.
%% Initialize parameters sampled using LHS
K = 1000000;
mortz_max = 0.5;
zetag_max = 2;
runs=1000000;
threatened_thresh = 500;
endangered_thresh = 125;

mortz_dist = makedist('Uniform','lower',0.1,'upper',mortz_max); % basic road-crossing mortality
a_dist = makedist('Uniform','lower',0.053,'upper',0.24); % alpha
r_dist = makedist('Uniform','lower',3000,'upper',7000); % r
d_dist = makedist('Uniform','lower',0.0012,'upper',0.0021); % delta
dx_dist = makedist('Uniform','lower',0.00211,'upper',0.0333); % mu_x
dy_dist = makedist('Uniform','lower',0.0041,'upper',0.0211); % mu_y
dz_dist = makedist('Uniform','lower',0.0016,'upper',0.0041); % mu_z
zetag_dist = makedist('Uniform', 'lower', 1.0, 'upper', zetag_max); % zeta_g
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

%% Run Simulations
for i=1:runs     
    params = [mortzs(i), as(i), ds(i), rs(i), zetaxs(i), dxs(i), dys(i), dzs(i)]; 

    if ( (1-zetags(i)*mortzs(i))<0 )
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
end
%% Find simulations such that the final adult population is less than the required threshold
Epoor_zetag = zetags(popMin <= endangered_thresh);
Tpoor_zetag = zetags((popMin < threatened_thresh) & (popMin > endangered_thresh));
good_zetag = zetags(popMin >= threatened_thresh);
zetag = 1:0.02:zetag_max;
figure;
hold on;
gh = histogram(good_zetag,1:0.02:zetag_max,'Facecolor','green', 'FaceAlpha',0.6);
eh = histogram(Epoor_zetag,1:0.02:zetag_max,'Facecolor','red', 'FaceAlpha', 0.6); 
th = histogram(Tpoor_zetag,1:0.02:zetag_max,'Facecolor','yellow', 'FaceAlpha', 1); 
grid
yt = get(gca, 'YTick'); 
top = max([eh.Values,th.Values,gh.Values]);
axis([1 zetag_max 0 top])
ytix = linspace(0, top, 11);  % Rebecca's code
set(gca, 'YTick',ytix, 'YTickLabel',round(ytix/top,2))
set(gca,'FontSize',14)
set(findall(gcf,'type','text'),'FontSize',14)
xlabel('gravid mortality factor ($\zeta_g$)', 'interpreter', 'latex')
ylabel('Probability of Endangerment', 'interpreter', 'latex')
title(['$K = $\ ', num2str(K), ', $\max(m_z)=$\ ', num2str(mortz_max)], 'interpreter', 'latex')

% Now make line plots where I'm sure that the sume of the three
% probabilities, at each x value, is 1
hgh = gh.BinCounts;
heh = eh.BinCounts;
hth = th.BinCounts;
hgh_norm = hgh./(hgh+heh+hth);
heh_norm = heh./(hgh+heh+hth);
hth_norm = hth./(hgh+heh+hth);
ebins = gh.BinEdges;
ebins_left = ebins(1:end-1);
figure;
plot(ebins_left, hgh_norm, 'g', ebins_left, hth_norm, 'y', ...
    ebins_left, heh_norm, 'r', 'LineWidth', 3);
xlabel('gravid mortality factor ($\zeta_g$)', 'interpreter', 'latex')
ylabel('Probability of Endangerment', 'interpreter', 'latex')

% Now make line plots where I'm sure that the sume of the three
% probabilities, at each x value, is 1
% add a curve that is the sum of the probabilities of being threatened or
% endangered
hgh = gh.BinCounts;
heh = eh.BinCounts;
hth = th.BinCounts;
hgh_norm = hgh./(hgh+heh+hth);
heh_norm = heh./(hgh+heh+hth);
hth_norm = hth./(hgh+heh+hth);
ebins = gh.BinEdges;
ebins_left = ebins(1:end-1);
figure;
plot(ebins_left, hgh_norm, 'g', ebins_left, hth_norm, 'y', ...
    ebins_left, heh_norm, 'r', 'LineWidth', 3);
hold on
plot(ebins_left, hth_norm+heh_norm, 'k', 'LineWidth', 1);
hold off
xlabel('gravid mortality factor ($\zeta_g$)', 'interpreter', 'latex')
ylabel('Probability of Status', 'interpreter', 'latex')
grid;

% experimenting with colours
% try these:
%    instead of "red", use "#A2142F" (dark red)
%    instead of "green", use "#4DBEE" (light blue)
%    instead of "yellow", use "#EDB129" (gold)

% stacked histogram with colourblind-safe palette
bstep = 0.02;
bvec = 1:bstep:zetag_max;
[NEpoor,edges_Epoor] = histcounts(Epoor_zetag, bvec);
[NTpoor,edges_Tpoor] = histcounts(Tpoor_zetag, bvec);
[Ngood,edges_good] = histcounts(good_zetag, bvec);
Ntot = NEpoor + NTpoor + Ngood;
figure;
bEndg = bar(1+bstep/2:bstep:zetag_max-bstep/2, [NEpoor./Ntot; NTpoor./Ntot; Ngood./Ntot], ...
    'stacked', 'BarWidth', 1);
bEndg(2).FaceColor = "#ffd700"; 
bEndg(1).FaceColor = "#ea5f94";
bEndg(3).FaceColor = "#9d02d7";
legend('Endangered', 'Threatened', 'Healthy', 'interpreter', 'latex', ...
    'FontSize', 16);
xlabel('gravid mortality factor ($\zeta_g$)', 'interpreter', 'latex', ...
    'FontSize', 18)
ylabel('Probability of Status', 'interpreter', 'latex', ...
    'FontSize', 18)
axis([1 zetag_max 0 1]);

