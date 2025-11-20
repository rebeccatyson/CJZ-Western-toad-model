% Created by Noah Marshall - Summer 2020

%{
Plots a histogram of calculated values of the extinction threshold. Extinction threshold is the value of road survivorship that will lead to extinction.

Values are estimated through repeated sampling from a triangular dist. with values described in 'Supplementary Material' of Toad Draft
%}

clear; clc;
%% Use the model to find the average sustainable level of highway mortality
% Idea of this is to run simulations over the parameter space to determine
% the average minimal sustainable survivorship level (maximal sustainable
% mortality) in the long run. Thus the model should be ran for a long time
% and some populatin level should be chosen that represents a unsustainable
% population level.
%% Initialize parameters sampled using LHS
runs=100000;
years = 50;
zetag_max = 3;

% sz_dist = makedist('Triangular','a',0.837,'b',0.899,'c',0.966);
% g_dist = makedist('Triangular','a',0.5,'b',0.75,'c',1);
% alpha_dist = makedist('Triangular','a',0.053,'b',0.085,'c',0.24);
% r_dist = makedist('Triangular','a',2900,'b',4000,'c',5900);
% delta_dist = makedist('Triangular','a',0.0014,'b',0.00165,'c',0.0021);
% mux_dist = makedist('Triangular','a',0.00211,'b',0.0177,'c',0.0333);
% muy_dist = makedist('Triangular','a',0.0041,'b',0.0126,'c',0.0211);
% muz_dist = makedist('Triangular','a',0.0016,'b',0.00285,'c',0.0041);

mortz_dist = makedist('Uniform', 'lower', 0, 'upper', 1);
%g_dist = makedist('Uniform','lower',0.5,'upper',1);
alpha_dist = makedist('Uniform','lower',0.053,'upper',0.24);
r_dist = makedist('Uniform','lower',3000,'upper',7000);
delta_dist = makedist('Uniform','lower',0.0014,'upper',0.0021);
mux_dist = makedist('Uniform','lower',0.00211,'upper',0.0333); % mu_x
muy_dist = makedist('Uniform','lower',0.0041,'upper',0.0211);  % mu_y
muz_dist = makedist('Uniform','lower',0.0016,'upper',0.0041);  % mu_z
zetag_dist = makedist('Uniform', 'lower', 1.0, 'upper', zetag_max); % zeta_g
zetax_dist = makedist('Uniform', 'lower', 1.0, 'upper', 2.0); % zeta_x

mortzs = random(mortz_dist, runs, 1);
alphas = random(alpha_dist,runs,1);
deltas = random(delta_dist,runs,1);
rs = random(r_dist,runs,1);
%gs = random(g_dist,runs,1);
muxs = random(mux_dist,runs,1);
muys = random(muy_dist,runs,1);
muzs = random(muz_dist,runs,1);
zetags = random(zetag_dist, runs,1);
zetaxs = random(zetax_dist, runs,1);

% open the data file for debugging
fileID = fopen('dataExtinctionThreshold_RTzeta_NEW.txt', 'w');
fprintf(fileID, '%+8s, %+8s, %+8s, %+8s, %+8s, %+8s, %+8s, %+8s, %+8s\n', ...
    'mzc', 'alpha', 'r', 'delta', 'zetax', 'zetag', 'mux', 'muy', 'muz')

for i = 1:runs

    mortz = mortzs(i);
%    gamma = gs(i);
    alpha = alphas(i);
    r = rs(i);
    delta = deltas(i);
    mux = muxs(i);
    muy = muys(i);
    muz = muzs(i);
    zetag = zetags(i);
    zetax = zetaxs(i);    
  
    threshold(i) = getExtinctThreshold_RTzeta_NEW(delta, mux, muy, muz, alpha, zetag, zetax, r);
    ext(i) = threshold(i) <= mortz;

    % print to the debugging data file
    if ((zetag*threshold(i) > 1) || (zetax*threshold(i) > 1))
         fprintf(fileID, '%8.4f, %8.4f, %8.4f, %8.4f, %8.4f, %8.4f, %8.4f, %8.4f, %8.4f\n', ...
            threshold(i), alpha, r, delta, zetax, zetag, mux, muy, muz);
    end

end

% close the data file
fclose(fileID)

% statistics
% first remove all the NaN entries, as if there are too many they mess up
% the mean and median
thresh = sort(threshold);
[MT,iMT] = max(thresh);
threshReal = thresh(1:iMT);

SEM = std(threshReal)/sqrt(length(threshReal));   % Standard Error
ts = tinv([0.025  0.975],length(threshReal)-1);      % T-Score
CI = mean(threshReal) + ts*SEM;                      % Confidence Intervals

h_fig = figure('Name', 'Mortality')
histogram(threshold,'Normalization','probability')
grid
set(gca,'FontSize',14)
set(findall(gcf,'type','text'),'FontSize',14)
xlabel('Adult Road-Crossing Mortality Extirpation Threshold, $m^c_z$', 'interpreter', 'latex')
ylabel('Relative Probability', 'interpreter', 'latex')
me = mean(threshReal)
med = median(threshReal)
