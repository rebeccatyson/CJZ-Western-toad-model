% Plot zmin as a function of gravid female mortality factor
% This code is replacing the maple code so that I can make figures with a
% consistent look for the Toad paper

% In this version, the code returns NaN if (1-zetag*mz)<0.

% clear all
clear all;

% Parameters
mz_vec = [0.1, 0.3];
alpha = 0.11;
r = 6000;
K_vec = [1000000, 2000000, 4000000];
delta = 0.0017;
mux = 0.027;
muy = 0.01;
muz = 0.0029;
zetax = 1.3;
zetag_vec = 1.0:0.01:3;
L = 365;
T = 90;

% composite parameter
Qright1 = exp(-muz*(L-T))*(1 - exp(-(delta+muy)*L));
Qright2 = exp(-(delta+muy)*(L-T))*(exp(-muz*L)-1);
Qfrac = exp(-mux*T)/(1-exp(-(delta+muy)*L));
Q = Qfrac*(Qright1 + Qright2);

% now compute zmin
zmin = zeros(length(K_vec), length(mz_vec), length(zetag_vec)); 
ik = 0;
iz = 0;
im = 0;
for K = K_vec
    ik = ik+1;
    for mz = mz_vec
        im = im + 1;
        for zetag = zetag_vec
            iz = iz + 1;

            if ((1-zetag*mz)<0 || (1-zetax*mz)<0)
                disp(['zetag=', num2str(zetag), ', mz=', num2str(mz), ', zetax=', num2str(zetax)]);
                zmin(ik, im, iz) = NaN;
            else
                zminfrac1 = alpha*(1-mz)*K/(2*(r-1));
                zminfrac2 = Q*r*delta/(delta+muy-muz);
                zminfrac3 = (1-zetax*mz)*(1-zetag*mz)/(1-alpha*(1-zetag*mz)*mz*exp(-muz*L));
                zmin(ik, im, iz) = max(zminfrac1*(zminfrac2*zminfrac3 - 1),0);
            end

        end
        iz = 0;
    end
    im = 0;
end

% plots
% zmin as a function of zetag
hmort = figure('Name', 'zmin vs mortality');

zmin11 = reshape(zmin(1,1,:),1,[]);
zmin21 = reshape(zmin(2,1,:),1,[]);
zmin31 = reshape(zmin(3,1,:),1,[]);
zmin12 = reshape(zmin(1,2,:),1,[]);
zmin22 = reshape(zmin(2,2,:),1,[]);
zmin32 = reshape(zmin(3,2,:),1,[]);

plot(zetag_vec, zmin11, 'r--', zetag_vec, zmin21, 'r-', zetag_vec, zmin31, 'r:', ...
    zetag_vec, zmin12, 'b--', zetag_vec, zmin22, 'b-', zetag_vec, zmin32, 'b:', ...
    'LineWidth', 3);
legend('$K=1,000,000$ \& $m_z=0.1$', '$K=2,000,000$ \& $m_z=0.1$', '$K=4,000,000$ \& $m_z=0.1$', ...
    '$K=1,000,000$ \& $m_z=0.3$', '$K=2,000,000$ \& $m_z=0.3$', '$K=4,000,000$ \& $m_z=0.3$', ...
    'interpreter', 'latex', 'FontSize', 16);
xlabel('Additional gravid female mortality factor $\zeta_g$', 'Interpreter', 'Latex', 'FontSize', 16);
ylabel('Minimum breeding population $z_{min}$', 'Interpreter', 'Latex', 'FontSize', 16);
hold on
plot(zetag_vec, 500*ones(size(zmin11)), 'k:', 'LineWidth', 2, 'HandleVisibility', 'Off');
plot(zetag_vec, 125*ones(size(zmin11)), 'k:', 'LineWidth', 2, 'HandleVisibility', 'Off');
hold off

hmort2 = figure('Name', 'zmin vs mortality, better colours');
p11 = plot(zetag_vec, zmin11, '--', 'Color', '#ea5f94', 'LineWidth', 3);
hold on
p12 = plot(zetag_vec, zmin21, '-', 'Color', '#ea5f94', 'LineWidth', 3);
p13 = plot(zetag_vec, zmin31, ':', 'Color', '#ea5f94', 'LineWidth', 3);
p21 = plot(zetag_vec, zmin12, '--', 'Color', '#9d02d7', 'LineWidth', 3);
p22 = plot(zetag_vec, zmin22, '-', 'Color', '#9d02d7', 'LineWidth', 3);
p23 = plot(zetag_vec, zmin32, ':', 'Color', '#9d02d7', 'LineWidth', 3);
legend('$K=1,000,000$ \& $m_z=10$\%', '$K=2,000,000$ \& $m_z=10$\%', '$K=4,000,000$ \& $m_z=10$\%', ...
    '$K=1,000,000$ \& $m_z=30$\%', '$K=2,000,000$ \& $m_z=30$\%', '$K=4,000,000$ \& $m_z=30$\%', ...
    'interpreter', 'latex', 'FontSize', 16);
xlabel('Additional gravid female mortality factor $\zeta_g$', 'Interpreter', 'Latex', 'FontSize', 16);
ylabel('Minimum breeding population $z_{min}$', 'Interpreter', 'Latex', 'FontSize', 18);
plot(zetag_vec, 500*ones(size(zmin11)), 'k:', 'LineWidth', 2, 'HandleVisibility', 'Off');
plot(zetag_vec, 125*ones(size(zmin11)), 'k:', 'LineWidth', 2, 'HandleVisibility', 'Off');
hold off

hmort3 = figure('Name', 'zmin vs mortality, better colours', 'Position', [680 415 830 470]);
p11 = plot(zetag_vec, zmin11, '--', 'Color', '#A2142F', 'LineWidth', 3);
hold on
p12 = plot(zetag_vec, zmin21, '-', 'Color', '#A2142F', 'LineWidth', 3);
p13 = plot(zetag_vec, zmin31, ':', 'Color', '#A2142F', 'LineWidth', 3);
p21 = plot(zetag_vec, zmin12, '--', 'Color', '#4DBEEE', 'LineWidth', 3);
p22 = plot(zetag_vec, zmin22, '-', 'Color', '#4DBEEE', 'LineWidth', 3);
p23 = plot(zetag_vec, zmin32, ':', 'Color', '#4DBEEE', 'LineWidth', 3);
legend('$K=1,000,000$ \& $m_z=10$\%', '$K=2,000,000$ \& $m_z=10$\%', '$K=4,000,000$ \& $m_z=10$\%', ...
    '$K=1,000,000$ \& $m_z=30$\%', '$K=2,000,000$ \& $m_z=30$\%', '$K=4,000,000$ \& $m_z=30$\%', ...
    'interpreter', 'latex', 'FontSize', 16);
xlabel('Gravid female increased mortality factor $\zeta_g$', 'Interpreter', 'Latex', 'FontSize', 16);
ylabel('Minimum breeding population $z_{min}$', 'Interpreter', 'Latex', 'FontSize', 18);
plot(zetag_vec, 500*ones(size(zmin11)), 'k:', 'LineWidth', 2, 'HandleVisibility', 'Off');
plot(zetag_vec, 125*ones(size(zmin11)), 'k:', 'LineWidth', 2, 'HandleVisibility', 'Off');
hold off
