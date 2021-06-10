set(0,'DefaultFigureWindowStyle','docked')
set(gca,'fontname','arial')

clear all;
close all;
%% Parameters from literature
beta1 = 1/80;           %basic trancription rate - nM/min (Semsey et al. 2012) -  original value  1/80
beta2 = 1/80;           %basic trancription rate - nM/min (Semsey et al. 2012) -  original value  1/80

% 
% delta_p = 0.0274;        % Protein degradation rate - 1/min - Boda et al. 2016 - original value 0.0274
% delta_CepR = 0.002 ;     %  degradation rate of CepR, unit min-1 - Weber and Buceta 2013 - original value 0 .002 min^−1
% delta_CepI = 0.01;       %  degradation rate of CepI, unit min-1 - Weber and Buceta 2013 - original value 0.01 min-1
% delta_CepRstar = 0.002;  %  degradation rate of CepR*, unit min-1 - Weber and Buceta 2013 - original value 0 .002 min^−1
% delta_CciRstar = 0.002;  %  degradation rate of CciR*, unit min-1 - Weber and Buceta 2013 - original value 0 .002 min^−1
% delta_CciI = 0.01;       %  degradation rate of Ccii, unit min-1 - Weber and Buceta 2013 - original value 0.01 min-1
% delta_CciR = 0.002 ;     %  degradation rate of CciR, unit min-1 - Weber and Buceta 2013 - original value 0 .002 min^−1
% delta_HSL = 0.001 ;    % Abiotic degradation rate of C8 and C6 AHL - Weber and Buceta 2013  original value = 0.001 min^-1 - OR	- https://doi.org/10.1111/j.1574-6941.2009.00828.x - original value: 0.005545h^-1
% delta_OACP = 0.001 ;    % Abiotic degradation rate of C8 and C6 AHL - Weber and Buceta 2013  original value = 0.001 min^-1 - OR	- https://doi.org/10.1111/j.1574-6941.2009.00828.x - original value: 0.005545h^-1
delta = 0.01; % Combined protein and metabolite degradation rate.

alpha_C8 =10 ;       % C8HSL diffusion rate - min^-1  - Weber and Buceta 2013 - original value = 10
alpha_C6 = 10;        % C6HSL diffusion rate - min^-1  - Weber and Buceta 2013
%alpha_O = 0.0001;    % OACP diffusion rate
alpha_OACP = .1;  %OACP production rate, from induced production rate of AHL in Fekete et al. 2010

gamma_C8colony = .015;      % rate of C8HSL production - Weber and Buceta 2013 - original value: 0.015 ??? CONTROLLARE ???
gamma_C6colony = .015;      % rate of C6HSL production - Weber and Buceta 2013 - original value: 0.015 ??? CONTROLLARE ???

k_on_CCIR_star = 0.1;  % CciR* - 1/min*nM - C6HSL association constant - Weber and Buceta 2013 - original value: 0.1
k_off_CCIR_star = 10; % CciR* - 1/min - C6HSL dissociation rate - Weber and Buceta 2013
v_max_CepI = 0.0041;   % CepI masimal rate - 1/min - Buroni et al. 2018
km_CepI = 0.068e-6;    % CepI MM constant - nM - Buroni et al. 2018 - original value: 0.068e-6
v_max_CciI = 0.0041;   % CciI masimal rate - 1/min - Buroni et al. 2018
km_CccI = 0.0068e-6;   % CciI MM constant - nM - Buroni et al. 2018
k_on_CepR_star = 0.1;  % CepR* C8 association constant - 1/min*nM - Weber and Buceta 2013 - original value 0.1
k_off_CepR_star = 10;  % CepR* C8 dissociation rate - 1/min - Weber and Buceta 2013  - original value 10

%       K_cat_CepI = 4.5;      %original value 0.075; %CepI catalytic rate
%       constant, unit s-1 - From Scoffone et al.
%       K_m_CepI = 75000; %original value 0.075 ;%CepI Michaelis  constant,
%       unit mM - From Scoffone et al.
n = 1.7;               % CepR Hill coefficient
m = 1.7;               % CepR Hill coefficient
%mu = .01;             % Cell growth rate, unit min-1
w = 1.7;               % CepR* Hill coefficient

%% Parameters without literature knowledge
k1 = 10;             % CepR*-CepI activation coefficient  http://2018.igem.org/Team:Tsinghua/Model?1
k2 = 10;             % CciR*-CepI/R repression coefficient
k3 = 10;             % CepR* - CepR repression coefficient
k4 = 10;
k5 = 10;
%% Initial concentrations
CEPI = 0;     % CepI protein
CEPR=0;      % CepR protein
CEPRstar = 0; % CepR* activated form of CepR protein
OACP = 0;    % Oxooctanoyl acyl carrierp protein , reuìquired for C8/C6-HSL biosynthesis
C8_I =0;  % N-Octanoyl-L-homoserine lactone, intracellular
C8_E =0;  % N-Octanoyl-L-homoserine lactone, extracellular

C6_I = 0;  % N-Hexanoyl-L-homoserine lactone, intracellular
C6_E = 0;  % N-Hexanoyl-L-homoserine lactone, extracellular
CCII = 0;     % CciI protein
CCIR = 0;     % CciR protein
CCIRstar = 0; % CciR* activated form of CciR protein

%% Simulation time (minutes)

tspan = [0 1200];
t_last = [ 50000, 50000, 50000, 50000, 50000, 50000]
l = [ 0.0001; 0.0005; 0.001; 0.005; 0.01; .1];
t_last = [ 50000]
l = [ 0.0001];


for i = 1:10:50
    
    for j = 1:1:length(l)
        
    tspan = [0 t_last(j)];
    delta= l(j);
    
%% core model (M1) - Simulation
fprintf('\n\nrunning simulation\n');

C8_E = i;

%parameters_M1 = [v_max_CepI; km_CepI; w; k1; alpha_OACP; k_on_CepR_star; k_off_CepR_star; beta; k3; n; alpha_C8; delta_HSL; delta_CepR; delta_CepI; delta_CepRstar; delta_OACP; mu;gamma];
parameters_M1 = [v_max_CepI; km_CepI; w; k1; alpha_OACP; k_on_CepR_star; k_off_CepR_star; beta1; k3; n; alpha_C8; delta; gamma_C8colony];
parameters_M1_names = ["v_max_CepI"; "km_CepI"; "w"; "k1"; "alpha_OACP"; "k_on_CepR_star"; "k_off_CepR_star"; "beta1"; "k3"; "n"; "alpha_C8"; "delta"; "gamma_C8colony"];
species_names_M1 = ["CEPI";"C8_E";"C8_I";"CEPR";"CEPRstar"; "OACP"];

Y_0 = [CEPI;C8_E;C8_I;CEPR;CEPRstar; OACP];

[T, Ys] =  ode23(@(t,y) dgdt_M1(t,y,parameters_M1), tspan, Y_0);

[T_sens,Ys_sens,DYDP] = sens_sys('dgdt_M1_sens',tspan,Y_0, [] ,parameters_M1);

%plot the result of sensitivity analysis 
% for j=1:1:length(parameters_M1)
% 
%     for i=1:1:length(Y_0)
%        
%     figure(100)
%     subplot(6,3,j)
%     plot(T_sens, DYDP(:,i,j), 'LineWidth', 2)
%     parameter_plot_M1 = parameters_M1_names(j)
%     title(parameter_plot_M1)
%     hold on
%     
%     figure(101)
%     subplot(6,4,j)
%     plot(T_sens, DYDP(:,i,j).*(parameters_M1(j)./Ys_sens(:,i)), 'LineWidth', 2)
%     parameter_plot_M1 = parameters_M1_names(j)
%     title(parameter_plot_M1)
%     hold on
% 
%     end
% 
%         
% end
% legend(species_names_M1)

figure(1)
subplot(2,3,1)
plot(T,Ys(:,1));
title('CEPI')
hold on
subplot(2,3,2)
plot(T,Ys(:,2));
title('C8_E;')
hold on
subplot(2,3,3)
plot(T,Ys(:,3));
title('C8_I')
hold on
subplot(2,3,4)
plot(T,Ys(:,4));
title('CEPR')
hold on
subplot(2,3,5)
plot(T,Ys(:,5));
title('CEPRstar')
hold on
subplot(2,3,6)
plot(T,Ys(:,6));
title('OACP')
hold on

sgtitle('Concentration of the different species in the core model')



%% complete model (M2) - Simulation

%close all

parameters_M2 = [v_max_CepI; km_CepI; w; k1; alpha_OACP; k_on_CepR_star; k_off_CepR_star; beta1; beta2; k3; n; alpha_C8; delta;gamma_C8colony;gamma_C6colony; k_on_CCIR_star; k_off_CCIR_star; k2; m; v_max_CciI; km_CccI; alpha_C6; k4; k5 ];
parameters_M2_names = ["v_max_CepI"; "km_CepI"; "w"; "k1"; "alpha_OACP"; "k_on_CepR_star"; "k_off_CepR_star"; "beta1"; "beta2"; "k3"; "n"; "alpha_C8"; "delta";"gamma_C8colony";"gamma_C6colony"; "k_on_CCIR_star"; "k_off_CCIR_star"; "k2"; "m"; "v_max_CciI"; "km_CccI"; "alpha_C6"; "k4"; "k5"];;
species_names_M2 = ["CEPI";"C8_E";"C8_I";"CEPR";"CEPRstar"; "OACP";"CCII"; "C6_E"; "C6_I"; "CCIR"; "CCIRstar"];

Y2_0 = [CEPI;C8_E;C8_I;CEPR;CEPRstar; OACP; CCII; C6_E; C6_I; CCIR; CCIRstar; ];

[T2, Y2s] =  ode23(@(t,y) dgdt_M2(t,y,parameters_M2), tspan, Y2_0);

[T2_sens,Y2s_sens,DYDP2] = sens_sys('dgdt_M2_sens',tspan,Y2_0, [] ,parameters_M2);

%plot the result of sensitivity analysis 
% for j=1:1:length(parameters_M2_names)
% 
%     for i=1:1:length(species_names_M2)
%        
%     figure(102)
%     subplot(6,4,j)
%     plot(T2_sens, DYDP2(:,i,j), 'LineWidth', 2)
%     parameter_plot_M2 = parameters_M2_names(j)
%     title(parameter_plot_M2)
%     hold on
%     
%     figure(103)
%     subplot(6,4,j)
%     plot(T2_sens, DYDP2(:,i,j).*(parameters_M2(j)./Y2s_sens(:,i)), 'LineWidth', 2)
%     parameter_plot_M2 = parameters_M2_names(j)
%     title(parameter_plot_M2)
%     hold on
% 
%     end
% 
%     
% end
% legend(species_names_M2)
% print -painters -depsc Sensitivity_analysis_M2.pdf

figure(2)
subplot(2,3,1)
plot(T2,Y2s(:,1));
title('CEPI')
hold on
subplot(2,3,2)
plot(T2,Y2s(:,2));
title('C8_E;')
hold on
subplot(2,3,3)
plot(T2,Y2s(:,3));
title('C8_I')
hold on
subplot(2,3,4)
plot(T2,Y2s(:,4));
title('CEPR')
hold on
subplot(2,3,5)
plot(T2,Y2s(:,5));
title('CEPRstar')
hold on
subplot(2,3,6)
plot(T2,Y2s(:,6));
title('OACP')
hold on

sgtitle('Concentration of the different species in the complete model')


figure(3)
subplot(6,2,1)
plot(T2,Y2s(:,1));
title('CEPI (complete model)')
hold on
subplot(6,2,2)
plot(T2,Y2s(:,2));
title('C8_E (complete model)')
hold on
subplot(6,2,3)
plot(T2,Y2s(:,3));
title('C8_I (complete model)')
hold on
subplot(6,2,4)
plot(T2,Y2s(:,4));
title('CEPR (complete model)')
hold on
subplot(6,2,5)
plot(T2,Y2s(:,5));
title('CEPRstar (complete model)')
hold on
subplot(6,2,6)
plot(T2,Y2s(:,6));
title('OACP (complete model)')
hold on
subplot(6,2,7)
plot(T2,Y2s(:,7));
title('CCII (complete model)')
hold on
subplot(6,2,8)
plot(T2,Y2s(:,8));
title('CCIR (complete model)')
hold on
subplot(6,2,9)
plot(T2,Y2s(:,9));
title('C6_E (complete model)')
hold on
subplot(6,2,10)
plot(T2,Y2s(:,10));
title('C6_I (complete model)')
hold on
subplot(6,2,11)
plot(T2,Y2s(:,11));
title('CCIRstar (complete model)')
hold on

sgtitle('Concentration of the different species in the complete model')





%%%% plot the comparison between the two models %%%%

figure(4)
subplot(2,3,6)
plot(T,Ys(:,1), '*');
title('CepI;')
hold on
subplot(2,3,2)
plot(T,Ys(:,2), '*');
title('C8_E;')
hold on
subplot(2,3,3)
plot(T,Ys(:,3), '*');
title('C8_I')
hold on
subplot(2,3,4)
plot(T,Ys(:,5)./Ys(:,4), '*');
ylabel('Concentration (a.u.)');
xlabel('Time (minutes)');
title('CEPR')
hold on
subplot(2,3,5)
plot(T,Ys(:,5), '*');
title('CEPRstar')
hold on
subplot(2,3,1)
plot(T,Ys(:,6), '*');
title('OACP')
hold on
subplot(2,3,1)
plot(T2,Y2s(:,6), 'o');
title('OACP')
ylabel('Concentration (nM)');
xlabel('Time (minutes)');
hold on
subplot(2,3,2)
plot(T2,Y2s(:,2), 'o');
legend('core','complete')
hold on
subplot(2,3,3)
plot(T2,Y2s(:,3), 'o');
legend('core','complete')
hold on
subplot(2,3,4)
plot(T2,Y2s(:,5)./Y2s(:,4), 'o');
legend('core','complete')
hold on
subplot(2,3,5)
plot(T2,Y2s(:,5), 'o');
legend('core','complete')
hold on
subplot(2,3,6)
plot(T2,Y2s(:,1), 'o');
title('CEPI')
legend('core','complete')


print -painters -depsc Steady_state.pdf


sgtitle('Comparison between the concentration of the different species in the core and in the complete model')






% plot values relative to steady state values


figure(5)
subplot(3,2,1)
plot(T,Ys(:,1)./Ys(end,1), 'LineWidth', 2);
title('CEPI (core model)')
xlabel('Time (min)','FontSize',16)
ylabel('Relative concentration','FontSize',16)
hold on
subplot(3,2,3)
plot(T,Ys(:,4)./Ys(end,4), 'LineWidth', 2);
title('CEPR (core model)')
xlabel('Time (min)','FontSize',16)
ylabel('Relative concentration','FontSize',16)
hold on
subplot(3,2,5)
plot(T,Ys(:,5)./Ys(end,5), 'LineWidth', 2);
title('CEPRstar (core model)')
xlabel('Time (min)','FontSize',16)
ylabel('Relative concentration','FontSize',16)
hold on
subplot(3,2,2)
plot(T2,Y2s(:,1)./Y2s(end,1), 'LineWidth', 2);
title('CEPI (complete model)')
xlabel('Time (min)','FontSize',16)
ylabel('Relative concentration','FontSize',16)
hold on
subplot(3,2,4)
plot(T2,Y2s(:,4)./Y2s(end,4), 'LineWidth', 2);
title('CEPR (complete model)')
xlabel('Time (min)','FontSize',16)
ylabel('Relative concentration','FontSize',16)
hold on
subplot(3,2,6)
plot(T2,Y2s(:,5)./Y2s(end,5), 'LineWidth', 2);
title('CEPRstar (complete model)', 'LineWidth', 2)
xlabel('Time (min)','FontSize',16)
ylabel('Relative concentration','FontSize',16)
hold on
leg = legend(num2str(C8_I));
hold on

%store sd and mean

std_cepI_core = std(Ys(:,1)./Ys(end,1))
std_cepI_complete = std(Y2s(:,1)./Y2s(end,1))


% subplot(4,2,7)
% plot(T2,Y2s(:,7)./Y2s(end,7), 'LineWidth', 2);
% title('CCII (complete model)', 'LineWidth', 2)
% hold on
% subplot(4,2,8)
% plot(T2,Y2s(:,10)./Y2s(end,10), 'LineWidth', 2);
% title('CCIR (complete model)', 'LineWidth', 2)
% hold on
% subplot(4,2,6)
% plot(T2,Y2s(:,11)./Y2s(end,11), 'LineWidth', 2);
% title('CCIR* (complete model)', 'LineWidth', 2)
% hold on

sgtitle('Relative concentration of the main species in the model')

% plot only cepi, cepr and cepr* concentrations

figure(6)
subplot(1,3,1)
plot(T,Ys(:,1)./Ys(end,1), 'LineWidth', 2);
title('\fontsize{22}CepI')
xlabel('Time (min)','FontSize',22)
ylabel('Relative concentration','FontSize',22)
hold on
subplot(1,3,2)
plot(T,Ys(:,4)./Ys(end,4), 'LineWidth', 2);
title('\fontsize{22}CepR')
xlabel('Time (min)','FontSize',22)
ylabel('Relative concentration','FontSize',22)
hold on
subplot(1,3,3)
plot(T,Ys(:,5)./Ys(end,5), 'LineWidth', 2);
title('\fontsize{22}CepR^*')
xlabel('Time (min)', 'FontSize',22)
ylabel('Relative concentration','FontSize',22)
hold on

figure(7)
subplot(1,3,1)
plot(T2,Y2s(:,1)./Y2s(end,1), 'LineWidth', 2);
title('\fontsize{22}CepI')
xlabel('Time (min)','FontSize',22)
ylabel('Relative concentration','FontSize',22)
hold on
subplot(1,3,2)
plot(T2,Y2s(:,4)./Y2s(end,4), 'LineWidth', 2);
title('\fontsize{22}CepR')
xlabel('Time (min)','FontSize',22)
ylabel('Relative concentration','FontSize',22)
hold on
subplot(1,3,3)
plot(T2,Y2s(:,5)./Y2s(end,5), 'LineWidth', 2);
title('\fontsize{22}CepR^*')
xlabel('Time (min)','FontSize',22)
ylabel('Relative concentration','FontSize',22)
hold on



%
% 
% [r,c,v] = find(Ys(:,1)./Ys(end,1)<0.51 & Ys(:,1)./Ys(end,1) > 0.5);
% RT_1 = T(r(1))
% [r2,c2,v2] = find(Y2s(:,1)./Y2s(end,1)<0.51 & Y2s(:,1)./Y2s(end,1) > 0.5)  ;
% RT_2 = T2(r2(1))
% RT_2./RT_1

end

end



