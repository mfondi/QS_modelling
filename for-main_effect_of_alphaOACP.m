%% in this codemwe test the effect of different values of parameters k (k1 to k5) on the resulting dynamics of the syste,

set(0,'DefaultFigureWindowStyle','docked')


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

gamma_C8colony = .015;      % rate of C8HSL production - Weber and Buceta 2013 - original value: 0.015
gamma_C6colony = .015;      % rate of C6HSL production - Weber and Buceta 2013 - original value: 0.015

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
tspan = [0 1200];

%% core model (M1) - Simulation

% l = [ 0.00001; .0001; .0005; .001; 0.002; 0.003; 0.004; .005; .01; .05; 0.1;];

l = [ .01; .015; .020; .025; 1;];

colors=jet(length(l));
RT_ratio_delta = [];
RT_1_array = [];
RT_2_array = []; 
RT_C8_1_array = [];
RT_C8_2_array = [];
RT_C8_ratio_delta = [];
RT_C6_1_array = [];
RT_C6_2_array = [];
RT_C6_ratio_delta = [];
steady_state_cepI_core = [];
steady_state_cepI_complete = [];
steady_state_C8_core = [];
steady_state_C8_complete = [];
AI_ratio_C6 = [];
AI_ratio_C8 = [];

%% k1
fprintf('\n\nrunning simulation for different values of k_1\n');

for i = 1:1:length(l)
    
alpha_OACP= l(i);
parameters_M1 = [v_max_CepI; km_CepI; w; k1; alpha_OACP; k_on_CepR_star; k_off_CepR_star; beta1; k3; n; alpha_C8; delta; gamma_C8colony];
parameters_M1_names = ["v_max_CepI"; "km_CepI"; "w"; "k1"; "alpha_OACP"; "k_on_CepR_star"; "k_off_CepR_star"; "beta1"; "k3"; "n"; "alpha_C8"; "delta"; "gamma_C8colony"];
species_names_M1 = ["CEPI";"C8_E";"C8_I";"CEPR";"CEPRstar"; "OACP"];

Y_0 = [CEPI;C8_E;C8_I;CEPR;CEPRstar; OACP];

[T, Ys] =  ode23(@(t,y) dgdt_M1(t,y,parameters_M1), tspan, Y_0);





%% full model %%


parameters_M2 = [v_max_CepI; km_CepI; w; k1; alpha_OACP; k_on_CepR_star; k_off_CepR_star; beta1; beta2; k3; n; alpha_C8; delta;gamma_C8colony;gamma_C6colony; k_on_CCIR_star; k_off_CCIR_star; k2; m; v_max_CciI; km_CccI; alpha_C6; k4; k5 ];
parameters_M2_names = ["v_max_CepI"; "km_CepI"; "w"; "k1"; "alpha_OACP"; "k_on_CepR_star"; "k_off_CepR_star"; "beta1"; "beta2"; "k3"; "n"; "alpha_C8"; "delta";"gamma_C8colony";"gamma_C6colony"; "k_on_CCIR_star"; "k_off_CCIR_star"; "k2"; "m"; "v_max_CciI"; "km_CccI"; "alpha_C6"; "k4"; "k5"];;
species_names_M2 = ["CEPI";"C8_E";"C8_I";"CEPR";"CEPRstar"; "OACP";"CCII"; "C6_E"; "C6_I"; "CCIR"; "CCIRstar"];

Y2_0 = [CEPI;C8_E;C8_I;CEPR;CEPRstar; OACP; CCII; C6_E; C6_I; CCIR; CCIRstar; ];

[T2, Y2s] =  ode23(@(t,y) dgdt_M2(t,y,parameters_M2), tspan, Y2_0);

% plot only cepi, cepr and cepr* concentrations



figure(3)
subplot(2,3,6)
plot(T,Ys(:,1));
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
plot(T,Ys(:,5)./Ys(:,4));
ylabel('Concentration (a.u.)');
xlabel('Time (minutes)');
title('CEPR')
hold on
subplot(2,3,5)
plot(T,Ys(:,5));
title('CEPRstar')
hold on
subplot(2,3,1)
plot(T,Ys(:,6));
title('OACP')
hold on
subplot(2,3,1)
plot(T2,Y2s(:,6));
title('OACP')
ylabel('Concentration (nM)');
xlabel('Time (minutes)');
hold on
subplot(2,3,2)
plot(T2,Y2s(:,2));
legend('core','complete')
hold on
subplot(2,3,3)
plot(T2,Y2s(:,3));
legend('core','complete')
hold on
subplot(2,3,4)
plot(T2,Y2s(:,5)./Y2s(:,4));
legend('core','complete')
hold on
subplot(2,3,5)
plot(T2,Y2s(:,5));
legend('core','complete')
hold on
subplot(2,3,6)
plot(T2,Y2s(:,1));
title('CEPI')
legend('core','complete')

figure(8)
subplot(1,3,1)
a = zeros(2);
plot(T2,Y2s(:,1)./Y2s(end,1), 'LineWidth', 2, 'Color', colors(i,:));
hold on
plot(T,Ys(:,1)./Ys(end,1),'--', 'LineWidth', 2, 'Color', colors(i,:),'HandleVisibility','off');
xlabel('Time (min)')
ylabel('Relative concentration')
hold on
subplot(1,3,2)
plot(T2,Y2s(:,4)./Y2s(end,4), 'LineWidth', 2, 'Color', colors(i,:));
hold on
plot(T,Ys(:,4)./Ys(end,4),'--', 'LineWidth', 2, 'Color', colors(i,:),'HandleVisibility','off');
xlabel('Time (min)')
ylabel('Relative concentration')
hold on
subplot(1,3,3)
plot(T2,Y2s(:,5)./Y2s(end,5), 'LineWidth', 2, 'Color', colors(i,:));
plot(T,Ys(:,5)./Ys(end,5),'--', 'LineWidth', 2, 'Color', colors(i,:),'HandleVisibility','off');
xlabel('Time (min)')
ylabel('Relative concentration')
hold on


figure(16)
a = zeros(2);
subplot(1,2,1)
plot(T2,Y2s(:,1)./Y2s(end,1), 'LineWidth', 2, 'Color', colors(i,:));
hold on
xlabel('Time (min)', 'FontSize',18)
ylabel('Relative concentration', 'FontSize',18);
set(gca, 'FontSize',14);
plot(T,Ys(:,1)./Ys(end,1),'--', 'LineWidth', 2, 'Color', colors(i,:),'HandleVisibility','off');
hold on


[r,c,v] = find(Ys(:,1)./Ys(end,1)<0.51 & Ys(:,1)./Ys(end,1) > 0.5)  ;
RT_1 = T(r(1));
RT_1_array(end+1) = RT_1;
[r2,c2,v2] = find(Y2s(:,1)./Y2s(end,1)<0.51 & Y2s(:,1)./Y2s(end,1) > 0.5)  
RT_2 = T2(r2(1));
RT_2_array(end+1) = RT_2;
RT_ratio_delta(end+1) = RT_1-RT_2;


% [r,c,v] = find(Ys(:,1)./Ys(end,1)<0.51 & Ys(:,1)./Ys(end,1) > 0.5);
% RT_1 = T(r(1))
% [r2,c2,v2] = find(Y2s(:,1)./Y2s(end,1)<0.51 & Y2s(:,1)./Y2s(end,1) > 0.5)  ;
% RT_2 = T(r2(1))

figure(10)
plot(T2,Y2s(:,2)./Y2s(end,2), 'LineWidth', 2, 'Color', colors(i,:)); 
hold on
plot(T,Ys(:,2)./Ys(end,2),'--', 'LineWidth', 2, 'Color', colors(i,:)); 

 [r,c,v] = find(Ys(:,2)./Ys(end,2)<0.51 & Ys(:,2)./Ys(end,2) > 0.5)  ;
 RT_C8_1 = T(r(1));
 RT_C8_1_array(end+1) = RT_C8_1;
 [r2,c2,v2] = find(Y2s(:,2)./Y2s(end,2)<0.51 & Y2s(:,2)./Y2s(end,2) > 0.5)  
 RT_C8_2 = T2(r2(1));
 RT_C8_2_array(end+1) = RT_C8_2;
 RT_C8_ratio_delta(end+1) = RT_C8_1 - RT_C8_2;

 



 [r,c,v] = find(Ys(:,3)./Ys(end,3)<0.51 & Ys(:,3)./Ys(end,3) > 0.5)  ;
 RT_C6_1 = T(r(1));
 RT_C6_1_array(end+1) = RT_C6_1;
 [r2,c2,v2] = find(Y2s(:,3)./Y2s(end,3)<0.51 & Y2s(:,3)./Y2s(end,3) > 0.5)  
 RT_C6_2 = T2(r2(1));
 RT_C6_2_array(end+1) = RT_C6_2;
 RT_C6_ratio_delta(end+1) = RT_C6_1 - RT_C6_2;

 
 steady_state_cepI_core(end+1)=Ys(end,1);
 steady_state_cepI_complete(end+1)=Y2s(end,1);
 steady_state_C8_core(end+1)=Ys(end,2);
 steady_state_C8_complete(end+1)=Y2s(end,2);
 
 figure(14)
 plot(T2,Y2s(:,2), 'LineWidth', 2, 'Color', colors(i,:)); 
 hold on
 figure(15)
 plot(T,Ys(:,2), 'LineWidth', 2, 'Color', colors(i,:)); 
 hold on
 
 %compute the ratio of steady-stae levels of AI in the two models 
 
 AI_ratio_C8(end+1) = Y2s(end,2)./Ys(end,2);
 AI_ratio_C6(end+1) = Y2s(end,3)./Ys(end,3);


end


leg = legend(num2str(l));

figure(16)
subplot(1,2,2);
plot(1:1:length(l), RT_ratio_delta, 'o-', 'LineWidth', 2, 'MarkerSize', 10);
hold on
plot(1:1:length(l), RT_C8_ratio_delta, 'o-', 'LineWidth', 2, 'MarkerSize', 10);
% hold on
% plot(1:1:length(l), RT_C6_ratio_delta, '*-', 'LineWidth', 2, 'MarkerSize', 10)
xlabel('Growth rate (h^-^1)', 'FontSize',18);
ylabel('Delta RT (min)', 'FontSize',18);
set(gca,'XTickLabel',l, 'FontSize',14);

title(leg, 'Delta values (h^-1)');

saveas(figure(16),'../../text/delta_m1vsm2.pdf');
 
figure(17)
plot(1:1:length(l), AI_ratio_C8, 'o-', 'LineWidth', 2, 'MarkerSize', 10);
xlabel('Growth rate (h^-^1)', 'FontSize',18);
ylabel('Steady state C8-HSL concentration ratio', 'FontSize',18);
set(gca,'XTickLabel',l, 'FontSize',14);

figure(18)
plot(1:1:length(l), AI_ratio_C6, 'o-', 'LineWidth', 2, 'MarkerSize', 10);
xlabel('Growth rate (h^-^1)', 'FontSize',18);
ylabel('Steady state C6-HSL concentration ratio', 'FontSize',18);
set(gca,'XTickLabel',l, 'FontSize',14);


RT_1
RT_2