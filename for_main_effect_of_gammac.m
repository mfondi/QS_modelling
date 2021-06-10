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
tspan = [0 400];

%% core model (M1) - Simulation
fprintf('\n\nrunning simulation for different values of gamma_colony\n');

%l = [0; .0001; .001; .01; .02; .03; 1; 2; 3; 10; 100];
l = [ 0.00015;.0015; .0075; .015; .03; .15; 1; 2; 3;10;];


colors=jet(length(l));

for i = 1:1:length(l)
    
gamma_C6colony = l(i);
%gamma_C8colony = l(i);

%parameters_M1 = [v_max_CepI; km_CepI; w; k1; alpha_OACP; k_on_CepR_star; k_off_CepR_star; beta; k3; n; alpha_C8; delta_HSL; delta_CepR; delta_CepI; delta_CepRstar; delta_OACP; mu;gamma];
parameters_M1 = [v_max_CepI; km_CepI; w; k1; alpha_OACP; k_on_CepR_star; k_off_CepR_star; beta1; k3; n; alpha_C8; delta; gamma_C8colony];

Y_0 = [CEPI;C8_E;C8_I;CEPR;CEPRstar; OACP];

[T, Ys] =  ode23(@(t,y) dgdt_M1(t,y,parameters_M1), tspan, Y_0);
% 


sgtitle('Concentration of the different species in the core model')

%% full model %%


parameters_M2 = [v_max_CepI; km_CepI; w; k1; alpha_OACP; k_on_CepR_star; k_off_CepR_star; beta1; beta2; k3; n; alpha_C8; delta;gamma_C8colony;gamma_C6colony; k_on_CCIR_star; k_off_CCIR_star; k2; m; v_max_CciI; km_CccI; alpha_C6; k4; k5 ];

Y2_0 = [CEPI;C8_E;C8_I;CEPR;CEPRstar; OACP; CCII; C6_E; C6_I; CCIR; CCIRstar; ];

[T2, Y2s] =  ode45(@(t,y) dgdt_M2(t,y,parameters_M2), tspan, Y2_0);


% plot only cepi, cepr and cepr* concentrations


figure(1)
subplot(1,3,1)
a = zeros(2);
plot(T2,Y2s(:,1)./Y2s(end,1), 'LineWidth', 2, 'Color', colors(i,:));
hold on
xlabel('Time (min)')
ylabel('Relative concentration')
plot(T,Ys(:,1)./Ys(end,1),'--', 'LineWidth', 2, 'Color', colors(i,:),'HandleVisibility','off');
%leg = legend(num2str(l));
%title(leg, 'k_1 (nM)');
hold on

figure(1)
subplot(1,3,2)
plot(T2,Y2s(:,4)./Y2s(end,4), 'LineWidth', 2, 'Color', colors(i,:));
hold on
plot(T,Ys(:,4)./Ys(end,4),'--', 'LineWidth', 2, 'Color', colors(i,:),'HandleVisibility','off');
xlabel('Time (min)')
ylabel('Relative concentration')
%leg = legend(num2str(l));
%title(leg, 'k_1 (nM)');
hold on

figure(1)
subplot(1,3,3)
plot(T2,Y2s(:,5)./Y2s(end,5), 'LineWidth', 2, 'Color', colors(i,:));
hold on
plot(T,Ys(:,5)./Ys(end,5),'--', 'LineWidth', 2, 'Color', colors(i,:),'HandleVisibility','off');
xlabel('Time (min)')
ylabel('Relative concentration')
%leg = legend(num2str(l));
%title(leg, 'k_1 (nM)');
hold on


end

saveas(figure(1),'../../text/gamma6colony_values_.pdf');


for i = 1:1:length(l)
    
%gamma_C6colony = l(i);
gamma_C8colony = l(i);

%parameters_M1 = [v_max_CepI; km_CepI; w; k1; alpha_OACP; k_on_CepR_star; k_off_CepR_star; beta; k3; n; alpha_C8; delta_HSL; delta_CepR; delta_CepI; delta_CepRstar; delta_OACP; mu;gamma];
parameters_M1 = [v_max_CepI; km_CepI; w; k1; alpha_OACP; k_on_CepR_star; k_off_CepR_star; beta1; k3; n; alpha_C8; delta; gamma_C8colony];

Y_0 = [CEPI;C8_E;C8_I;CEPR;CEPRstar; OACP];

[T, Ys] =  ode23(@(t,y) dgdt_M1(t,y,parameters_M1), tspan, Y_0);
% 


sgtitle('Concentration of the different species in the core model')

%% full model %%


parameters_M2 = [v_max_CepI; km_CepI; w; k1; alpha_OACP; k_on_CepR_star; k_off_CepR_star; beta1; beta2; k3; n; alpha_C8; delta;gamma_C8colony;gamma_C6colony; k_on_CCIR_star; k_off_CCIR_star; k2; m; v_max_CciI; km_CccI; alpha_C6; k4; k5 ];

Y2_0 = [CEPI;C8_E;C8_I;CEPR;CEPRstar; OACP; CCII; C6_E; C6_I; CCIR; CCIRstar; ];

[T2, Y2s] =  ode45(@(t,y) dgdt_M2(t,y,parameters_M2), tspan, Y2_0);


% plot only cepi, cepr and cepr* concentrations


figure(2)
subplot(1,3,1)
a = zeros(2);
plot(T2,Y2s(:,1)./Y2s(end,1), 'LineWidth', 2, 'Color', colors(i,:));
hold on
xlabel('Time (min)')
ylabel('Relative concentration')

figure(3)
subplot(1,3,1)
plot(T,Ys(:,1)./Ys(end,1),'--', 'LineWidth', 2, 'Color', colors(i,:),'HandleVisibility','off');
%leg = legend(num2str(l));
%title(leg, 'k_1 (nM)');
hold on

figure(2)
subplot(1,3,2)
plot(T2,Y2s(:,4)./Y2s(end,4), 'LineWidth', 2, 'Color', colors(i,:));
hold on

figure(3)
subplot(1,3,2)
plot(T,Ys(:,4)./Ys(end,4),'--', 'LineWidth', 2, 'Color', colors(i,:),'HandleVisibility','off');
xlabel('Time (min)')
ylabel('Relative concentration')
%leg = legend(num2str(l));
%title(leg, 'k_1 (nM)');
hold on

figure(2)
subplot(1,3,3)
plot(T2,Y2s(:,5)./Y2s(end,5), 'LineWidth', 2, 'Color', colors(i,:));
hold on

figure(3)
subplot(1,3,3)
plot(T,Ys(:,5)./Ys(end,5),'--', 'LineWidth', 2, 'Color', colors(i,:),'HandleVisibility','off');
xlabel('Time (min)')
ylabel('Relative concentration')
%leg = legend(num2str(l));
%title(leg, 'k_1 (nM)');
hold on

figure(9)
plot(gamma_C8colony, Y2s(end,1), 'LineWidth', 2);
hold on

end
leg = legend(num2str(l));

saveas(figure(2),'../../text/gamma8colony_values_.pdf');
