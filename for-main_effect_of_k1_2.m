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
tspan = [0 1000];

%% core model (M1) - Simulation

l = [  01; .1; 1; 10; 100;1000 ];


colors=jet(length(l));

%% k1
fprintf('\n\nrunning simulation for different values of k_1\n');
RT_ratio_k1 = [];
RT_ratio_k2 = [];
RT_ratio_k3 = [];
RT_ratio_k4 = [];
RT_ratio_k5 = [];

k1_value = [];
k2_value = [];


for i = 1:1:length(l)
    
    for j = 1:1:length(l)

k3= l(i)
k4= l(j)
%parameters_M1 = [v_max_CepI; km_CepI; w; k1; alpha_OACP; k_on_CepR_star; k_off_CepR_star; beta; k3; n; alpha_C8; delta_HSL; delta_CepR; delta_CepI; delta_CepRstar; delta_OACP; mu;gamma];
parameters_M1 = [v_max_CepI; km_CepI; w; k1; alpha_OACP; k_on_CepR_star; k_off_CepR_star; beta1; k3; n; alpha_C8; delta; gamma_C8colony];

Y_0 = [CEPI;C8_E;C8_I;CEPR;CEPRstar; OACP];

[T, Ys] =  ode23(@(t,y) dgdt_M1(t,y,parameters_M1), tspan, Y_0);



sgtitle('Concentration of the different species in the core model')

%% full model %%


parameters_M2 = [v_max_CepI; km_CepI; w; k1; alpha_OACP; k_on_CepR_star; k_off_CepR_star; beta1; beta2; k3; n; alpha_C8; delta;gamma_C8colony;gamma_C6colony; k_on_CCIR_star; k_off_CCIR_star; k2; m; v_max_CciI; km_CccI; alpha_C6; k4; k5 ];

Y2_0 = [CEPI;C8_E;C8_I;CEPR;CEPRstar; OACP; CCII; C6_E; C6_I; CCIR; CCIRstar; ];

[T2, Y2s] =  ode45(@(t,y) dgdt_M2(t,y,parameters_M2), tspan, Y2_0);

% plot only cepi, cepr and cepr* concentrations


figure(8)
subplot(5,3,1)
a = zeros(2);
plot(T2,Y2s(:,1)./Y2s(end,1), 'LineWidth', 2, 'Color', colors(i,:));
hold on
plot(T,Ys(:,1)./Ys(end,1),'--', 'LineWidth', 2, 'Color', colors(i,:),'HandleVisibility','off');
xlabel('Time (min)')
ylabel('Relative concentration')

hold on

figure(8)
subplot(5,3,2)
plot(T2,Y2s(:,4)./Y2s(end,4), 'LineWidth', 2, 'Color', colors(i,:));
hold on
plot(T,Ys(:,4)./Ys(end,4),'--', 'LineWidth', 2, 'Color', colors(i,:),'HandleVisibility','off');
xlabel('Time (min)')
ylabel('Relative concentration')

hold on

figure(8)
subplot(5,3,3)
plot(T2,Y2s(:,5)./Y2s(end,5), 'LineWidth', 2, 'Color', colors(i,:));
plot(T,Ys(:,5)./Ys(end,5),'--', 'LineWidth', 2, 'Color', colors(i,:),'HandleVisibility','off');
xlabel('Time (min)')
ylabel('Relative concentration')

hold on


figure(9)
subplot(1,5,1)
a = zeros(2);
plot(T2,Y2s(:,1)./Y2s(end,1), 'LineWidth', 2, 'Color', colors(i,:));
hold on

xlabel('Time (min)')
ylabel('Relative concentration')
plot(T,Ys(:,1)./Ys(end,1),'--', 'LineWidth', 2, 'Color', colors(i,:),'HandleVisibility','off');


hold on

% compute response time

[r,c,v] = find(Ys(:,1)./Ys(end,1)<0.51 & Ys(:,1)./Ys(end,1) > 0.5)  
RT_1 = T(r(1))
[r2,c2,v2] = find(Y2s(:,1)./Y2s(end,1)<0.51 & Y2s(:,1)./Y2s(end,1) > 0.5)  
RT_2 = T(r2(1))
RT_ratio_k1(end+1) = RT_2./RT_1;
k2_value(end+1) = k3;
k1_value(end+1) = k4;

    end
    
end
figure
scatter_data = [k1_value' k2_value' RT_ratio_k1']
scatter(k1_value', k2_value',[], RT_ratio_k1)

%plot response time
figure(100)
plot(l,RT_ratio_k1, 'o');
hold on


%% k2
fprintf('\n\nrunning simulation for different values of k_2\n');

for i = 1:1:length(l)
    
k2= l(i);
%parameters_M1 = [v_max_CepI; km_CepI; w; k1; alpha_OACP; k_on_CepR_star; k_off_CepR_star; beta; k3; n; alpha_C8; delta_HSL; delta_CepR; delta_CepI; delta_CepRstar; delta_OACP; mu;gamma];
parameters_M1 = [v_max_CepI; km_CepI; w; k1; alpha_OACP; k_on_CepR_star; k_off_CepR_star; beta1; k3; n; alpha_C8; delta; gamma_C8colony];

Y_0 = [CEPI;C8_E;C8_I;CEPR;CEPRstar;OACP];

[T, Ys] =  ode23(@(t,y) dgdt_M1(t,y,parameters_M1), tspan, Y_0);


sgtitle('Concentration of the different species in the core model')

%% full model %%


parameters_M2 = [v_max_CepI; km_CepI; w; k1; alpha_OACP; k_on_CepR_star; k_off_CepR_star; beta1; beta2; k3; n; alpha_C8; delta;gamma_C8colony;gamma_C6colony; k_on_CCIR_star; k_off_CCIR_star; k2; m; v_max_CciI; km_CccI; alpha_C6; k4; k5 ];

Y2_0 = [CEPI;C8_E;C8_I;CEPR;CEPRstar; OACP; CCII; C6_E; C6_I; CCIR; CCIRstar; ];

[T2, Y2s] =  ode45(@(t,y) dgdt_M2(t,y,parameters_M2), tspan, Y2_0);


% plot only cepi, cepr and cepr* concentrations


figure(8)
subplot(5,3,4)
plot(T2,Y2s(:,1)./Y2s(end,1), 'LineWidth', 2, 'Color', colors(i,:));
hold on
xlabel('Time (min)')
ylabel('Relative concentration')
plot(T,Ys(:,1)./Ys(end,1),'--', 'LineWidth', 2, 'Color', colors(i,:),'HandleVisibility','off');


hold on

figure(8)
subplot(5,3,5)
plot(T2,Y2s(:,4)./Y2s(end,4), 'LineWidth', 2, 'Color', colors(i,:));
hold on

plot(T,Ys(:,4)./Ys(end,4),'--', 'LineWidth', 2, 'Color', colors(i,:),'HandleVisibility','off');
xlabel('Time (min)')
ylabel('Relative concentration')


hold on

figure(8)
subplot(5,3,6)
plot(T2,Y2s(:,5)./Y2s(end,5), 'LineWidth', 2, 'Color', colors(i,:));
hold on

plot(T,Ys(:,5)./Ys(end,5),'--', 'LineWidth', 2, 'Color', colors(i,:),'HandleVisibility','off');
xlabel('Time (min)')
ylabel('Relative concentration')


hold on


figure(9)
subplot(1,5,2)
a = zeros(2);
plot(T2,Y2s(:,1)./Y2s(end,1), 'LineWidth', 2, 'Color', colors(i,:));
hold on

xlabel('Time (min)')
ylabel('Relative concentration')
plot(T,Ys(:,1)./Ys(end,1),'--', 'LineWidth', 2, 'Color', colors(i,:),'HandleVisibility','off');


hold on


% compute response time

[r,c,v] = find(Ys(:,1)./Ys(end,1)<0.51 & Ys(:,1)./Ys(end,1) > 0.5)  
RT_1 = T(r(1))
[r2,c2,v2] = find(Y2s(:,1)./Y2s(end,1)<0.51 & Y2s(:,1)./Y2s(end,1) > 0.5)  
RT_2 = T(r2(1))
RT_ratio_k2(end+1) = RT_2./RT_1;


end

%plot response time
figure(100)
plot(l,RT_ratio_k2, 'o');
hold on


%% k3
fprintf('\n\nrunning simulation for different values of k_3\n');

for i = 1:1:length(l)
    
k3= l(i);
%parameters_M1 = [v_max_CepI; km_CepI; w; k1; alpha_OACP; k_on_CepR_star; k_off_CepR_star; beta; k3; n; alpha_C8; delta_HSL; delta_CepR; delta_CepI; delta_CepRstar; delta_OACP; mu;gamma];
parameters_M1 = [v_max_CepI; km_CepI; w; k1; alpha_OACP; k_on_CepR_star; k_off_CepR_star; beta1; k3; n; alpha_C8; delta; gamma_C8colony];

Y_0 = [CEPI;C8_E;C8_I;CEPR;CEPRstar; OACP];

[T, Ys] =  ode23(@(t,y) dgdt_M1(t,y,parameters_M1), tspan, Y_0);


sgtitle('Concentration of the different species in the core model')

%% full model %%


parameters_M2 = [v_max_CepI; km_CepI; w; k1; alpha_OACP; k_on_CepR_star; k_off_CepR_star; beta1; beta2; k3; n; alpha_C8; delta;gamma_C8colony;gamma_C6colony; k_on_CCIR_star; k_off_CCIR_star; k2; m; v_max_CciI; km_CccI; alpha_C6; k4; k5 ];

Y2_0 = [CEPI;C8_E;C8_I;CEPR;CEPRstar; OACP; CCII; C6_E; C6_I; CCIR; CCIRstar; ];

[T2, Y2s] =  ode45(@(t,y) dgdt_M2(t,y,parameters_M2), tspan, Y2_0);


% plot only cepi, cepr and cepr* concentrations


figure(8)
subplot(5,3,7)
plot(T2,Y2s(:,1)./Y2s(end,1), 'LineWidth', 2, 'Color', colors(i,:));
hold on

xlabel('Time (min)')
ylabel('Relative concentration')
plot(T,Ys(:,1)./Ys(end,1),'--', 'LineWidth', 2, 'Color', colors(i,:),'HandleVisibility','off');


hold on

figure(8)
subplot(5,3,8)
plot(T2,Y2s(:,4)./Y2s(end,4), 'LineWidth', 2, 'Color', colors(i,:));
hold on

plot(T,Ys(:,4)./Ys(end,4),'--', 'LineWidth', 2, 'Color', colors(i,:),'HandleVisibility','off');
xlabel('Time (min)')
ylabel('Relative concentration')


hold on

figure(8)
subplot(5,3,9)
plot(T2,Y2s(:,5)./Y2s(end,5), 'LineWidth', 2, 'Color', colors(i,:));
hold on

plot(T,Ys(:,5)./Ys(end,5),'--', 'LineWidth', 2, 'Color', colors(i,:),'HandleVisibility','off');
xlabel('Time (min)')
ylabel('Relative concentration')


hold on


figure(9)
subplot(1,5,3)
a = zeros(2);
plot(T2,Y2s(:,1)./Y2s(end,1), 'LineWidth', 2, 'Color', colors(i,:));
hold on

xlabel('Time (min)')
ylabel('Relative concentration')
plot(T,Ys(:,1)./Ys(end,1),'--', 'LineWidth', 2, 'Color', colors(i,:),'HandleVisibility','off');


hold on

% compute response time

[r,c,v] = find(Ys(:,1)./Ys(end,1)<0.51 & Ys(:,1)./Ys(end,1) > 0.5)  
RT_1 = T(r(1))
[r2,c2,v2] = find(Y2s(:,1)./Y2s(end,1)<0.51 & Y2s(:,1)./Y2s(end,1) > 0.5)  
RT_2 = T(r2(1))
RT_ratio_k3(end+1) = RT_2./RT_1;


end

%plot response time
figure(100)
plot(l,RT_ratio_k2, 'o');
hold on


%% k4
fprintf('\n\nrunning simulation for different values of k_4\n');

for i = 1:1:length(l)
    
k4= l(i);
%parameters_M1 = [v_max_CepI; km_CepI; w; k1; alpha_OACP; k_on_CepR_star; k_off_CepR_star; beta; k3; n; alpha_C8; delta_HSL; delta_CepR; delta_CepI; delta_CepRstar; delta_OACP; mu;gamma];
parameters_M1 = [v_max_CepI; km_CepI; w; k1; alpha_OACP; k_on_CepR_star; k_off_CepR_star; beta1; k3; n; alpha_C8; delta; gamma_C8colony];

Y_0 = [CEPI;C8_E;C8_I;CEPR;CEPRstar; OACP];

[T, Ys] =  ode23(@(t,y) dgdt_M1(t,y,parameters_M1), tspan, Y_0);


sgtitle('Concentration of the different species in the core model')

%% full model %%


parameters_M2 = [v_max_CepI; km_CepI; w; k1; alpha_OACP; k_on_CepR_star; k_off_CepR_star; beta1; beta2; k3; n; alpha_C8; delta;gamma_C8colony;gamma_C6colony; k_on_CCIR_star; k_off_CCIR_star; k2; m; v_max_CciI; km_CccI; alpha_C6; k4; k5 ];

Y2_0 = [CEPI;C8_E;C8_I;CEPR;CEPRstar; OACP; CCII; C6_E; C6_I; CCIR; CCIRstar; ];

[T2, Y2s] =  ode45(@(t,y) dgdt_M2(t,y,parameters_M2), tspan, Y2_0);


% plot only cepi, cepr and cepr* concentrations


figure(8)
subplot(5,3,10)
plot(T2,Y2s(:,1)./Y2s(end,1), 'LineWidth', 2, 'Color', colors(i,:));
hold on

xlabel('Time (min)')
ylabel('Relative concentration')
plot(T,Ys(:,1)./Ys(end,1),'--', 'LineWidth', 2, 'Color', colors(i,:),'HandleVisibility','off');


hold on

figure(8)
subplot(5,3,11)
plot(T2,Y2s(:,4)./Y2s(end,4), 'LineWidth', 2, 'Color', colors(i,:));
hold on

plot(T,Ys(:,4)./Ys(end,4),'--', 'LineWidth', 2, 'Color', colors(i,:),'HandleVisibility','off');
xlabel('Time (min)')
ylabel('Relative concentration')


hold on

figure(8)
subplot(5,3,12)
plot(T2,Y2s(:,5)./Y2s(end,5), 'LineWidth', 2, 'Color', colors(i,:));
hold on

plot(T,Ys(:,5)./Ys(end,5),'--', 'LineWidth', 2, 'Color', colors(i,:),'HandleVisibility','off');
xlabel('Time (min)')
ylabel('Relative concentration')


hold on


figure(9)
subplot(1,5,4)
a = zeros(2);
plot(T2,Y2s(:,1)./Y2s(end,1), 'LineWidth', 2, 'Color', colors(i,:));
hold on

xlabel('Time (min)')
ylabel('Relative concentration')
plot(T,Ys(:,1)./Ys(end,1),'--', 'LineWidth', 2, 'Color', colors(i,:),'HandleVisibility','off');


hold on


% compute response time

[r,c,v] = find(Ys(:,1)./Ys(end,1)<0.51 & Ys(:,1)./Ys(end,1) > 0.5)  
RT_1 = T(r(1))
[r2,c2,v2] = find(Y2s(:,1)./Y2s(end,1)<0.51 & Y2s(:,1)./Y2s(end,1) > 0.5)  
RT_2 = T(r2(1))
RT_ratio_k4(end+1) = RT_2./RT_1;


end

%plot response time
figure(100)
plot(l,RT_ratio_k4, 'o');
hold on


%% k5
fprintf('\n\nrunning simulation for different values of k_5\n');

for i = 1:1:length(l)
    
k5= l(i);
%parameters_M1 = [v_max_CepI; km_CepI; w; k1; alpha_OACP; k_on_CepR_star; k_off_CepR_star; beta; k3; n; alpha_C8; delta_HSL; delta_CepR; delta_CepI; delta_CepRstar; delta_OACP; mu;gamma];
parameters_M1 = [v_max_CepI; km_CepI; w; k1; alpha_OACP; k_on_CepR_star; k_off_CepR_star; beta1; k3; n; alpha_C8; delta; gamma_C8colony];

Y_0 = [CEPI;C8_E;C8_I;CEPR;CEPRstar; OACP];

[T, Ys] =  ode23(@(t,y) dgdt_M1(t,y,parameters_M1), tspan, Y_0);


sgtitle('Concentration of the different species in the core model')

%% full model %%


parameters_M2 = [v_max_CepI; km_CepI; w; k1; alpha_OACP; k_on_CepR_star; k_off_CepR_star; beta1; beta2; k3; n; alpha_C8; delta;gamma_C8colony;gamma_C6colony; k_on_CCIR_star; k_off_CCIR_star; k2; m; v_max_CciI; km_CccI; alpha_C6; k4; k5 ];

Y2_0 = [CEPI;C8_E;C8_I;CEPR;CEPRstar; OACP; CCII; C6_E; C6_I; CCIR; CCIRstar; ];

[T2, Y2s] =  ode45(@(t,y) dgdt_M2(t,y,parameters_M2), tspan, Y2_0);


% plot only cepi, cepr and cepr* concentrations


figure(8)
subplot(5,3,13)
plot(T2,Y2s(:,1)./Y2s(end,1), 'LineWidth', 2, 'Color', colors(i,:));
hold on

xlabel('Time (min)')
ylabel('Relative concentration')
plot(T,Ys(:,1)./Ys(end,1),'--', 'LineWidth', 2, 'Color', colors(i,:),'HandleVisibility','off');


hold on

figure(8)
subplot(5,3,14)
plot(T2,Y2s(:,4)./Y2s(end,4), 'LineWidth', 2, 'Color', colors(i,:));
hold on

plot(T,Ys(:,4)./Ys(end,4),'--', 'LineWidth', 2, 'Color', colors(i,:),'HandleVisibility','off');
xlabel('Time (min)')
ylabel('Relative concentration')


hold on

figure(8)
subplot(5,3,15)
plot(T2,Y2s(:,5)./Y2s(end,5), 'LineWidth', 2, 'Color', colors(i,:));
hold on

plot(T,Ys(:,5)./Ys(end,5),'--', 'LineWidth', 2, 'Color', colors(i,:),'HandleVisibility','off');
xlabel('Time (min)')
ylabel('Relative concentration')


hold on


figure(9)
subplot(1,5,5)
a = zeros(2);
plot(T2,Y2s(:,1)./Y2s(end,1), 'LineWidth', 2, 'Color', colors(i,:));
hold on

xlabel('Time (min)')
ylabel('Relative concentration')
plot(T,Ys(:,1)./Ys(end,1),'--', 'LineWidth', 2, 'Color', colors(i,:),'HandleVisibility','off');

hold on


% compute response time

[r,c,v] = find(Ys(:,1)./Ys(end,1)<0.51 & Ys(:,1)./Ys(end,1) > 0.5)  
RT_1 = T(r(1))
[r2,c2,v2] = find(Y2s(:,1)./Y2s(end,1)<0.51 & Y2s(:,1)./Y2s(end,1) > 0.5)  
RT_2 = T(r2(1))
RT_ratio_k5(end+1) = RT_2./RT_1;


end

%plot response time
figure(100)
plot(l,RT_ratio_k5, 'o');
hold on

%leg = legend(num2str(l));
leg = legend( '01', '.1', '1', '10', '100', '1000')

title(leg, 'k_{1-5} values (nM)');

saveas(figure(8),'../../text/CepIRvsK_values.pdf');
saveas(figure(9),'../../text/CepIvsK_values.pdf');




