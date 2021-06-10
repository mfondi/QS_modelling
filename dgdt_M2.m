function Ys = dgdt_M2(t, Y, parameters)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

CEPI = Y(1);
C8_E = Y(2);
C8_I = Y(3);
CEPR = Y(4);
CEPRstar = Y(5);
OACP = Y(6);
CCII = Y(7);
C6_E = Y(8);
C6_I = Y(9);
CCIR = Y(10);
CCIRstar = Y(11);



v_max_CepI = parameters(1);
km_CepI = parameters(2);
w = parameters(3);
k1 = parameters(4);
alpha_OACP = parameters(5);
k_on_CepR_star = parameters(6);
k_off_CepR_star = parameters(7);
beta1 = parameters(8);
beta2 = parameters(9);
k3 = parameters(10);
n = parameters(11);
alpha_C8 = parameters(12);
delta = parameters(13);
gamma_c8colony = parameters(14);  
gamma_c6colony = parameters(15);  
k_on_CCIR_star = parameters(16);
k_off_CCIR_star = parameters(17);
k2 = parameters(18);
m = parameters(19);
v_max_CciI = parameters(20);
km_CccI = parameters(21);
alpha_C6 = parameters(22);
k4 = parameters(23);
k5 = parameters(24);



%% M1 (core model) part
CEPI_dot = -delta*CEPI  + beta2/(1+ (CCIRstar/k2)^m)+ beta1*(CEPRstar)^w/(k1 + CEPRstar^w)  ;

OACP_dot = alpha_OACP - delta*OACP - ((v_max_CepI*OACP)/(km_CepI+OACP))*(OACP * CEPI) - ((v_max_CciI*OACP)/(km_CccI+OACP))*(OACP*CCII) ;
% 
CepR_dot = -delta*CEPR - k_on_CepR_star*(CEPR * C8_I) + k_off_CepR_star*CEPRstar + beta1/(1 + (CEPRstar/k3)^n) ;
%
%CepR_dot = 0;

C8_E_dot =  - delta*C8_E + alpha_C8*C8_I + gamma_c8colony -  alpha_C8*C8_E  ;

CEPRstar_dot = k_on_CepR_star*(CEPR * C8_I) - k_off_CepR_star*(CEPRstar) - delta*CEPRstar ;

%CEPRstar_dot = 0

C8_I_dot = -delta*C8_I -k_on_CepR_star*(CEPR * C8_I) + k_off_CepR_star*(CEPRstar) + ((v_max_CepI*OACP)/(km_CepI+OACP))*(OACP * CEPI)  - alpha_C8*C8_I + alpha_C8*C8_E;

%% M2 (complete model) part

CCII_dot = beta2 - delta*CCII ;

CCIR_dot = -delta*CCIR - k_on_CCIR_star*(CCIR * C6_I) + k_off_CCIR_star*CCIRstar + beta2/(1+ (CCIRstar/k4)^m) + beta2*(CEPRstar)^w/(k5 + CEPRstar^w) ; 

C6_I_dot = -delta*C6_I + ((v_max_CciI*OACP)/(km_CccI+OACP))*(OACP * CCII)  - alpha_C6*C6_I + alpha_C6*C6_E; ;

CCIRstar_dot = k_on_CCIR_star*(CCIR * C6_I) - k_off_CCIR_star*(CCIRstar) - delta*CCIRstar ;

C6_E_dot =  - delta*C6_E  + gamma_c6colony  -  alpha_C6*C6_E + alpha_C6*C6_I   ;

% figure(100)
% plot(C8_E_dot,C8_E )
% hold on

Ys = [CEPI_dot; C8_E_dot; C8_I_dot; CepR_dot; CEPRstar_dot; OACP_dot; CCII_dot; C6_E_dot;  C6_I_dot; CCIR_dot;  CCIRstar_dot;] ;


end

