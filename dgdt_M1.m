function Ys = dgdt_M1(t, Y, parameters)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

CEPI = Y(1);
C8_E = Y(2);
C8_I = Y(3);
CEPR = Y(4);
CEPRstar = Y(5);
OACP = Y(6);


v_max_CepI = parameters(1);
km_CepI = parameters(2);
w = parameters(3);
k1 = parameters(4);
alpha_OACP = parameters(5);
k_on_CepR_star = parameters(6);
k_off_CepR_star = parameters(7);
beta1 = parameters(8);
k3 = parameters(9);
n = parameters(10);
alpha_C8 = parameters(11);
delta = parameters(12);
gamma_C8colony = parameters(13);





CEPI_dot = -delta*CEPI  + beta1*(CEPRstar)^w/(k1 + CEPRstar^w) ;

OACP_dot = alpha_OACP - delta*OACP - ((v_max_CepI*OACP)/(km_CepI+OACP))*(OACP * CEPI) ;
 
CepR_dot = -delta*CEPR - k_on_CepR_star*(CEPR * C8_I) + k_off_CepR_star*CEPRstar + beta1/(1 + (CEPRstar/k3)^n) ;

C8_E_dot =  -delta*C8_E + alpha_C8*C8_I + gamma_C8colony -  alpha_C8*C8_E   ;

CEPRstar_dot = k_on_CepR_star*(CEPR * C8_I) - k_off_CepR_star*(CEPRstar) - delta*CEPRstar;

C8_I_dot = -delta*C8_I -k_on_CepR_star*(CEPR * C8_I) + k_off_CepR_star*(CEPRstar) + ((v_max_CepI*OACP)/(km_CepI+OACP))*(OACP * CEPI) - alpha_C8*C8_I + alpha_C8*C8_E;

Ys = [CEPI_dot; C8_E_dot; C8_I_dot; CepR_dot;CEPRstar_dot;OACP_dot] ;

end

