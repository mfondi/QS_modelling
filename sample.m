%Program to simulate the erythroid-myeloid switch
%parameters in a structure
struct.alpha1=0.001;
struct.alpha2=0.25;
struct.beta1=0.001;
struct.beta2=0.25;
struct.beta3=1;
struct.gamma1=0.01;
struct.delta1=5;
struct.delta2=0.25;
struct.epsilon1=5;
struct.epsilon2=0.25;
struct.epsilon3=1;
struct.gamma2=0.01;
%set maximum tspan
tmax=2000;
%call ode45 to solve the system of differential equations
%first gene is GATA-1 (initial condition = 0)
%second gene is PU.1 (initial condition = 0)
%set the options for ode45
optionsODE=odeset('AbsTol',[1e-8 1e-8],'RelTol',1e-8);
[T Y]=ode45('bistableSwitch',[0 tmax],[0 0],optionsODE,struct);
GATA1=Y(:,1);
PU1=Y(:,2);
%plot the gene profiles versus time
figure(1)
plot(T,GATA1,'-b',T,PU1,'-r');
xlabel('time');
ylabel('expression level');
legend('GATA-1','PU.1');