beta=1; gamma=1; delta=0.2; Kmfeed=20; kfeed=2; sigma=1; Lext=2;

dL/dt = beta*Lext*Y-gamma*L;
dY/dt = delta+kfeed*L^4/(Kmfeed+L^4)-sigma*Y

nullclineL=@(L)gamma.*L./(beta*Lext);
nullclineY=@(L)delta./sigma+kfeed.*L.^4./(sigma*(Kmfeed+L.^4));
L=1:.1:5;
figure(3)
plot(L,nullclineL(L),'b','LineWidth',2); hold on;
plot(L,nullclineY(L),'r','LineWidth',2);
