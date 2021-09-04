% ModelIS   [ S, R_A^A,R_A^B,R_B^A,R_B^B,X ]

 y0=[rand(6,1);0]
 y0=y0./sum(y0);
% y0=[ swap*(yOut(end,1:6)') ;0];
% y0= [(yOut(end,1:6)');0];

T=250;
mu=1/15; %(immgration/emigration rate)
beta= [1,0.99,0.99]; %Infection rates;

chi=@(t) 0.5;

ts=linspace(0,10*T,1000);

smallImport=0.01;


filterArrivals=   [1;0;1;0;1;1];%DumbMixing
filterInfections= [1;1;0;1;0;1];
importation= [0.5,smallImport,smallImport,smallImport,smallImport,1]';
importation=  importation*mu/(filterArrivals'*importation);

betaMatrix= [beta(1),0,0,0,0,0;
             0,beta(2),beta(2),0,0,0;
             0,beta(2),beta(2),0,0,0;
             0,0,0,beta(3),beta(3),0;
             0,0,0,beta(3),beta(3),0;
             -beta(1),-beta(2),-beta(2),-beta(3),-beta(3),0];
         
g=1/10; %Ten day recovery without meds;
q=1/2.5; %treated recovery.
recovery= diag(-[q,g,q,q,g,0]);
recovery(6,:)=[q,g,q,q,g,0];
         
Deriv =@(t,V) [(filterArrivals.*importation + filterInfections.*betaMatrix*V(1:6)*V(6) + recovery*V(1:6)-mu*V(1:6)); V(6)];

[tOut,yOut] = ode45(Deriv,[0,200],y0);

% figure(1)
% plot(tOut,yOut);
% hold on
% plot(tOut,chi(tOut));

figure(1)
subplot(2,2,1)
smoosh= [1,0,0,0,0,0,0;
         0,1,1,0,0,0,0;
         0,0,0,1,1,0,0;
         0,0,0,0,0,1,0];
     
plot(tOut,yOut*smoosh(1:3,:)');
hold on
plot(tOut,yOut*smoosh(4,:)','lineWidth',2);
xlabel('t');
ylabel('population');
ylim([0,0.8])
xlim([0,200])
meanX=yOut(end,end)./tOut(end)
title('Antibiotic Mixing');
% dumbApproxA= importation(2)/(q-g) *exp(-q*tOut);
% dumbApproxB= (q/g)*importation(2)/(q-g) *(1-exp(-g*tOut));
% 
% plot(tOut,dumbApproxA)
% plot(tOut,dumbApproxB)
