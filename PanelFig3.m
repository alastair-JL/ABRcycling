% ModelIS   [ S, R_A,R_B,X ]

 y0=[rand(4,1);0]
% y0=[ swap*(yOut(end,1:6)') ;0];
% y0= [(yOut(end,1:6)');0];

Tlist= logspace(-0.5,2.5,1250);
Tlist=[0,Tlist];
Result=0*Tlist;
mu=1/15; %(immgration/emigration rate)
beta= [1,0.99,0.99]; %Infection rates;
smallImport=0.01;

importation= [0.5,smallImport,smallImport,1]';
importation=  importation*mu/([1,1,1,1]*importation);

g=1/10; %Ten day recovery without meds;
q=1/2.5; %treated recovery.


tau= (q-g);
MixEquil= (g + (tau*g)/(q+g))/beta(2)


Xequil= (g+mu)/beta(2);
Sequil= importation(1)/(q+mu-beta(1)*Xequil)
Bequil= importation(3)/(q+mu-beta(1)*Xequil)
Aequil= 1-Xequil-Sequil-Bequil
tSat= (log(Aequil)-log(Bequil)-1)*2/tau;

y0=[Sequil;Aequil;0;0;Bequil;Xequil;0];

extra=8;
ExitNumbers=-ones(length(Tlist),4+extra);
Maxs=-ones(length(Tlist),1);
Mins=-ones(length(Tlist),1);
Means=-ones(length(Tlist),1);
meanDownSwing=-ones(length(Tlist),1);
meanUpSwing=-ones(length(Tlist),1);


TA=22;
TB=22;

chi=@(t) 1*(mod(t,TA+TB)>TA);
ts=linspace(0,5*(TA+TB),1000);


importation= [0.5,smallImport,smallImport,1]';
importation=  importation*mu/([1,1,1,1]*importation);

b=0.5;
betaMatrix= [beta(1),0,0,0;
             0,beta(2),0,0;
             0,0,beta(3),0;
             -beta(1),-beta(2),-beta(3),0];
         
recoveryA= diag(-[q,g,q,0]);
recoveryA(end,:)=[q,g,q,0];
         
recoveryB= diag(-[q,q,g,0]);
recoveryB(end,:)=[q,q,g,0];
         
recovery= @(t)  recoveryB+chi(t)*(recoveryA-recoveryB);

Deriv =@(t,V) [(importation + betaMatrix*V(1:4)*V(4) + recovery(t)*V(1:4)-mu*V(1:4)); V(1:4); importation(2)./(V(2)); importation(3)./(V(3)); V(3)*chi(t)+(V(2)*(1-chi(t))); V(2)*V(3)];

 y0=[ones(4,1)/4;zeros(extra,1)];
 
warmupTime= 4*(TA+TB);

[tOut,yOut] = ode45(Deriv,[0,125*(TA+TB)+200],y0);

Select=find(tOut>125*(TA+TB));

figure(1)
subplot(2,2,3)

plot(tOut(Select)-tOut(Select(1)),yOut(Select,1:3));
hold on
plot(tOut(Select)-tOut(Select(1)),yOut(Select,4),'lineWidth',2);
ylim([0,0.8])
xlim([0,125])
xlabel('t');
ylabel('population');
title('Fast Cycling');