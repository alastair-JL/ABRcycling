% ModelIS   [ S, R_A,R_B,X ]

 y0=[rand(4,1);0]
% y0=[ swap*(yOut(end,1:6)') ;0];
% y0= [(yOut(end,1:6)');0];

Tlist= logspace(-0.5,2.5,250);
Tlist=[0,Tlist];
Result=0*Tlist;
mu=1/15; %(immgration/emigration rate)
beta= [1,0.99,0.99]; %Infection rates;
smallImport=0.001;

importation= [0.5,smallImport,smallImport,1]';
importation=  importation*mu/([1,1,1,1]*importation);

g=1/10; %Ten day recovery without meds;
q=1/4; %treated recovery.


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

phaseOffset=0

for(lll=length(Tlist):-1:2)

T=Tlist(lll);

chi=@(t) ceil( mod((t+2*T*phaseOffset)/T,2)+10^-8)-1;
ts=linspace(0,10*T,1000);


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
 
warmupTime= ceil(500/(2*T))*(2*T);
[tOut,yOut] = ode45(Deriv,[T/2,warmupTime],y0);
y0=[yOut(end,1:4)';zeros(extra,1)];
[tOut,yOut] = ode45(Deriv,[0,warmupTime+T],y0);
y0=[(y0(1:4)+yOut(end,[1,3,2,4])')/2; zeros(extra,1)]

[tOut,yOut] = ode45(Deriv,[0,20*T],y0);

Maxs(lll)= max(yOut(:,2)+yOut(:,3));
Mins(lll)= min(yOut(:,2)+yOut(:,3));
Means(lll)=(yOut(end,8)+yOut(end,9))/(tOut(end)-tOut(1));
indexSwitch=find(abs(diff(sign(diff(yOut(:,2)+(yOut(:,3))))))==2); %This gives index of each moment we switch directions.

switchTimes=diff(tOut(indexSwitch));
meanDownSwing(lll)=mean(switchTimes(1:2:length(switchTimes)));
meanUpSwing(lll)=mean(switchTimes(2:2:length(switchTimes)));


Result(lll)= (yOut(end,7)-yOut(1,7))/(tOut(end)-tOut(1));

ExitNumbers(lll,1:4)= yOut(end,1:4);
ExitNumbers(lll,5:end)= yOut(end,5:end)/(tOut(end)-tOut(1));

end

lll=1;
chi=@(t) 0.5;

recovery= @(t)  recoveryB+chi(t)*(recoveryA-recoveryB);

Deriv =@(t,V) [(importation + betaMatrix*V(1:4)*V(4) + recovery(t)*V(1:4)-mu*V(1:4)); V(1:4); importation(2)./(V(2)); importation(3)./(V(3)); V(3)*chi(t)+(V(2)*(1-chi(t))); V(2)*V(3)];

 y0=[ones(4,1)/4;zeros(extra,1)];
 
warmupTime= 1500;
[tOut,yOut] = ode45(Deriv,[0,warmupTime],y0);

Maxs(lll)= max(yOut(end,2)+yOut(end,3));
Mins(lll)= min(yOut(end,2)+yOut(end,3));
Means(lll)=(yOut(end,2)+yOut(end,3))/(tOut(end)-tOut(1));
meanDownSwing(lll)=nan;
meanUpSwing(lll)=nan;
Result(lll)= yOut(end,3);

ExitNumbers(lll,1:4)= yOut(end,1:4);
ExitNumbers(lll,5:8)= yOut(end,1:4);

ExitNumbers(lll,9)= importation(2)./yOut(end,2);
ExitNumbers(lll,10)= importation(3)./yOut(end,3);
ExitNumbers(lll,11)= 0.5*yOut(end,3)+0.5*yOut(end,2);
ExitNumbers(lll,12)= yOut(end,3)*yOut(end,2);


figure(23);

plot(Tlist,((ExitNumbers(:,8)-DoubleStateEquilibrium)./max((ExitNumbers(:,8)-DoubleStateEquilibrium))),'LineWidth',2);
hold on
plot(Tlist,((ExitNumbers(:,8)-DoubleStateEquilibrium)./(ExitNumbers(:,6)+ExitNumbers(:,7)))./max(((ExitNumbers(:,8)-DoubleStateEquilibrium)./(ExitNumbers(:,6)+ExitNumbers(:,7)))),'LineWidth',2);
plot(Tlist,(((ExitNumbers(:,8)-DoubleStateEquilibrium)./ExitNumbers(:,11)))./max((((ExitNumbers(:,8)-DoubleStateEquilibrium)./ExitNumbers(:,11)))),'LineWidth',2);
plot(Tlist,(((ExitNumbers(:,8)-DoubleStateEquilibrium)./ExitNumbers(:,12)))./max((((ExitNumbers(:,8)-DoubleStateEquilibrium)./ExitNumbers(:,12)))),'LineWidth',2);
plot([tSat,tSat],[0,1.1],'k');

legend({'M_{import}','M_{base}','M_{select}','M_{HGT}'})
xlabel('T')
ylabel('$\int_0^{t_c} X - X_0 dt $','Interpreter','latex')
title('4 State System')
