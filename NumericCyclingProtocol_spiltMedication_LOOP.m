% ModelIS   [ S, R_A^A,R_A^B,R_B^A,R_B^B,X ]

 y0=[rand(6,1);0]
% y0=[ swap*(yOut(end,1:6)') ;0];
% y0= [(yOut(end,1:6)');0];

Tlist= logspace(-0.5,2.5,50);
Tlist=[0,Tlist];

Result=0*Tlist;
mu=1/15; %(immgration/emigration rate)
beta= [1,0.99,0.99]; %Infection rates;
smallImport=0.001;

importation= [0.5,smallImport,smallImport,smallImport,smallImport,1]';
importation=  importation*mu/([1,1,0,1,0,1]*importation);

g=1/10; %Ten day recovery without meds;
q=1/4; %treated recovery. NOT the same thing a "q" in the manuscript!!!


tau= (q-g);

Xequil= (g+mu)/beta(2);
Sequil= importation(1)/(q+mu-beta(1)*Xequil)
Bequil= importation(3)/(q+mu-beta(1)*Xequil)
Aequil= 1-Xequil-Sequil-Bequil
tSat= (log(Aequil)-log(Bequil)-1)*2/tau;


MixEquil= (g + (tau*g)/(q+g))/beta(2)


Xequil= (g+mu)/beta(2);
Sequil= importation(1)/(q+mu-beta(1)*Xequil)
Bequil= importation(4)/(q+mu-beta(1)*Xequil)
Aequil= 1-Xequil-Sequil-Bequil

y0=[Sequil;Aequil;0;0;Bequil;Xequil;0];

ExitNumbers=-ones(length(Tlist),16);
Maxs=-ones(length(Tlist),1);
Mins=-ones(length(Tlist),1);
Means=-ones(length(Tlist),1);
meanDownSwing=-ones(length(Tlist),1);
meanUpSwing=-ones(length(Tlist),1);


for(lll=2:length(Tlist))

T=Tlist(lll);


chi=@(t) ceil( mod(t/T,2)+10^-8)-1;

ts=linspace(0,10*T,1000);



filter= @(t) [1;chi(t);1-chi(t);chi(t);1-chi(t);1];

importation= [0.5,smallImport,smallImport,smallImport,smallImport,1]';
importation=  importation*mu/(filter(1)'*importation);

b=0.5;
betaMatrix= [beta(1),0,0,0,0,0;
             0,beta(2),beta(2),0,0,0;
             0,beta(2),beta(2),0,0,0;
             0,0,0,beta(3),beta(3),0;
             0,0,0,beta(3),beta(3),0;
             -beta(1),-beta(2),-beta(2),-beta(3),-beta(3),0];
         
recovery= diag(-[q,g,q,q,g,0]);
recovery(6,:)=[q,g,q,q,g,0];
         
Deriv =@(t,V) [(filter(t).*importation + filter(t).*betaMatrix*V(1:6)*V(6) + recovery*V(1:6)-mu*V(1:6)); V(1:6); importation(2)./(V(2)+V(3)); importation(4)./(V(4)+V(5)); V(3)+V(4);(V(2)+V(3))*(V(4)+V(5))  ];

 y0=[ones(6,1)/6;zeros(10,1)];
 
warmupTime= ceil(500/(2*T))*(2*T);
[tOut,yOut] = ode45(Deriv,[T/2,warmupTime],y0);
y0=[yOut(end,1:6)';zeros(10,1)];
[tOut,yOut] = ode45(Deriv,[0,warmupTime+T],y0);
y0=[(y0(1:6)+yOut(end,[1,5,4,3,2,6])')/2; zeros(10,1)]

[tOut,yOut] = ode45(Deriv,[0,20*T],y0);

Maxs(lll)= max(yOut(:,2)+yOut(:,3));
Mins(lll)= min(yOut(:,2)+yOut(:,3));
Means(lll)=(yOut(end,8)+yOut(end,9))/(tOut(end)-tOut(1));
indexSwitch=find(abs(diff(sign(diff(yOut(:,2)+(yOut(:,3))))))==2); %This gives index of each moment we switch directions.

switchTimes=diff(tOut(indexSwitch));
meanDownSwing(lll)=mean(switchTimes(1:2:length(switchTimes)));
meanUpSwing(lll)=mean(switchTimes(2:2:length(switchTimes)));




Result(lll)= (yOut(end,7)-yOut(1,7))/(tOut(end)-tOut(1));

ExitNumbers(lll,1:6)= yOut(end,1:6);
ExitNumbers(lll,7:end)= yOut(end,7:end)/(tOut(end)-tOut(1));

end

lll=1;
chi=@(t) 0.5;

filter= @(t) [1;chi(t);1-chi(t);chi(t);1-chi(t);1];

importation= [0.5,smallImport,smallImport,smallImport,smallImport,1]';
importation=  importation*mu/(filter(1)'*importation);

b=0.5;
betaMatrix= [beta(1),0,0,0,0,0;
             0,beta(2),beta(2),0,0,0;
             0,beta(2),beta(2),0,0,0;
             0,0,0,beta(3),beta(3),0;
             0,0,0,beta(3),beta(3),0;
             -beta(1),-beta(2),-beta(2),-beta(3),-beta(3),0];
         
recovery= diag(-[q,g,q,q,g,0]);
recovery(6,:)=[q,g,q,q,g,0];
         
Deriv =@(t,V) [(filter(t).*importation + filter(t).*betaMatrix*V(1:6)*V(6) + recovery*V(1:6)-mu*V(1:6)); V(1:6); importation(2)./(V(2)+V(3)); importation(4)./(V(4)+V(5)); tau*(V(3))./(V(2)+V(3)) ;tau*(V(4))./(V(4)+V(5))  ];

 y0=[ones(6,1)/6;zeros(10,1)];
  
warmupTime= 1500;
[tOut,yOut] = ode45(Deriv,[0,warmupTime],y0);

Maxs(lll)= max(yOut(end,2)+yOut(end,3));
Mins(lll)= min(yOut(end,2)+yOut(end,3));
Means(lll)=(yOut(end,2)+yOut(end,3))/(tOut(end)-tOut(1));
meanDownSwing(lll)=nan;
meanUpSwing(lll)=nan;
Result(lll)= yOut(end,3);

ExitNumbers(lll,1:6)= yOut(end,1:6);
ExitNumbers(lll,7:12)= yOut(end,1:6);

ExitNumbers(lll,13)= importation(2)./yOut(end,2);
ExitNumbers(lll,14)= importation(3)./yOut(end,3);
ExitNumbers(lll,15)= yOut(end,3)+yOut(end,4);
ExitNumbers(lll,16)= (yOut(end,4)+yOut(end,5))*(yOut(end,2)+yOut(end,3));


figure(3)
subplot(1,2,2)
% plot(Tlist,(sum(ExitNumbers(:,8:11),2))./max(sum(ExitNumbers(:,8:11),2)),'LineWidth',2);
% hold on
% plot(Tlist,(ExitNumbers(:,15))./max(ExitNumbers(:,15)),'LineWidth',2);
% plot(Tlist,(ExitNumbers(:,16))./max(ExitNumbers(:,16)),'LineWidth',2);
% plot(Tlist,1+0*Tlist,'k');
plot(Tlist,1+0*Tlist,'LineWidth',2);
hold on
plot(Tlist,(sum(ExitNumbers(:,8:11),2))./(sum(ExitNumbers(1,8:11),2)),'LineWidth',2);
plot(Tlist,(ExitNumbers(:,15))./(ExitNumbers(1,15)),'LineWidth',2);
plot(Tlist,(ExitNumbers(:,16))./(ExitNumbers(1,16)),'LineWidth',2);

legend({'M_{import}','M_{base}','M_{select}','M_{HGT}'})
xlabel('T')
ylabel('M')
title('6 State System')

DoubleStateEquilibrium= (mu+g)./(beta(2)*beta(2));

figure(23);
subplot(1,2,2)

plot(Tlist,((ExitNumbers(:,12)-DoubleStateEquilibrium))./max(((ExitNumbers(:,12)-DoubleStateEquilibrium) )),'LineWidth',2);
hold on
plot(Tlist,((ExitNumbers(:,12)-DoubleStateEquilibrium)./(sum(ExitNumbers(:,8:11),2)))./max(((ExitNumbers(:,12)-DoubleStateEquilibrium)./(sum(ExitNumbers(:,8:11),2) ))),'LineWidth',2);
plot(Tlist,(((ExitNumbers(:,12)-DoubleStateEquilibrium)./ExitNumbers(:,15)))./max((((ExitNumbers(:,12)-DoubleStateEquilibrium)./ExitNumbers(:,15)))),'LineWidth',2);
plot(Tlist,(((ExitNumbers(:,12)-DoubleStateEquilibrium)./ExitNumbers(:,16)))./max((((ExitNumbers(:,12)-DoubleStateEquilibrium)./ExitNumbers(:,16)))),'LineWidth',2);

%plot(Tlist,1+0*Tlist,'k');
% plot(Tlist,((ExitNumbers(:,12)-DoubleStateEquilibrium)./(sum(ExitNumbers(:,8:11),2)))./((ExitNumbers(1,12)-DoubleStateEquilibrium)./(sum(ExitNumbers(1,8:11),2) )),'LineWidth',2);
% hold on
% plot(Tlist,(((ExitNumbers(:,12)-DoubleStateEquilibrium)./ExitNumbers(:,15)))./(((ExitNumbers(1,12)-DoubleStateEquilibrium)./ExitNumbers(1,15))),'LineWidth',2);
% plot(Tlist,(((ExitNumbers(:,12)-DoubleStateEquilibrium)./ExitNumbers(:,16)))./(((ExitNumbers(1,12)-DoubleStateEquilibrium)./ExitNumbers(1,16))),'LineWidth',2);
% plot(Tlist,1+0*Tlist,'k');
legend({'M_{import}','M_{base}','M_{select}','M_{HGT}'})
xlabel('T')
ylabel('$\int_0^{t_c} X - X_0 dt $','Interpreter','latex')
title('6 State System')

figure(73);
subplot(1,2,2)

plot(Tlist,(ExitNumbers(:,12)./max(ExitNumbers(:,12))),'LineWidth',2);
hold on
plot(Tlist,((ExitNumbers(:,12))./(sum(ExitNumbers(:,8:11),2)))./max(((ExitNumbers(:,12))./(sum(ExitNumbers(:,8:11),2)))),'LineWidth',2);
plot(Tlist,(((ExitNumbers(:,12))./ExitNumbers(:,15)))./max((((ExitNumbers(:,12))./ExitNumbers(:,15)))),'LineWidth',2);
plot(Tlist,(((ExitNumbers(:,12))./ExitNumbers(:,16)))./max((((ExitNumbers(:,12))./ExitNumbers(:,16)))),'LineWidth',2);

% plot(Tlist,((ExitNumbers(:,12))./(sum(ExitNumbers(:,8:11),2)))./((ExitNumbers(1,8))./(sum(ExitNumbers(1,8:11),2))),'LineWidth',2);
% hold on
% plot(Tlist,(((ExitNumbers(:,12))./ExitNumbers(:,15)))./(((ExitNumbers(1,12))./ExitNumbers(1,15))),'LineWidth',2);
% plot(Tlist,(((ExitNumbers(:,12))./ExitNumbers(:,16)))./(((ExitNumbers(1,12))./ExitNumbers(1,16))),'LineWidth',2);
% plot(Tlist,1+0*Tlist,'k');
legend({'M_{import}','M_{base}','M_{select}','M_{HGT}'})
xlabel('T')
ylabel('$\int_0^{t_c} X dt $','Interpreter','latex')
title('6 State System')




figure(29);

BetaList=linspace(0.6,1.2,1000);
maxTimes=zeros(length(BetaList),4);

for(bbb=1:length(BetaList))
    DoubleStateEquilibrium= (mu+g)./BetaList(bbb);
[~,indx]=max(((ExitNumbers(:,12)-DoubleStateEquilibrium)));
maxTimes(bbb,1)=Tlist(indx);

[~,indx]=max(((ExitNumbers(:,12)-DoubleStateEquilibrium)./(sum(ExitNumbers(:,8:11),2))));
maxTimes(bbb,2)=Tlist(indx);

  [~,indx]=max(((ExitNumbers(:,12)-DoubleStateEquilibrium)./ExitNumbers(:,15)));
maxTimes(bbb,3)=Tlist(indx);

  [~,indx]=max(((ExitNumbers(:,12)-DoubleStateEquilibrium)./ExitNumbers(:,16)));
maxTimes(bbb,4)=Tlist(indx);
end
 

subplot(1,2,2)

plot(BetaList,maxTimes,'LineWidth',2);
hold on
plot([0,1.5],[tSat,tSat],'k');
plot([beta(2),beta(2)],[0,Tlist(end)],'--k');
xlim([0.65,1])
ylim([0,300])
legend({'M_{import}','M_{base}','M_{select}','M_{HGT}'})
xlabel('$\beta_{AB}$','Interpreter','latex')
ylabel('Optimal cycle Time','Interpreter','latex')
title('6 State System')
