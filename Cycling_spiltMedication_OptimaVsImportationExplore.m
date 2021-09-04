% ModelIS   [ S, R_A^A,R_A^B,R_B^A,R_B^B,X ]

 y0=[rand(6,1);0]
% y0=[ swap*(yOut(end,1:6)') ;0];
% y0= [(yOut(end,1:6)');0];

importationRates=logspace(-5,-1,120);

maxTimes=zeros(length(importationRates),6);
Tlist=zeros(1,450);
%Tlist(1)=0.05;

for(bbb=1:length(importationRates))
    bbb
    
mu=1/15; %(immgration/emigration rate)
beta= [1,0.99,0.99]; %Infection rates;
smallImport=importationRates(bbb);

importation= [0.5,smallImport,smallImport,smallImport,smallImport,1]';
importation=  importation*mu/([1,1,0,1,0,1]*importation);

g=1/10; %Ten day recovery without meds;
q=1/4; %treated recovery.

tau= (q-g);

Xequil= (g+mu)/beta(2);
Sequil= importation(1)/(q+mu-beta(1)*Xequil)
Bequil= importation(4)/(q+mu-beta(1)*Xequil)
Aequil= 1-Xequil-Sequil-Bequil
tSat= (log(Aequil)-log(Bequil)-1)*2/tau;

MixEquil= (g + (tau*g)/(q+g))/beta(2);



Tlist= tSat*logspace(-0.5, 0.5, length(Tlist));


Result=0*Tlist;


y0=[Sequil;Aequil;0;0;Bequil;Xequil;0];

ExitNumbers=-ones(length(Tlist),16);
Maxs=-ones(length(Tlist),1);
Mins=-ones(length(Tlist),1);
Means=-ones(length(Tlist),1);
meanDownSwing=-ones(length(Tlist),1);
meanUpSwing=-ones(length(Tlist),1);


 y0=[ones(6,1)/6;zeros(10,1)];
 
for(lll=1:length(Tlist))
[bbb,lll]
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

 
warmupTime= ceil(500/(2*T))*(2*T);
[tOut,yOut] = ode45(Deriv,[T/2,warmupTime],y0);
y0=[yOut(end,1:6)';zeros(10,1)];
[tOut,yOut] = ode45(Deriv,[0,warmupTime+T],y0);
y0=[(y0(1:6)+yOut(end,[1,5,4,3,2,6])')/2; zeros(10,1)];

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




 DoubleStateEquilibrium= (mu+g)./beta(3)^2;%%Seems a reasonable beta...
    
[~,indx]=max(((ExitNumbers(:,12)-DoubleStateEquilibrium)./(sum(ExitNumbers(:,8:11),2))));
maxTimes(bbb,1)=Tlist(indx);

  [~,indx]=max(((ExitNumbers(:,12)-DoubleStateEquilibrium)./ExitNumbers(:,15)));
maxTimes(bbb,2)=Tlist(indx);

  [~,indx]=max(((ExitNumbers(:,12)-DoubleStateEquilibrium)./ExitNumbers(:,16)));
maxTimes(bbb,3)=Tlist(indx);

maxTimes(bbb,4)=tSat;

maxTimes(bbb,5)=max(Tlist);
maxTimes(bbb,6)=min(Tlist);
end


figure(18);
subplot(1,2,2);
plot(importationRates,maxTimes(:,1:3),'lineWidth',2);
hold on
plot(importationRates,maxTimes(:,4),'k');
plot(importationRates,maxTimes(:,5),'k:');
plot(importationRates,maxTimes(:,6),'k:');

title("6 box model")
xlabel("import rate m_A=m_B")
ylabel("Optimal T")
legend({"M_{base}","M_{select}","M_{HGT}","t_{sat}"})
 set(gca, 'XScale', 'log')
 
 
 
figure(19);
subplot(1,2,2);
plot(importationRates,maxTimes(:,1:3)./maxTimes(:,4),'lineWidth',2);
hold on
plot(importationRates,maxTimes(:,4)./maxTimes(:,4),'k');
plot(importationRates,maxTimes(:,5)./maxTimes(:,4),'k:');
plot(importationRates,maxTimes(:,6)./maxTimes(:,4),'k:');


title("6 box model")
xlabel("import rate $m_A=m_B$")
ylabel("Optimal T / t_{sat}")
legend({"M_{base}","M_{select}","M_{HGT}","t_{sat}"})
set(gca, 'XScale', 'log')


figure(29);
subplot(1,2,2);
plot(importationRates,maxTimes(:,1:3)-maxTimes(:,4),'lineWidth',2);
hold on
plot(importationRates,maxTimes(:,4)-maxTimes(:,4),'k');
plot(importationRates,maxTimes(:,5)-maxTimes(:,4),'k:');
plot(importationRates,maxTimes(:,6)-maxTimes(:,4),'k:');


title("6 box model")
xlabel("import rate $m_A=m_B$")
ylabel("Optimal T - t_{sat}")
legend({"M_{base}","M_{select}","M_{HGT}","t_{sat}"})
set(gca, 'XScale', 'log')