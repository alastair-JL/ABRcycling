
importationRates=logspace(-5,-1,120);

maxTimes=zeros(length(importationRates),6);
Tlist=zeros(1,450);
extra=8;

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



Tlist(1:end)= tSat*logspace(-0.5, 0.5, length(Tlist));


Result=0*Tlist;


y0=[Sequil;Aequil;Bequil;Xequil];

ExitNumbers=-ones(length(Tlist),12);
Maxs=-ones(length(Tlist),1);
Mins=-ones(length(Tlist),1);
Means=-ones(length(Tlist),1);
meanDownSwing=-ones(length(Tlist),1);
meanUpSwing=-ones(length(Tlist),1);


 y0=[ones(4,1)/4;zeros(extra,1)];

for(lll=1:length(Tlist))
[bbb,lll]
T=Tlist(lll);

chi=@(t) ceil( mod(t/T,2)+10^-8)-1;
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
y0=[(y0(1:4)+yOut(end,[1,3,2,4])')/2; zeros(extra,1)];

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


 DoubleStateEquilibrium= (mu+g)./beta(3)^2;%%Seems a reasonable beta...
    
[~,indx]=max(((ExitNumbers(:,8)-DoubleStateEquilibrium)./(sum(ExitNumbers(:,6:7),2))));
maxTimes(bbb,1)=Tlist(indx);

  [~,indx]=max(((ExitNumbers(:,8)-DoubleStateEquilibrium)./ExitNumbers(:,11)));
maxTimes(bbb,2)=Tlist(indx);

  [~,indx]=max(((ExitNumbers(:,8)-DoubleStateEquilibrium)./ExitNumbers(:,12)));
maxTimes(bbb,3)=Tlist(indx);

maxTimes(bbb,4)=tSat;
maxTimes(bbb,5)=max(Tlist);
maxTimes(bbb,6)=min(Tlist);
end

 figure(18)
subplot(1,2,1);
plot(importationRates,maxTimes(:,1:3),'lineWidth',2);
hold on
plot(importationRates,maxTimes(:,4),'k');
plot(importationRates,maxTimes(:,5),'k:');
plot(importationRates,maxTimes(:,6),'k:');

title("4 box model")
xlabel("import rate m_A=m_B")
ylabel("Optimal T")
legend({"M_{base}","M_{select}","M_{HGT}","t_{sat}"})
 set(gca, 'XScale', 'log')
 
 
figure(19);
subplot(1,2,1);
plot(importationRates,maxTimes(:,1:3)./maxTimes(:,4),'lineWidth',2);
hold on
plot(importationRates,maxTimes(:,4)./maxTimes(:,4),'k');
plot(importationRates,maxTimes(:,5)./maxTimes(:,4),'k:');
plot(importationRates,maxTimes(:,6)./maxTimes(:,4),'k:');


title("4 box model")
xlabel("import rate $m_A=m_B$")
ylabel("Optimal T / t_{sat}")
legend({"M_{base}","M_{select}","M_{HGT}","t_{sat}"})
set(gca, 'XScale', 'log')


figure(29);
subplot(1,2,1);
plot(importationRates,maxTimes(:,1:3)-maxTimes(:,4),'lineWidth',2);
hold on
plot(importationRates,maxTimes(:,4)-maxTimes(:,4),'k');
plot(importationRates,maxTimes(:,5)-maxTimes(:,4),'k:');
plot(importationRates,maxTimes(:,6)-maxTimes(:,4),'k:');


title("4 box model")
xlabel("import rate $m_A=m_B$")
ylabel("Optimal T - t_{sat}")
legend({"M_{base}","M_{select}","M_{HGT}","t_{sat}"})
set(gca, 'XScale', 'log')