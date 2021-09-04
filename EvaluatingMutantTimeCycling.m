% ModelIS   [ S, R_A^A,R_A^B,R_B^A,R_B^B,X ,R_AB]

TList= linspace(1,300,550);

Result=0*TList;
mu=1/15; %(immgration/emigration rate)
beta= [1,0.99,0.99,0.98]; %Infection rates;
smallImport=0.001;

mutationRate=0.1;

importation= [0.5,smallImport,smallImport,smallImport,smallImport,1,0]';
importation=  importation*mu/([1,1,0,1,0,1,1]*importation);

g=1/10; %Ten day recovery without meds;
q=1/2.5; %treated recovery.
durgCorrection=0.01; %Rate at which incorrect drug use is replaced by correct drug use.

tau= (q-g);

MutantArrivalRates= zeros(length(TList),4);
M_IMPORT=1;
M_BASE=2;
M_SELECT=3;
M_HGT=4;

X365=zeros(length(TList),4);
X365_lowMutant=zeros(length(TList),4);
DoesThisIntegrateToOne=zeros(length(TList),4);
Expected_T_half=zeros(length(TList),4);
X_T=zeros(length(TList),4);
X_Tstar=zeros(length(TList),4);
X_bar=zeros(length(TList),1);

epsilon=10^-2.5;

barXm=(g+mu)/beta(4)

RhoObserved=zeros(length(TList),2);
RhoPredicted=zeros(length(TList),2);

importation= [0.5,smallImport,smallImport,smallImport,smallImport,1,0]';
importation=  importation*mu/(filter(1)'*importation);

b=0.5;
betaMatrix= [beta(1),0,0,0,0,0,0;
             0,beta(2),beta(2),0,0,0,0;
             0,beta(2),beta(2),0,0,0,0;
             0,0,0,beta(3),beta(3),0,0;
             0,0,0,beta(3),beta(3),0,0;
             -beta(1),-beta(2),-beta(2),-beta(3),-beta(3),0,-beta(4);
             0,0,0,0,0,0,beta(4)];
      
recovery= diag(-[q,g,q,q,g,0,g]);
recovery(6,:)=[q,g,q,q,g,0,g];


flipr= eye(6);
flipr(2:5,:)= flipr(5:-1:2,:);

filterFake= [1;0;1;0;1;1;1];
DerivFake =@(t,V) [(filterFake.*importation + filterFake.*betaMatrix*V(1:7)*V(6) + recovery*V(1:7)-mu*V(1:7)); V(1:7); (V(2)+V(3))*(V(4)+V(5)) ];


for(ccc=1:length(TList))
    
    T=TList(ccc);
    
    chi= @(t) 1*(mod(t,2*T)>T);
    
filter= @(t) [1;chi(t);1-chi(t);chi(t);1-chi(t);1;1];

Deriv =@(t,V) [(filter(t).*importation + filter(t).*betaMatrix*V(1:7)*V(6) + recovery*V(1:7)-mu*V(1:7)); V(1:7); (V(2)+V(3))*(V(4)+V(5)) ];

 y0=[ones(6,1)/6;zeros(9,1)];
  
warmupTime= ceil(500/(2*T))*2*T;

[tOut,yOut] = ode45(Deriv,[0,warmupTime],y0); %This is a very derp plan to get the equilibrium values. 

y0=[yOut(end,(1:6))';zeros(9,1)];

for(qmp=1:50)
    [qmp,ccc]
    [tOut,yOut]=ode45(DerivFake,[0,T],y0);
    y0(1:6)= (y0(1:6)+flipr*yOut(end,1:6)')/2;
end

%There are FAR smarter plans avaliable. Oh well.
 
% [tOut,yOut] = ode45(Deriv,[0,warmupTime],y0); %This is a very derp plan to get the equilibrium values. 

yFinal=yOut(end,7+(1:8));

MutantArrivalRates(ccc,M_IMPORT)=1;
MutantArrivalRates(ccc,M_BASE)=sum(yFinal(2:5))/T;
MutantArrivalRates(ccc,M_SELECT)=sum(yFinal(3:4))/T;
MutantArrivalRates(ccc,M_HGT)=sum(yFinal(8))/T;
end

figure(15)
subplot(2,2,1)
plot(TList,MutantArrivalRates(end,1)./MutantArrivalRates(:,1),'k','LineWidth',2)
xlabel('T')
title('M_{import}')

subplot(2,2,2)
plot(TList,MutantArrivalRates(end,2)./MutantArrivalRates(:,2),'k','LineWidth',2)
xlabel('T')
title('M_{base}')

subplot(2,2,3)
plot(TList,MutantArrivalRates(end,3)./MutantArrivalRates(:,3),'k','LineWidth',2)
xlabel('T')
title('M_{select}')

subplot(2,2,4)
plot(TList,MutantArrivalRates(end,4)./MutantArrivalRates(:,4),'k','LineWidth',2)
xlabel('T')
title('M_{HGT}')




function [value, isterminal, direction] = myEvent(T, Y)
value      =  (sum(Y(1:5))-Y(7))*(Y(7)- 10^-3);
isterminal = 1;   % Stop the integration
direction  = 0;
end
