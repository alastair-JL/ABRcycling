% ModelIS   [ S, R_A^A,R_A^B,R_B^A,R_B^B,X ,R_AB]

ChiList= linspace(0,1,500);

Result=0*ChiList;
mu=1/15; %(immgration/emigration rate)
beta= [1,0.99,0.99,0.98]; %Infection rates;
smallImport=0.001;

mutationRate=0.1;

importation= [0.5,smallImport,smallImport,smallImport*2,smallImport*2,1,0]';
importation=  importation*mu/([1,1,0,1,0,1,1]*importation);

g=1/10; %Ten day recovery without meds;
q=1/2.5; %treated recovery.
durgCorrection=1/6; %Rate at which incorrect drug use is replaced by correct drug use.

tau= (q-g);

MutantArrivalRates= zeros(length(ChiList),4);
M_IMPORT=1;
M_BASE=2;
M_SELECT=3;
M_HGT=4;

X365=zeros(length(ChiList),4);
X365_lowMutant=zeros(length(ChiList),4);
DoesThisIntegrateToOne=zeros(length(ChiList),4);
Expected_T_half=zeros(length(ChiList),4);
X_T=zeros(length(ChiList),4);
X_Tstar=zeros(length(ChiList),4);
X_bar=zeros(length(ChiList),1);

epsilon=10^-2.5;

barXm=(g+mu)/beta(4)

RhoObserved=zeros(length(ChiList),2);
RhoPredicted=zeros(length(ChiList),2);

for(ccc=1:length(ChiList))
    
    chi=ChiList(ccc);
    
filter= [1;chi;1-chi;chi;1-chi;1;1];

importation= [0.5,smallImport,smallImport,smallImport*2,smallImport*2,1,0]';
importation=  importation*mu/(filter'*importation);

b=0.5;
betaMatrix= [beta(1),0,0,0,0,0,0;
             0,beta(2),beta(2),0,0,0,0;
             0,beta(2),beta(2),0,0,0,0;
             0,0,0,beta(3),beta(3),0,0;
             0,0,0,beta(3),beta(3),0,0;
             -beta(1),-beta(2),-beta(2),-beta(3),-beta(3),0,-beta(4);
             0,0,0,0,0,0,beta(4)];
      
correctionMatrix=zeros(7,7);         
correctionMatrix(2,2)=-durgCorrection;
correctionMatrix(3,2)=durgCorrection;
correctionMatrix(5,5)=-durgCorrection;
correctionMatrix(4,5)=durgCorrection;



recovery= diag(-[q,g,q,q,g,0,g]);
recovery(6,:)=[q,g,q,q,g,0,g];



Deriv =@(t,V) [(filter.*importation + filter.*betaMatrix*V(1:7)*V(6) + (recovery+correctionMatrix)*V(1:7)-mu*V(1:7)); V(1:7) ];

 y0=[ones(6,1)/6;zeros(8,1)];
  
warmupTime= 1500;
[tOut,yOut] = ode45(Deriv,[0,warmupTime],y0); %This is a very derp plan to get the equilibrium values. 
                                            %There are FAR smarter plans avaliable. Oh well.
           
yFinal=yOut(end,1:7);

RhoObserved(ccc,1)= yFinal(3)/(yFinal(2)+yFinal(3));
RhoObserved(ccc,2)= yFinal(4)/(yFinal(4)+yFinal(5));
RhoPredicted(ccc,1)= ((mu+g)*chi + durgCorrection)/(tau*(1-chi)+mu+g+durgCorrection);
RhoPredicted(ccc,2)= ((mu+g)*(1-chi) + durgCorrection)/(tau*chi+mu+g+durgCorrection);

MutantArrivalRates(ccc,M_IMPORT)=1;
MutantArrivalRates(ccc,M_BASE)=sum(yFinal(2:5));
MutantArrivalRates(ccc,M_SELECT)=sum(yFinal(3:4));
MutantArrivalRates(ccc,M_HGT)=sum(yFinal(2:3))*sum(yFinal(4:5));
X_bar(ccc)=yFinal(6);

yFinal(7)=epsilon;

tIn= 365-fliplr([0,logspace(-9,0,2000)*365]);

[tOut,yOut] = ode45(Deriv,tIn,[yFinal',zeros(7,1)]); %This is a very derp plan to get the equilibrium values. 

LogRfinal=log(yOut(:,7));
logRInitial= log(epsilon);

PostArrivalIntegral= barXm*tOut +LogRfinal-logRInitial;
%The behavior of the system post R_AB arrival can be found using the Log
%forumal in the appendix.

PreArrivalIntegral= (tOut(end)-tOut)*yFinal(6);
%The integral prior to T_epsilon is just 

TotalIntegral=PostArrivalIntegral+PreArrivalIntegral;


for(mmm=1:size(MutantArrivalRates,2))
M=MutantArrivalRates(ccc,mmm)
PDF=  exp(-M*(tOut(end)-tOut))*M;
%The probability density function of arrival times is exponential, with
%decay rate proptional to M_i
widths= [tOut(2)-tOut(1);tOut(3:end)-tOut((3:end)-2);tOut(end)-tOut(end-1)];

if(M<=10^-7)
   PDF(1:end)=0; 
end
X365(ccc,mmm)=  sum( widths.*TotalIntegral.*PDF )/2+ exp(-M*tOut(end))*yFinal(6)*365;
  %Here we count all the probabilities for if T_eps<365, along with the
  %alternative possibility that T_eps>365, in which case we use the
  %pre-mutant equilibrium across the whole integral.
 DoesThisIntegrateToOne(ccc,mmm)=sum( widths.*PDF )/2+ exp(-M*tOut(end));

end
  
  
for(mmm=1:size(MutantArrivalRates,2))
M=MutantArrivalRates(ccc,mmm)*mutationRate;
PDF=  exp(-M*(tOut(end)-tOut))*M;
%The probability density function of arrival times is exponential, with
%decay rate proptional to M_i
widths= [tOut(2)-tOut(1);tOut(3:end)-tOut((3:end)-2);tOut(end)-tOut(end-1)];

if(M<=10^-7)
   PDF(1:end)=0; 
end
X365_lowMutant(ccc,mmm)=  sum( widths.*TotalIntegral.*PDF )/2+ exp(-M*tOut(end))*yFinal(6)*365;
  %Here we count all the probabilities for if T_eps<365, along with the
  
end
  
  
if(yFinal(6)<= barXm)
    Expected_T_half(ccc,:)=inf;
else
    
    Opt    = odeset('Events', @myEvent);
    
    warmupTime= 3650;
    [tOut,yOut] = ode45(Deriv,[0,warmupTime],[yFinal',zeros(7,1)],Opt); %This is a very derp plan to get the equilibrium values. 
    
    if( (tOut(end)>=3640) || (yOut(end,7)<sum(yOut(end,1:5))) )
        Expected_T_half(ccc,:)=inf;
        X_T(ccc,:)=inf;
    else
        Expected_T_half(ccc,:)=tOut(end)+1./MutantArrivalRates(ccc,:);
        X_T(ccc,:)=yOut(end,13)+yFinal(6)./MutantArrivalRates(ccc,:);
    end
    
end

X_Tstar(ccc,:)=  (yFinal(6)- barXm)./MutantArrivalRates(ccc,:);



end

figure()
plot(ChiList,DoesThisIntegrateToOne)

figure(11)
plot(ChiList,RhoObserved);
hold on
plot(ChiList,RhoPredicted,':');
xlabel('\chi_A')
ylabel('\rho')

figure(5)
subplot(2,1,1)
plot(ChiList,MutantArrivalRates)
xlabel('\chi_A')
ylabel('M')
legend({'import','base','select','HGT'});

subplot(2,1,2)
plot(ChiList,X_bar)
hold on
plot(ChiList,0*X_bar+g/beta(4),'k:');
xlabel('\chi_A')
ylabel('\bar X')


figure(1)
subplot(3,2,1)
set(gca, 'ColorOrder', StandardColours,'NextPlot', 'replacechildren');
hold on
plot(ChiList,X365,'--','lineWidth',2 )
xlabel('\chi_A')
ylabel('X_{365}')
title('High mutation rate')

subplot(3,2,2)
set(gca, 'ColorOrder', StandardColours,'NextPlot', 'replacechildren');
hold on
plot(ChiList,X365_lowMutant,'--','lineWidth',2 )
xlabel('\chi_A')
ylabel('X_{365}')
title('Low mutation rate')

subplot(3,2,3)
set(gca, 'ColorOrder', StandardColours,'NextPlot', 'replacechildren');
hold on
plot(ChiList,Expected_T_half,'--','lineWidth',2 )
xlabel('\chi_A')
ylabel('T_{1/2}')

subplot(3,2,4)
set(gca, 'ColorOrder', StandardColours,'NextPlot', 'replacechildren');
hold on
plot(ChiList,X_T,'--','lineWidth',2 )
xlabel('\chi_A')
ylabel('X_{T}')

subplot(3,2,5)
set(gca, 'ColorOrder', StandardColours,'NextPlot', 'replacechildren');
hold on
plot(ChiList,X_Tstar,'--','lineWidth',2 )
xlabel('\chi_A')
ylabel('X_{T*}')
title('Raw Calculation')

subplot(3,2,6)
set(gca, 'ColorOrder', StandardColours,'NextPlot', 'replacechildren');
hold on
plot(ChiList,X_Tstar./max(X_TstarSave),'--','lineWidth',2 )
xlabel('\chi_A')
ylabel('X_{T*}')
legend({'import','base','select','HGT'});
title('Rescaled')
ylim([-0.1,1.1])

function [value, isterminal, direction] = myEvent(T, Y)
value      =  (sum(Y(1:5))-Y(7))*(Y(7)- 10^-3);
isterminal = 1;   % Stop the integration
direction  = 0;
end
