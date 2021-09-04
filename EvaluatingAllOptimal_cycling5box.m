% ModelIS   [ S, R_A,R_B,X ]

 y0=[rand(4,1);0]
% y0=[ swap*(yOut(end,1:6)') ;0];
% y0= [(yOut(end,1:6)');0];

Tlist= logspace(0,2,120);
Tlist=[0,Tlist];
Result=0*Tlist;
mu=1/15; %(immgration/emigration rate)
beta= [1,0.99,0.99,0.98]; %Infection rates;
smallImport=0.001;

importation= [0.5,smallImport,smallImport,1,0]';
importation=  importation*mu/([1,1,1,1,1]*importation);

gamma=1/10; %Ten day recovery without meds;
tau=1/4-1/10; %treated recovery.



MutantArrivalRates= zeros(length(Tlist),4);
M_IMPORT=1;
M_BASE=2;
M_SELECT=3;
M_HGT=4;


Xequil= (gamma+mu)/beta(2);
Sequil= importation(1)/(tau+mu-beta(1)*Xequil)
Bequil= importation(3)/(tau+mu-beta(1)*Xequil)
Aequil= 1-Xequil-Sequil-Bequil
tSat= (log(Aequil)-log(Bequil)-1)*2/tau;

yFinal=[Sequil;Aequil;Bequil;Xequil;0];

extra=8;
ExitNumbers=-ones(length(Tlist),4+extra);
Maxs=-ones(length(Tlist),1);
Mins=-ones(length(Tlist),1);
Means=-ones(length(Tlist),1);
meanDownSwing=-ones(length(Tlist),1);
meanUpSwing=-ones(length(Tlist),1);

epsilon=10^-2.5;

fineness=120;

nu=0.1*[1,2,7,100]; %These are the scaling factors on our various 
%mutation rates. and are basically picked so as to keep the various
%mutation channels RELATIVELY comparable; IE   M_base is often around
%0.5ish, M_HGT is often closer to 1/100
nuSmall=0.1*nu;

X365=ones(length(Tlist),4);
X365small=ones(length(Tlist),4);
Thalf=ones(length(Tlist),4);
X_T=ones(length(Tlist),4);
X_Tstar=ones(length(Tlist),4);


for(lll=length(Tlist):-1:2)

T=Tlist(lll);

chi=@(t) ceil( mod(t/T,2)+10^-8)-1;

b=0.5;
betaMatrix= [beta(1),0,0,0,0;
             0,beta(2),0,0,0;
             0,0,beta(3),0,0;
             -beta(1),-beta(2),-beta(3),0,-beta(4);
             0,0,0,0,beta(4)];
         
recoveryA= diag(-[tau+gamma,gamma,tau+gamma,0,gamma]);
recoveryA(4,:)=[tau+gamma,gamma,tau+gamma,0,gamma];
         
recoveryB= diag(-[tau+gamma,tau+gamma,gamma,0,gamma]);
recoveryB(4,:)=[tau+gamma,tau+gamma,gamma,0,gamma];
         
recovery= @(t)  recoveryB+chi(t)*(recoveryA-recoveryB);

Deriv =@(t,V) [(importation + betaMatrix*V(1:5)*V(4) + recovery(t)*V(1:5)-mu*V(1:5)); V(1:5); 1; V(3)+V(2) ; V(3)*chi(t)+(V(2)*(1-chi(t))); V(2)*V(3)];

 y0=[ones(4,1)/4;zeros(extra,1)];
 
warmupTime= ceil(50/(2*T))*(2*T);

[tOut,yOut] = ode45(Deriv,[0,warmupTime],[yFinal;zeros(9,1)]);

[tOutReference,yOutReference] = ode45(Deriv,[linspace(0,T-epsilon,fineness+1),linspace(T+epsilon,2*T,fineness+1)],[yOut(end,1:5)';zeros(9,1)]);
%This is going to double count t=0=2*T. That's fine.

%%Now, we run into a bit of a problem regarding the initial conditions and
%%initial "phase"... and also regarding the effects of the exact moment
%%that mutatant arrives. Which leads to... complexity....

%%Okay, so lets say the initial phase (Phi) is choosen uniformly at random
phis=tOut(1:end-1);


%So, I want a "time/offset to PDF" fucntion.
%Now... its' not TOO hard to make a "Time/offest to CDF function... I
%think"

yOutReferenceMutant=zeros(size(yOutReference,1),4);
yOutReferenceMutant(:,1)=1;
yOutReferenceMutant(:,2)=yOutReference(:,3)+yOutReference(:,2);
yOutReferenceMutant(:,3)=chi(tOutReference).*yOutReference(:,3) + (1-chi(tOutReference)).*yOutReference(:,2);
yOutReferenceMutant(:,4)=yOutReference(:,3).*yOutReference(:,2);

CDF= @(time,phi,Mtype) exp(-nu(Mtype)*(yOutReference(end,Mtype+10)*(floor((time+phi)/(2*T))) +  interp1(tOutReference,yOutReference(:,Mtype+10),((time+phi) - floor((time+phi)./(2*T)).*(2*T)),'linear','extrap')  - interp1(tOutReference,yOutReference(:,Mtype+10),phi))) ;
PDF= @(time,phi,Mtype) nu(Mtype).*interp1(tOutReference,yOutReferenceMutant(:,Mtype),((time+phi) - floor((time+phi)./(2*T)).*(2*T)),'linear','extrap').*CDF(time,phi,Mtype);
XDF=@(time,phi) interp1(tOutReference,yOutReference(:,4),((time+phi) - floor((time+phi)./(2*T)).*(2*T)),'linear','extrap');


CDFsmall= @(time,phi,Mtype) exp(-nuSmall(Mtype)*(yOutReference(end,Mtype+10)*floor((time+phi)/(2*T)) +  interp1(tOutReference,yOutReference(:,Mtype+10),((time+phi) - floor((time+phi)./(2*T)).*(2*T)),'linear','extrap')  - interp1(tOutReference,yOutReference(:,Mtype+10),phi))) ;
PDFsmall= @(time,phi,Mtype) nuSmall(Mtype).*interp1(tOutReference,yOutReferenceMutant(:,Mtype),((time+phi) - floor((time+phi)./(2*T)).*(2*T)),'linear','extrap').*CDFsmall(time,phi,Mtype);



PriorXintegral= @(time,phi) (yOutReference(end,9)*floor((time+phi)/(2*T)) +  interp1(tOutReference,yOutReference(:,9),mod(time+phi,2*T))  - interp1(tOutReference,yOutReference(:,9),phi)) ;
%This gives the Cumulitive distribution function and probability density
%function for like... time till arrival, assuming a fixed phase offset phi
%(which is in seconds, not fracctions of a cycle. Cool.

Mcycle= exp(-nu.*yOutReference(end,11:14));
Tcycle=2*T;
Xcycle=yOutReference(end,9);
meancycles=Mcycle./(1-Mcycle);
McycleRaw=nu.*yOutReference(end,11:14);


ExpectedValsX365=zeros(fineness,4);
ExpectedValsX365small=zeros(fineness,4);

DoesThisSumToOne=zeros(fineness,4);

DelayTime=zeros(fineness,1);
DelayX=zeros(fineness,1);



for(tktktk= [1:fineness,(1:fineness)+1+fineness]  )
    %%Assume that mutant occurs at time %tOut(tktktk)% During the loop, and
    %%Solve from there to 365.
    [tktktk,lll]
    yStart=yOutReference(tktktk,1:5);
    yStart(5)=epsilon;
    
    time= @(t) 365-t;%If our mutant arrives, with t days to go before 365,  this is how many days happened BEFORE it arrived.
    phi= @(t) mod(t+tOutReference(tktktk)-365,2*T); %If our mutant arrives, with t days to go before 365, this is the phase at time 0.
    
    
DerivWeave =@(t,V) [(importation + betaMatrix*V(1:5)*V(4) + recovery(t+tOut(tktktk) )*V(1:5)-mu*V(1:5));...
   V(4).*(1-CDF(365-t,phi(t),1))+ XDF(phi-t)*CDF(t,phi(t),1);...
   V(4).*(1-CDF(365-t,phi(t),2))+ XDF(phi-t)*CDF(t,phi(t),2);...
   V(4).*(1-CDF(365-t,phi(t),3))+ XDF(phi-t)*CDF(t,phi(t),3);...
   V(4).*(1-CDF(365-t,phi(t),4))+ XDF(phi-t)*CDF(t,phi(t),4);...
   V(4).*(1-CDFsmall(365-t,phi(t),1))+ XDF(phi-t)*CDFsmall(t,phi(t),1);...
   V(4).*(1-CDFsmall(365-t,phi(t),2))+ XDF(phi-t)*CDFsmall(t,phi(t),2);...
   V(4).*(1-CDFsmall(365-t,phi(t),3))+ XDF(phi-t)*CDFsmall(t,phi(t),3);...
   V(4).*(1-CDFsmall(365-t,phi(t),4))+ XDF(phi-t)*CDFsmall(t,phi(t),4);
];
%%Okay, so this function takes the Probability density function for the
%%probability of arrival at a particular time, then multiplies it by the
%%total X_365 integral ASSUMING we arrived at the time, and with phase
%%tOut(tktktk).  phi (phase at time zero) is backtracked out from this. 
%% It's... kind of a mess, but I'm 80% sure it works. :|
%%This DOESN'T account for the probability associated with NOT arriving
%%throughout the full 360, that is accounted for below.

%%This ALSO doesn't deal with the possiblities of different phases at
%%arrival time. This is dealt with elsewhere. 

    Opt    = odeset('Events', @ThalfEvent);
    
[tOut,yOut,tEvent,yEvent,~] = ode45(DerivWeave,linspace(0,365,4000),[yStart';zeros(16,1)],Opt);

ExpectedValsX365(tktktk,:)=yOut(end,6:9)+[CDF(time(0),phi(0),1)*PriorXintegral(time(0),phi(0));CDF(time(0),phi(0),2)*PriorXintegral(time(0),phi(0));CDF(time(0),phi(0),3)*PriorXintegral(time(0),phi(0));CDF(time(0),phi(0),4)*PriorXintegral(time(0),phi(0))]';

ExpectedValsX365small(tktktk,:)=yOut(end,14:17)+[CDFsmall(time(0),phi(0),1)*PriorXintegral(time(0),phi(0));CDFsmall(time(0),phi(0),2)*PriorXintegral(time(0),phi(0));CDFsmall(time(0),phi(0),3)*PriorXintegral(time(0),phi(0));CDFsmall(time(0),phi(0),4)*PriorXintegral(time(0),phi(0))]';

DoesThisSumToOne(tktktk,:)=yOut(end,10:13)+[CDF(time(0),phi(0),1);CDF(time(0),phi(0),2);CDF(time(0),phi(0),3);CDF(time(0),phi(0),4)]';

DoesThisSumToOnesmall(tktktk,:)=yOut(end,18:21)+[CDFsmall(time(0),phi(0),1);CDFsmall(time(0),phi(0),2);CDFsmall(time(0),phi(0),3);CDFsmall(time(0),phi(0),4)]';

if(isempty(tEvent))
    DelayTime(tktktk)=inf;
    DelayX(tktktk)=inf;
else
    DelayTime(tktktk)=tEvent(1);
    DelayX(tktktk)= (log(yEvent(5))-log(epsilon)+ tEvent(1)*(gamma+mu))/beta(4);
end


end

delayTFromTime= @(t) interp1(tOutReference,[DelayTime;DelayTime(1)]',t);
delayXFromTime=@(t) interp1(tOutReference ,[DelayX;DelayX(1)]',t);

%%There are multiple parts of the cycle where we could "start counting",
%%and there are multiple places where the mutant might appear.
%%In order to calculate T_{1/2}, we need to figure out the expected 
%%NUMBER of cycles before the mutant appears (which depends on exp(-nu*M)
%%at the end of a cycle), the and then also account for the one incomplete
%%cycle, and the DELAY in time between T_epsilon and T_{1/2}. THEN we need
%%to account for probabilities, which we rescale to "1" over the length of
%%one cycle (multiple cycles already accounted for.

%Also, there's probably a normalization constant unaccounted for, because
%we have both the phi and t direction. I'm pretty sure I've normalized well
%in the t direction but no in the phi direction, which is probably an extra
%factor of 2*T.

for(Mtype=1:4)
    ThalfIntegral= @(time,phi) (PDF(time,phi,Mtype)./ (1-Mcycle(Mtype)) ).*(meancycles(Mtype)*Tcycle  + time +delayTFromTime(mod(phi+time,2*T)) );
    X_T_Integral= @(time,phi) (PDF(time,phi,Mtype)./ (1-Mcycle(Mtype)) ).*(meancycles(Mtype)*Xcycle  + PriorXintegral(time,phi) +delayXFromTime(mod(phi+time,2*T)) );
    Base= @(time,phi) (PDF(time,phi,Mtype)./ (1-Mcycle(Mtype)) );

        %%Because M_{select} is discontinuous, (due to chi multiple), we
        %%need to use a integration method more robust to discontinuities.
        Scale = integral2(Base,0,2*T,0,2*T,'method','iterated');
        Thalf(lll,Mtype)=integral2(ThalfIntegral,0,2*T,0,2*T,'method','iterated')./Scale;
        X_T(lll,Mtype)=integral2(X_T_Integral,0,2*T,0,2*T,'method','iterated')./Scale;
        X_Tstar= (Xcycle./Tcycle - (gamma+mu)/beta(4))./((McycleRaw(Mtype))./Tcycle);
end

X365(lll,:)=mean(ExpectedValsX365)./mean(DoesThisSumToOne);
X365small(lll,:)=mean(ExpectedValsX365small)./mean(DoesThisSumToOnesmall);
lll
end


figure(7)
subplot(3,2,1)
plot(Tlist,X365,'lineWidth',1.5 )
xlabel('\chi_A')
ylabel('X_{365}')
title('High mutation rate')

subplot(3,2,2)
plot(Tlist,X365small,'lineWidth',1.5 )
xlabel('\chi_A')
ylabel('X_{365}')
title('Low mutation rate')

subplot(3,2,3)
plot(Tlist,Thalf,'lineWidth',1.5 )
xlabel('\chi_A')
ylabel('T_{1/2}')

subplot(3,2,4)
plot(Tlist,X_T,'lineWidth',1.5 )
xlabel('\chi_A')
ylabel('X_{T}')

subplot(3,2,5)
plot(Tlist,X_Tstar,'lineWidth',1.5 )
xlabel('\chi_A')
ylabel('X_{T*}')
title('Raw Calculation')

subplot(3,2,6)
plot(Tlist,X_Tstar./max(X_Tstar),'lineWidth',1.5 )
xlabel('\chi_A')
ylabel('X_{T*}')
legend({'import','base','select','HGT'});
title('Rescaled')
ylim([-0.1,1.1])






function [position,isterminal,direction] = ThalfEvent(t,y)
position = sum(y(1:3))-y(5); % The value that we want to be zero
isterminal = 0;  % Halt integration 
direction = 0;   % The zero can be approached from either direction
end