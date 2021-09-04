% ModelIS   [ S, R_A,R_B,X ]

yFinal=[rand(4,1);0]
% y0=[ swap*(yOut(end,1:6)') ;0];
% y0= [(yOut(end,1:6)');0];

Tlist= logspace(-1,2.3,50);
%Tlist=[0,Tlist];
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
Sequil= importation(1)/(tau+gamma+mu-beta(1)*Xequil)
Bequil= importation(3)/(tau+gamma+mu-beta(1)*Xequil)
Aequil= 1-Xequil-Sequil-Bequil

tSat= (log(Aequil)-log(Bequil)-1)*2/tau;

CycleVector=zeros(length(Tlist),5);

yFinal=[Sequil;Aequil;Bequil;Xequil;0];

extra=8;
ExitNumbers=-ones(length(Tlist),4+extra);
Maxs=-ones(length(Tlist),1);
Mins=-ones(length(Tlist),1);
Means=-ones(length(Tlist),1);
meanDownSwing=-ones(length(Tlist),1);
meanUpSwing=-ones(length(Tlist),1);

epsilon=10^-2.5;

fineness=250;

nu=0.1*[1;2;7;500]; %These are the scaling factors on our various 
%mutation rates. and are basically picked so as to keep the various
%mutation channels RELATIVELY comparable; IE   M_base is often around
%0.5ish, M_HGT is often closer to 1/100
nuSmall=0.2*nu;
nuSmall=0.01*nu;

X365=ones(length(Tlist),4);
X365small=ones(length(Tlist),4);
Thalf=ones(length(Tlist),4);
X_T=ones(length(Tlist),4);
X_Tstar=ones(length(Tlist),4);
X_Tstar_Bastard=ones(length(Tlist),4);
X_Tstar_small=ones(length(Tlist),4);

XbarRecord=ones(length(Tlist),1);
MbarRecord=ones(length(Tlist),4);
expectedMutationTimesRecord=ones(length(Tlist),4);
expectedMutationTimesRecordSmall=ones(length(Tlist),4);
for(lll=length(Tlist):-1:(1+(Tlist(1)==0)))

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

M= @(V,t) [1;(V(3)+V(2));(V(3)*chi(t)+(V(2)*(1-chi(t))));(V(2)*V(3))];

Deriv =@(t,V) [(importation + betaMatrix*V(1:5)*V(4) + recovery(t)*V(1:5)-mu*V(1:5));...
        1-nu.*M(V,t).*V(6:9); 1-nuSmall.*M(V,t).*V(10:13);....
        V(4).*V(6:9)- nu.*M(V,t).*V(14:17);V(4).*V(10:13)- nuSmall.*M(V,t).*V(18:21); nu.*M(V,t);V(4);nuSmall.*M(V,t)];

warmupTime= ceil(500/(2*T))*(2*T);

[tOut,yOut] = ode45(Deriv,linspace(0,warmupTime,40000),[yFinal; 10*ones(4,1)/nu(1); 10*ones(4,1)/nuSmall(1); zeros(17,1)]);
[tOut,yOut] = ode45(Deriv,linspace(0,warmupTime,40000),yOut(end,:)');
%plot(tOut,yOut(6:9))
%plot(tOut,yOut(10:13))
yOut(end,5+4+(1:4))= (tOut(end)-tOut(1))./((yOut(end,end-4+(1:4)))-(yOut(1,end-4+(1:4))));
warmupTime= ceil(3500/(2*T))*(2*T);
[tOut,yOut] = ode45(Deriv,linspace(0,warmupTime,40000),yOut(end,:)');

[tOutWoop,yOutWoop] = ode45(Deriv,[0,T],[[yOut(end,1:4),0]';yOut(end,6:end)']);
Mixer=eye(5);
Mixer(2:3,2:3)=[0,1;1,0];
Mixer(5,5)=0;

for(iii=1:20)
[tOutWoop,yOutWoop] = ode45(Deriv,[0,T],[0.5*(Mixer*yOutWoop(end,1:5)'+yOutWoop(1,1:5)');yOutWoop(end,6:end)']);
end
if(max(abs(Mixer*[yOutWoop(end,1:5)]'-yOutWoop(1,1:5)'))>0.0001)
    for(iii=1:200)
        [tOutWoop,yOutWoop] = ode45(Deriv,[0,T],[0.5*(Mixer*yOutWoop(end,1:5)'+yOutWoop(1,1:5)');yOutWoop(end,6:end)']);
    end
    if(max(abs(Mixer*[yOutWoop(end,1:5)]'-yOutWoop(1,1:5)'))>0.0001)
        warning('have not converged to symetry!');
    end
end

yFinal=Mixer*[yOutWoop(end,1:5)]';

    chi=@(t) ceil( mod(t/T,2)+10^-8)-1;

Deriv =@(t,V) [(importation + betaMatrix*V(1:5)*V(4) + recovery(t)*V(1:5)-mu*V(1:5));...
        1-nu.*M(V,t).*V(6:9); 1-nuSmall.*M(V,t).*V(10:13);....
        V(4).*V(6:9)- nu.*M(V,t).*V(14:17);V(4).*V(10:13)- nuSmall.*M(V,t).*V(18:21);...
        V(1:13);nu.*M(V,t).*V(14:17);nuSmall.*M(V,t).*V(18:21); M(V,t)];

[tOutReference,yOutReference] = ode45(Deriv,linspace(0,2*T,fineness+1),[yOutWoop(1,1:21)';zeros(25,1)]);
%This is going to double count t=0=T. That's fine.

expectedMutationTimes= yOutReference(end,(21+5+(1:4)))/(2*T)
expectedMutationTimesSmall= yOutReference(end,(21+9+(1:4)))/(2*T)

expectedMutationXintegral= yOutReference(end,(21+8+5+(1:4)))/(2*T)
expectedMutationXintegralSmall= yOutReference(end,(21+8+9+(1:4)))/(2*T)

Xcycle=yOutReference(end,(21+4));

CycleVector(lll,:)=yOutReference(end,(21+(1:5)));

if(   abs(log(CycleVector(lll,2))-log(CycleVector(lll,3)))>log(4)  )
    error("well that ain't balanced");
end

% Mcycle=yOutReference(end,22:25)' 
% Mcycle(1)=2*T;
% Mcycle(2)=yOutReference(end,(21+2))+yOutReference(end,(21+3)); 
% Mcycle(3)=yOutReference(end,(21+2)); 
% Mcycle(4)=yOutReference(end,end); 
Mcycle=yOutReference(end,42+(1:4))'; 

TestThing=expectedMutationTimesSmall'.*nuSmall.*Mcycle/(2*T);

yOutReferenceMutant=zeros(size(yOutReference,1),4);
yOutReferenceMutant(:,1)=1;
yOutReferenceMutant(:,2)=yOutReference(:,3)+yOutReference(:,2);
%yOutReferenceMutant(:,3)=chi(tOutReference).*yOutReference(:,3) + (1-chi(tOutReference)).*yOutReference(:,2);
yOutReferenceMutant(:,3)=yOutReference(:,2);
yOutReferenceMutant(:,4)=yOutReference(:,3).*yOutReference(:,2);

ThreadDensity= @(time,Mtype) interp1(tOutReference,yOutReference(:,Mtype+5),time) ;
ThreadDensity_small= @(time,Mtype) interp1(tOutReference,yOutReference(:,Mtype+9),time) ;
PDF= @(time,Mtype) nu(Mtype).*interp1(tOutReference,yOutReferenceMutant(:,Mtype),time,'linear','extrap').*ThreadDensity(time,Mtype);
PDF_small= @(time,Mtype) nuSmall(Mtype).*interp1(tOutReference,yOutReferenceMutant(:,Mtype),time,'linear','extrap').*ThreadDensity_small(time,Mtype);

Zfunction= @(time,Mtype) interp1(tOutReference,yOutReference(:,Mtype+13),time) ;
Xintegral= @(time) interp1(tOut,yOut(:,26),time) ;

 Mtype=3;
% 
% bigM= @(t) interp1(tOut,yOut(:,21+Mtype),t)
% 
% tau=T/3;
% toIntegrate= @(phi) exp(-bigM(tau+ warmupTime-T)+bigM(tau+ warmupTime-T-phi)) ;
% DirectIntegral=integral(toIntegrate,0,(tau+ warmupTime-T))/T
% indirectMethod=ThreadDensity(tau,Mtype)
% phis=linspace(0,(tau+ warmupTime-T),1000);
% GO TO "Scratchpad 3.
%% % % This section pretty much shows that our analytic approximations are... reasonablly appropriate. Cool.



yInitResists=[0.05*(1-(gamma+mu)/beta(4));0.05*(1-(gamma+mu)/beta(4));0.05*(1-(gamma+mu)/beta(4));(gamma+mu)/beta(4);0.85*(1-(gamma+mu)/beta(4));0];

Deriv =@(t,V) [(importation + betaMatrix*V(1:5)*V(4) + recovery(t)*V(1:5)-mu*V(1:5)); log(V(5))];
[tOutReferenceResists,yOutReferenceResists] = ode45(Deriv,[0,warmupTime],yInitResists');
yOutReferenceResists(end,end)=0;
[tOutReferenceResists,yOutReferenceResists] = ode45(Deriv,[0,T],yOutReferenceResists(end,:)');

X365LogTerm=(yOutReferenceResists(end,6)/T-log(epsilon))/beta(4);



DelayTime=zeros(fineness,1);
DelayX=zeros(fineness,1);

chi=@(t) ceil( mod(t/T,2)+10^-8)-1;


X365terms=-ones(fineness,4);
X365termsSmall=-ones(fineness,4);
    
for(tktktk= [1:fineness]  )
    %%Assume that mutant occurs at time %tOut(tktktk)% During the loop, and
    %%Solve from there to 365.
    [tktktk,lll]
    yStart=yOutReference(tktktk,1:5);
    yStart(5)=epsilon;
    
    
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

    
[~,yOut] = ode45(DerivWeave,linspace(0,365,4000),[yStart';zeros(8,1)]);

X365terms(tktktk,:)=yOut(end,6:9);
X365termsSmall(tktktk,:)=yOut(end,10:13);
    
if(Xcycle./(2*T)< (mu+gamma)/beta(4))
%if(true)
     DelayTime(:)=inf;
     DelayX(:)=inf;
else
    
DerivWeave =@(t,V) [(importation + betaMatrix*V(1:5)*V(4) + recovery(t+tOutReference(tktktk) )*V(1:5)-mu*V(1:5))];

    Opt    = odeset('Events', @ThalfEvent);
    
[tOut,yOut,tEvent,yEvent,~] = ode45(DerivWeave,linspace(0,365,4000),[yStart'],Opt);

 DelayTime(tktktk)=tEvent(1);
 DelayX(tktktk)= (log(yEvent(1,5))-log(epsilon)+ tEvent(1)*(gamma+mu))/beta(4);
end%%This is the end of the "If-Else statement that determines if Delay is infinite.


end %%This is the end of the tktktk for loops

delayTFromTime= @(t) interp1(tOutReference,[DelayTime;DelayTime(1)]',t);
delayXFromTime=@(t) interp1(tOutReference ,[DelayX;DelayX(1)]',t);

for(Mtype=1:4)
    
    if(expectedMutationTimes(Mtype)>200)
       warning('this could be a problem'); 
    end
%             
% expectedMutationTimes= yOutReference(end,(21+5+(1:4)))
% expectedMutationTimesSmall= yOutReference(end,(21+9+(1:4)))
% 
% expectedMutationXintegral= yOutReference(end,(21+8+5+(1:4)))/T
% expectedMutationXintegralSmall= yOutReference(end,(21+8+9+(1:4)))/T
            
%These are pretty detailed NUMERICAL approximations.
            if(max(DelayTime)>9999999)    
                Thalf(lll,Mtype)= inf;
                X_T(lll,Mtype)= inf;
                %For these two, I may be assuming/approximating too much.
                %Consider revising...
                X365(lll,Mtype)= Xcycle./(2*T)*365;
                X365small(lll,Mtype)=Xcycle./(2*T)*365;
            else
                    toIntegrate= @(t) delayTFromTime(t).*PDF(t,Mtype);
                    Thalf(lll,Mtype)= expectedMutationTimes(Mtype) + integral(toIntegrate,0,2*T);
            
                    toIntegrate= @(t) delayXFromTime(t).*PDF(t,Mtype);
                    X_T(lll,Mtype)= expectedMutationXintegral(Mtype) + integral(toIntegrate,0,2*T);
    
                    
                    toIntegrate= @(t) interp1(tOutReference,[X365terms(:,Mtype);X365terms(1,Mtype)]',t);
                    X365(lll,Mtype)= integral(toIntegrate,2*T);
                    toIntegrate= @(t) interp1(tOutReference,[X365terms(:,Mtype);X365termsSmall(1,Mtype)]',t);
                    X365small(lll,Mtype)=integral(toIntegrate,2*T);
            end
            
            if(X365(lll,Mtype)<0)
                warning('lul,wut');
            end
            
            %%This is exact.
            X_Tstar(lll,Mtype)= (Xcycle./(2*T) - (gamma+mu)/beta(4)).*(expectedMutationTimes(Mtype)); 
            X_Tstar_Bastard(lll,Mtype)= (Xcycle-2*T*(gamma+mu)/beta(4))./Mcycle(Mtype);
            X_Tstar_small(lll,Mtype)= (Xcycle./(2*T) - (gamma+mu)/beta(4)).*(expectedMutationTimesSmall(Mtype)); 
            
            XbarRecord(lll)=Xcycle/(2*T);
            MbarRecord(lll,Mtype)=Mcycle(Mtype)/(2*T);
            expectedMutationTimesRecord(lll,:)=expectedMutationTimes;
            expectedMutationTimesRecordSmall(lll,:)=expectedMutationTimesSmall;
end

lll
end














if((Tlist(1)==0))
lll=1;

T=5;

chi=@(t) 0.5;

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

M= @(V,t) [1;(V(3)+V(2));(V(3)*chi(t)+(V(2)*(1-chi(t))));(V(2)*V(3))];

Deriv =@(t,V) [(importation + betaMatrix*V(1:5)*V(4) + recovery(t)*V(1:5)-mu*V(1:5));...
        1-nu.*M(V,t).*V(6:9); 1-nuSmall(1)*M(V,t).*V(10:13);....
        V(4).*V(6:9)- nu.*M(V,t).*V(14:17);V(4).*V(10:13)- nuSmall.*M(V,t).*V(18:21); nu.*M(V,t);V(4)];


warmupTime= ceil(1500/(2*T))*(2*T);

[tOut,yOut] = ode45(Deriv,linspace(0,warmupTime,40000),[yFinal;zeros(21,1)]);

[tOutWoop,yOutWoop] = ode45(Deriv,[0,T],[[yOut(end,1:4),0]';yOut(end,6:end)']);
Mixer=eye(5);
Mixer(2:3,2:3)=[0,1;1,0];
Mixer(5,5)=0;

for(iii=1:20)
[tOutWoop,yOutWoop] = ode45(Deriv,[0,T],[0.5*(Mixer*yOutWoop(end,1:5)'+yOutWoop(1,1:5)');yOutWoop(end,6:end)']);
end

yFinal=Mixer*[yOutWoop(end,1:5)]';

chi=@(t) 0.5

Deriv =@(t,V) [(importation + betaMatrix*V(1:5)*V(4) + recovery(t)*V(1:5)-mu*V(1:5));...
        1-nu.*M(V,t).*V(6:9); 1-nuSmall.*M(V,t).*V(10:13);....
        V(4).*V(6:9)- nu.*M(V,t).*V(14:17);V(4).*V(10:13)- nuSmall.*M(V,t).*V(18:21);V(1:13);nu.*M(V,t).*V(14:17);nuSmall.*M(V,t).*V(18:21)];

[tOutReference,yOutReference] = ode45(Deriv,linspace(0,2*T,fineness+1),[yOut(end,1:21)';zeros(21,1)]);
%This is going to double count t=0=T. That's fine.

expectedMutationTimes= yOutReference(end,(21+5+(1:4)))/(2*T);
expectedMutationTimesSmall= yOutReference(end,(21+9+(1:4)))/(2*T);

expectedMutationXintegral= yOutReference(end,(21+8+5+(1:4)))/(2*T);
expectedMutationXintegralSmall= yOutReference(end,(21+8+9+(1:4)))/(2*T);

Xcycle=yOutReference(end,(21+4));


yOutReferenceMutant=zeros(size(yOutReference,1),4);
yOutReferenceMutant(:,1)=1;
yOutReferenceMutant(:,2)=yOutReference(:,3)+yOutReference(:,2);
%yOutReferenceMutant(:,3)=chi(tOutReference).*yOutReference(:,3) + (1-chi(tOutReference)).*yOutReference(:,2);
yOutReferenceMutant(:,3)=yOutReference(:,2);
yOutReferenceMutant(:,4)=yOutReference(:,3).*yOutReference(:,2);

ThreadDensity= @(time,Mtype) interp1(tOutReference,yOutReference(:,Mtype+5),time) ;
ThreadDensity_small= @(time,Mtype) interp1(tOutReference,yOutReference(:,Mtype+9),time) ;
PDF= @(time,Mtype) nu(Mtype).*interp1(tOutReference,yOutReferenceMutant(:,Mtype),time,'linear','extrap').*ThreadDensity(time,Mtype);
PDF_small= @(time,Mtype) nuSmall(Mtype).*interp1(tOutReference,yOutReferenceMutant(:,Mtype),time,'linear','extrap').*ThreadDensity_small(time,Mtype);

Zfunction= @(time,Mtype) interp1(tOutReference,yOutReference(:,Mtype+13),time) ;
Xintegral= @(time) interp1(tOut,yOut(:,26),time) ;

 Mtype=3;
% 
% bigM= @(t) interp1(tOut,yOut(:,21+Mtype),t)
% 
% tau=T/3;
% toIntegrate= @(phi) exp(-bigM(tau+ warmupTime-T)+bigM(tau+ warmupTime-T-phi)) ;
% DirectIntegral=integral(toIntegrate,0,(tau+ warmupTime-T))/T
% indirectMethod=ThreadDensity(tau,Mtype)
% phis=linspace(0,(tau+ warmupTime-T),1000);
% GO TO "Scratchpad 3.
%% % % This section pretty much shows that our analytic approximations are... reasonablly appropriate. Cool.



yInitResists=[0.05*(1-(gamma+mu)/beta(4));0.05*(1-(gamma+mu)/beta(4));0.05*(1-(gamma+mu)/beta(4));(gamma+mu)/beta(4);0.85*(1-(gamma+mu)/beta(4));0];

Deriv =@(t,V) [(importation + betaMatrix*V(1:5)*V(4) + recovery(t)*V(1:5)-mu*V(1:5)); log(V(5))];
[tOut,yOut] = ode45(Deriv,[0,warmupTime],yInitResists');
[tOutReferenceResists,yOutReferenceResists] = ode45(Deriv,[0,T],yOut(end,:)');

X365LogTerm=(yOutReferenceResists(end,6)/T-log(epsilon))/beta(4);



DelayTime=zeros(fineness,1);
DelayX=zeros(fineness,1);


chi=@(t) 0.5;

    
for(tktktk= [1:fineness]  )
    %%Assume that mutant occurs at time %tOut(tktktk)% During the loop, and
    %%Solve from there to 365.
    [tktktk,lll]
    yStart=yOutReference(tktktk,1:5);
    yStart(5)=epsilon;

    
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

    
[~,yOut] = ode45(DerivWeave,linspace(0,365,4000),[yStart';zeros(8,1)]);

X365terms(tktktk,:)=yOut(end,6:9);
X365termsSmall(tktktk,:)=yOut(end,10:13);



if(Xcycle./(2*T)< (mu+gamma)/beta(4))
%if(true)
    DelayTime(:)=inf;
     DelayX(:)=inf;
else

DerivWeave =@(t,V) [(importation + betaMatrix*V(1:5)*V(4) + recovery(t+tOutReference(tktktk) )*V(1:5)-mu*V(1:5))];

    Opt    = odeset('Events', @ThalfEvent);
    
[tOut,yOut,tEvent,yEvent,~] = ode45(DerivWeave,linspace(0,365,4000),[yStart'],Opt);

 DelayTime(tktktk)=tEvent(1);
 DelayX(tktktk)= (log(yEvent(1,5))-log(epsilon)+ tEvent(1)*(gamma+mu))/beta(4);
end%%This is the end of the "If-Else statement that determines if Delay is infinite.

end 

delayTFromTime= @(t) interp1(tOutReference,[DelayTime;DelayTime(1)]',t);
delayXFromTime=@(t) interp1(tOutReference ,[DelayX;DelayX(1)]',t);

for(Mtype=1:4)
    
    if(expectedMutationTimes(Mtype)>200)
       warning('this could be a problem'); 
    end
%             
% expectedMutationTimes= yOutReference(end,(21+5+(1:4)))
% expectedMutationTimesSmall= yOutReference(end,(21+9+(1:4)))
% 
% expectedMutationXintegral= yOutReference(end,(21+8+5+(1:4)))/T
% expectedMutationXintegralSmall= yOutReference(end,(21+8+9+(1:4)))/T
            
%These are pretty detailed NUMERICAL approximations.
            toIntegrate= @(t) delayTFromTime(t).*PDF(t,Mtype);
            Thalf(lll,Mtype)= expectedMutationTimes(Mtype) + integral(toIntegrate,0,T);
            
            toIntegrate= @(t) delayXFromTime(t).*PDF(t,Mtype);
            X_T(lll,Mtype)= expectedMutationXintegral(Mtype) + integral(toIntegrate,0,T);
            
            %For these two, I may be assuming/approximating too much.
            %Consider revising...
            
                    toIntegrate= @(t) interp1(tOutReference,[X365terms(:,Mtype);X365terms(1,Mtype)]',t);
                    X365(lll,Mtype)= integral(toIntegrate,2*T);
                    toIntegrate= @(t) interp1(tOutReference,[X365terms(:,Mtype);X365termsSmall(1,Mtype)]',t);
                    X365small(lll,Mtype)=integral(toIntegrate,2*T);
                    
                    
            %%This is exact.
            X_Tstar(lll,Mtype)= (Xcycle./(2*T) - (gamma+mu)/beta(4)).*(expectedMutationTimes(Mtype)); 
            X_Tstar_small(lll,Mtype)= (Xcycle./(2*T) - (gamma+mu)/beta(4)).*(expectedMutationTimesSmall(Mtype)); 
            X_Tstar_Bastard(lll,Mtype)= (Xcycle-2*T*(gamma+mu)/beta(4))./Mcycle(Mtype);
            
end
            expectedMutationTimesRecord(lll,:)=expectedMutationTimes;
            expectedMutationTimesRecordSmall(lll,:)=expectedMutationTimesSmall;
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

% subplot(3,2,5)
% plot(Tlist,X_Tstar,'lineWidth',1.5 )
% xlabel('\chi_A')
% ylabel('X_{T*}')
% title('Raw Calculation')
% hold on
% plot([tSat,tSat],[0,1.1],'k');
% 
subplot(3,2,5)
plot(Tlist,X_Tstar_Bastard,'lineWidth',1.5 )
xlabel('\chi_A')
ylabel('X_{T*}')
title('Bastard simplification')
hold on
plot([tSat,tSat],[0,1.1],'k');

subplot(3,2,6)
plot(Tlist,X_Tstar./max(X_Tstar),'lineWidth',1.5 )
xlabel('\chi_A')
ylabel('X_{T*}')
legend({'import','base','select','HGT'});
title('Rescaled')
ylim([-0.1,1.1])
hold on
plot([tSat,tSat],[0,1.1],'k');


figure(8)
subplot(1,3,1)
plot(Tlist,X_Tstar_Bastard./max(X_Tstar_Bastard),'lineWidth',1.5 )
xlabel('\chi_A')
ylabel('X_{T*}')
title('Bastard simplification')
hold on
plot([tSat,tSat],[0,1.1],'k');

subplot(1,3,2)
plot(Tlist,X_Tstar_small./max(X_Tstar_small),'lineWidth',1.5 )
xlabel('T')
ylabel('X_{T*} SMALL')
legend({'import','base','select','HGT'});
title('Rescaled')
ylim([-0.1,10.1])
hold on
plot([tSat,tSat],[0,1.1],'k');


subplot(1,3,3)
plot(Tlist,X_Tstar./max(X_Tstar),'lineWidth',1.5 )
xlabel('\chi_A')
ylabel('X_{T*}')
legend({'import','base','select','HGT'});
title('Rescaled')
ylim([-0.1,10.1])
hold on
plot([tSat,tSat],[0,1.1],'k');

figure(12)
plot(Tlist,CycleVector)

figure(57)
subplot(1,3,1)
plot(Tlist,MbarRecord(1,:)./MbarRecord,'lineWidth',1.5 )
subplot(1,3,2)
plot(Tlist,expectedMutationTimesRecord./expectedMutationTimesRecord(1,:),'lineWidth',1.5 )
subplot(1,3,3)
plot(Tlist,expectedMutationTimesRecordSmall./expectedMutationTimesRecordSmall(1,:),'lineWidth',1.5 )


function [position,isterminal,direction] = ThalfEvent(t,y)
position = sum(y(1:3))-y(5); % The value that we want to be zero
isterminal = 1;  % Halt integration 
direction = 0;   % The zero can be approached from either direction
end