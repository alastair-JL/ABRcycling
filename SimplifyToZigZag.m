
P=rand(5,1);

lockedParams=[0.2,1,0.3,0.00,0.04,0.04,0.17,2,0.1,0,0,0];

mu= lockedParams(1); %discharge rate.
m=[0.1,0.001,0.001,1]; %InfluxRate S,A,B,X

beta= lockedParams(3); % 0.3 %Rate of infection
sigma=  lockedParams(4); %0.00; %relative rate of replacement infection

cA=  lockedParams(5);% 0.04;
cB=  lockedParams(6);% 0.04;
cAB=  lockedParams(7);%0.17; %Cost of antibiotic resistance.
betaS=beta;
betaA=beta*(1-cA);
betaB=beta*(1-cB);

tau=  lockedParams(8);%1.2; %rate of recovery due to treatment
gam=  lockedParams(9);% 0.1;  %Spontaneous recovery.

nuA=  lockedParams(10);%10^-6; %Spotaneous rate of resistance mutation thing.
nuB=  lockedParams(11);% 10^-6;
nuAB=  lockedParams(12);%10^-13;

sigma=sum(m)/mu;
SatX = (gam+mu+tau/2)/betaB;
z= m(1)/(gam+tau+mu-betaS*SatX);
Fulcrum= (sigma-SatX-z)/2;
slope= tau*Fulcrum/2;


T=5;

RaFunct= @(t) -slope*(t-T/2)+Fulcrum  ;
RbFunct= @(t) +slope*(t-T/2)+Fulcrum;


RaFunct= @(t) exp(-slope*(t-T/2)/Fulcrum)*Fulcrum  ;
RbFunct= @(t) exp(+slope*(t-T/2)/Fulcrum)*Fulcrum;


tvect=[linspace(0,T/10),linspace(T/10,T)];
RAvect=RaFunct(tvect);
RBvect=RbFunct(tvect);

plot([tvect,T+tvect],[RAvect,RBvect]);
hold on
plot([tvect,T+tvect],[RBvect,RAvect]);

m=[m(1:3),0,m(4)]'; %InfluxRate S,A,B,AB, X

chi=[0,0,1,0];

 ABRsim= @(t,y) ABRsimulationZZ(chi,y,lockedParams,t,T,m);
 
 y0=rand(1,5);
y0(1)=z;
y0(2)=RAvect(1);
y0(3)=RBvect(1);
y0(5)=SatX;
 y0(4)=0;
 
 [t,yQQ]=ode15s(ABRsim, linspace(0,T*20,10000), y0, odeset('NonNegative',1))
 figure(44);
 plot(t,yQQ,'LineWidth',2);
 hold on

 t_cycle=[tvect,T+tvect,2*T+tvect,3*T+tvect,4*T+tvect,5*T+tvect,6*T+tvect,7*T+tvect,8*T+tvect,9*T+tvect];
 t_cycle=[t_cycle,t_cycle+10*T]
 cyclePredict=[RAvect,RBvect,RAvect,RBvect,RAvect,RBvect,RAvect,RBvect,RAvect,RBvect];
 cyclePredict=[cyclePredict,cyclePredict]
 plot(t_cycle,cyclePredict,'o');
 
% 
% plot(t_cycle(5:5:end),cyclePredict(:,5:5:end),'o');

lowA= m(2)/(mu+gam-betaA*SatX);
lowB= m(3)/(tau+mu+gam-betaB*SatX);
figure(69);
plot(t,log(yQQ(:,2)-lowA),'LineWidth',2);
hold on
plot(t,log(yQQ(:,3)-lowB),'LineWidth',2);

 