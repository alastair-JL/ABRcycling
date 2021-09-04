
P=rand(5,1);

lockedParams=[0.2,1,0.3,0.00,0.04,0.04,0.17,2,0.1,0,0,0];

mu= lockedParams(1); %discharge rate.
m=[0.1,0.0000001,0.0000001,1]; %InfluxRate S,A,B,X

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
SatX = (gam+mu)/betaB;
z= m(1)/(gam+tau+mu-betaS*SatX);
bottomA= m(2)/(gam+tau+mu-betaA*SatX);
bottomB= m(3)/(gam+tau+mu-betaB*SatX);
topA= sigma-bottomB-z;
topB= sigma-bottomA-z;


kA= (gam+mu+tau)- betaA*(sigma-z-m(3)/(gam+tau+mu-betaB*SatX) );
kB= (gam+mu)- betaB*(sigma-z-m(2)/(gam+tau+mu-betaA*SatX) );

%kA= (gam+mu+tau)- betaA*(sigma-z-topB );
%kB= (gam+mu)- betaB*(sigma-z-bottomA );

rA = roots([-betaA;-kA;m(2)]);
rA_Approx= [m(2)/kA,-kA/betaA];

rB = roots([-betaB;-kB;m(3)]);
rB_Approx= [m(3)/kB,-kB/betaB];

%%%%%%%%%%%%%5
%%Evil cheat hacking.

sigma=sum(m)/mu;
SatX = (gam+mu)/betaB;
z= m(1)/(gam+tau+mu-betaS*SatX);
bottomA= m(2)/(gam+tau+mu-betaA*SatX);
bottomB= m(3)/(gam+tau+mu-betaB*SatX);
topA= sigma-bottomB-z-SatX;
topB= sigma-bottomA-z-SatX;

rA(2)=bottomA;
rB(1)=topB;


kA= (gam+mu+tau)- betaA*(sigma-z-m(3)/(gam+tau+mu-betaB*SatX) );
kB= (gam+mu)- betaB*(sigma-z-bottomA );

kB=@(t) ((gam+mu)- betaB*(sigma-z-bottomA))*t - kA*betaB*(topA-bottomA)*exp(-kA*t)+kA*betaB*(topA-bottomA);

%%%%%%%%


T=20;

RaFunct= @(t,C) rA(1) +(rA(2)-rA(1))./(1+C*exp(-kA*t));
RbFunct= @(t,C) rB(1) +(rB(2)-rB(1))./(1+C*exp(-kB(t)));

RaFunct_C=@(t,C)  -(rA(2)-rA(1)).*exp(-kA*t)./(1+C*exp(-kA*t)).^2;
RbFunct_C= @(t,C) -(rB(2)-rB(1)).*exp(-kB(t))./(1+C*exp(-kB(t))).^2;

Cb=(rB(2)-rB(1))/(rA(2)-rB(1))-1; %Some initial approximations,
Ca= (rA(2)-rA(1))./(rB(1)-rA(1))-1;

C=[Ca;Cb]
 
for(iii=1:10)
    Vect=[RaFunct(T,C(1))-RbFunct(0,C(2));RaFunct(0,C(1))-RbFunct(T,C(2))]
    M=[RaFunct_C(T,C(1)),-RbFunct_C(0,C(2));RaFunct_C(0,C(1)),-RbFunct_C(T,C(2))]
    
    C=C - M\Vect
    
end

Vect

C

tvect=[linspace(0,T/10),linspace(T/10,T)];
RAvect=RaFunct(tvect,C(1));
RBvect=RbFunct(tvect,C(2));

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
y0(5)=sigma-sum(y0(1:3));
 y0(4)=0;
 
 [t,yQQ]=ode15s(ABRsim, linspace(0,T*20,10000), y0', odeset('NonNegative',1))
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