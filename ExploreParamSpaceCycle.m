
P=rand(5,1);

lockedParams=[0.2,1,0.3,0.00,0.04,0.04,0.17,1,0.1,0,0,0];

mu= lockedParams(1); %discharge rate.
%m=[0.1,0.0001,0.0001,1]; %InfluxRate S,A,B,X
m=[0.1,0.03,0.03,1]; %InfluxRate S,A,B,X


beta= lockedParams(3); % 0.3 %Rate of infection
betaA=lockedParams(3)*(1-lockedParams(5))
betaB=lockedParams(3)*(1-lockedParams(6))
betaS=lockedParams(3)
sigma=  lockedParams(4); %0.00; %relative rate of replacement infection

tau=  lockedParams(8);%0.4; %rate of recovery due to treatment
gam=  lockedParams(9);% 0.02;  %Spontaneous recovery.


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
SatX = (gam+mu+tau/2)/betaB;
z= m(1)/(gam+tau+mu-betaS*SatX);
Fulcrum= (sigma-SatX-z)/2;

y0=[z,Fulcrum,Fulcrum,0,SatX, 0,0,0,0,0]';
m=[m(1:3),0,m(4)]';

Ts= logspace(-0.5,2,200);
Xmeans=-ones(1,length(Ts));

for(tttt=1:length(Ts))

    %T=[Ts(tttt),Ts(tttt)*1]
    T=Ts(tttt);
chi=[0,0,1,0];

 ABRsim= @(t,y) ABRsimulationParamExplore(chi,y,lockedParams,t,T,m);
 
  [t,yQQ]=ode15s(ABRsim, [-12.5*sum(T),sum(T)*10], y0', odeset('NonNegative',1,'RelTol',1e-8))
 

  tttt
  
  
  Xmeans(tttt)=yQQ(end,end)/(sum(T)*10);
  if(Xmeans(tttt)==0)
     warning('wut'); 
  end
  
end

subplot(1,2,2)
plot(t,yQQ(:,[1,2,3,5]),'LineWidth',2)

figure(9)
subplot(3,1,1)
plot(t,yQQ(:,3),'LineWidth',2)
hold on
plot([0,0],[0,5],'Color',[0.7,0.7,0.7])
plot([25,25],[0,5],'Color',[0.7,0.7,0.7])

xlim([-10,50])
xlabel('t');
ylabel('R_A');

subplot(3,1,2)
plot(t,m(3)./yQQ(:,3),'LineWidth',2)
xlim([-10,50])
xlabel('t');
ylabel('m_A/R_A');
hold on
plot([0,0],[0,tau],'Color',[0.7,0.7,0.7])
plot([25,25],[0,tau],'Color',[0.7,0.7,0.7])

plot([0,0],[m(3)/tau,5],'Color',[0.7,0.7,0.7])

plot(linspace(0,10),exp(-linspace(0,10)),'k')

figure(12)
plot(Ts,Xmeans,'LineWidth',2)
alpha=10;

spiralB= (log(topA)-log(bottomA))-1;
hold on
%plot(Ts, -log(exp(-alpha*(mu+gam+tau/2))*ones(size(Ts)) + exp(-alpha*(mu+gam + spiralB./Ts)))./alpha );
plot(Ts, min( (mu+gam+tau/2),mu+gam + spiralB./Ts)/betaB,'LineWidth',2);
