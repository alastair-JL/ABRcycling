
lockedParams=[0.2,1,0.3,0.00,0.04,0.04,0.17,1,0.1,0,0,0];

mu= lockedParams(1); %discharge rate.
%m=[0.1,0.001,0.001,1]; %InfluxRate S,A,B,X
tau=  lockedParams(8);%1.2; %rate of recovery due to treatment
gam=  lockedParams(9);% 0.1;  %Spontaneous recovery.


betaS= lockedParams(3); % 0.3 %Rate of infection
betaA=lockedParams(3)*(1-lockedParams(5))
betaB=lockedParams(3)*(1-lockedParams(6))

sigma=sum(m)/mu;
SatX = (gam+mu)/betaB;
z= m(1)/(gam+tau+mu-betaS*SatX);
bottomA= m(2)/(gam+tau+mu-betaA*SatX);
BeneathA= -(gam+tau+mu-betaA*(sigma-z))/betaA;
bottomB= m(3)/(gam+tau+mu-betaB*SatX);
topA= sigma-bottomB-z;
topB= sigma-bottomA-z;


FlowDownA= @(t) bottomA + (topA-bottomA)*exp(-(tau)*t);
FlowDownB= @(t) bottomA + (topA-bottomA)*exp(-(gam+tau+mu-betaA*sigma)*t);
FlowDownC= @(t) bottomA + (topA*(gam+tau+mu-betaA*(sigma-z))/tau-bottomA)*(1./((gam+tau+mu-betaA*(sigma-z))/tau-1+exp((gam+tau+mu-betaA*(sigma-z))*t)) );
FlowUpA= @(t) bottomA*exp(tau*t);


figure(5)
plot(t,yQQ(:,2:3))
hold on
timeVect=linspace(-1,25);
plot(timeVect,FlowDownA(timeVect),'o');
plot(timeVect,FlowDownB(timeVect),'o');
plot(timeVect,FlowDownC(timeVect),'o');

plot(timeVect,FlowUpA(timeVect),'x');

figure(7)
plot(t,m(2)./yQQ(:,2:3))
hold on
timeVect=linspace(-1,25);
plot(timeVect,m(2)./FlowDownA(timeVect),'o');
plot(timeVect,m(2)./FlowDownB(timeVect),'o');
plot(timeVect,m(2)./FlowDownC(timeVect),'o');
plot(timeVect,m(2)./FlowUpA(timeVect),'x');

figure(13)
plot(t,log(yQQ(:,2:3)))
hold on
timeVect=linspace(-1,25);
plot(timeVect,log(FlowDownA(timeVect)),'o');
plot(timeVect,log(FlowDownB(timeVect)),'o');
plot(timeVect,log(FlowDownC(timeVect)),'o');
plot(timeVect,log(FlowUpA(timeVect)),'x');


IntralThing_mOverRA= cumsum((m(2)./yQQ(2:end,2)+m(2)./yQQ(1:(end-1),2)).*diff(t))/2;
IntralThing_X= cumsum((yQQ(2:end,5)+yQQ(1:(end-1),5)).*diff(t))/2;

plot(t(2:end),IntralThing_X*betaA)
hold on
plot(t(2:end),IntralThing_mOverRA)
plot(t(2:end),(t(2:end)-t(1))*(tau/2+mu+gam))

plot(t(2:end),IntralThing_X*betaA+IntralThing_mOverRA) 

%%Okay, so this little plotting scheme is good confirmation that what we
%%are doing *works*. That doesn't mean we have the algebra to make it work
%%well, and it doesn't mean that this is an approximation that can be
%%justified by the real world situation, but it is a thing.

figure(14)
plot(t,betaA*yQQ(:,end)+m(2)*yQQ(:,7))
hold on
plot(t,yQQ(:,6)+log(yQQ(:,2)))
