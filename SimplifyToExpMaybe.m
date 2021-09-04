
P=rand(5,1);

lockedParams=[0.02,1,0.1,0.02,0.04,0.04,0.17,0.75,0.01,10^-6,10^-6,10^-13];
lockedParams=[0.02,1,0.1,0.02,0.04,0.04,0.17,1.2,0.04,10^-6,10^-6,0];

chi=[0,0,1,0];
P=rand(5,1);
P(4)=0.5;
P=P/sum(P);

for(iii=1:20)
[d,J] = ABRsimulation(chi,P,lockedParams);
P=P-J\d;
end

y0=P;

s=y0(2);
y0(2)=y0(3);
y0(3)=s;

tspan=[0,300];

ABRsim= @(t,y) ABRsimulation(chi,y,lockedParams)

[t,y]=ode45(ABRsim, tspan, y0)
 

%%Okay, make the local approximation y'=A+By+Cy^2
A= lockedParams(2)*lockedParams(1)/3;
B= -lockedParams(8)-lockedParams(9) - lockedParams(1)+ lockedParams(3)*(1-lockedParams(5))*(lockedParams(2)-y(1,1)-y(1,3)-y(1,4));
C= lockedParams(3)*(1-lockedParams(5));

rA = roots([C;B;A]);
rA_Approx= [C/B,B/A]

mid=mean(rA);
Derive= C*mid*mid + B*mid+A;
qA= -4*Derive/(rA(2)-rA(1));
FA= (y(1,2)-rA(1))/(rA(2)-y(1,2));

PredictionA= @(t) rA(1) + (rA(2)-rA(1))*FA./(FA + exp(qA*t));

hold on
plot(t, PredictionA(t),'r:' ,'lineWidth',2);



%figure(4)
% plot(t,y(:,3));

%% Simplified equation  --->    y' = Ay^2 + By +C
%% Step one: Find the params!
%% 
%%%%%Professor Neil Johnston   
%Mathematical Biologist imperial college 
%

%%Okay, make the local approximation y'=A+By+Cy^2
A= lockedParams(2)*lockedParams(1)/3;
B= -lockedParams(9) - lockedParams(1)+ lockedParams(3)*(1-lockedParams(5))*(lockedParams(2)-y(end,1)-y(end,2)-y(end,4));
C= -lockedParams(3)*(1-lockedParams(5));

rB = roots([C;B;A]);

mid=mean(rB);
Derive= C*mid*mid + B*mid+A;
qB= -4*Derive/(rB(2)-rB(1));
FB= (y(1,3)-rB(1))/(rB(2)-y(1,3));

PredictionB= @(t) rB(1) + (rB(2)-rB(1))*FB./(FB + exp(qB*t));

%hold on
%plot(t, PredictionB(t),'r:','lineWidth',2)


%%This gives a NOT BAD approximation of the actual curve we care about. 
%%Not great... but not bad either.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

%figure(5)
%plot(t,y(:,5));
%hold on
%plot(t, (lockedParams(2)-y(end,1)-y(end,4)-PredictionB(t)-PredictionA(t) ),'r:','lineWidth',2);


%figure(6)

Integrator=[t(2)-t(1); t(3:end)-t(1:(end-2)); t(end)-t(end-1)]/2

%plot(t, cumsum(Integrator.*y(:,5)) )
integralApproximation= @(T) T*(lockedParams(2)-y(end,1)-y(end,4) -rA(2)-rB(2)) + (rB(2)-rB(1))/qB*(log(FB + exp(qB*T))-log(FB + 1)) + (rA(2)-rA(1))/qA*(log(FA + exp(qA*T))-log(FA + 1));

%hold on
%plot(t, integralApproximation(t),'r:','lineWidth',2);

%%%%%%%%%%%%%%%%%%%%%

%Okay, and the next step, if I REALLY want to make use of this
%approximation, is to make it circular- IE, pick a time, and then ensure
%that the Left boundary of one function matches the Right hand boundary of
%the other.

%% Probably easiest just to start with the approximation of a collapsed 
%fast exp, [A(inf)=B(0)] then feed that into the slow growth of the resistant strain
%And then use [B(T)= A(0)], at which point A(T) is approx A(inf) anyway.

%Or, alternatively, I can just try to solve exactly. That's a thing.

%
% rA(2) - (rA(2)-rA(1))*exp(qA*T)/(FA+exp(qA*T))= rB(2) -(rB(2)-rB(1))/(FB+1)
% rB(2) - (rB(2)-rB(1))*exp(qB*T)/(FB+exp(qB*T))= rA(2) -(rA(2)-rA(1))/(FA+1)
%
%  ==> Algebra  ==>
%
% 0= FA*(exp(qB*T)-1)*(rB(1)-rA(2)) +FB*(exp(qA*T)-1)*(rA(1)-rB(2)) + (exp(qB*T)-exp(qA*T))*(rB(1)-rA(1))
% 0= FA*FB*(rB(2)-rA(2)) + FA*exp(qB*T)*(rB(1)-rA(2)) - FB*(rA(1)-rB(2))+exp(qB*T)*(rB(1)-rA(1))

clear FB
clear FA

T=100;

FAfunc= @(FB,T) -(FB.*(exp(qA*T)-1)*(rA(1)-rB(2)) + (exp(qB*T)-exp(qA*T))*(rB(1)-rA(1)) )/( (exp(qB*T)-1)*(rB(1)-rA(2)) );
Combine= @(FA,FB,T) FA.*FB.*(rB(2)-rA(2)) + FA.*exp(qB*T).*(rB(1)-rA(2)) - FB.*(rA(1)-rB(2))+exp(qB*T)*(rB(1)-rA(1));


FB= linspace(0,10);
plot(FB,FAfunc(FB,T))
plot(FB,Combine(FAfunc(FB,T),FB,T))

%Okay, so the algebra checks out, I now just need to find roots. The
%Co-efficients are horrible to calculate, but fundementally, all we have
%here is quadratic with very very ugly co-efficients. So...

coEfficients= [-(rB(2)-rA(2))*(exp(qA*T)-1)*(rA(1)-rB(2))/((exp(qB*T)-1)*(rB(1)-rA(2)));
               -(rB(2)-rA(2))*(exp(qB*T)-exp(qA*T))*(rB(1)-rA(1))/( (exp(qB*T)-1)*(rB(1)-rA(2)))-exp(qB*T)*(rB(1)-rA(2))*((exp(qA*T)-1)*(rA(1)-rB(2)))/( (exp(qB*T)-1)*(rB(1)-rA(2)))-(rA(1)-rB(2));
               exp(qB*T)*(rB(1)-rA(1))-exp(qB*T)*(rB(1)-rA(2))*((exp(qB*T)-exp(qA*T))*(rB(1)-rA(1)))/((exp(qB*T)-1)*(rB(1)-rA(2)))];
           
    FBvals= roots(coEfficients)      
  
    FB=max(FBvals);
    FA=FAfunc(FB,T);

    if(T>10)
        t=[linspace(0,4,50),linspace(4,T,30)];
    else
        t=[linspace(0,T,ceil(15*T))];
    end
    
PredictionB= @(t) rB(1) + (rB(2)-rB(1))*FB./(FB + exp(qB*t));
PredictionA= @(t) rA(1) + (rA(2)-rA(1))*FA./(FA + exp(qA*t));

cyclePredict= [PredictionB(t),PredictionA(t),PredictionB(t),PredictionA(t),PredictionB(t),PredictionA(t),PredictionB(t),PredictionA(t),PredictionB(t),PredictionA(t);
               PredictionA(t),PredictionB(t),PredictionA(t),PredictionB(t),PredictionA(t),PredictionB(t),PredictionA(t),PredictionB(t),PredictionA(t),PredictionB(t)];
t_cycle= [t,t+T,t+2*T,t+3*T,t+4*T,t+5*T,t+6*T,t+7*T,t+8*T,t+9*T];

%figure(42)
%plot(t_cycle,cyclePredict);

%Next steps: Compare to the ODEs
%NOTE: I suspect I am doing the initial boundary conditions wrong, in that
%I am using the original ODE to produce boundary conditions, and then doing
%then switch/matching my FA and FB parameters, but not others. This may be
%problem?

%This is PROBABLY a problem.

%Also, remember to step back and look at system from a wider angle.
%Currently playing to much in the details.


ABRsim= @(t,y) ABRsimulationCycle(chi,y,lockedParams,t,T)

y0=[min(cyclePredict(:,1));cyclePredict(:,1);0;1-sum([min(cyclePredict(:,1));cyclePredict(:,1)]) ]

[t,yQQ]=ode45(ABRsim, [0,T*10], y0, odeset('NonNegative',1))
figure(44);
plot(t,yQQ,'LineWidth',2);
hold on
%plot(t_cycle,cyclePredict,'*');


plot(t_cycle(5:5:end),cyclePredict(:,5:5:end),'o');