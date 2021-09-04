%%
%%
%%
%% Here I attempt to implement the ODE's as used in Uecker's review paper:
%    "Antibiotic treatment protocols revisited: The challenges
%     of a conclusive assessment by mathematical modeling"
%
% All the models originally came from elsewhere. I don't care.

%
%Things to note: It seems silly to measure all strategies based on all
%things at once. We can Determine which strategy will be best AFTER ABR has
%taken hold. 
%Then, using the value of that strategy, we can look at how much PRE-ABR
%strategies improve on that (The "Delta-utility" as it were). We then look
%at how long we expect them to last - (IE, Delta-utility * time till ABR).
%
% We might also be persuaded to look at the immediate aftermath of ABR
% arrival.
%
% But ranking based on time of arrival doesn't make sense unless we know
% how valuable each day of ABR is going to "cost". 
%
% WE might also hypothesize that there is just plain dimishing returns in

% With all that in mind....


S=1;
RA=2;
RB=3;
RAB=4;
X=5;
numCompartments=5;

mu= 0.02;  %turnover rate. Rate at which folks leave hospital. Completly uniform for no damn reason.

InTotal=1;
m=ones(numCompartments,1); %fractional rate of arrival;
m(RAB)=0;
m(X)=0;
m=m/sum(m);

beta=0.1; %Rate of infection
sigma= 0.02; %relative rate of replacement infection

cA=0.04;
cB=0.04;
cAB=0.17; %Cost of antibiotic resistance.


tau=0.4; %rate of recovery due to treatment
gam=0.02;  %Spontaneous recovery.

nuA= 10^-6; %Spotaneous rate of resistance mutation thing.
nuB= 10^-6;
nuAB= 10^-13;


P=ones(1,numCompartments);
chi= [0;0.5;0.5;0]; %Rate of application of antibiotics. This is the "mixing" paradigm, for example.
chi= [0;0.0;0.0;1]; %Rate of application of antibiotics. This is the "mixing" paradigm, for example.
chi= [1;0.0;0.0;0]; 

turnOvers=[logspace(0,1.9,30),logspace(1.9,2.1,60),logspace(2.1,2.5,30)]

Cstor=0*turnOvers;
Rstor=0*turnOvers;
Vstor=0*turnOvers; %Average "VALUE"

for(qqq=1:length(turnOvers))

TurnOver=turnOvers(qqq);
chiA= @(t) 1*(mod(t,2*TurnOver)< TurnOver);
chiB= @(t) 1*(mod(t,2*TurnOver)>= TurnOver);
chiAB= @(t) 0;

%%Now we apply "Double push" thing.
% chiA= @(t) 0;
% chiB= @(t) 0;
% chiAB= @(t) 1;


%%Now we apply "Do Nothing".
% chiA= @(t) 0;
% chiB= @(t) 0;
% chiAB= @(t) 0;


%%Now we apply "Mixing" thing.
% chiA= @(t) 0.5;
% chiB= @(t) 0.5;
% chiAB= @(t) 0;

recoveryMatrix = @(t) [-(gam+tau),0,0,0,0;
                  0,-(gam+tau*(chiB(t)+chiAB(t))),0,0,0;
                  0,0,-gam-tau*(chiA(t)+chiAB(t)),0,0;
                  0,0,0,-gam,0;
                  gam+tau,gam+tau*(chiB(t)+chiAB(t)),gam+tau*(chiA(t)+chiAB(t)),gam,0];

              
              
evolutionMatrix=  @(t) [-(chiA(t)*nuA+chiB(t)*(nuB)+nuAB),0,0,0,0;
                        chiA(t)*nuA,-(chiB(t)+chiAB(t) )*nuB,0,0,0;
                        chiB(t)*(nuB),0,-(chiA(t)+chiAB(t))*nuA,0,0;
                        nuAB,(chiB(t)+chiAB(t))*nuB,(chiA(t)+chiAB(t))*nuA,0,0;
                        0,0,0,0,0];
              
%Still need to do infection matrix and such.

infectionRates= diag([beta;beta*(1-cA);beta*(1-cB);beta*(1-cAB);0])
infectionRates(end,:)=-sum(infectionRates);



ode= @(t,P) InTotal*m*mu -mu*P + recoveryMatrix(t)*P + evolutionMatrix(t)*P + infectionRates*P*P(end)+...
        [-sigma*beta*(1-cA)*chiA(t)*P(RA)*P(S)-sigma*beta*(1-cA)*chiB(t)*P(RB)*P(S)-sigma*beta*(1-cAB)*P(RAB)*P(S);
        sigma*beta*(1-cA)*chiA(t)*P(RA)*P(S)+sigma*beta*(1-cA)*chiA(t)*P(RA)*P(RB)-sigma*beta*(1-cB)*chiB(t)*P(RB)*P(RA)-sigma*beta*(1-cAB)*(chiB(t)+chiAB(t))*P(RAB)*P(RA);
        sigma*beta*(1-cB)*chiB(t)*P(RB)*P(S)-sigma*beta*(1-cA)*chiA(t)*P(RA)*P(RB)+sigma*beta*(1-cB)*chiB(t)*P(RB)*P(RA)-sigma*beta*(1-cAB)*(chiA(t)+chiAB(t))*P(RAB)*P(RB);
        sigma*beta*(1-cAB)*P(RAB)*P(S)+sigma*beta*(1-cAB)*(chiB(t)+chiAB(t))*P(RAB)*P(RA)+sigma*beta*(1-cAB)*(chiA(t)+chiAB(t))*P(RAB)*P(RB)
        0]; %This line is co-infections.

integroODE= @(t,y)  [ode(t,y(1:X));y(X)];
    
tspan = [0 16000];
y0 = rand(numCompartments+1,1);
%y0(RAB)=0;
y0(end)=0;
y0=y0/sum(y0);


[tClunk,yClunk] = ode45(integroODE, (tspan(1):TurnOver:tspan(2)), y0);
    
Vstor(qqq)= (yClunk(end,end)-yClunk(end-4,end))/(tClunk(end)-tClunk(end-4));


Q=yClunk(:,5)-yClunk(end,5)

Mid= floor(length(tClunk)/2);

Mid=min(Mid,12);
Mid=max(Mid,4);


%Assume  yClunk(:,5) = A *r^k +C.s
%Find remaining values.
 rsqr=((yClunk(Mid+3,4)-yClunk(Mid-1,4))./(yClunk(Mid+1,4)-yClunk(Mid-3,4)));  %(A r^13 -A r^9)/(A r^11 -A r^7) = r^2
 r=sqrt(max(rsqr,0));
 if(r==0)
    A=0;
    C= mean(yClunk( (Mid-2):(Mid+2),4)); 
 else
     A= (yClunk(Mid+2,4)-yClunk(Mid-2,4))/(r^(Mid+2)-r^(Mid-2));
    C= yClunk(Mid,4)-A*r^Mid;
 end
    Cstor(qqq)=C;
    Rstor(qqq)=r;
 
end

plot(turnOvers,Cstor)
hold on
plot(turnOvers,Rstor)
plot(turnOvers,Vstor)

