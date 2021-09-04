
function [deriv,Jacobian] = ABRsimulationCycle_withSignFlipper(chi,P,lockedParams,t,cycleTime)

S=1;
RA=2;
RB=3;
RAB=4;
X=5;
numCompartments=5;

P=P(1:5);


if(mod(t,2*cycleTime)<cycleTime)
 chiA=chi(2);
 chiB=chi(3);
else
 chiA=chi(3);
 chiB=chi(2);    
end

chiAB=chi(4);


if(~exist('lockedParams','var'))
    lockedParams=[0.02,1,0.1,0.02,0.04,0.04,0.17,0.4,0.02,10^-6,10^-6,10^-13];
end

mu= lockedParams(1);% 0.02;  %turnover rate. Rate at which folks leave hospital. Completly uniform for no damn reason.

InTotal= lockedParams(2);
m=ones(numCompartments,1); %fractional rate of arrival;
m(RAB)=0;
m(X)=0;
m=m/sum(m);

beta= lockedParams(3); % 0.1 %Rate of infection
sigma=  lockedParams(4);0.02; %relative rate of replacement infection

cA=  lockedParams(5);% 0.04;
cB=  lockedParams(6);% 0.04;
cAB=  lockedParams(7);%0.17; %Cost of antibiotic resistance.


tau=  lockedParams(8);%0.4; %rate of recovery due to treatment
gam=  lockedParams(9);% 0.02;  %Spontaneous recovery.

nuA=  lockedParams(10);%10^-6; %Spotaneous rate of resistance mutation thing.
nuB=  lockedParams(11);% 10^-6;
nuAB=  lockedParams(12);%10^-13;


recoveryMatrix = [-(gam+tau),0,0,0,0;
                  0,-(gam+tau*(chiB+chiAB)),0,0,0;
                  0,0,-gam-tau*(chiA+chiAB),0,0;
                  0,0,0,-gam,0;
                  gam+tau,gam+tau*(chiB+chiAB),gam+tau*(chiA+chiAB),gam,0];

              
              
evolutionMatrix= [-(chiA*nuA+chiB*(nuB)+nuAB),0,0,0,0;
                        chiA*nuA,-(chiB+chiAB )*nuB,0,0,0;
                        chiB*(nuB),0,-(chiA+chiAB)*nuA,0,0;
                        nuAB,(chiB+chiAB)*nuB,(chiA+chiAB)*nuA,0,0;
                        0,0,0,0,0];
              
%Still need to do infection matrix and such.

infectionRates= diag([beta;beta*(1-cA);beta*(1-cB);beta*(1-cAB);0]);
infectionRates(end,:)=-sum(infectionRates);

deriv=  InTotal*m*mu -mu*P + recoveryMatrix*P + evolutionMatrix*P + infectionRates*P*P(end)+...
        [-sigma*beta*(1-cA)*chiA*P(RA)*P(S)-sigma*beta*(1-cB)*chiB*P(RB)*P(S)-sigma*beta*(1-cAB)*P(RAB)*P(S);
        sigma*beta*(1-cA)*chiA*P(RA)*P(S)+sigma*beta*(1-cA)*chiA*P(RA)*P(RB)-sigma*beta*(1-cB)*chiB*P(RB)*P(RA)-sigma*beta*(1-cAB)*(chiB+chiAB)*P(RAB)*P(RA);
        sigma*beta*(1-cB)*chiB*P(RB)*P(S)-sigma*beta*(1-cA)*chiA*P(RA)*P(RB)+sigma*beta*(1-cB)*chiB*P(RB)*P(RA)-sigma*beta*(1-cAB)*(chiA+chiAB)*P(RAB)*P(RB);
        sigma*beta*(1-cAB)*P(RAB)*P(S)+sigma*beta*(1-cAB)*(chiB+chiAB)*P(RAB)*P(RA)+sigma*beta*(1-cAB)*(chiA+chiAB)*P(RAB)*P(RB)
        0]; %This line is co-infections.

    
Jacobian= -mu*eye(numCompartments) + recoveryMatrix+ evolutionMatrix+infectionRates*P(end)+ [zeros(numCompartments,numCompartments-1),infectionRates*P]...
   + [-sigma*beta*(1-cA)*chiA*P(RA)-sigma*beta*(1-cA)*chiB*P(RB)-sigma*beta*(1-cAB)*P(RAB), -sigma*beta*(1-cA)*chiA*P(S), -sigma*beta*(1-cA)*chiB*P(S),-sigma*beta*(1-cAB)*P(S),0;
        sigma*beta*(1-cA)*chiA*P(RA), sigma*beta*(1-cA)*chiA*P(S)+sigma*beta*(1-cA)*chiA*P(RB)-sigma*beta*(1-cB)*chiB*P(RB)-sigma*beta*(1-cAB)*(chiB+chiAB)*P(RAB), +sigma*beta*(1-cA)*chiA*P(RA)-sigma*beta*(1-cB)*chiB*P(RA), -sigma*beta*(1-cAB)*(chiB+chiAB)*P(RA),0;
        sigma*beta*(1-cB)*chiB*P(RB),-sigma*beta*(1-cA)*chiA*P(RB)+sigma*beta*(1-cB)*chiB*P(RB), sigma*beta*(1-cB)*chiB*P(S)-sigma*beta*(1-cA)*chiA*P(RA)+sigma*beta*(1-cB)*chiB*P(RA)-sigma*beta*(1-cAB)*(chiA+chiAB)*P(RAB), -sigma*beta*(1-cAB)*(chiA+chiAB)*P(RB),0;
        sigma*beta*(1-cAB)*P(RAB), sigma*beta*(1-cAB)*(chiB+chiAB)*P(RAB), sigma*beta*(1-cAB)*(chiA+chiAB)*P(RAB),sigma*beta*(1-cAB)*P(S)+sigma*beta*(1-cAB)*(chiB+chiAB)*P(RA)+sigma*beta*(1-cAB)*(chiA+chiAB)*P(RB),0;
        0,0,0,0,0];
    
    
   deriv=[deriv;sign((mod(t,2*cycleTime)-cycleTime))];
    
end