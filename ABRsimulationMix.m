
function [deriv,Jacobian] = ABRsimulationMix(chi,P,lockedParams)

S=1;
RA=2;
RB=3;
RAB=4;
X=5;
numCompartments=5;

 chiA=chi(2);
 chiB=chi(3);
chiAB=chi(4);


if(~exist('lockedParams','var'))
    lockedParams=[0.02,1,0.1,0.02,0.04,0.04,0.17,0.4,0.02,10^-6,10^-6,10^-13];
end

mu= lockedParams(1);% 0.02;  %turnover rate. Rate at which folks leave hospital. Completly uniform for no damn reason.

InTotal= lockedParams(2);
m=[1/3;chiA/3;(1-chiA)/3;(1-chiB)/3;chiB/3;0;0]; %fractional rate of arrival;

beta= lockedParams(3); % 0.1 %Rate of infection
sigma=  lockedParams(4); %0.02; %relative rate of replacement infection

cA=  lockedParams(5);% 0.04;
cB=  lockedParams(6);% 0.04;
cAB=  lockedParams(7);%0.17; %Cost of antibiotic resistance.


tau=  lockedParams(8);%0.4; %rate of recovery due to treatment
gam=  lockedParams(9);% 0.02;  %Spontaneous recovery.

nuA=  lockedParams(10);%10^-6; %Spotaneous rate of resistance mutation thing.
nuB=  lockedParams(11);% 10^-6;
nuAB=  lockedParams(12);%10^-13;

% S, RA_A, RA_B,RB_A,RB_B,RAB,X

RecDiag=[(gam+tau),gam,(gam+tau),(gam+tau),gam,gam,0];

recoveryMatrix = -diag(RecDiag);
recoveryMatrix(end,:)=RecDiag;
             

infectionRates= diag([beta;beta*(1-cA)*chiA;beta*(1-cA)*(1-chiA);beta*(1-cB)*(1-chiB);beta*(1-cB)*chiB;beta*(1-cAB);0]);
infectionRates(end,:)=-sum(infectionRates);

deriv=  InTotal*m*mu -mu*P + recoveryMatrix*P + infectionRates*P*P(end);

Jacobian= -mu*eye(size(P,1)) + recoveryMatrix+infectionRates*P(end)+ [zeros(size(P,1),size(P,1)-1),infectionRates*P];
    
end