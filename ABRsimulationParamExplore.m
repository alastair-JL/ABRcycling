
function [deriv] = ABRsimulationParamExplore(chi,P,lockedParams,t,cycleTime,m)

S=1;
RA=2;
RB=3;
RAB=4;
X=5;
numCompartments=5;

if(length(cycleTime)==1)
if(mod(t,cycleTime*2)<cycleTime)
 chiA=chi(2);
 chiB=chi(3);
else
 chiA=chi(3);
 chiB=chi(2);    
end
elseif(length(cycleTime)==2)
    if(mod(t,sum(cycleTime))<cycleTime(1))
        chiA=chi(2);
        chiB=chi(3);
    else
        chiA=chi(3);
        chiB=chi(2);    
    end
else
   error('Cycle Time vector does not have one or two entries');
end
chiAB=chi(4);


if(~exist('lockedParams','var'))
    lockedParams=[0.02,1,0.1,0.02,0.04,0.04,0.17,0.4,0.02,10^-6,10^-6,10^-13];
end

mu= lockedParams(1);% 0.02;  %turnover rate. Rate at which folks leave hospital. Completly uniform for no damn reason.

InTotal= lockedParams(2);

if(~exist('m','var'))
m=ones(numCompartments,1); %fractional rate of arrival;
m(RAB)=0;
m(X)=0;
m=m/sum(m);
end

beta= lockedParams(3); % 0.1 %Rate of infection

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
P=P(1:5);
P=P.*(P>0);

deriv=  m -mu*P + recoveryMatrix*P + evolutionMatrix*P + infectionRates*P*P(end);

P(2:3)=1./P(2:3);
P(1)=(gam+tau*(chiB+chiAB))+mu;
deriv=[deriv; (t>=0).*P];

end