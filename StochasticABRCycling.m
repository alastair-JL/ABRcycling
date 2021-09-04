mu=1/10; %(exit rate)
m= [2,0.1,0.1,20]; %entrance rate [S,R_A,R_B,X];
gam= 1/15; %daseline recovery.
tau= 1/1.5; %Recovery with appropriate treatment.
recovVect=[gam+tau,gam,gam+tau,gam+tau,gam,0];

 %expected population,
R0=2.6; %this is my target R0 for an untreated infection in a fully succeptible population.

beta= [1,0.95,0.95]; %Infection rates [S,R_A,R_B];
betaVect=R0*(gam+mu)*[beta(1),beta(2),beta(2),beta(3),beta(3),0]/(sum(m)/mu); %rescalled infection rate.

chi=0.5; %Probablity of drug A.

CurrentState= [m(1),m(2)*chi,m(2)*(1-chi),m(3)*chi,m(3)*(1-chi),m(4)]/mu; % [S,R_A^A,R_A^B,R_B^A,R_B^B,X];
CurrentState=poissrnd(CurrentState);


%%Some Constants
IMPORT=0;
EXPORT=1;
INFECT=2;
RECOV=3;

ChiEr=[-1,chi,chi,chi,chi,2];

integralBase=0;
integralSelect=0;
integralHGT=0;

t=0;

cycleTime=20;

SelectRolladex=rand(8000,1);
ChiRolladex=rand(length(SelectRolladex),1);
TimeRolladex=-log(rand(length(SelectRolladex),1));
integrals=zeros(4,length(SelectRolladex));
stateRecord=zeros(6,length(SelectRolladex));

for(iii=(1:length(SelectRolladex)))
    actionList=[[m(1),m(2)*chi,m(2)*(1-chi),m(3)*chi,m(3)*(1-chi),m(4)],...
                CurrentState*mu,...
                betaVect.*CurrentState*CurrentState(6),...
                recovVect.*CurrentState];
     
     dt= TimeRolladex(iii)/sum(actionList); 
     t=t+dt;
     integralSelect=integralSelect+dt*(CurrentState*[0;0;1;1;0;0]);
     integralBase=integralBase+dt*(CurrentState*[0;1;1;1;1;0]);
     integralHGT=integralHGT+dt*(CurrentState*[0;1;1;0;0;0])*(CurrentState*[0;0;0;1;1;0]);
     integrals(1,iii)=t;
     integrals(2,iii)=integralBase;
     integrals(3,iii)=integralSelect;
     integrals(4,iii)=integralHGT;
     Select=sum(cumsum(actionList)<SelectRolladex(iii)*sum(actionList));
    
     sanctionedAction= floor(Select/6);
     Cut=mod(Select,6)+1;
     if(sanctionedAction==IMPORT)
         CurrentState(Cut)=CurrentState(Cut)+1;
     end
     if(sanctionedAction==EXPORT)
         CurrentState(Cut)=CurrentState(Cut)-1;
     end
     if(sanctionedAction==INFECT)
         
     chiTab= mod(t,2*cycleTime)<cycleTime;
     ChiEr=[1,chiTab,chiTab,chiTab,chiTab,0];
         Cut=Cut-mod(Cut,2)+ (ChiEr(Cut));
         CurrentState(Cut)=CurrentState(Cut)+1;
         CurrentState(end)=CurrentState(end)-1;
     end
     if(sanctionedAction==RECOV)
         CurrentState(Cut)=CurrentState(Cut)-1;
         CurrentState(end)=CurrentState(end)+1;
     end
     
     
     stateRecord(:,iii)=CurrentState';
end

integralBase=0;
integralSelect=0;
integralHGT=0;

t=0;

SelectRolladex=rand(80000,1);
ChiRolladex=rand(length(SelectRolladex),1);
TimeRolladex=-log(rand(length(SelectRolladex),1));
integrals=zeros(4,length(SelectRolladex));
stateRecord=zeros(6,length(SelectRolladex));

for(iii=(1:length(SelectRolladex)))
    actionList=[[m(1),m(2)*chi,m(2)*(1-chi),m(3)*chi,m(3)*(1-chi),m(4)],...
                CurrentState*mu,...
                betaVect.*CurrentState*CurrentState(6),...
                recovVect.*CurrentState];
     
     dt= TimeRolladex(iii)/sum(actionList); 
     t=t+dt;
     integralSelect=integralSelect+dt*(CurrentState*[0;0;1;1;0;0]);
     integralBase=integralBase+dt*(CurrentState*[0;1;1;1;1;0]);
     integralHGT=integralHGT+dt*(CurrentState*[0;1;1;0;0;0])*(CurrentState*[0;0;0;1;1;0]);
     integrals(1,iii)=t;
     integrals(2,iii)=integralBase;
     integrals(3,iii)=integralSelect;
     integrals(4,iii)=integralHGT;
     Select=sum(cumsum(actionList)<SelectRolladex(iii)*sum(actionList));
    
     sanctionedAction= floor(Select/6);
     Cut=mod(Select,6)+1;
     
     if(sanctionedAction==IMPORT)
         CurrentState(Cut)=CurrentState(Cut)+1;
     end
     if(sanctionedAction==EXPORT)
         CurrentState(Cut)=CurrentState(Cut)-1;
     end
     if(sanctionedAction==INFECT)
         
     chiTab= mod(t,2*cycleTime)<cycleTime;
     ChiEr=[1,chiTab,chiTab,chiTab,chiTab,0];
         Cut=Cut-mod(Cut,2)+ (ChiEr(Cut));
         CurrentState(Cut)=CurrentState(Cut)+1;
         CurrentState(end)=CurrentState(end)-1;
     end
     if(sanctionedAction==RECOV)
         CurrentState(Cut)=CurrentState(Cut)-1;
         CurrentState(end)=CurrentState(end)+1;
     end
     
     
     stateRecord(:,iii)=CurrentState';
end



plot(integrals(1,:),stateRecord);
legend({"S","R_A^A","R_A^B","R_B^A","R_B^B","X"})

