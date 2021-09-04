mu=1/10; %(exit rate)
m= [5,0.03,0.03,10]; %entrance rate [S,R_A,R_B,X];
gam= 1/15; %daseline recovery.
tau= 1/1.5; %Recovery with appropriate treatment.
recovVect=[gam+tau,gam,gam+tau,gam+tau,gam,0];

 %expected population,
R0=2.6; %this is my target R0 for an untreated infection in a fully succeptible population.

beta= [1,0.95,0.95,0.9]; %Infection rates [S,R_A,R_B];
betaVect=R0*(gam+mu)*[beta(1),beta(2),beta(2),beta(3),beta(3),0]/(sum(m)/mu); %rescalled infection rate.

betaVectAB=R0*(gam+mu)*[beta(4)]/(sum(m)/mu);
MultiresistanceBase= (gam+mu)/betaVectAB; %The hypothetical "How bad would things be with multiresistance.

chi=0.5; %Probablity of drug A.

CurrentState= [m(1),m(2)*chi,m(2)*(1-chi),m(3)*chi,m(3)*(1-chi),m(4)]/mu; % [S,R_A^A,R_A^B,R_B^A,R_B^B,X];
CurrentState=poissrnd(CurrentState);


%%Some Constants
IMPORT=0;
EXPORT=1;
INFECT=2;
RECOV=3;

ChiEr=[-1,chi,chi,chi,chi,2];

t=0;


CycleTimes=[0.0001,0.1,linspace(0.1,50,60),linspace(60,300,30)];

integralFinals=zeros(4,length(CycleTimes))



for(ccc=1:length(CycleTimes))
iii=1;
ccc

integralBase=0;
integralSelect=0;
integralHGT=0;
integralX=0;
integralX_compensated=0;

CycleTime=CycleTimes(ccc);
warmUpTime= ceil(250/(2*CycleTime))*(2*CycleTime);
CompletionTime= ceil(150000/(2*CycleTime))*(2*CycleTime);

SelectRolladex=rand(8000,1);
TimeRolladex=-log(rand(length(SelectRolladex),1));

while(t<warmUpTime)
    actionList=[[m(1),m(2)*chi,m(2)*(1-chi),m(3)*chi,m(3)*(1-chi),m(4)],...
                CurrentState*mu,...
                betaVect.*CurrentState*CurrentState(6),...
                recovVect.*CurrentState];
     
     dt= TimeRolladex(iii)/sum(actionList); 
     t=t+dt;
   %  integralSelect=integralSelect+dt*(CurrentState*[0;0;1;1;0;0]);
   %  integralBase=integralBase+dt*(CurrentState*[0;1;1;1;1;0]);
   %  integralHGT=integralHGT+dt*(CurrentState*[0;1;1;0;0;0])*(CurrentState*[0;0;0;1;1;0]);

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
         
     chiTab= mod(t,2*CycleTime)<CycleTime;
     ChiEr=[1,chiTab,chiTab,chiTab,chiTab,0];
         Cut=Cut-mod(Cut,2)+ (ChiEr(Cut));
         CurrentState(Cut)=CurrentState(Cut)+1;
         CurrentState(end)=CurrentState(end)-1;
     end
     if(sanctionedAction==RECOV)
         CurrentState(Cut)=CurrentState(Cut)-1;
         CurrentState(end)=CurrentState(end)+1;
     end
     
     
    if(iii==length(SelectRolladex))
        [ccc,t]
        SelectRolladex=rand(28000,1);
        TimeRolladex=-log(rand(length(SelectRolladex),1));
        iii=1;
    else
        iii=iii+1;
    end
end

t=0;
while(t<CompletionTime)
    actionList=[[m(1),m(2)*chi,m(2)*(1-chi),m(3)*chi,m(3)*(1-chi),m(4)],...
                CurrentState*mu,...
                betaVect.*CurrentState*CurrentState(6),...
                recovVect.*CurrentState];
     
     dt= TimeRolladex(iii)/sum(actionList); 
     t=t+dt;
     integralSelect=integralSelect+dt*(CurrentState*[0;0;1;1;0;0]);
     integralBase=integralBase+dt*(CurrentState*[0;1;1;1;1;0]);
     integralHGT=integralHGT+dt*(CurrentState*[0;1;1;0;0;0])*(CurrentState*[0;0;0;1;1;0]);
        
     integralX=integralX+dt*(CurrentState*[0;0;0;0;0;1]);
     integralX_compensated=integralX_compensated+dt*(CurrentState*[0;0;0;0;0;1]-MultiresistanceBase);
        
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
         
     chiTab= mod(t,2*CycleTime)<CycleTime;
     ChiEr=[1,chiTab,chiTab,chiTab,chiTab,0];
         Cut=Cut-mod(Cut,2)+ (ChiEr(Cut));
         CurrentState(Cut)=CurrentState(Cut)+1;
         CurrentState(end)=CurrentState(end)-1;
     end
     if(sanctionedAction==RECOV)
         CurrentState(Cut)=CurrentState(Cut)-1;
         CurrentState(end)=CurrentState(end)+1;
     end
     
     
    if(iii==length(SelectRolladex))
        [ccc,t]
        SelectRolladex=rand(28000,1);
        TimeRolladex=-log(rand(length(SelectRolladex),1));
        iii=1;
    else
        iii=iii+1;
    end
end

integralFinals(1,ccc)=t;
integralFinals(2,ccc)=integralSelect;
integralFinals(3,ccc)=integralBase;
integralFinals(4,ccc)=integralHGT;
integralFinals(5,ccc)=integralX;
integralFinals(6,ccc)=integralX_compensated;

save('StochCycleData3.mat');
end

Rup= sum(m)./mu - (mu+gam)./betaVect(2)- m(2)./tau- m(1)./( (betaVect(2)-betaVect(1))*(gam+mu)/betaVect(1) +tau)
Rdown= m(2)./tau

tSat=2/tau*(log(Rup)-log(Rdown)-1)

thingZ=integralFinals([3,2,4],:)./integralFinals(1,:)./integralFinals([3,2,4],1);
plot(CycleTimes,thingZ,'lineWidth',2);
hold on
plot([tSat,tSat],[0,max(max(thingZ))],'k')
legend({"M_{base}","M_{select}","M_{HGT}"});

figure()

thingA=integralFinals(5,:)./integralFinals([3,2,4],:);
thingA=thingA./max(thingA,[],2);
plot(CycleTimes,thingA,'lineWidth',2);
hold on
plot([tSat,tSat],[0,1],'k')
legend({"M_{base}","M_{select}","M_{HGT}"});

figure();
thingB=integralFinals(6,:)./integralFinals([3,2,4],:);
thingB=thingB./max(thingB,[],2);
plot(CycleTimes,thingB,'lineWidth',2);
hold on
plot([tSat,tSat],[0,1],'k')
legend({"M_{base}","M_{select}","M_{HGT}"});