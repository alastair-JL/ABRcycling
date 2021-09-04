
if(~exist('mu','var')) %%Load a bunch of variables, but only if not already loaded. Take mu as an "indicator" variable.

mu=1/10; %(exit rate)
m= [2,0.1,0.1,20]; %entrance rate [S,R_A,R_B,X];
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
end

ChiEr=[-1,chi,chi,chi,chi,2];

t=0;

integralFinalsSmashHighest=zeros(6,1);
integralFinalsSmashLowestNonZero=zeros(6,1);

iii=1
ccc=1;

integralBase=0;
integralSelect=0;
integralHGT=0;
integralX=0;
integralX_compensated=0;

CycleTime=CycleTimes(ccc);
warmUpTime= 250;
CompletionTime= 250000;

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
         
     chiTab= ((CurrentState(5)+CurrentState(4)) <(CurrentState(2)+CurrentState(3)));
     
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
         
      chiTab= ((CurrentState(5)+CurrentState(4)) <(CurrentState(2)+CurrentState(3)));
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

integralFinalsSmashHighest(1,ccc)=t;
integralFinalsSmashHighest(2,ccc)=integralSelect;
integralFinalsSmashHighest(3,ccc)=integralBase;
integralFinalsSmashHighest(4,ccc)=integralHGT;
integralFinalsSmashHighest(5,ccc)=integralX;
integralFinalsSmashHighest(6,ccc)=integralX_compensated;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


iii=1
ccc=2;

integralBase=0;
integralSelect=0;
integralHGT=0;
integralX=0;
integralX_compensated=0;

CycleTime=CycleTimes(ccc);
warmUpTime= 250;
CompletionTime= 250000;

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
         
         if((CurrentState(5)+CurrentState(4))==0)
              chiTab=1;
         elseif((CurrentState(2)+CurrentState(3))==0)
              chiTab=0;
         else
              chiTab= ((CurrentState(5)+CurrentState(4)) >(CurrentState(2)+CurrentState(3)));
         end
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
         if((CurrentState(5)+CurrentState(4))==0)
              chiTab=1;
         elseif((CurrentState(2)+CurrentState(3))==0)
              chiTab=0;
         else
              chiTab= ((CurrentState(5)+CurrentState(4)) >(CurrentState(2)+CurrentState(3)));
         end
         
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

integralFinalsSmashLowestNonZero(1,1)=t;
integralFinalsSmashLowestNonZero(2,1)=integralSelect;
integralFinalsSmashLowestNonZero(3,1)=integralBase;
integralFinalsSmashLowestNonZero(4,1)=integralHGT;
integralFinalsSmashLowestNonZero(5,1)=integralX;
integralFinalsSmashLowestNonZero(6,1)=integralX_compensated;

%save('StochCycleData.mat');






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Rup= sum(m)./mu - (mu+gam)./betaVect(2)- m(2)./tau- m(1)./( (betaVect(2)-betaVect(1))*(gam+mu)/betaVect(1) +tau)
Rdown= m(2)./tau

tSat=2/tau*(log(Rup)-log(Rdown)-1)

figure()
thingZ=integralFinals([3,2,4],:)./integralFinals(1,:)./integralFinals([3,2,4],1);
plot(CycleTimes,thingZ,'lineWidth',2);
hold on
plot([tSat,tSat],[0,max(max(thingZ))],'k')
legend({"M_{base}","M_{select}","M_{HGT}"});
xlabel("T, cycle time");
ylabel("M, mutant arrival rate");

figure()

thingA=integralFinals(5,:)./integralFinals([3,2,4],:);
thingA=thingA./max(thingA,[],2);
plot(CycleTimes,thingA,'lineWidth',2);
hold on
plot([tSat,tSat],[0,1],'k')
legend({"M_{base}","M_{select}","M_{HGT}"});
xlabel("T, cycle time");
ylabel("\int X dt");



figure();
thingB=integralFinals(6,:)./integralFinals([3,2,4],:);
thingB=thingB./max(thingB,[],2);
plot(CycleTimes,thingB,'lineWidth',2);
hold on
plot([tSat,tSat],[0,1],'k')
legend({"M_{base}","M_{select}","M_{HGT}"});
xlabel("T, cycle time");
ylabel("\int X -X_M dt");

% 

figure(9)
subplot(1,3,1)
thingB=integralFinals(6,:)./integralFinals(3,:);
thingB_LNZ=integralFinalsSmashLowestNonZero(6)./integralFinalsSmashLowestNonZero(3);
thingB_High=integralFinalsSmashHighest(6)./integralFinalsSmashHighest(3);

plot([CycleTimes(1),CycleTimes(end)],[1,1]*thingB_High./max(thingB,[],2),'m:','LineWidth',2);
hold on
plot([CycleTimes(1),CycleTimes(end)],[1,1]*thingB_LNZ./max(thingB,[],2),'g--','LineWidth',2);
thingB=thingB./max(thingB,[],2);
plot(CycleTimes,thingB,'lineWidth',2);
hold on
plot([tSat,tSat],[0,1],'k')
title('M_{base}');
xlabel('T');
ylabel('\int X-X_0 dt');

subplot(1,3,2)
thingB=integralFinals(6,:)./integralFinals(2,:);
thingB_LNZ=integralFinalsSmashLowestNonZero(6)./integralFinalsSmashLowestNonZero(2);
thingB_High=integralFinalsSmashHighest(6)./integralFinalsSmashHighest(2);

plot([CycleTimes(1),CycleTimes(end)],[1,1]*thingB_High./max(thingB,[],2),'m:','LineWidth',2);
hold on
plot([CycleTimes(1),CycleTimes(end)],[1,1]*thingB_LNZ./max(thingB,[],2),'g--','LineWidth',2);
thingB=thingB./max(thingB,[],2);
plot(CycleTimes,thingB,'lineWidth',2);
hold on
plot([tSat,tSat],[0,1],'k')
title('M_{select}');
xlabel('T');
ylabel('\int X-X_0 dt');


subplot(1,3,3)
thingB=integralFinals(6,:)./integralFinals(4,:);
thingB_LNZ=integralFinalsSmashLowestNonZero(6)./integralFinalsSmashLowestNonZero(4);
thingB_High=integralFinalsSmashHighest(6)./integralFinalsSmashHighest(4);

plot([CycleTimes(1),CycleTimes(end)],[1,1]*thingB_High./max(thingB,[],2),'m:','LineWidth',2);
hold on
plot([CycleTimes(1),CycleTimes(end)],[1,1]*thingB_LNZ./max(thingB,[],2),'g--','LineWidth',2);
thingB=thingB./max(thingB,[],2);
plot(CycleTimes,thingB,'lineWidth',2);
hold on
plot([tSat,tSat],[0,1],'k')
title('M_{HGT}');
xlabel('T');
ylabel('\int X-X_0 dt');

