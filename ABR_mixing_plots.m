beta=[1,0.95,0.95]*0.01;
m= [5,3,0.1,0.1]; %import rates X,S,RA RB

p=0.5; %%Probability of drug A
q=1/3; %Three day time to correct AB use. 
tau=2; %half day recover time with antibiotitics 
mu=1/10; % ten day recovery time without.

gamma=1/20; %Discharge rate, average of ten days.
%
qVals=logspace(-3,0.5,100);
InsQ=zeros(6,length(qVals));
OutsQ=zeros(6,length(qVals));

for(qqq=1:length(qVals))
    q=qVals(qqq);
    [InsQ(:,qqq),OutsQ(:,qqq)] = ABR_mixingProtocolFunction(beta, m,p,q,tau,mu,gamma);   
    
end

    rhoB=((mu+gamma)*p + qVals)./(tau.*(1-p) + mu + gamma + qVals);
    rhoA=((mu+gamma)*(1-p) + qVals)./(tau.*p + mu + gamma + qVals);
    
    XAnalytic=min( [(mu+gamma + tau *rhoB)./(beta(3));(mu+gamma + tau *rhoA)./(beta(2))]);
    
figure(1)
subplot(1,2,1)
plot(qVals,InsQ(1,:),'LineWidth',2)
hold on
%plot(qVals,XAnalytic,'k:','LineWidth',2)

%set(gca,'XScale','log');
xlabel({'q';'Drug correction rate'})
ylabel('X(\infty)');
hold on
sigma= sum(m)/mu;

plot(qVals,XAnalytic,'k--','LineWidth',1)
plot(qVals, sigma-(mu+gamma+tau-sigma*beta(1)-sqrt((mu+gamma + tau -sigma*beta(1))^2+4*beta(1)*m(2)))/(-2*beta(1))+0*qVals ,'k:')
plot(qVals, sigma-sum(m(2:4))/(tau+gamma) +0*qVals ,'k--')
ylim([1.9*(gamma+mu)/beta(2),sigma])
ylim([0,sigma])
xlim([0,max(qVals)])


pVals=linspace(0.01,0.99,100);
InsP=zeros(6,length(pVals));
OutsP=zeros(6,length(pVals));
q=1/3;
for(qqq=1:length(pVals))
    p=pVals(qqq);
    [InsP(:,qqq),OutsP(:,qqq)] = ABR_mixingProtocolFunction(beta, m,p,q,tau,mu,gamma);
    
end

figure(1)
subplot(1,2,2)
%plot(pVals,InsP(1,:),'LineWidth',2)


beta=[1,0.99,0.8]*0.01;

pVals=linspace(0.01,0.99);
InsP=zeros(6,length(pVals));
OutsP=zeros(6,length(pVals));
q=1/3;
for(qqq=1:length(pVals))
    p=pVals(qqq);
    [InsP(:,qqq),OutsP(:,qqq)] = ABR_mixingProtocolFunction(beta, m,p,q,tau,mu,gamma);
    
end
figure(1)
hold on
plot(pVals,InsP(1,:),'LineWidth',2)
xlabel({'p';'Drug A probability'})
ylabel('X(\infty)');


  rhoB=((mu+gamma)*pVals + q)./(tau.*(1-pVals) + mu + gamma + q);
    rhoA=((mu+gamma)*(1-pVals) + q)./(tau.*pVals + mu + gamma + q);
    
    XAnalytic=min( [(mu+gamma + tau *rhoB)./(beta(3));(mu+gamma + tau *rhoA)./(beta(2))]);
    
plot(pVals,(mu+gamma + tau *rhoB)./(beta(3)),'k--')
plot(pVals,(mu+gamma + tau *rhoA)./(beta(2)),'k--')
ylim([0,sigma])
%plot(pVals,(mu+gamma + tau *rhoB)./(0.95*0.01),'k:')
%plot(pVals,(mu+gamma + tau *rhoA)./(0.95*0.01),'k:')
