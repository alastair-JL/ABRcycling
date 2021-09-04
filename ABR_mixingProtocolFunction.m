function [InVect,OutVect] = ABR_mixingProtocolFunction(beta, m,p,q,tau,mu,gamma)
% 
% beta=[1,0.99,0.99];
% m= [10,1,0.1,0.1]; %import rates X,S,RA RB
% 
% p=0.5; %%Probability of drug A
% q=1/3; %Three day time to correct AB use. 
% tau=2; %half day recover time with antibiotitics 
% mu=1/10; % ten day recovery time without.
% 
% gamma=1/10; %Discharge rate, average of ten days.
% %

FullImport= [m(1);m(2); p*m(3);(1-p)*m(3); (1-p)*m(4);p*m(4)];

 InVect=FullImport/mu; %Default   %% X,S RAA RAB RBB RBA
    InVect=InVect.* (1.1.^randn(6,1)); %Blebu.

    CorrectiveMedicine= zeros(6,6);
    CorrectiveMedicine(3,3)=-q;
    CorrectiveMedicine(4,3)=q;
    CorrectiveMedicine(5,5)=-q;
    CorrectiveMedicine(6,5)=q;
    DerivDischarge= -mu*eye(6);
    DerivRecover= -diag([0;tau+gamma;gamma;tau+gamma;gamma;tau+gamma]);
    DerivRecover(1,:)=[0;tau+gamma;gamma;tau+gamma;gamma;tau+gamma]';
    
    InfectGlobulizer= [0,0,0,0,0,0;
                    0,beta(1),0,0,0,0;
                    0,0,beta(2)*p,beta(2)*p,0,0;
                    0,0,beta(2)*(1-p),beta(2)*(1-p),0,0;
                    0,0,0,0,beta(3)*(1-p),beta(3)*(1-p);
                    0,0,0,0,beta(3)*p,beta(3)*p;
                    ];
    InfectGlobulizer(1,:)= -sum(InfectGlobulizer);
            sigma= -sum(FullImport)./DerivDischarge(1,1);
             
             
     for(iii=1:3450)

FullInfect= (InfectGlobulizer*InVect)*(InVect(1));

OutVect=FullImport+DerivDischarge*InVect+DerivRecover*InVect+FullInfect + CorrectiveMedicine*InVect;

if(any(isnan(OutVect)) )
    error('what')
end

InVect=InVect+OutVect*0.001;
InVect=max(InVect,0*InVect);
InVect=InVect*sigma./sum(InVect);

     end
     
 for(iii=1:145)


FullInfect= (InfectGlobulizer*InVect)*(InVect(1));

OutVect=FullImport+DerivDischarge*InVect+DerivRecover*InVect+FullInfect + CorrectiveMedicine*InVect;

DerivInfect=InfectGlobulizer*(InVect(1));
DerivInfect(:,1)= InfectGlobulizer*InVect;

TotalDeriv=DerivDischarge+DerivRecover+DerivInfect+CorrectiveMedicine;

dir=TotalDeriv\OutVect;


 if(any(InVect<dir))

     alpha=1;
     while (any(InVect+OutVect*alpha<0) )
         alpha=alpha/2; 
     end     
     
     InVect=InVect+OutVect*alpha;
 else
     InVect=InVect-dir;
 end

 if(any(InVect<0))
     error('wut');
    InVect=FullImport/gamma; %Default   %% X,S RAA RAB RBB RBA
    InVect=InVect.* (1.1.^randn(6,1)); %Blebu.
 end


 end
 
InVect;
 
 
end