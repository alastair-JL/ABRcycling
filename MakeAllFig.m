beta=[1,0.99,0.99];
m= [10,1,0.1,0.1]; %import rates X,S,RA RB

p=1; %%Probability of drug A
q=1/3; %Three day time to correct AB use. 
tau=2; %half day recover time with antibiotitics 
mu=1/10; % ten day recovery time without.

gamma=1/10; %Discharge rate, average of ten days.
% %

FullImport= [m(1);m(2); p*m(3);(1-p)*m(3); (1-p)*m(4);p*m(4)];

 InVect=FullImport/gamma; %Default   %% X,S RAA RAB RBB RBA
    InVect=InVect.* (1.1.^randn(6,1)); %Blebu.

     for(iii=1:45)

FullImport= [m(1);m(2); p*m(3);(1-p)*m(3); (1-p)*m(4);p*m(4)];
FullDischarge= -gamma*InVect;
FullRecover= -InVect.*[0;tau+mu;mu;tau+mu;mu;tau+mu];
FullRecover(1)=-sum(FullRecover);

FullInfect= InVect.*[0;beta(1);beta(2)*p;(1-p)*beta(2);beta(3)*p;(1-p)*beta(3)]*(InVect(1));
FullInfect(1)=-sum(FullInfect);

OutVect=FullImport+FullDischarge+FullRecover+FullInfect;

InVect=InVect+OutVect*0.001
     end
     
 for(iii=1:45)

FullImport= [m(1);m(2); p*m(3);(1-p)*m(3); (1-p)*m(4);p*m(4)];
FullDischarge= -gamma*InVect;
FullRecover= -InVect.*[0;tau+mu;mu;tau+mu;mu;tau+mu];
FullRecover(1)=-sum(FullRecover);

FullInfect= InVect.*[0;beta(1);beta(2)*p;(1-p)*beta(2);beta(3)*p;(1-p)*beta(3)]*(InVect(1));
FullInfect(1)=-sum(FullInfect);

OutVect=FullImport+FullDischarge+FullRecover+FullInfect;

DerivDischarge= -gamma*eye(6);
DerivRecover= -diag([0;tau+mu;mu;tau+mu;mu;tau+mu]);
DerivRecover(1,:)=[0;tau+mu;mu;tau+mu;mu;tau+mu]';
DerivInfect=diag([0;beta(1);beta(2)*p;(1-p)*beta(2);beta(3)*p;(1-p)*beta(3)]*(InVect(1)));
DerivInfect(1,:)= -[0;beta(1);beta(2)*p;(1-p)*beta(2);beta(3)*p;(1-p)*beta(3)]'*(InVect(1));
DerivInfect(:,1)= InVect.*[0;beta(1);beta(2)*p;(1-p)*beta(2);beta(3)*p;(1-p)*beta(3)];
DerivInfect(1,1)= -InVect'*[0;beta(1);beta(2)*p;(1-p)*beta(2);beta(3)*p;(1-p)*beta(3)];

TotalDeriv=DerivDischarge+DerivRecover+DerivInfect;




dir=TotalDeriv\OutVect;


oldIn=InVect;
oldOut=OutVect;

 if(any(InVect<dir))

     alpha=1;
     while (any(InVect+OutVect*alpha<0) )
         alpha=alpha/2; 
     end     
     
     InVect=InVect+OutVect*alpha
 else
     InVect=InVect-dir
 end

end

OutVect
InVect

plot()