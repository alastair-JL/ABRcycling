% ModelIS   [ S, R_A,R_B,X ]

 y0=[rand(4,1);0]
% y0=[ swap*(yOut(end,1:6)') ;0];
% y0= [(yOut(end,1:6)');0];

Tlist= logspace(-0.5,2.5,1250);
Tlist=[0,Tlist];
Result=0*Tlist;
mu=1/15; %(immgration/emigration rate)
beta= [1,0.99,0.99,0.97]; %Infection rates;
smallImport=0.01;

importation= [0.5,smallImport,smallImport,1]';
importation=  importation*mu/([1,1,1,1]*importation);

g=1/10; %Ten day recovery without meds;
q=1/2.5; %treated recovery.


tau= (q-g);
MixEquil= (g + (tau*g)/(q+g))/beta(2)


Xequil= (g+mu)/beta(2);
Sequil= importation(1)/(q+mu-beta(1)*Xequil)
Bequil= importation(3)/(q+mu-beta(1)*Xequil)
Aequil= 1-Xequil-Sequil-Bequil
tSat= (log(Aequil)-log(Bequil)-1)*2/tau;

y0=[Sequil;Aequil;0;0;Bequil;Xequil;0];

extra=8;
ExitNumbers=-ones(length(Tlist),4+extra);
Maxs=-ones(length(Tlist),1);
Mins=-ones(length(Tlist),1);
Means=-ones(length(Tlist),1);
meanDownSwing=-ones(length(Tlist),1);
meanUpSwing=-ones(length(Tlist),1);


TA=50;
TB=50;

chi=@(t) 1*(mod(t,TA+TB)>TA);
ts=linspace(0,5*(TA+TB),1000);


importation= [0.5,smallImport,smallImport,1,0]';
importation=  importation*mu/([1,1,1,1,1]*importation);

b=0.5;
betaMatrix= [beta(1),0,0,0,0;
             0,beta(2),0,0,0;
             0,0,beta(3),0,0;
             -beta(1),-beta(2),-beta(3),0,-beta(4);
             0,0,0,0,beta(4)];
         
recoveryA= diag(-[q,g,q,0,g]);
recoveryA(4,:)=[q,g,q,0,g];
         
recoveryB= diag(-[q,q,g,0,g]);
recoveryB(4,:)=[q,q,g,0,g];
         
recovery= @(t)  recoveryB+chi(t)*(recoveryA-recoveryB);

Deriv =@(t,V) [1;1;1;1;1*t>240].*[(importation + betaMatrix*V(1:5)*V(4) + recovery(t)*V(1:5)-mu*V(1:5))];

 y0=[ones(4,1)/4;0];
 
warmupTime= 4*(TA+TB);

[tOut,yOut] = ode45(Deriv,[-5*(TA+TB),240],y0);
y0=yOut(end,:)';
y0(5)=0.003;
[tOut2,yOut2] = ode45(Deriv,[240,400],y0);

tOut=[tOut;tOut2];

yOut(:,end)=NaN
yOut=[yOut;yOut2]

Select=find(tOut>0);

figure(1)

area(tOut(Select)-tOut(Select(1)),yOut(Select,4),'facecolor',[0.9,0.6,0.9]);
hold on
plot(tOut(Select)-tOut(Select(1)),yOut(Select,[1,2,3]));
plot(tOut(Select)-tOut(Select(1)),yOut(Select,4),'lineWidth',2);
plot(tOut(Select)-tOut(Select(1)),yOut(Select,[5]));
xlim([0,400])
xlabel('t');
ylabel('population');

title('Slow Cycling');

ylabel('population');