% ModelIS   [ S, R_A^A,R_A^B,R_B^A,R_B^B,X ]

 y0=[rand(6,1);0]
% y0=[ swap*(yOut(end,1:6)') ;0];
% y0= [(yOut(end,1:6)');0];

T=250;
mu=1/15; %(immgration/emigration rate)
beta= [1,0.99,0.99]; %Infection rates;

chi=@(t) ceil( mod(t/T,2)+10^-8)-1;

ts=linspace(0,10*T,1000);

plot(ts,chi(ts));

smallImport=0.00001;

filter= @(t) [1;chi(t);1-chi(t);chi(t);1-chi(t);1];
importation= [0.5,smallImport,smallImport,smallImport,smallImport,1]';
importation=  importation*mu/(filter(1)'*importation);

b=0.5;
betaMatrix= [b,0,0,0,0,0;
             0,b,b,0,0,0;
             0,b,b,0,0,0;
             0,0,0,b,b,0;
             0,0,0,b,b,0;
             -b,-b,-b,-b,-b,0];
         
g=1/10; %Ten day recovery without meds;
q=1/4; %treated recovery.
recovery= diag(-[q,g,q,q,g,0]);
recovery(6,:)=[q,g,q,q,g,0];
         
Deriv =@(t,V) [(filter(t).*importation + filter(t).*betaMatrix*V(1:6)*V(6) + recovery*V(1:6)-mu*V(1:6)); V(6)];

[tOut,yOut] = ode45(Deriv,[0,20*T],y0);

figure(1)
plot(tOut,yOut);
hold on
plot(tOut,chi(tOut));

figure(2)
smoosh= [1,0,0,0,0,0,0;
         0,1,1,0,0,0,0;
         0,0,0,1,1,0,0;
         0,0,0,0,0,1,0];
     
plot(tOut,yOut*smoosh');
hold on

meanX=yOut(end,end)./tOut(end)

% dumbApproxA= importation(2)/(q-g) *exp(-q*tOut);
% dumbApproxB= (q/g)*importation(2)/(q-g) *(1-exp(-g*tOut));
% 
% plot(tOut,dumbApproxA)
% plot(tOut,dumbApproxB)
