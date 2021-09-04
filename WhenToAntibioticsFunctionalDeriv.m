


alph=1;
bigbet=0.35;
bet = @(t) 0.6+0*t;
betPert = @(t) 0.5*exp(-(t-12)^2);
epsilon=10^-4;

deriv = @(x,t) [-alph*x(1)*x(2);alph*x(1)*x(2)-bet(t)*x(2);bet(t)*x(2)];
derivB = @(x,t) [-alph*x(1)*x(2);alph*x(1)*x(2)-(epsilon*betPert(t)+bet(t))*x(2);(epsilon*betPert(t)+bet(t))*x(2)];

xVals= zeros(3,46000);
xValsPert=xVals;
xVals(:,1)=[1-10^-4;10^-4;0];
xValsPert(:,1)=[1-10^-4;10^-4;0];
dt=10^-3;


tim=0;
timSpace= (1:length(xVals))*dt;

marker= 0*timSpace;

for(qqq=2:length(xVals))
    tim=tim+dt;
  
    xVals(:,qqq)=xVals(:,qqq-1)+deriv(xVals(:,qqq-1),tim)*dt;
    xValsPert(:,qqq)=xValsPert(:,qqq-1)+derivB(xValsPert(:,qqq-1),tim)*dt;
    
end

figure(1)
subplot(1,2,1)
axis([0,max(timSpace),0,1])
hold on
plot(timSpace,xVals(:,:),'LineWidth',2)
xlabel('time');

subplot(1,2,2)
axis([0,max(timSpace),0,1])
hold on
plot(timSpace,xValsPert(:,:),'LineWidth',2)
xlabel('time');



figure(2)
axis([0,max(timSpace),-1*epsilon,1*epsilon])
hold on
plot(timSpace,xValsPert(:,:)-xVals(:,:))
xlabel('time');


integrals= cumsum(xVals,2)*dt;
integralBet= cumsum(bet(timSpace).*xVals(2,:))*dt;

bigIntegral= 0*integrals;
bigIntegral=bigIntegral(1:2,:);
predictedPert= bigIntegral;

predictedPertB= bigIntegral;


for(qqq=2:length(bigIntegral))
    MatrixExpNow=expm([-alph*(integrals(2,qqq)-integrals(1,qqq)), +alph*integrals(1,qqq) ;-integralBet(qqq),-integralBet(qqq)]);
    bigIntegral(:,qqq)=MatrixExpNow\[0;betPert(timSpace(qqq))];
    predictedPert([1,3],qqq)=MatrixExpNow*sum(bigIntegral(:,1:qqq),2)*dt;
    predictedPert(2,qqq)= -predictedPert(1,qqq)-predictedPert(3,qqq);
    nowM=[-alph*(xVals(2,qqq)-xVals(1,qqq)), alph*xVals(1,qqq) ;-bet(timSpace(qqq)),-bet(timSpace(qqq))];
    predictedPertB(:,qqq) = [0;betPert(timSpace(qqq-1))*xVals(2,qqq-1)*dt] + expm(nowM*dt)*predictedPertB(:,qqq-1);
    
end


predictedPertB(3,:)= -predictedPertB(1,:)-predictedPertB(2,:)
plot(timSpace,predictedPert*epsilon,'--','LineWidth',2)
predictedPertB=flipud(predictedPertB)

plot(timSpace,predictedPertB*epsilon,':','LineWidth',3)