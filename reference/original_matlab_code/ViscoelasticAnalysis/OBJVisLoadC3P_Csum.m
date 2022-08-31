function errvect=OBJVisLoadC3P_Csum(X,expdata,Xguess,riseTime,rate,a,m, weights)
    numParam = (length(Xguess)-1)/2;
    Csum =X(1)*Xguess(1); % enforce G0 to be possitive
    C = X(2:numParam + 1).*Xguess(2:numParam + 1);
    C0 = Csum + sum(C);
    T = X(numParam + 2:length(Xguess)).*Xguess(numParam + 2:length(Xguess));
    fit=[];
    for k=1:length(expdata(:,1))
        if expdata(k,1)<riseTime
            fit(k)=(1/a)^(1/m)*(C0*rate*expdata(k,1)-sum(C.*T.*rate.*(1-exp(-expdata(k,1)./T))))^(1/m);
        end
        if expdata(k,1)>=riseTime
            fit(k)=(1/a)^(1/m)*(C0*rate*riseTime-sum(C.*T.*rate.*exp(-expdata(k,1)./T).*(exp(riseTime./T)-1)))^(1/m);
        end
    end
  
  errvect=(fit'-expdata(:,2))/mean(expdata(:,2));
  errvect = errvect.*weights;
 
       
%      figure(1)
%       subplot(2,1,1);plot(expdata(:,1),expdata(:,2),'k.',expdata(:,1),fit(:),'b-');
%        subplot(2,1,2);bar(expdata(:,1),errvect);
