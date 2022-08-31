function errvect=OBJVisDispC3P(X,expdata,Xguess,riseTime,hmax,a,m)
    numParam = (length(Xguess)-1)/2;
    C0 =X(1)*Xguess(1);
    C = X(2:numParam + 1).*Xguess(2:numParam + 1);
    T = X(numParam + 2:length(Xguess)).*Xguess(numParam + 2:length(Xguess));
    fit=[];
    for k=1:length(expdata(:,1))
        fit(k)=C0*a*hmax^m+sum(C.*a*hmax^m.*T./riseTime.*(exp(riseTime./T)-1).*exp(-expdata(k,1)./T));
    end
    errvect=(fit'-expdata(:,2))/mean(expdata(:,2));

%    figure(2)
%    subplot(2,1,1);plot(expdata(:,1),expdata(:,2),'k.',expdata(:,1),fit,'b-');
%    subplot(2,1,2);bar(expdata(:,1),errvect);
 