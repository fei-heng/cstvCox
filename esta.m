function  [betahat]=esta(nsamp,ncovt,time,covar,covart,deltacs,wkernk,beta,wipw,rhohat,delta)
% estimate the time-varying coefficient in the cox model (aipw)
error = 1;
iter = 1;
while error>0.001 && iter<15
    U=zeros(ncovt, 1);
    F=zeros(ncovt, ncovt);
    
    % bz = covar * beta;
    bzt = zeros(nsamp, nsamp);
    for i=1:nsamp
        bzt(:,i) = covart(:,:,i) * beta;
    end
    
    wt = wkernk.*delta;
    s0 = zeros(nsamp,1);
    s1 = zeros(ncovt, nsamp);
    s2 = zeros(ncovt, ncovt, nsamp);
    
    for i=1:nsamp
        if wt(i)>0
            for l=1:nsamp
                if time(l) >= time(i)
                    s0(i) = s0(i) + exp(bzt(l,i));
                    for k=1:ncovt
                        s1(k,i) = s1(k,i) + exp(bzt(l,i))*covart(l,k,i);
                        for j=1:ncovt
                            s2(j,k,i) = s2(j,k,i) + exp(bzt(l,i))*covart(l,j,i)*covart(l,k,i);
                        end
                    end
                end
            end
            
            for k=1:ncovt
                U(k) = U(k)+wt(i)*(covar(i,k)-s1(k,i)/s0(i))*(wipw(i)*deltacs(i)+(1-wipw(i))*rhohat(i));
                for j=1:ncovt
                    F(j,k) = F(j,k)+wt(i)*(s2(j,k,i)/s0(i)-s1(j,i)*s1(k,i)/s0(i)^2)...
                        *(wipw(i)*deltacs(i)+(1-wipw(i))*rhohat(i));
                end
            end
        end
    end
    try
        change = pinv(F)*U;
        beta = beta + change;
        error = sum(abs(change));
        betahat = beta;
    catch
        betahat = nan(ncovt,1);
        break;
    end
    iter = iter + 1;
end

if iter >= 15
    betahat = nan(ncovt,1);
end

if abs(betahat(1,1))>5
    betahat = nan(ncovt,1);
end















