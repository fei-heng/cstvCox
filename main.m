
data = csvread('simdata.csv');

%% INPUT
time=data(:,1);
delta=data(:,2);
z1=data(:,3:4);
R=data(:,5);
A=data(:,6);
cause=data(:,7);
missingmodel=[z1,A]; % model for the missing probability

nsamp = length(time);
ncs=length(unique(cause(~isnan(cause)&cause~=0)));
ncov=size(z1,2); % number of covariates
ncovt=ncov*2;
npr=1+size(missingmodel,2); % number of parameters for predicting missing probability r(W, \psi)
% nqr=1; % number of parameters for h(a|t,k,z)

%% Parameter settings
h=0.3; % bandwidth
b=h; % bandwidth for lambda_0(t)

% set grid points
tau=2;
t1=h;
t2=2-h;
tstep=0.025; % grid point width
tgrid = tstep:tstep:tau;
ngrid=fix(tau/tstep); % number of grids



%% covariate matrix
z2 = zeros(nsamp,ncov);
covar = [z1, z2];
covar0 = [ones(nsamp, 1), missingmodel]; % covar3
covar2 = z1;

%% initialization
sbeta_f = zeros(ncov,ngrid,ncs);
sstd_f = zeros(ncov,ngrid,ncs);

sbeta_c = zeros(ncov,ngrid,ncs);
sstd_c = zeros(ncov,ngrid,ncs);

sbeta_ic = zeros(ncov,ngrid,ncs);
sstd_ic = zeros(ncov,ngrid,ncs);

sbeta_acc = zeros(ncov,ngrid,ncs);
sstd_acc = zeros(ncov,ngrid,ncs);
sinv_acc = zeros(ncov,ncov,ngrid,ncs);

%% complete case only:
nsampc = sum(R);
indexc = find(R == 1);
timec = time(indexc);
covarc = z1(indexc,:);
covar2c = z1(indexc,:);
causec = cause(indexc);

%% ipw and aipw only:
% correctly specified model of r(W)
% ul = score function
% fl = negative second derivative
% tz = theta*Z

[psi,dev,stats]=glmfit(missingmodel(delta==1,:),R(delta==1),'binomial');

% psi = zeros(npr, 1);
% error = 1;
% while error > 0.001
%     ul = zeros(npr, 1);
%     fl = zeros(npr, npr);
%     tz = covar0 * psi;
%     p = exp(tz)./(1+exp(tz));
%     % ul = covar0' * (R.*delta.*(1-p)-(1-R).*p);
%     for i=1:nsamp
%         for k=1:npr
%             ul(k)=ul(k)+(R(i)*delta(i)*(1-p(i))-(1-R(i))*p(i))*covar0(i,k);
%             for j=1:npr
%                 fl(j,k)=fl(j,k)+(R(i)*delta(i)+(1-R(i)))*p(i)*(1-p(i))*covar0(i,k)*covar0(i,j);
%             end
%         end
%     end
%     
%     %         if (rcond(fl)<0.00000001 )
%     %             break;
%     %         end
%     psi = psi+pinv(fl)*ul;
%     error = sum(abs(pinv(fl)*ul));
% end
ph = 1./(1+exp(-covar0*psi)); % r(W, \hat{\psi})
wipw = R./(delta.*ph+1-delta);% R / \pi(Q, \hat{\psi})

%% aipw only:
% only use complete data to estimate \hat{\rho}
hhat = zeros(nsamp,2);
%% \hat{\rho} for cause 1
% logistic fit model (correctly specified)
nsampc1 = sum(R.*(cause==1));
indexc1 = find(R.*(cause==1));
A1 = A(indexc1);
theta1 = log(sum(A1)/(nsampc1-sum(A1)));
hhat(:,1) = (1/(1+exp(-theta1)))*A + (1/(1+exp(theta1)))*(1-A);

%% \hat{\rho} for cause 2
% logistic fit model (correctly specified)
nsampc2 = sum(R.*(cause==2));
indexc2 = find(R.*(cause==2));
A2 = A(indexc2);
theta2 = log(sum(A2)/(nsampc2-sum(A2)));
hhat(:,2) = (1/(1+exp(-theta2)))*A + (1/(1+exp(theta2)))*(1-A);


%% estimation at grid points for each cause
covart = zeros(nsamp, ncovt, nsamp);
covartc = zeros(nsampc, ncovt, nsampc);
wkernk = zeros(nsamp, 1);
wkernkc = zeros(nsampc, 1);

%% for full, complete case and ipw
ipwbeta_c = zeros(ncov,ngrid,ncs);
ipwbeta_c_all = zeros(ncov,nsamp,ncs);
beta_acc_all = zeros(ncov,nsamp,ncs);

ipwbeta_m = zeros(ncov,ngrid,ncs);
ipwbeta_m_all = zeros(ncov,nsamp,ncs);

for ics=1:ncs
    deltacs = cause == ics; % deltacs: indicator of non-censoring due to each cause
    deltacsc = causec == ics;
    tvalue = tstep;
    for igrid=1:ngrid
        % calculate covar and covart: Z and Z_tilda
        for j=1:ncov
            ncovj = ncov + j;
            for i=1:nsamp
                covar(i, ncovj) = covar(i, j)*(time(i)-tvalue);
                for k=1:nsamp
                    covart(k,j,i) = covar(k, j);
                    covart(k,ncovj,i) = covar(k, j)*(time(i)-tvalue);
                end
            end
        end
        % complete case only: calculate covarc and covartc
        for j=1:ncov
            ncovj = ncov + j;
            for i=1:nsampc
                covarc(i, ncovj) = covarc(i, j)*(timec(i)-tvalue);
                for k=1:nsampc
                    covartc(k,j,i) = covarc(k, j);
                    covartc(k,ncovj,i) = covarc(k, j)*(timec(i)-tvalue);
                end
            end
        end
        
        
        % calculate wkernk = K_h(u-t)
        for i=1:nsamp
            wkernk(i) = epanker(time(i), tvalue, h);
        end
        % complete case only: calculate wkernkc = K_h(u-t)
        for i=1:nsampc
            wkernkc(i) = epanker(timec(i), tvalue, h);
        end
        
                %% full estimation
                beta0_f = zeros(ncovt,1);
                beta_sig_f = zeros(ncov,1);
                beta_f = estf(nsamp,ncovt,time,covar,covart,deltacs,wkernk,beta0_f);
        
                if ~isnan(beta_f(1,1))
                    [beta_f, var_f] = varf(nsamp,ncov,time,covar2,deltacs,wkernk,beta_f);
                    if ~isnan(var_f(1,1))
                        for j=1:ncov
                            beta_sig_f(j) = (var_f(j,j)/h)^(1/2);
                        end
                    else
                        beta_sig_f = nan(ncov,1);
                    end
                else
                    beta_sig_f = nan(ncov,1);
                end
        
        %% complete estimation
        beta0_c = zeros(ncovt,1);
        beta_sig_c = zeros(ncov,1);
        beta_c = estf(nsampc,ncovt,timec,covarc,covartc,deltacsc,wkernkc,beta0_c);
        
        if ~isnan(beta_c(1,1))
            [beta_c, var_c] = varf(nsampc,ncov,timec,covar2c,deltacsc,wkernkc,beta_c);
            if ~isnan(var_c(1,1))
                for j=1:ncov
                    beta_sig_c(j) = (var_c(j,j)/h)^(1/2);
                end
            else
                beta_sig_c = nan(ncov,1);
            end
        else
            beta_sig_c = nan(ncov,1);
        end
        
        %% ipw-c estimation
        beta0_ic = zeros(ncovt,1);
        beta_sig_ic = zeros(ncov,1);
        beta_ic = esti(nsamp,ncovt,time,covar,covart,deltacs,wkernk,beta0_ic,wipw);
        
        if ~isnan(beta_ic(1,1))
            [beta_ic, var_ic] = vari(nsamp,ncovt,time,covar,covart,deltacs,wkernk,beta_ic,wipw);
            if ~isnan(var_ic(1,1))
                for j=1:ncov
                    beta_sig_ic(j) = var_ic(j,j)^(1/2);
                end
            else
                beta_sig_ic = nan(ncov,1);
            end
        else
            beta_sig_ic = nan(ncov,1);
        end
        
        %% save the results
        for j=1:ncov
            %% full result
            %             sbeta_f(j,igrid,ics) = beta_f(j);
            %             sstd_f(j,igrid,ics) = beta_sig_f(j);
            %% complete result
            sbeta_c(j,igrid,ics) = beta_c(j);
            sstd_c(j,igrid,ics) = beta_sig_c(j);
            %% ipw-c result
            ipwbeta_c(j,igrid,ics) = beta_ic(j);
            sbeta_ic(j,igrid,ics) = beta_ic(j);
            sstd_ic(j,igrid,ics) = beta_sig_ic(j);
        end
        tvalue = tvalue + tstep;
    end
end

%% for aipw
% stage 1: ipw estimator of lambda_0(t) at each failure time point?
[tsort, tindex] = sort(time);
% ipw-c
lamt0_c = zeros(nsamp,ncs);
% lamt0s_c = zeros(nsamp,ncs);
lamt_c = zeros(nsamp,ncs);
for ics=1:ncs
    for j=1:ncov
        temp0 = ipwbeta_c(j,:,ics);
        nan0 = isnan(temp0);
        temp=interp1(tgrid(~nan0),temp0(~nan0),tsort,'linear','extrap');
        ipwbeta_c_all(j,tindex,ics)=temp;
    end
end
for ics=1:ncs
    deltacs = cause == ics; % deltacs: indicator of non-censoring due to each cause
    % lambda_0(t)
    SI_c = zeros(nsamp,1);
    for i=1:nsamp
        for k=1:nsamp
            if time(k) >= time(i) && time(i) < tau
                zb = 0;
                for j=1:ncov
                    zb = zb + ipwbeta_c_all(j,i,ics)*covar2(k,j);
                end
                SI_c(i) = SI_c(i) + exp(zb)*wipw(k);
            end
        end
    end
    for it=1:nsamp
        tvalue=time(it);
        % if tvalue <= tau
        for i=1:nsamp
            wkernk(i) = epanker(time(i), tvalue, b);
        end
        wt = wkernk.*wipw.*deltacs;
        for m=1:nsamp
            if wt(m) > 0
                lamt0_c(it,ics) = lamt0_c(it,ics) + wt(m)/SI_c(m);
            end
        end
        % lamt0s_c(it,ics) = lamt0s_c(it,ics) + lamt0_c(it,ics);
        % end
    end
    % lambda(t)
    for i=1:nsamp
        ezz = 0;
        for j=1:ncov
            ezz = ezz + ipwbeta_c_all(j,i,ics)*covar2(i,j);
        end
        lamt_c(i,ics) = lamt0_c(i,ics)*exp(ezz);
    end
end

rhohat_cc = zeros(nsamp,ncs);
for ics=1:ncs
    for i=1:nsamp
        % rhohat_cc
        if lamt_c(i,1)~=0 || lamt_c(i,2)~=0
            rhohat_cc(i,ics)=(lamt_c(i,ics)*hhat(i,ics))/(lamt_c(i,1)*hhat(i,1)+lamt_c(i,2)*hhat(i,2));
        end
    end
end

% stage 2 aipw estimation
for ics=1:ncs
    deltacs = cause == ics; % deltacs: indicator of non-censoring due to each cause
    tvalue = tstep;
    for igrid=1:ngrid
        % calculate covar and covart
        for j=1:ncov
            ncovj = ncov + j;
            for i=1:nsamp
                covar(i, ncovj) = covar(i, j)*(time(i)-tvalue);
                for k=1:nsamp
                    covart(k,j,i) = covar(k, j);
                    covart(k,ncovj,i) = covar(k, j)*(time(i)-tvalue);
                end
            end
        end
        
        % calculate wkernk = K_h(u-t)
        for i=1:nsamp
            wkernk(i) = epanker(time(i), tvalue, h);
        end
        
        %% aipw-cc estimation
        beta0_acc = zeros(ncovt,1);
        beta_sig_acc = zeros(ncov,1);
        beta_acc = esta(nsamp,ncovt,time,covar,covart,deltacs,wkernk,beta0_acc,wipw,rhohat_cc(:,ics),delta);
        inv_acc = zeros(ncov,ncov);
        
        if sum(isnan(beta_acc))==0
            [beta_acc, var_acc, inv_acc] = vara(nsamp,ncovt,time,covar,covart,deltacs,wkernk,beta_acc,wipw,rhohat_cc(:,ics),delta);
            if sum(isnan(var_acc))==0
                for j=1:ncov
                    beta_sig_acc(j) = var_acc(j,j)^(1/2);
                end
            else
                beta_sig_acc = nan(ncov,1);
            end
        else
            beta_sig_acc = nan(ncov,1);
        end
        
        %% save aipw results
        for j=1:ncov
            % aipw_cc result
            sbeta_acc(j,igrid,ics) = beta_acc(j);
            sstd_acc(j,igrid,ics) = beta_sig_acc(j);
            for jj=1:ncov
                sinv_acc(j,jj,igrid,ics) = inv_acc(j,jj);
            end
        end
        % igrid
        tvalue = tvalue + tstep;
    end
    % ics
end

%% calculate the cumulative baseline function

[tsort, tindex] = sort(time);
LamtA0_c = zeros(ngrid,ncs);

for ics=1:ncs
    for j=1:ncov
        temp0 = sbeta_acc(j,:,ics);
        nan0 = isnan(temp0);
        temp=interp1(tgrid(~nan0),temp0(~nan0),tsort,'linear','extrap');
        beta_acc_all(j,tindex,ics)=temp;
    end
end

for ics=1:ncs
    deltacs = cause == ics; % deltacs: indicator of non-censoring due to each cause
    % lambda_0(t)
    Sf_acc = zeros(nsamp,1);
    for i=1:nsamp
        for k=1:nsamp
            if time(k) >= time(i) && time(i) < tau
                zb = 0;
                for j=1:ncov
                    zb = zb + beta_acc_all(j,i,ics)*covar2(k,j);
                end
                Sf_acc(i) = Sf_acc(i) + exp(zb);
            end
        end
    end
    
    for it=1:ngrid
        tvalue = tgrid(it);
        wt = (wipw.*deltacs+(1-wipw).*rhohat_cc(:,ics)).*delta;
        for m=1:nsamp
            if wt(m) > 0 && time(m) < tvalue
                LamtA0_c(it,ics) = LamtA0_c(it,ics) + wt(m)/SI_c(m);
            end
        end
        % lamt0s_c(it,ics) = lamt0s_c(it,ics) + lamt0_c(it,ics);
        % end
    end
end

%% hypothesis testing
% calculate lambda0
% sbeta_acc(1,:,ics,isim)
% only for one covariate, ncov = 1

% nt1 = 1;
% t1 = tstep*nt1;
% ngap = 1;
% nt11 = nt1 + ngap;
% t11 = tstep*nt11;
% nt2 = 20;
% t2 = tstep*nt2;

nt1=max(fix(t1/tstep),1);
t1 = tstep*nt1;
ngap = 1;
nt11 = nt1 + ngap;
t11 = tstep*nt11; % t1* for test 2
nt2=min(fix(t2/tstep),ngrid);
t2 = tstep*nt2;



D = 1000; % Gaussian resampling


pv_test1_a1=NaN([ncov ,ncs]);
pv_test1_a2=NaN([ncov ,ncs]);
pv_test1_m1=NaN([ncov ,ncs]);
pv_test1_m2=NaN([ncov ,ncs]);
pv_test2_a1=NaN([ncov ,ncs]);
pv_test2_a2=NaN([ncov ,ncs]);
pv_test2_m1=NaN([ncov ,ncs]);
pv_test2_m2=NaN([ncov ,ncs]);

for ics=1:ncs
    % beta_all = zeros(nsamp);
    beta_hat = sbeta_acc(:,:,ics);
    inv_hat = sinv_acc(:,:,:,ics);
    beta_all=zeros(nsamp, ncov);
    inv_all=zeros(ncov, ncov,nsamp);
    for icov=1:ncov
        beta_all(:,icov)=interp1(tgrid,beta_hat(icov,:),time,'linear','extrap');
        for jcov=1:ncov
            inv_all(icov,jcov,:)=interp1(tgrid,reshape(inv_hat(icov,jcov,:),ngrid,1),time,'linear','extrap');
        end
    end
    
    % bz = covar * beta;
    bzt = sum(covar2 .* beta_all,2);
    
    U=zeros(ncov, 1);
    F=zeros(ncov, ncov);
    % FU=zeros(ncov, ncov);
    
    for icov=1:ncov
        H_hat = zeros(nsamp,nt2,ncov);
        % s2 = zeros(ncov, ncov, nsamp);
        
        for igrid=1:nt2
            tvalue = igrid*tstep;
            deltacs = cause == ics;
            % calculate covar and covart
            %             for j=1:ncov
            %                 ncovj = ncov + j;
            %                 for i=1:nsamp
            %                     covar(i, ncovj) = covar(i, j)*(time(i)-tvalue);
            %                     for k=1:nsamp
            %                         covart(k,j,i) = covar(k, j);
            %                         covart(k,ncovj,i) = covar(k, j)*(time(i)-tvalue);
            %                     end
            %                 end
            %             end
            
            
            for i=1:nsamp
                wkernk(i) = epanker(time(i), tvalue, h);
            end
            
            %wt = wkernk.*delta.*(time>=t1).*(time<=t2);
            wt = wkernk.*delta.*(time<=t2);
            s0 = zeros(nsamp,1);
            s1 = zeros(nsamp,ncov);
            %H_hat_0 = zeros(nsamp,1);
            
            for i=1:nsamp
                if wt(i)>0
                    s0(i)=sum((time >= time(i)).*exp(bzt(i)));
                    for iicov=1:ncov
                        s1(i,iicov)=sum((time >= time(i)).*exp(bzt(i)).*covar2(:,iicov));
                    end
                    % also try use covart and ncov*2
                    H_hat(i,igrid,:) = wt(i)*nsamp*(covar2(i,:)-s1(i,:)/s0(i))*inv_hat(:,:,igrid)*(wipw(i)*deltacs(i)+(1-wipw(i))*rhohat_cc(i,ics));
                    %H_hat(i,:) = (tgrid >= time(i)).*(time(i)>=t1)*H_hat_0(i);
                end
            end
        end
        
        Beta_hat = cumsum(beta_hat(icov,nt1:nt2-1))*tstep;
        H_hat = cumsum(H_hat(:,nt1:nt2-1,icov),2)*tstep;
        
        %test 1
        p_t1_a1=0; %zeros(1,ncov);
        p_t1_a2=p_t1_a1;
        p_t1_m1=p_t1_a1;
        p_t1_m2=p_t1_a1;
        
        G=normrnd(0,1,D,nsamp);
        
        W_hat=sqrt(nsamp)*Beta_hat;
        if ics==1
            W_hat1_cs_1 = W_hat;
        else
            W_hat1_cs_2 = W_hat;
        end
        
        V_hat=H_hat;
        W_tuta=G*V_hat/sqrt(nsamp);
        if ics==1
            W_tuta1_cs_1 = W_tuta;
        else
            W_tuta1_cs_2 = W_tuta;
        end
        
        F1_hat=max(abs(W_hat));
        F2_hat=nansum(W_hat.^2)*tstep;
        F3_hat=min(W_hat);
        F4_hat=nansum(W_hat)*tstep;
        
        F1_tuta=max(abs(W_tuta),[],2);
        F2_tuta=sum(W_tuta.^2,2)*tstep;
        F3_tuta=min(W_tuta,[],2);
        F4_tuta=sum(W_tuta,2)*tstep;
        
        
        p_t1_a1=mean(F1_tuta>F1_hat);
        p_t1_a2=mean(F2_tuta>F2_hat);
        p_t1_m1=mean(F3_tuta<F3_hat);
        p_t1_m2=mean(F4_tuta<F4_hat);
        
        %test 2
        p_t2_a1=0;
        p_t2_a2=p_t2_a1;
        p_t2_m1=p_t2_a1;
        p_t2_m2=p_t2_a1;
        
        G=normrnd(0,1,D,nsamp);
        
        V_hat=((H_hat(:,ngap:nt2-nt1))./(ngap:nt2-nt1)-...
            (H_hat(:,nt2-nt1))./(nt2-nt1))/tstep;
        W_hat=sqrt(nsamp)*((Beta_hat(ngap:nt2-nt1))./(ngap:nt2-nt1)-...
            (Beta_hat(nt2-nt1))./(nt2-nt1))/tstep;
        if ics==1
            W_hat2_cs_1 = W_hat;
        else
            W_hat2_cs_2 = W_hat;
        end
        W_tuta=G*V_hat/sqrt(nsamp);
        if ics==1
            W_tuta2_cs_1 = W_tuta;
        else
            W_tuta2_cs_2 = W_tuta;
        end
        
        
        F1_hat=max(abs(W_hat));
        F2_hat=nansum(W_hat.^2)*tstep;
        F3_hat=min(W_hat);
        F4_hat=nansum(W_hat)*tstep;
        
        F1_tuta=max(abs(W_tuta),[],2);
        F2_tuta=sum(W_tuta.^2,2)*tstep;
        F3_tuta=min(W_tuta,[],2);
        F4_tuta=sum(W_tuta,2)*tstep;
        
        
        p_t2_a1=mean(F1_tuta>F1_hat);
        p_t2_a2=mean(F2_tuta>F2_hat);
        p_t2_m1=mean(F3_tuta<F3_hat);
        p_t2_m2=mean(F4_tuta<F4_hat);
        
        
        pv_test1_a1(icov,ics)=p_t1_a1;
        pv_test1_a2(icov,ics)=p_t1_a2;
        pv_test1_m1(icov,ics)=p_t1_m1;
        pv_test1_m2(icov,ics)=p_t1_m2;
        
        pv_test2_a1(icov,ics)=p_t2_a1;
        pv_test2_a2(icov,ics)=p_t2_a2;
        pv_test2_m1(icov,ics)=p_t2_m1;
        pv_test2_m2(icov,ics)=p_t2_m2;
    end
end

save simres.mat


