% Code for evaluation of second-order SPSA (2SPSA) versus first-order SPSA
% (1SPSA).  Function measurements include added Gaussian noise.
% Code includes the capability for initializing 2SPSA by running 1SPSA for
% N measurements.  This code allows for averaging of the
% SP gradients and Hessian estimates at EACH iteration after the initial (N)
% measurements where only 1SPSA is used for estimating theta.  We use "theta2"
% for the 2SPSA recursion and "theta" for the 1SPSA recursion.

% Because of the 1SPSA component to this code and the comparative aspects, there is
% material here that can be eliminated for basic 2SPSA implementation.  This is test
% code, not production code, so it includes more material than necessary.  A user can
% easily strip out this "extra" material if desired.

global z sigma p;       %declaration of random var. used for normal noise
                        %generation in loss fn. calls given seed above
p=10;
%value of numerator in a_k sequence for all iterations of 1SPSA
%and first N-measurement-based iterations in the initialization of 2SPSA
a1=.01;
%value of numerator in a_k in 2SPSA part
a2=.5;
A1=100;           %stability constant for 1SPSA
A2=100;           %stability constant for 2SPSA
c1=.19;         	%numerator in c_k for 1SPSA
c2=2*c1;         	%numerator in c_k for 2SPSA
ctilda=2*c2;     	%numerator in ctilda_k for 2SPSA;
alpha1=.602;      %a_k decay rate for 1SPSA
alpha2=.602;      %a_k decay rate for 2SPSA
gamma1=.101;      %c_k decay rate for 1SPSA
gamma2=.101;      %c_k decay rate for 2SPSA
n=16000;	    		%total no. of function measurements
N=0;					%no. of function meas. for 1SPSA
                                        %initialization

sigma=.05;
lossfinaleval=quartic_loss_factory(p); 	%loss function for "true" evaluation of algorithm (no noise)
loss=noisy_function_factory(lossfinaleval, sigma);	  				%loss function for use in algorithm (usually with noise)
cases=5;         			%number of cases (replications) of 2SPSA and 1SPSA
gH_avg=1;               %no. of averaged gradients/Hessian in 2SPSA
toleranceloss=0;			%tolerance in loss-based blocking step
avg=0;						%no. of loss evals. per loss-based blocking step (1SPSA&2SPSA)
tolerancetheta=1000;		%max. allowable change in elements of theta
rand('seed',31415297)
randn('seed',3111113)
%
%the loop 1:cases below is for doing multiple cases for use in averaging to
%evaluate the relative performance.
%
%the first loop in 2SPSA below uses the standard 1SPSA form to initialize 2SPSA
%
%the second loop does 2SPSA following guidelines in Spall, 1999 ASP
%
%lines below initialize various recuresions for the gradient/Hess. averaging
%and for final error reporting based on the average of the solutions for
%"cases" replications.
%
meanHbar=0;
errtheta=0;				%cum. sum of theta errors
errthetaIA=0;
errtheta2=0;
losstheta=0;				%cum. sum of loss values
lossthetaIA=0;
losstheta2=0;
lossthetasq=0;				%cum. sum of loss squared values
lossthetaIAsq=0;
losstheta2sq=0;
truetheta=ones(p,1);
theta_0=-.52276*ones(p,1);
%dummy statement for setting dimension of Hhat (avoids occasional error message)
Hhat=eye(p);
for j=1:cases
%INITIALIZATION OF PARAMETER AND HESSIAN ESTIMATES
  theta=theta_0;
  theta2=theta;
  Hbar=500*eye(p);
  lossold=0;	%lossold calculation is for use in loss-based blocking step below
  for i=1:avg
    lossold=lossold+feval(loss,theta2);
  end
%*********2SPSA*********
%Initial iterations of 1SPSA based on N measurements (prior to 2SPSA iterations)
  for k=1:(N-avg)/(2+avg)    %use of N-avg is to account for avg used in setting lossold
    a_k=a1/(k+A1)^alpha1;
    c_k=c1/k^gamma1;
    delta=2*round(rand(p,1))-1;
    thetaplus=theta2+c_k*delta;
    thetaminus=theta2-c_k*delta;
    yplus=feval(loss,thetaplus);
    yminus=feval(loss,thetaminus);
    ghat=(yplus-yminus)./(2*c_k*delta);
%   theta update
    theta2lag=theta2;
    theta2=theta2-a_k*ghat;
    % Steps below perform "blocking" step with "avg" no. of loss evaluations
    % This blocking is based on extra loss evaluation(s)
    lossnew=0;
    for i=1:avg
      lossnew=lossnew+feval(loss,theta2);
    end
    if lossnew > lossold-avg*toleranceloss; %if avg=0, this statement is always false
      theta2=theta2lag;
    else                                    %statements to follow are harmless when avg=0
      lossold1=lossold;
      lossold=lossnew;
    end
    % Blocking below is based on magnitude of theta change (no loss evals.)
    if max(abs(theta2lag-theta2)) > tolerancetheta;
      theta2=theta2lag;
      lossold=lossold1;     %only relevant if also using loss-based blocking
    end
 end
 caseiter=j   %print-out of iteration number (for monitoring progress)
%
% START 2SPSA ITERATIONS FOLLOWING INITIALIZATION
%
for k=(N-avg)/(2+avg)+1:(N-avg)/(2+avg)+(n-N)/(gH_avg*4+avg)
    a_k=a2/(k+A2)^alpha2;
    c_k=c2/k^gamma2;
    ctilda_k=ctilda/k^gamma2;
    ghatinput=0;
    Hhatinput=0;
% GENERATION OF AVERAGED GRADIENT AND HESSIAN (NO AVERAGING IF gH_avg=1)
    for m=1:gH_avg
      delta=2*round(rand(p,1))-1;
      thetaplus=theta2+c_k*delta;
      thetaminus=theta2-c_k*delta;
      yplus=feval(loss,thetaplus);
      yminus=feval(loss,thetaminus);
      ghat=(yplus-yminus)./(2*c_k*delta);
% GENERATE THE HESSIAN UPDATE
      deltatilda=2*round(rand(p,1))-1;
      thetaplustilda=thetaplus+ctilda_k*deltatilda;
      thetaminustilda=thetaminus+ctilda_k*deltatilda;
% LOSS FUNCTION CALLS
      yplustilda=feval(loss,thetaplustilda);
      yminustilda=feval(loss,thetaminustilda);
      ghatplus=(yplustilda-yplus)./(ctilda_k*deltatilda);
      ghatminus=(yminustilda-yminus)./(ctilda_k*deltatilda);
% STATEMENT PROVIDING AN AVERAGE OF SP GRAD. APPROXS. PER ITERATION
      ghatinput=((m-1)/m)*ghatinput+ghat/m;
      deltaghat=ghatplus-ghatminus;
      for i=1:p
        Hhat(:,i)=deltaghat(i)./(2*c_k*delta);
      end
      Hhat=.5*(Hhat+Hhat');
      Hhatinput=((m-1)/m)*Hhatinput+Hhat/m;
    end
    Hbar=((k-(N-avg)/(2+avg))/(k-(N-avg)/(2+avg)+1))*Hbar+Hhatinput/(k-(N-avg)/(2+avg)+1);
%   THE THETA UPDATE (FORM BELOW USES GAUSSIAN ELIMINATION TO AVOID DIRECT
%   COMPUTATION OF HESSIAN INVERSE)
    Hbarbar=sqrtm(Hbar*Hbar)+.000001*eye(p)/k;
    theta2lag=theta2;
    theta2=theta2-a_k*(Hbarbar\ghatinput);
%   Steps below perform "blocking" step with "avg" no. of loss evaluations
    lossnew=0;
    for i=1:avg
      lossnew=lossnew+feval(loss,theta2);
    end
    if lossnew > lossold-avg*toleranceloss;
      theta2=theta2lag;
    else
      lossold1=lossold;
      lossold=lossnew;
    end
    if max(abs(theta2lag-theta2)) > tolerancetheta;
      theta2=theta2lag;
      lossold=lossold1;
   end
 end
theta2;
%
%********1SPSA*************
% The iterations below are the 1SPSA iterations.  Uses the same gain sequences
% as the 1SPSA loop above (where 2SPSA is not fully engaged).  Uses same number
% of loss function measurements.  The overall loop is broken into two parts to
% accomodate a sliding window of the last IA iterates for an iterate averaging
% solution.
%
 IA=200;		%no. of sliding window iterations in iterate averaging
 lossold=0;	%lossold calculation is for use in loss-based blocking step below
 for i=1:avg
   lossold=lossold+feval(loss,theta);
 end
 for k=1:(n-avg)/(2+avg)-IA
    a_k=a1/(k+A1)^alpha1;
    c_k=c1/k^gamma1;
    delta=2*round(rand(p,1))-1;
    thetaplus=theta+c_k*delta;
    thetaminus=theta-c_k*delta;
    yplus=feval(loss,thetaplus);
    yminus=feval(loss,thetaminus);
    ghat=(yplus-yminus)./(2*c_k*delta);
%   theta update
    thetalag=theta;
    theta=theta-a_k*ghat;
%   Steps below perform "blocking" steps as in 1SPSA part of 2SPSA above
    lossnew=0;
    for i=1:avg
      lossnew=lossnew+feval(loss,theta);
    end
    if lossnew > lossold-avg*toleranceloss; %statement always false if avg=0
      theta=thetalag;
    else
      lossold1=lossold;
      lossold=lossnew;
    end
    if max(abs(thetalag-theta)) > tolerancetheta;
      theta=thetalag;
      lossold=lossold1;
   end
  end
  thetabar=0;
  for k=(n-avg)/(2+avg)+1-IA:(n-avg)/(2+avg)
    a_k=a1/(k+A1)^alpha1;
    c_k=c1/k^gamma1;
    delta=2*round(rand(p,1))-1;
    thetaplus=theta+c_k*delta;
    thetaminus=theta-c_k*delta;
%   THE NEXT TWO LINES SHOULD BE CHANGED AS LOSS CHANGES
    yplus=feval(loss,thetaplus);
    yminus=feval(loss,thetaminus);
    ghat=(yplus-yminus)./(2*c_k*delta);
%   theta update
    thetalag=theta;
    theta=theta-a_k*ghat;
%   Steps below perform "blocking" steps as in 1SPSA part of 2SPSA above
    lossnew=0;
    for i=1:avg
      lossnew=lossnew+feval(loss,theta);
    end
    if lossnew > lossold-avg*toleranceloss;
      theta=thetalag;
    else
      lossold1=lossold;
      lossold=lossnew;
    end
    if max(abs(thetalag-theta)) > tolerancetheta;
      theta=thetalag;
      lossold=lossold1;
    end
    thetabar=thetabar+theta;
  end
  %theta
  thetabar=thetabar/IA;
  meanHbar=meanHbar+Hbar;
  errtheta=errtheta+(theta-truetheta)'*(theta-truetheta);
  errthetaIA=errthetaIA+(thetabar-truetheta)'*(thetabar-truetheta);
  errtheta2=errtheta2+(theta2-truetheta)'*(theta2-truetheta);
  lossthetasq=lossthetasq+feval(lossfinaleval,theta)^2;
  lossthetaIAsq=lossthetaIAsq+feval(lossfinaleval,thetabar)^2;
  losstheta2sq=losstheta2sq+feval(lossfinaleval,theta2)^2;
  losstheta=losstheta+feval(lossfinaleval,theta);
  lossthetaIA=lossthetaIA+feval(lossfinaleval,thetabar);
  losstheta2=losstheta2+feval(lossfinaleval,theta2);
end
'meanHar/cases'; meanHbar/cases;
% normalized results of 1SPSA and 2SPSA
if norm(theta_0-truetheta)~= 0
  disp('norm(theta_0-truetheta)~= 0; compute normalized results of 1SPSA and 2SPSA');
  ((errtheta/cases)^.5)/norm(theta_0-truetheta)
  ((errthetaIA/cases)^.5)/norm(theta_0-truetheta)
  ((errtheta2/cases)^.5)/norm(theta_0-truetheta)
end
% standard dev. of mean of normalized loss values; these are by multiplied by
% (cases/(cases-1))^.5 to account for loss of degree of freedom in standard
% deviation calculation before using with t-test
if cases > 1
    disp('cases > 1');
  (cases^(-.5))*((cases/(cases-1))^.5)*(lossthetasq/(cases*feval(lossfinaleval,theta_0)^2)-(losstheta/(cases*feval(lossfinaleval,theta_0)))^2)^.5
  (cases^(-.5))*((cases/(cases-1))^.5)*(lossthetaIAsq/(cases*feval(lossfinaleval,theta_0)^2)-(lossthetaIA/(cases*feval(lossfinaleval,theta_0)))^2)^.5
  (cases^(-.5))*((cases/(cases-1))^.5)*(losstheta2sq/(cases*feval(lossfinaleval,theta_0)^2)-(losstheta2/(cases*feval(lossfinaleval,theta_0)))^2)^.5
end
disp('normalized loss values');
losstheta/(cases*feval(lossfinaleval,theta_0))
lossthetaIA/(cases*feval(lossfinaleval,theta_0))
losstheta2/(cases)%*feval(lossfinaleval,theta_0))
exit;