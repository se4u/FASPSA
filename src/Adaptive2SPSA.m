function [k, theta, condition_number, time_taken, ls, md] = ...
    Adaptive2SPSA(budget, Y, InitTheta, p, a_num, c_num, L, OptTheta, OptL)

assert(mod(budget, 5)==0);
MAD=@(m) mean(abs(m));
del = @() 2*round(rand(p,1))-1;
max_iterations=budget/5;
A=max(100, floor(max_iterations/10)); 
a = @(k) a_num /(k+1+A)^1;
c = @(k) c_num /(k+1)^0.1667;
c_tilda = @(k) (2*c_num)/(k+1)^.1667;
symmetric=@(m) (m+m')/2;
Hbar=0;
theta=InitTheta;
budget_cntr=0;

k=0;
ls=NaN(1, max_iterations);
md=NaN(1, max_iterations);
start_loss=L(theta);
ls(1)=1;
md(1)=1;
tic;
cur_theta_loss=Y(theta);
budget_cntr=budget_cntr+1;
while budget_cntr <= budget-5
    del_k=del();
    del_tilda_k=del();
    c_k=c(k);
    c_tilda_k=c_tilda(k);
    theta_plus=theta+c_k*del_k;
    theta_minus=theta-c_k*del_k;
    theta_plus_tilda=theta_plus+c_tilda_k*del_tilda_k;
    theta_minus_tilda=theta_minus+c_tilda_k*del_tilda_k;
    % Create gradient
    Del_Y = (Y(theta_plus)-Y(theta_minus));
    Del_tilda_Y=(Y(theta_plus_tilda)-Y(theta_minus_tilda));
    
    % Create Hessian
    Hhat_by_kp1=((Del_tilda_Y-Del_Y)/((2*((k+1)*c_k))*c_tilda_k))*...
        symmetric(del_tilda_k*del_k');
         
    Hbar = (k/(k+1))*Hbar + Hhat_by_kp1;
    Hbarbar = Hbar + ((1e-6)/(k+1))*eye(size(Hbar));
    
    % Update Theta
    theta_delta=(a(k)/c_k)*(Del_Y/2)*(Hbarbar\del_k);
    theta=theta-theta_delta;
    
    % A dash of greediness, and algorithm B
    tmp=Y(theta);
    if tmp > cur_theta_loss + 2
        theta=theta+theta_delta;
    else
        cur_theta_loss=tmp;
    end
    % We Made 5 loss measurements and one iteration
    budget_cntr = budget_cntr+5;
    k=k+1;
    ls(k+1)=abs(OptL-L(theta))/abs(start_loss-OptL);
    md(k+1)=MAD(theta-OptTheta)/MAD(InitTheta-OptTheta);
end
time_taken=toc;
condition_number=cond(inv(Hbar));
assert(k==max_iterations-1);