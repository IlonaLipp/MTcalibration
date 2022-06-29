pars0 = [0, 0]; %%% starting point
newfun = @(cpars,iv)(cpars(1) + iv * cpars(1) * cpars(2)); %%% cpars(1) = delta_corr, cpars(2) = A
[BETA,R,J,COVB,MSE] = nlinfit(iv,dv,newfun,pars0); 
prediction = newfun(BETA, iv);
[r p] = corrcoef(dv,prediction);

se = sqrt(diag(COVB)); %%% coefficient standard errors according to https://uk.mathworks.com/help/stats/nlinfit.html
C_se = se(2) * nominal_mt_pulse_to_calibrate_to_in_rad; %%% se of C, needs to be multiplied by pulse to convert from A
%%% The diagonal elements of COVB are the (estimated) variances of
%%% each coefficient distribution, and their square roots are the standard
%%% deviations. So, assuming by "uncertainty" you mean one standard
%%% deviation, then you are correct. (https://uk.mathworks.com/matlabcentral/answers/41738-fitting-with-nlinfit)

MAPS.R2(x,y,z) = r(1,2)^2;
MAPS.A(x,y,z) = BETA(2);
MAPS.C(x,y,z) = MAPS.A(x,y,z) * nominal_mt_pulse_to_calibrate_to_in_rad; %%% like in new theory section
ci = nlparci(BETA, R,'Jacobian',J); %%% confidence intervals of parmaeters
MAPS.C_CI_low(x,y,z) = ci(2,1) * nominal_mt_pulse_to_calibrate_to_in_rad;
MAPS.C_CI_high(x,y,z) = ci(2,2) * nominal_mt_pulse_to_calibrate_to_in_rad;
MAPS.C_SE(x,y,z) = C_se;

prediction4d(x,y,z,idx) = prediction ./ (local_ft_brain^2); %%% because it predictes the *local_ft_brain^2
