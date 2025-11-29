var 
    x      $x_t$                  (long_name = 'Output gap')
    piH    $\pi_t^H$              (long_name = 'Domestic inflation')
    r      $i_t$                  (long_name = 'Nominal interest rate')
    y      $y_t$                  (long_name = 'Output')
    c      $c_t$                  (long_name = 'Consumption')
    s      $s_t$                  (long_name = 'Terms of trade')
    n      $n_t$                  (long_name = 'Employment')
    a      $a_t$                  (long_name = 'Technology')
    ystar  $y_t^*$                (long_name = 'Foreign output')
    rnat   $r_t^{nat}$            (long_name = 'Natural real rate')
    pi     $\pi_t$                (long_name = 'CPI inflation')
    q      $q_t$                  (long_name = 'Real exchange rate')
    e      $e_t$                  (long_name = 'Nominal exchange rate (log)')
;

varexo 
    ea      $\varepsilon_t^a$     (long_name = 'Technology shock')
    eystar  $\varepsilon_t^{y^*}$ (long_name = 'Foreign demand shock')
    em      $\varepsilon_t^m$     (long_name = 'Monetary policy shock')
;

parameters
    beta          ${\beta}$               (long_name = 'Discount factor')
    sigma         ${\sigma}$              (long_name = 'Inverse IES')
    varphi        ${\varphi}$             (long_name = 'Inverse Frisch elasticity')
    lambda        ${\lambda}$             (long_name = 'Home bias')
    eta           ${\eta}$                (long_name = 'Elasticity of substitution')
    epsilon       ${\epsilon}$            (long_name = 'Price markup parameter')
    phi_p         ${\phi_p}$              (long_name = 'Rotemberg adjustment cost')
    phi_pi        ${\phi_{\pi}}$          (long_name = 'Taylor rule inflation coef.')
    rho_a         ${\rho_a}$              (long_name = 'Persistence technology shock')
    rho_star      ${\rho_*}$              (long_name = 'Persistence foreign demand')
    omega         ${\omega}$              (long_name = 'Composite parameter \omega')
    sigma_lambda  ${\sigma_{\lambda}}$    (long_name = 'Effective risk aversion')
    psi_rot       ${\psi}$                (long_name = 'Rotemberg slope parameter')
    kappa_lambda  ${\kappa_{\lambda}}$    (long_name = 'Slope of NKPC')
    chi           ${\chi}$                (long_name = 'Share of foreign goods')
;

beta     = 0.99;
sigma    = 1.00;
varphi   = 1.00;
lambda   = 0.30;
eta      = 1.50;
epsilon  = 6.00;
phi_p    = 60.00;
phi_pi   = 1.50;
rho_a    = 0.90;
rho_star = 0.90;
chi      = 0.20;

omega        = eta*sigma + (1 - lambda)*(eta*sigma - 1);
sigma_lambda = sigma / ( (1 - lambda) + lambda*omega );
psi_rot      = (epsilon - 1)/phi_p;
kappa_lambda = psi_rot * (sigma_lambda + varphi);

model;
    // NKPC (Rotemberg)
    piH = beta*piH(+1) + kappa_lambda*x;

    // IS equation
    x   = x(+1) - (1/sigma_lambda)*( r - piH(+1) - rnat );

    // Taylor monetary rule
    r   = phi_pi*piH + em;

    // AR(1) with technology
    a      = rho_a*a(-1) + ea;

    // AR(1) foreign demand
    ystar  = rho_star*ystar(-1) + eystar;

    // real natural rate 
    rnat = -sigma_lambda*(1 - rho_a)*a + lambda*sigma_lambda*(omega + chi)*( ystar - ystar(-1) );

    // terms of trade
    s = (sigma/(1 - lambda))*(c - ystar);

    // goods mkt clearing condition
    y = c + (lambda*omega/sigma)*s;

    // production
    n = y - a;

    // output gap
    x = y - a;

    // CPI inflation
    pi = piH + lambda * ( s - s(-1) );

    // real exchange rate
    q = (1 - lambda) * s;

    // nominal exchange rate (log) 
    e = e(-1) + ( q - q(-1) ) + pi;

end;

initval;
  x      = 0;
  piH    = 0;
  r      = 0;
  y      = 0;
  c      = 0;
  s      = 0;
  n      = 0;
  a      = 0;
  ystar  = 0;
  rnat   = 0;
  pi     = 0;
  q      = 0;
  e      = 0;
end;

steady;

steady_state_model;
  x     = 0;
  piH   = 0;
  r     = 0;
  y     = 0;
  c     = 0;
  s     = 0;
  n     = 0;
  a     = 0;
  ystar = 0;
  rnat  = 0;
  pi     = 0;
  q      = 0;
  e      = 0;
end;

shocks;
    var ea;     stderr 0.01;   // productivity shock
    var eystar; stderr 0;
    var em;     stderr 0;
end;

stoch_simul(order=1, irf=20);

shocks;
  var ea;     stderr 0;   
  var eystar; stderr 0.01;  // shock to foreign demand
  var em;     stderr 0.0;
end;

stoch_simul(order=1, irf=20);

write_latex_definitions;
write_latex_parameter_table;
write_latex_dynamic_model;
write_latex_steady_state_model;
collect_latex_files;