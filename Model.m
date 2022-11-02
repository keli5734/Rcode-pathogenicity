function  var = Model(t,y,...
                       s_M, delta_MR,...
                       s_V, k1, k2, alpha1,D_50, k_m1, k_m2,V50_M,...
                       delta_M1,delta_M2,...
                       g_T,T_max,beta,phi,xi_R,delta_I,kappa_F,...
                       p_I,delta_V,...
                       q_FI,q_FM,delta_F,kappa_D,delta_D,...
                       kappa_E, kappa_AS, q_prime,l)
                   
            
var = zeros(10,1);

var(1) = s_M - delta_MR * y(1) - k1 *  (y(6)/(y(6) + V50_M))  * ( 1/(1 + alpha1 * y(3))) * y(1) - k2 * ( y(9)/(y(9) + D_50) ) * y(2) * y(1) + k_m1 * y(2) + k_m2 * y(3);
var(2) = s_V * y(5) + k1 *  (y(6)/(y(6) + V50_M))  * ( 1/(1 + alpha1 * y(3))) * y(1) - delta_M1 * y(2) - k_m1 * y(2);
var(3) = k2 * ( y(9)/(y(9) + D_50) ) * y(2) * y(1) - delta_M2 * y(3) - k_m2 * y(3);

var(4) = g_T*(y(4)+y(8))*(1-(y(4)+y(5)+y(8)+y(10))/T_max) - beta*y(4)*y(6) - phi*y(7)*y(4) + xi_R *y(8); % T
var(5) = l * y(10) - delta_I*y(5)  - kappa_F*y(5)*y(7) - kappa_E * (t^4 / (t^4 + 10^4)) * y(5); % I
var(6) = p_I*y(5) - delta_V*y(6) - q_prime * y(1) * y(6) - kappa_AS *  (t^4 / (t^4 + 10^4)) * y(6); % V 

var(7) = q_FI*y(5) + q_FM*y(2) - delta_F*y(7); % F
var(8) = phi*y(7)*y(4) - xi_R*y(8); % R
var(9) = delta_I*y(5)  + kappa_F*y(5)*y(7) + kappa_E * (t^4 / (t^4 + 10^4)) * y(5) - kappa_D*y(2)*y(9) - delta_D*y(9); % D
var(10) = beta*y(4)*y(6) - l * y(10);

end 

