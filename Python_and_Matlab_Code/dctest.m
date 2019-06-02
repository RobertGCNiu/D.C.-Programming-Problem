clear;clc;
H = [0.4310, 0.0002, 0.2605, 0.0039; 
    0.0002, 0.3018, 0.0008, 0.0054; 
    0.0129, 0.0005, 0.4266, 0.1007; 
    0.0011, 0.0031, 0.0099, 0.0634];
M = 4;
w = [1.0/6, 1.0/6, 1.0/3, 1.0/3];
p_max = [0.7, 0.8, 0.9, 1.0];
ri = 0;
sigma = 0.0001;
p_initial = p_max';
t_old = 0;
max_value_dif = 100;
iter = 0;
while(abs(max_value_dif)>0.001)
    
iter = iter + 1;
cvx_begin quiet
variables p(M)
 
 f_p = 0;
 g_p = 0;
 for i = 1 : M
     H_f = 0;
     H_g = 0;
     for j = 1:M
         H_f = H_f +  H(j,i) * p(j);
        if j ~= i
            H_g = H_g + H(j,i) * p_initial(j);
        end
     end
     f_p  = f_p + w(i) * log(sigma + H_f)/log(2);   % f(p)
     g_p  = g_p + w(i) * log(sigma + H_g)/log(2);  % g(p)
 end
    gra_g = g(p_initial);   % gradient of g
    t= f_p - g_p - gra_g'*(p-p_initial);
    maximize(t)
    subject to
            for i = 1:M
                He = 0;
                for j = 1:M
                   if j ~= i
                        He = He + H(j,i) * p(j);
                   end
                end
                 H(i,i) * p(i) + (1-2^ri)*(He + sigma) >= 0;
                p(i)<=p_max(i);
                p(i)>=0;
            end
cvx_end
    p_initial = p;
    max_value_dif = t-t_old;
    t_old = t
end