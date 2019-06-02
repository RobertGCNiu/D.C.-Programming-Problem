clear;clc;
H = [-58.8228781068589 + 198.892794640064i,-1.21833110147180 + 2.40402487088321i,0.0243868782191938 - 0.0111384468559801i,-2.14935564487486 - 0.719225409695810i;
    -0.245272044320914 - 0.660935714985225i,-9.38377971976335 - 53.4355097996488i,0.0425170545042271 - 0.0779284384818814i,0.926405105143966 - 0.306989287332009i;
    -0.0141033193981806 + 0.0157636408248593i,0.211613810724130 - 0.164038954460624i,163.618477353170 + 2.26885011246024i,1.50830907749530 + 1.44863962969990i;
    1.03197981197851 + 1.00090787309430i,1.37347315295496 + 1.92724539124343i,-0.0909444831489928 + 1.67891254896088i,-94.8748793444937 + 91.1404341726689i];
M = 4;
K =2;
m_k = M/K;
%w = [1.0/6, 1.0/6, 1.0/3, 1.0/3];
p_max = [0.7, 0.8, 0.9, 1.0];
ri = 0;
sigma = 0.001;
p_initial_v = rand(m_k,1) + 1j*rand(m_k,1);
p_initial = p_initial_v * p_initial_v';
t_old = 0;
r_old = 0;

iter = 0;
for u = 1:M
    Fbb(:, (u-1)*m_k+1:u*m_k) =  p_initial;
    Fbb_2(:, (u-1)*m_k+1:u*m_k) =  p_initial;
end
a = rand(m_k,1);
b = rand(m_k,1);
p = (a+1j*b)*(a+1j*b)'/1000;
for u = 1:1
    Fbb_2(:, (u-1)*m_k+1:u*m_k) =  p_initial+p;
end
conver = 100;

%while(abs(conver)>0.001)
u=1;
f_p = 0;
g_p = 0;
g_p_2=0;
% for i = 1 : M
%     H_f = 0;
%     H_g = 0;
%     H_g_2 = 0;
%     for j = 1:M
%         H_f = H_f +  real(H(i,:) * Fbb(:, (j-1)*M+1:j*M) * H(i,:)');
%         if j ~= i
%             H_g = H_g + real(H(i,:) * Fbb(:, (j-1)*M+1:j*M)  * H(i,:)');
%             H_g_2 = H_g_2 + real(H(i,:) * Fbb_2(:, (j-1)*M+1:j*M)  * H(i,:)');
%         end
%     end
% end


for i = 1 : M
     H_f = 0;
     H_g = 0;
     H_g_2 = 0;
    for k = 1:K
        for j = 1+(k-1)*m_k : k*m_k
            if j ~= i
                H_g = H_g + real(H(i,1+(k-1)*m_k : k*m_k) * Fbb(:, (j-1)*m_k+1:j*m_k) * H(i,1+(k-1)*m_k : k*m_k)');
                H_g_2 = H_g_2 + real(H(i,1+(k-1)*m_k : k*m_k) * Fbb_2(:, (j-1)*m_k+1:j*m_k) * H(i,1+(k-1)*m_k : k*m_k)');
            end
        end
    end
    
end
%f_p  = f_p +  log(sigma + H_f)/log(2);   % f(p)
g_p  = g_p +  log(sigma + H_g)/log(2);  % g(p)m
g_p_2  = g_p_2 +  log(sigma + H_g_2)/log(2);  % g(p)
gra_g = g_block(H, Fbb,u,sigma);   % gradient of g

t2 = g_p + real((trace(gra_g'*p)));
t3 = g_p_2 ;
t2>  t3

