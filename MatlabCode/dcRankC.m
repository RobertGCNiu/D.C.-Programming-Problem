clear;clc;
data= load('snr16.mat');
a_TX = data.a_TX;
a_RX =data.a_RX;
H_c=data.H_c;
%H_c=data.H;
M = 4;
SNR = 10;
sigma = 1.0/SNR;
ri = 0;
ITER = 1;
r_all = 0;

for iter = 1: ITER  
t_all = [];
  for u=1:1:M
        Frf(:,u)=a_TX(:,u);
        Wrf(:,u)=a_RX(:,u);
  end      
    
    % Constructin the effective channels
    for u=1:1:M
        Channel=zeros(16,64);
        Channel(:,:)= H_c(u,:,:);
        H(u,:)=Wrf(:,u)'*Channel*Frf ;    % Effective channels
    end
    
%p_initial_v = [0.9120; 0.9514; 0.3460; 0.2902];
%p_initial_v = rand(M,1);
p_initial_v =[1;1;1;1] + rand(4,1);
p_initial = p_initial_v * p_initial_v';  
%p_initial = rand(M,M);
t_old = 0;
r_old = 0;
%Fbb = zeros(M,M*M);
iter = 0;
%  for u = 1:M
%      p_initial_v = data.Fbb(:,u);
%      Fbb(:, (u-1)*M+1:u*M) =  abs(p_initial_v * p_initial_v');
%  end
Fbb = data.Fbb;
conver = 100;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       w_k = 0.1;
   t_s = 1 ;
while(abs(conver)>0.01)
    r = 0;
    if iter >=100
        cannotsolve = 1;
        break
    end
  %w_k = w_k*t_s; 
 for u = 1:4 
  
p_initial = Fbb(:, (u-1)*M+1:u*M);
max_value_dif = 100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while(abs(max_value_dif)>0.01)    

cvx_begin sdp quiet 
variable p_0(M,M) hermitian
f_p = 0;
 g_p = 0;
 for i = 1 : M
     H_f = 0;
     H_g = 0;
     for j = 1 : M
         if j == u
         H_f = H_f +  H(i,:) * p_0 * H(i,:)';
         else
          H_f = H_f +  real(H(i,:) * Fbb(:, (j-1)*M+1:j*M) * H(i,:)');   
         end
        if j ~= i
                H_g = H_g + real(H(i,:) * Fbb(:, (j-1)*M+1:j*M)  * H(i,:)');
        end
     end
     f_p  = f_p +  log(sigma + H_f)/log(2);   % f(p)
     g_p  = g_p +  log(sigma + H_g)/log(2);  % g(p)
 end
    gra_g = g(H, Fbb,u,sigma,M);   % gradient of g
    t= -f_p +( g_p + real(trace(gra_g'*(p_0-p_initial))));
    minimize(t)
    subject to
    p_0>=0;
     norm(Frf*p_0*Frf')<=1;
            for i = 1:M
                He = 0;
                for j = 1:M
                   if j ~= i
                       if j ==u
                           He = He + H(i,:) * p_0 * H(i,:)';
                       else
                        He = He + H(i,:) * Fbb(:, (j-1)*M+1:j*M) * H(i,:)';
                       end
                   end
                  % p_0(i,j)>=0;       
                end
                 H(i,:) * p_0 *H(i,:)' + (1-2^ri)*(He + sigma) >= 0;
            end 
cvx_end
p_initial = p_0;
max_value_dif = t-t_old;
t_old= t;
cvx_status;
end
 Fbb(:, (u-1)*M+1:u*M) = p_0;
[vec_x, val_x] = eigs(p_0);    
t_s = 1.1;
w_k = 0.2;
p_initial = p_0;
e_k_old =val_x(2,2); 
iter=0;
max_value_dif = 100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while(abs(max_value_dif)>0.001)    

    iter = iter+1;
cvx_begin sdp quiet
variable p(M,M) hermitian
variable e_k
 f_p = 0;
 g_p = 0;
 for i = 1 : M
     H_f = 0;
     H_g = 0;
     for j = 1 : M
         if j == u
         H_f = H_f +  H(i,:) * p * H(i,:)';
         else
          H_f = H_f +  real(H(i,:) * Fbb(:, (j-1)*M+1:j*M) * H(i,:)');   
         end
        if j ~= i
                H_g = H_g + real(H(i,:) * Fbb(:, (j-1)*M+1:j*M)  * H(i,:)');
        end
     end
     f_p  = f_p +  log(sigma + H_f)/log(2);   % f(p)
     g_p  = g_p +  log(sigma + H_g)/log(2);  % g(p)
 end
    gra_g = g(H, Fbb,u,sigma,M);   % gradient of g
    t= -f_p + g_p + real(trace(gra_g'*(p-p_initial))) +  w_k*e_k;
    minimize(t)
    subject to
        p>=0;
        e_k*eye(3)- vec_x(:,2:end)'*p*vec_x(:,2:end)>=0;
        e_k<=e_k_old;
        norm(Frf*p*Frf')<=1;
            for i = 1:M
                He = 0;
                for j = 1:M
                   if j ~= i
                       if j ==u
                           He = He + H(i,:) * p * H(i,:)';
                       else
                        He = He + H(i,:) * Fbb(:, (j-1)*M+1:j*M) * H(i,:)';
                       end
                   end
                end
                 H(i,:) * p *H(i,:)' + (1-2^ri)*(He + sigma) >= 0;
            end
cvx_end  
    e_k_old = e_k;
    p_initial = p;
    Fbb(:, (u-1)*M+1:u*M) = p;
    max_value_dif = t-t_old;
    t_old = t
    t_all = [t_all, -t+w_k*e_k];
    cvx_status
    [vec_x, val_x] = eigs(p); 
end

    r =  t;
 end
 
 conver = r - r_old;
 r_old = r;
end
  r_all = r_all + r/ITER;
end