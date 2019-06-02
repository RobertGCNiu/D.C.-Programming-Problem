function [Fbb_dc]= rank_dcp(G_cl, rho,sigma,K)
Num_users = size(G_cl,1);
M = Num_users/K;

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
    gra_g = g(H, Fbb,u,sigma);   % gradient of g
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
    gra_g = g(H, Fbb,u,sigma);   % gradient of g
    t= -f_p + g_p + real(trace(gra_g'*(p-p_initial))) +  w_k*e_k;
    minimize(t)
    subject to
        p>=0;
        e_k*eye(M-1)- vec_x(:,2:end)'*p*vec_x(:,2:end)>=0;
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
