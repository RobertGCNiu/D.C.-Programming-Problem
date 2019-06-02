function [t_all] = BlockDCP_by_blkdiag(Frf, H, M, rho, Fbb_initial)
% H = [0.4310, 0.02, 0.2605, 0.39; 
%     0.02, 0.3018, 0.08, 0.54; 
%     0.0129, 0.05, 0.4266, 0.1007; 
%     0.11, 0.31, 0.0099, 0.0634];
%w = [1.0/6, 1.0/6, 1.0/3, 1.0/3];
K=2;
m_k = M/K;
ri = 0;
ITER = 1;
r_all = 0;
t_all = [];
 sigma = 1/rho;

t_old = 0;
r_old = 0;
iter = 0;

 for u = 1:M
     %Fbb_t(:, (u-1)*M+1:u*M) = p_initial_v*p_initial_v'/1000;
     Fbb_i(:, (u-1)*M+1:u*M)  = Fbb_initial(u,:)'*Fbb_initial(u,:);
 end
 
 Fbb =Fbb_i;

conver = 100;

while(abs(conver)>0.01)
    r = 0;
    iter = iter + 1;
 for u = 1:4
p_initial = Fbb(:, (u-1)*M+1:u*M);
max_value_dif = 100;
while(abs(max_value_dif)>0.01)    

cvx_begin sdp quiet
variable p_blk(m_k,m_k) hermitian
if u == 1|| u==2
    p = blkdiag(p_blk, zeros(m_k,m_k));
elseif u == 3|| u==4
      p = blkdiag(zeros(m_k,m_k),p_blk);
end
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
    t= f_p - g_p - real(trace(gra_g'*(p-p_initial)));
    maximize(t)
    subject to
    p>=0; 
   norm(Frf*p*Frf')<=1;
            for i = 1:M
                He = 0;
                for j = 1:M 
                   if j ~= i
                       if j ==u
                           He = He + real(H(i,:) * p * H(i,:)'); 
                       else
                        He = He + real(H(i,:) * Fbb(:, (j-1)*M+1:j*M) * H(i,:)');
                       end
                   end
                 % abs(p(i,j))<=2;           
                end
                 real(H(i,:) * p *H(i,:)') + (1-2^ri)*(He + sigma) >= 0;
            end 
cvx_end  
    %if t>t_old
    p_initial = p;
    Fbb(:, (u-1)*M+1:u*M) = p;
    max_value_dif = t-t_old;
    t_old = t
    t_all = [t_all t];
    cvx_status
    %end
end

 %   Fbb(:, (u-1)*M+1:u*M) = p_initial;
    r =  t;
 end
 
 conver = r - r_old;
 r_old = r;
end
  r_all = r_all + r/ITER;
%end