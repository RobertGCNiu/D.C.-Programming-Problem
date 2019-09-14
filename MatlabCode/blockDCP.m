%function [t_all] = blockDCP(Frf, H, M, rho, Fbb_initial)

% H = [0.4310, 0.02, 0.2605, 0.39; 
%     0.02, 0.3018, 0.08, 0.54; 
%     0.0129, 0.05, 0.4266, 0.1007; 
%     0.11, 0.31, 0.0099, 0.0634];
clear;clc
 M = 4;
 rho = db2pow(10);
 K = 2;
 m_k = M/K;
 sigma = 1/rho;
%w = [1.0/6, 1.0/6, 1.0/3, 1.0/3];
ri = 0;
ITER = 1;
r_all = 0;

 [H_c,a_TX,a_RX]=generate_channels(M,4,4,2,2,1); 
% % data = load('testdata.mat');
% % a_TX = data.a_TX;
% % a_RX =data.a_RX;
% % H_c = data.H_c;
   for u=1:1:M
         Frf(:,u)=a_TX(:,u);
         Wrf(:,u)=a_RX(:,u);
   end      
%     
%     % Constructin the effective channels
     for u=1:1:M
         Channel=zeros(4,16);
         Channel(:,:)= H_c(u,:,:);
        H(u,:)=Wrf(:,u)'*Channel*Frf ;    % Effective channels
     end
    

%p_initial_v = [0.9120; 0.9514; 0.3460; 0.2902];
p_initial_v = ones(m_k,1);
b = rand(m_k,1);
p_initial = p_initial_v * p_initial_v' +1j * (b*b');  
t_old = 0;
r_old = 0;
%Fbb = zeros(m_k,m_k*M);
iter = 0;
for k = 1:M
     Fbb(:, (k-1)*m_k+1:k*m_k) =p_initial;
end


conver = 100;


while(abs(conver)>0.1)
    iter = iter + 1;
 for u = 1:M
p_initial = Fbb(:, (u-1)*m_k+1:u*m_k);
max_value_dif = 100;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while(abs(max_value_dif)>0.1)    

cvx_begin sdp quiet 
cvx_precision best
variable p(m_k,m_k) hermitian
f_p = 0;
 g_p = 0;
 for i = 1 : M
         H_f = 0;
         H_g = 0;
         k=1;
   % for k = 1:K
         for j = (k-1)*m_k+1:k*m_k
             if j == u
                 H_f = H_f +  H(i,1+(k-1)*m_k : k*m_k) * p * H(i,1+(k-1)*m_k : k*m_k)';
             else
                 H_f = H_f +  real(H(i,1+(k-1)*m_k : k*m_k) * Fbb(:, (j-1)*m_k+1:j*m_k) * H(i,1+(k-1)*m_k : k*m_k)');
             end
             if j ~= i
                 H_g = H_g + real(H(i,1+(k-1)*m_k : k*m_k) * Fbb(:, (j-1)*m_k+1:j*m_k)  * H(i,1+(k-1)*m_k : k*m_k)');
             end
         end
    % end
     f_p  = f_p +  log(sigma + H_f)/log(2);   % f(p)  
     g_p  = g_p +  log(sigma + H_g)/log(2);  % g(p)   
 end
    gra_g = g_block(H, Fbb,u,sigma,M,K);   % gradient of g
    t= f_p - g_p- real(trace(gra_g'*(p-p_initial)));
    maximize(t)
    subject to
    p>=0;
    if u == 1 || u == 2 
 %norm(Frf*p*Frf')<=1;
   norm(Frf(:,1:m_k)*p*Frf(:,1:m_k)')<=1;
    elseif u ==4 || u == 3
    norm(Frf(:,m_k+1:end)*p*Frf(:,m_k+1:end)')<=1;
    end
%             for i = 1:M
%                 He = 0;
%                 for k = 1:K
%                 for j = 1+(k-1)*m_k : k*m_k
%                    if j ~= i
%                        if j ==u
%                            He = He + H(i,1+(k-1)*m_k : k*m_k) * p_0 * H(i,1+(k-1)*m_k : k*m_k)'; 
%                        else
%                         He = He + H(i,1+(k-1)*m_k : k*m_k) * Fbb(:, (j-1)*m_k+1:j*m_k) * H(i,1+(k-1)*m_k : k*m_k)';
%                        end
%                    end
%                   % abs(p(i,j))<=2;     
%                 end           
%                 end
%                 H(i,1+(k-1)*m_k : k*m_k) * p_0 *H(i,1+(k-1)*m_k : k*m_k)' + (1-2^ri)*(He + sigma) >= 0;
%             end 
cvx_end

p_initial = p;
Fbb(:, (u-1)*m_k+1:u*m_k) = p;
max_value_dif = t-t_old;
t_old = t
cvx_status;
 
end
%Fbb(:, (u-1)*m_k+1:u*m_k) = p_0;
%[vec_x, val_x] = eigs(p_0);    
% t_s = 1.1;
% w_k = 0.2;
% p_initial = p_0;
% e_k_old =val_x(2,2); 


% max_value_dif = 100
u
% while(abs(max_value_dif)>0.1)    
%  w_k = w_k*t_s; 
% cvx_begin sdp quiet
% variable p(m_k,m_k) hermitian
% variable e_k
%  f_p = 0;
%  g_p = 0;
%  for i = 1 : M
%          H_f = 0;
%          H_g = 0;
%     for k = 1:K
%          for j = (k-1)*m_k+1:k*m_k
%              if j == u
%                  H_f = H_f +  H(i,1+(k-1)*m_k : k*m_k) * p * H(i,1+(k-1)*m_k : k*m_k)';
%              else
%                  H_f = H_f +  real(H(i,1+(k-1)*m_k : k*m_k) * Fbb(:, (j-1)*m_k+1:j*m_k) * H(i,1+(k-1)*m_k : k*m_k)');
%              end
%              if j ~= i
%                  H_g = H_g + real(H(i,1+(k-1)*m_k : k*m_k) * Fbb(:, (j-1)*m_k+1:j*m_k)  * H(i,1+(k-1)*m_k : k*m_k)');
%              end
%          end
%      end
%      f_p  = f_p +  log(sigma + H_f)/log(2);   % f(p)
%      g_p  = g_p +  log(sigma + H_g)/log(2);  % g(p)   
%  end
%     gra_g = g_block(H, Fbb,u,sigma);   % gradient of g
%     t= -f_p + g_p + real(trace(gra_g'*(p-p_initial)))+  w_k*e_k;
%     minimize(t)
%     subject to
%     p>=0;
%     if u == 1 || u == 2 || u == 3
%     norm(Frf(:,1:m_k)*p*Frf(:,1:m_k)')<=1;
%     elseif u ==4 || u == 5|| u == 6
%     norm(Frf(:,m_k+1:end)*p*Frf(:,m_k+1:end)')<=1;
%     end
%        e_k*eye(m_k-1)- vec_x(:,2:end)'*p*vec_x(:,2:end)>=0;
%         e_k<=e_k_old;
%             for i = 1:M
%                     He = 0;
%                 for k = 1:K
%                 for j = 1+(k-1)*m_k : k*m_k
%                    if j ~= i
%                        if j ==u
%                            He = He + H(i,1+(k-1)*m_k : k*m_k) * p * H(i,1+(k-1)*m_k : k*m_k)'; 
%                        else
%                         He = He + H(i,1+(k-1)*m_k : k*m_k) * Fbb(:, (j-1)*m_k+1:j*m_k) * H(i,1+(k-1)*m_k : k*m_k)';
%                        end
%                    end
%                   % abs(p(i,j))<=2;     
%                 end                    
%                 end
%                  H(i,1+(k-1)*m_k : k*m_k) * p *H(i,1+(k-1)*m_k : k*m_k)' + (1-2^ri)*(He + sigma) >= 0;
%      %          norm(Frf*p(:,i))<=1;
%            end 
% cvx_end  
%     e_k_old = e_k;
%     p_initial = p;
%     Fbb(:, (u-1)*m_k+1:u*m_k) = p;
%     max_value_dif = t-t_old;
%     t_old = t;
%     t_all = -t+w_k*e_k;
%     cvx_status;
% end

 %   Fbb(:, (u-1)*M+1:u*M) = p_initial;
    r =  t;
 end
 
 conver = r - r_old;
 r_old = r;
end
t_all = r
%end