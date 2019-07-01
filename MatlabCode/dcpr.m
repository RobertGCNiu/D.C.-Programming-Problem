clear;clc;

% H = [0.4310, 0.02, 0.2605, 0.39; 
%     0.02, 0.3018, 0.08, 0.54; 
%     0.0129, 0.05, 0.4266, 0.1007; 
%     0.11, 0.31, 0.0099, 0.0634];
 M = 4
 r_all_all=[];
 r_all = [];
%w = [1.0/6, 1.0/6, 1.0/3, 1.0/3];
p_max = [0.7, 0.8, 0.9, 1.0];
ri = 0;
ITER = 1;
r_all = 0;
r_p_opt_all = [];
R_opt_all = [];
for iter = 1: ITER  

[H_c,a_TX,a_RX]=generate_channels(M,8,8,4,4,1); 
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
    


for rho = 25:5:30
sigma =  1/db2pow(rho);

%%%%%%%%%%%%%%%%%%%%%%%%zero-forcing
    p_opt =  pinv(H);
for u = 1:M
p_opt(:,u) = p_opt(:,u)/norm(Frf*p_opt(:,u));
end
R_opt_i = det(eye(M)+1/sigma*(H*(p_opt*p_opt')*H'));
r_p_opt = log2(R_opt_i);

%R_opt_all =[R_opt_all,log2(R_opt_i(1,1)),log2(R_opt_i(2,2))]; 

r_p_opt_all =[r_p_opt_all r_p_opt];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%p_initial_v = [0.9120; 0.9514; 0.3460; 0.2902];
p_initial_v = ones(M,1);
p_initial = p_initial_v * p_initial_v';  
t_old = 0;
r_old = 0;
Fbb = zeros(M,M*M);
iter = 0;
for u = 1:M
    Fbb(:, (u-1)*M+1:u*M) =  p_initial;
end

conver = 100;
    r = 0;
while(abs(conver)>0.1)

 for u = 1:4
p_initial = Fbb(:, (u-1)*M+1:u*M);
max_value_dif = 100;
while(abs(max_value_dif)>0.1)    

cvx_begin sdp quiet
cvx_precision best
variable p(M,M) hermitian
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
    gra_g = g(H,Fbb,u,sigma,M);   % gradient of g
    t= f_p - g_p - real(trace(gra_g'*(p-p_initial)));
    maximize(t)
    subject to
  %  p>=0;
%     for i = 1:M
%         if i == u
%         all = all + norm(Frf*p*Frf');
%         else 
%             all = all + norm(Frf*Fbb(:, (i-1)*M+1:i*M)*Frf');
%         end
%     end
%   all<=M;
    norm(Frf*p*Frf')<=1;
%             for i = 1:M
%                 He = 0;
%                 for j = 1:M
%                    if j ~= i
%                        if j ==u
%                            He = He + H(i,:) * p * H(i,:)'; 
%                        else
%                         He = He + H(i,:) * Fbb(:, (j-1)*M+1:j*M) * H(i,:)';
%                        end
%                    end
%                   % abs(p(i,j))<=2;           
%                 end
%                %norm(Frf*p(:,i))<=1;
%                  H(i,:) * Fbb(:, (i-1)*M+1:i*M) *H(i,:)' + (1-2^ri)*(He + sigma) >= 0;
%             end 
cvx_end  
    p_initial = p;
    Fbb(:, (u-1)*M+1:u*M) = p;
    max_value_dif = t-t_old;
    t_old = t
    cvx_status
      r_all_all = [r_all_all, t];
end

 %   Fbb(:, (u-1)*M+1:u*M) = p_initial;
    r =  t;
 end
 
 conver = r - r_old;
 r_old = r;
end
  r_all = [r_all,r];
end
end