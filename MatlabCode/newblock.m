clear;clc;

% H = [0.4310, 0.02, 0.2605, 0.39; 
%     0.02, 0.3018, 0.08, 0.54; 
%     0.0129, 0.05, 0.4266, 0.1007; 
%     0.11, 0.31, 0.0099, 0.0634];
 M = 8;
 K = 2;
 %sigma_0 = 0.001;
 m_k = M/K
 R_all=[];
 R_opt_all=[];
w = ones(M,1);
p_max = [0.7, 0.8, 0.9, 1.0];
ri = 0;
ITER = 1;
r_all = 0;
t_all = [];
%for iter = 1: ITER  
% data = load('data.mat');
% a_TX = data.a_TX;
% a_RX =data.a_RX;
% H_c = data.H_c;
%Fbb_v = data.Fbb;
at_num =8;
rt_num =2;
R_dc = [];
for iter_realization = 1:ITER
%     data = load('block_2.mat');
%     H_c = data.H_c;
%     a_TX = data.a_TX;
%     a_RX = data.a_RX;
  [H_c,a_TX,a_RX]=generate_channels(M,at_num,at_num,rt_num,rt_num,1); 
  for u=1:1:M
        Frf(:,u)=a_TX(:,u);
        Wrf(:,u)=a_RX(:,u);
  end      
    r_p_opt_all = [];
    r_all_all = [];
    % Constructin the effective channels
    for u=1:1:M
        Channel=zeros(rt_num^2,at_num^2);
        Channel(:,:)= H_c(u,:,:);
        H(u,:)=Wrf(:,u)'*Channel*Frf ;    % Effective channels
    end
 for SNR_dB = -15:5:30
     

   %sigma =  1/db2pow(rho);
   

   

%p_initial_v = [0.9120; 0.9514; 0.3460; 0.2902];
%p_initial_v =  1j * rand(M,1);
%p_initial = p_initial_v * p_initial_v';  
%p_initial = p_initial/sqrt(p_initial_v'*(Frf'*Frf)*p_initial_v);
t_old = 0;
r_old = 0;
%Fbb_i = zeros(M,M*M);
iter = 0;
p_opt = [];
for k = 1:K
    p_opt_k =  pinv(H((k-1)*m_k+1:k*m_k,(k-1)*m_k+1:k*m_k));
    p_opt = blkdiag(p_opt,p_opt_k);
end

   for u = 1:M
        p_opt(:,u) = p_opt(:,u)/norm(Frf*p_opt(:,u));
   end
   
   rho=db2pow(SNR_dB);
  sigma = 1/rho;
     
   Heff_p = abs(H * p_opt).^2;
   [p_dc,~] =  dcpower(Heff_p, M, M,sigma);
   
R_opt_i =  RGH(H,p_opt, rho);
 [R_each_user, ~] = PGH(Heff_p,p_dc,sigma);
 R_dc = [R_dc R_each_user];
%R_opt_all =[R_opt_all,log2(R_opt_i(1,1)),log2(R_opt_i(2,2))]; 

r_p_opt_all =[r_p_opt_all R_opt_i];

% if r_p_opt<=13
%     continue
% end

a = ones(m_k,1);
b = rand(m_k,1);
p_initial_v = a + 1j*b;  
p_initial_v = p_initial_v/sqrt(p_initial_v'*(Frf(:,1:m_k)'*Frf(:,1:m_k))*p_initial_v);
 for u = 1:M
     %Fbb_t(:, (u-1)*M+1:u*M) = p_initial_v*p_initial_v'/1000;
     %Fbb_i(:, (u-1)*M+1:u*M)  =p_opt(:,u)*p_opt(:,u)';% a*a'+1j*(b*b');
     Fbb(:, (u-1)*m_k+1:u*m_k)  =(a*a'+1j*(b*b'))/10; 
 end
 

conver = 100;

while(abs(conver)>0.1)

    iter = iter + 1;

 for u = 1:M
p_initial = Fbb(:, (u-1)*m_k+1:u*m_k);
max_value_dif = 100;
while(abs(max_value_dif)>0.1)    

cvx_begin sdp quiet 
cvx_precision best
variable p(m_k,m_k) hermitian
 f_p = 0;
 g_p = 0;
 for i = 1 : M
     H_f = 0;
     H_g = 0;
     for k = 1:K
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
     end
     f_p  = f_p +  log(sigma + H_f)/log(2);   % f(p)
     g_p  = g_p +  log(sigma + H_g)/log(2);  % g(p)
 end
    gra_g = g_block(H, Fbb,u,sigma,M,K);   % gradient of g
    all=0;
    t= f_p - g_p - real(trace(gra_g'*(p-p_initial)));
    maximize(t)
    subject to
   p==semidefinite(m_k);      
   for k = 1:K
    for i =(k-1)*m_k+1:k*m_k
        if u<=m_k
        if i == u
        all = all + norm(Frf(:,1:m_k)*p*Frf(:,1:m_k)');
        else 
            all = all + norm(Frf(:,1:m_k)*Fbb(:, (i-1)*m_k+1:i*m_k)*Frf(:,1:m_k)');
        end
        end
        if u>m_k
        if i == u
        all = all + norm(Frf(:,m_k+1:end)*p*Frf(:,m_k+1:end)');
        else 
            all = all + norm(Frf(:,m_k+1:end)*Fbb(:, (i-1)*m_k+1:i*m_k)*Frf(:,m_k+1:end)');
        end
        end
    end
   end
    all<=M;
 %norm(Frf*p*Frf')<=1;
%     if u <=m_k
%   norm(Frf(:,1:m_k)*p*Frf(:,1:m_k)')<=1;
%    elseif u>m_k
%    norm(Frf(:,m_k+1:end)*p*Frf(:,m_k+1:end)')<=1;
%    end
 if t_old>=2*ri
            for i = u
                He = 0;
                for k = 1:K
                for j = (k-1)*m_k+1:k*m_k
                   if j ~= i
                       if j ==u
                           He = He + real(H(i,:) * p * H(i,:)'); 
                       else
                        He = He + real(H(i,1+(k-1)*m_k : k*m_k) * Fbb(:, (j-1)*m_k+1:j*m_k) * H(i,1+(k-1)*m_k : k*m_k)');
                       end
                   end
                 % abs(p(i,j))<=2;           
                end
                if i==u
                    real(H(i,1+(k-1)*m_k : k*m_k) * p * H(i,1+(k-1)*m_k : k*m_k)') + (1-2^ri)*(He + sigma) >= 0;
                elseif i ~= u
                    real(H(i,1+(k-1)*m_k : k*m_k) * Fbb(:, (i-1)*m_k+1:i*m_k) *H(i,1+(k-1)*m_k : k*m_k)') + (1-2^ri)*(He + sigma) >= 0;
                end
                end 
            end
end
cvx_end  

    p_initial = p;
    Fbb(:, (u-1)*m_k+1:u*m_k) = p;
    max_value_dif = t-t_old;
    t_old = t
    t_all = [t_all t];
    cvx_status;
     if t > R_opt_i
         max_value_dif = 0;
         conver=0;
     end
end

u
 %   Fbb(:, (u-1)*M+1:u*M) = p_initial;
    r =  t;
 end
 

 
 conver = r - r_old;
 r_old = r;
       if t > R_opt_i
         max_value_dif = 0;
         conver=0;
       end
end

%%%%%%%%%===CDF of capacity
%             for i = 1:M
%                 He = 0;
%                 for j = 1:M 
%                        if j ~= i
%                         He = He + real(H(i,:) * Fbb(:, (j-1)*M+1:j*M) * H(i,:)');
%                        end
%                  end
%                  % abs(p(i,j))<=2;           
%                 R_i = log2(1+real(H(i,:) * Fbb(:, (i-1)*M+1:i*M) *H(i,:)') / (He + sigma));
%                 R_all = [R_all, R_i];  
%             end 
%%%%%%%

 % r_all = r_all + r/ITER;
  r_all_all = [r_all_all r];
%end

 end
end
