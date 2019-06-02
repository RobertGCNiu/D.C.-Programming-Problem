clear;clc;

% H = [0.4310, 0.02, 0.2605, 0.39; 
%     0.02, 0.3018, 0.08, 0.54; 
%     0.0129, 0.05, 0.4266, 0.1007; 
%     0.11, 0.31, 0.0099, 0.0634];
 M = 3
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
at_num = 6;
rt_num = 3;
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
for ITER = 1
for rho = -30:10:10
   sigma =  1/db2pow(rho);
   
   

%p_initial_v = [0.9120; 0.9514; 0.3460; 0.2902];
%p_initial_v =  1j * rand(M,1);
%p_initial = p_initial_v * p_initial_v';  
%p_initial = p_initial/sqrt(p_initial_v'*(Frf'*Frf)*p_initial_v);
t_old = 0;
r_old = 0;
%Fbb_i = zeros(M,M*M);
iter = 0;
p_opt =  pinv(H);

for u = 1:M
p_opt(:,u) = p_opt(:,u)/norm(Frf*p_opt(:,u));
end

r_p_opt = log2(det(eye(M)+diag(w)*1/sigma*(H*(p_opt*p_opt')*H')));
r_p_opt_all =[r_p_opt_all r_p_opt];
a = rand(M,1);
b = rand(M,1);
 for u = 1:M
     %Fbb_t(:, (u-1)*M+1:u*M) = p_initial_v*p_initial_v'/1000;
     %Fbb_i(:, (u-1)*M+1:u*M)  =p_opt(:,u)*p_opt(:,u)';% a*a'+1j*(b*b');
     Fbb_i(:, (u-1)*M+1:u*M)  = a*a'+1j*(b*b');
 end
 
 Fbb =Fbb_i;%-Fbb_t;

conver = 100;

while(abs(conver)>0.01)
    
    iter = iter + 1;

 for u = 1:M
p_initial = Fbb(:, (u-1)*M+1:u*M);
max_value_dif = 100;
while(abs(max_value_dif)>0.01)    

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
     f_p  = f_p +  w(i)*log(sigma + H_f)/log(2);   % f(p)
     g_p  = g_p +  w(i)*log(sigma + H_g)/log(2);  % g(p)
 end
    gra_g = g(H, Fbb,u,sigma,M,w);   % gradient of g
    all=0;
    t= f_p - g_p - real(trace(gra_g'*(p-p_initial)));
    maximize(t)
    subject to
    p>=0;     
%     for i = 1:M
%         if i == u
%         all = all + norm(Frf*p*Frf');
%         else 
%             all = all + norm(Frf*Fbb(:, (i-1)*M+1:i*M)*Frf');
%         end
%     end
%     all<=M;
  norm(Frf*p*Frf')<=1;
%             for i = 1:M
%                 He = 0;
%                 for j = 1:M 
%                    if j ~= i
%                        if j ==u
%                            He = He + real(H(i,:) * p * H(i,:)'); 
%                        else
%                         He = He + real(H(i,:) * Fbb(:, (j-1)*M+1:j*M) * H(i,:)');
%                        end
%                    end
%                  % abs(p(i,j))<=2;           
%                 end
%                  real(H(i,:) * Fbb(:, (i-1)*M+1:i*M) *H(i,:)') + (1-2^ri)*(He + sigma) >= 0;
%             end 
cvx_end  

    p_initial = p;
    Fbb(:, (u-1)*M+1:u*M) = p;
    max_value_dif = t-t_old;
    t_old = t
    t_all = [t_all t];
    cvx_status

end

 %   Fbb(:, (u-1)*M+1:u*M) = p_initial;
    r =  t;
 end
 
 conver = r - r_old;
 r_old = r;
end
 % r_all = r_all + r/ITER;
  r_all_all = [r_all_all r];
%end

end
end