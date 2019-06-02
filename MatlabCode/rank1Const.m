clear;
data= load('data.mat');
a_TX = data.a_TX;
a_RX =data.a_RX;
Fbb = data.Fbb;
H_c=data.H_c;
M = 4;
sigma = 0.001;

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

Fbb_new = zeros(M,M*M);
for u = 1:M
x_0 = Fbb(:,(u-1)*M+1:u*M);
A = H(u,:);


cvx_begin sdp quiet 
variable x(M,M) symmetric
f = norm(A*x*A'-A*x_0*A');
minimize(f);
subject to
    norm(Frf*x*Frf') <=1;
        x>=0;
cvx_end
[vec_x, val_x] = eigs(x);

t = 0.9;
w_k = 1;
epsi = 100;
f_old = 10; 
k = 1;
e_k_old =val_x(3,3);
while(epsi>0.001)
k = k+1;
w_k = w_k*t;

cvx_begin sdp quiet
variable x_k(M,M) symmetric
variable e_k
f = norm(A*x_k*A'-A*x_0*A') + w_k*e_k; 
minimize(f);
subject to
        norm(Frf*x_k*Frf') <=1;
        x_k>=0;
        e_k*eye(3)- vec_x(:,2:end)'*x_k*vec_x(:,2:end)>=0;
        e_k<=e_k_old;
cvx_end
e_k_old = e_k;
[vec_x, val_x] = eigs(x_k); 
epsi = f_old -f;
f_old = f;
end
Fbb_new(:,(u-1)*M+1:u*M) = x_k;
end