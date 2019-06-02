function [g_gradient] = g_block(H,Fbb,u,sigma)



M = 4;
K=2;
m_k = M/K;
g_gradient = 0;



for i = 1 : M
    H_d = 0;
    for k = 1:K
        for j = 1+(k-1)*m_k : k*m_k
            if j ~= i
                H_d = H_d + real(H(i,1+(k-1)*m_k : k*m_k) * Fbb(:, (j-1)*m_k+1:j*m_k) * H(i,1+(k-1)*m_k : k*m_k)');
            end
        end
    end
    if i~=u
        if i<=m_k
            k = 1;
        elseif i>m_k
            k = 2;
        end
        g_gradient = g_gradient + real(H(i,1+(k-1)*m_k : k*m_k)' * H(i,1+(k-1)*m_k : k*m_k)/log(2)/(sigma+H_d));
    end
    
end



end

