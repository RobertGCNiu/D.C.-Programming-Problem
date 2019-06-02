function [g_gradient] = g(H,Fbb,u,sigma,M,w)



g_gradient = 0;


for i = 1 : M
    H_d = 0;
    for j = 1: M
        if j ~= i
            H_d = H_d + real(H(i,:) * Fbb(:, (j-1)*M+1:j*M) * H(i,:)');
        end
    end
    if i~=u
    g_gradient = g_gradient + w(i)*real(H(i,:)' * H(i,:))/log(2)/(sigma+H_d);
    end
end

end

