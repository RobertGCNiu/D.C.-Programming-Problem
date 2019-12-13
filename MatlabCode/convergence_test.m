data = r_all_all;
j = 1;
new_data = zeros(600,212);
row = 0;
for i = 1:(length(data)-1)
    if data(i) > (data(i+1)+0.7)
        row = row+1;
        new_data(j,row) = data(i);
        new_data(j,row+1:end) = data(i);
        j = j+1;
        row = 0;
    else
        row = row+1;
        new_data(j,row) = data(i);
    end
    
%     if data(i) <=data(i+1)
%              row = row+1;
%             new_data(j,row) = data(i);
%     elseif data(i) > data(i+1)
%         row = row+1;
%         new_data(j,row) = data(i);
%             j = j+1;
%             row = 1;
%             
%     end
end

hold on;
for i = 0:5
    data_rho = 0;
for j = i +1 : 6 : 120
    data_rho = new_data(j,:) + data_rho;
end
plot(data_rho/20)
end
legend('-10 dB', '-5 dB', '0 dB', '5 dB', '10 dB', '15 dB')