S1 = [4.7, 2.0];
S2 = [1.6, 1.6];
S3 = [3.0, 1.5];
S4 = [1.8, 1.0];
sigma = 1.6;
S5 = [3.0, sigma; 1.0, sigma; 2.5, sigma; 0.9, sigma];

PE = [];
PE(1,:,:) = [S1; S2; S3; S4; S5(1,:)];
PE(2,:,:) = [S1; S2; S3; S4; S5(2,:)];
PE(3,:,:) = [S1; S2; S3; S4; S5(3,:)];
PE(4,:,:) = [S1; S2; S3; S4; S5(4,:)];

res_dolev = Dolev(PE);
[bound_sensor,bound_hyp] = Sensor_fuse(PE);


%% Dolev
function res = Dolev(PE)
    tau = 1;
    res = [];
    for iter=1:4
        for i = 1:4
            val = PE(i,:,1);
            val = sort(val);
            val = val(tau + 1: end-tau);
            res(i) = mean(val);
        end 
    end
end


%% Mahaney and Schnelder's algorithm
% % true value is in [2.7,2.8], for PE1 delta=4.0, for PE2 delta=2.7 PE3
% % delta = 1.8 PE4 = 2.0;
% delta = 2.2;%[4.0, 2.7, 1.8, 2.0];
% epsilon = 0.1;
% while (max(res) - min(res) > epsilon)
%     iter = iter + 1;
%     for i = 1:4
%         val = PE(i,:,1);
%         val = sort(val);
%         if abs(val(end) - val(1)) > delta
%             if val(end - 1) - val(1) > delta
%                 a = [a 1];
%             else if val(end) - val(2) > delta
%                 a = [a size(val, 1)];
%                 end
%             end
%         end
%         mean = (sum(val) - sum(val(a))) / size(val, 1);
%         val(a) = mean;
%         res(i) = mean(val);
%     end 
% end
% res

%% Sensor fusion algorithm and Hybrid

function [bound,bound_mean] =  Sensor_fuse(PE)
    bound = [];
    bound_mean = [];
    for sr=1:4
        lowB = PE(sr, :,1) - PE(sr,:,2);
        upB = PE(sr,:,1) + PE(sr,:,2);

        B = [lowB upB];
        [B, I] = sort(B);
        flag = zeros(size(B));
        for i=1:size(B, 2)
            if I(i) > size(B, 2) / 2
                continue
            end
            k = i;
            while I(k) ~= I(i)+ 5
                flag(k) = flag(k) + 1;
                k = k + 1;
            end
        end
        idx = 1:size(B, 2);
        idx = idx(flag >= 4);
        weight = flag(flag >= 4);
        ranges = [idx;idx+1];
        mean_r = mean(B(ranges));
        bound_mean(sr) = (mean_r * weight') / sum(weight);
        bound = [bound; [B(idx(1)), B(idx(end) + 1)]];
    end
end