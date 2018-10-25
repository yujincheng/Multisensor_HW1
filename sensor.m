function main
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
end

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