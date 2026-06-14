function [peak, average, papr, tx_power] = papr_calculation(s, s_sample, M, N)
        peak = max(abs(s).^2);
        average = mean(abs(s).^2 + abs(s_sample).^2);
        papr = 10*log10(peak./average);
        tx_power = sum(abs(s).^2)/(M*N);
end