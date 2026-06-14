function [Y_tilda, Y, Y_tf] = otfs_demodulation(r1, N, M, Fn, Fm)
        Y_tilda=reshape(r1,M,N);
        Y = Y_tilda*Fn;% delay-Doppler received symbols
        Y_tf = Fm*Y_tilda;
end