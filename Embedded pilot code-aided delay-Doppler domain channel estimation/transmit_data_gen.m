function [X, data, trans_info_bit] = transmit_data_gen(N, M, N_bits_perfram, N_syms_perfram, M_bits, M_mod, data_grid)
        trans_info_bit = randi([0,1],N_bits_perfram,1);
        data=qammod(reshape(trans_info_bit,M_bits,N_syms_perfram), M_mod,'gray','InputType','bit','UnitAveragePower',true);
        X = Generate_2D_data_grid(N,M,data,data_grid);
end