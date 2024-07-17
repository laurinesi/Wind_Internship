function [Uopt] = FrmIni (func, parameters, arg)
        nstoch = size(parameters,1);
        Uopt = zeros(nstoch,1);
%         Uopt = ones(nstoch,1);
