%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Correct mapping of coefficients consindering the structure given on page
% 542: Supported Modulator Topologies: CIFB Structure Even Order
% Return
% a: Feedback coefficients
% g: Resonator coefficient
% b: Feed-in coefficients from input
% c: Integrator inter-stage coefficients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [a, g, b, c] = mapABCD_CIFB_2nd_LSE(ABCD)
    [Acs, Bcs, Ccs, Dcs] = partitionABCD(ABCD);
    % Tobi's mapping:
    a = [-flip(Bcs(:,2))];
    b = flip([Dcs(2); Bcs(:,1)]);
    c = [Ccs(2); Acs(2,1)];
    g = [-Acs(1,2)];
end


