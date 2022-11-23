%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Correct chosen frequency to the simulation time vector
% info: necessary for correct fft and calculation of SNR
% Return
% inSigFreq: frequency which is contained in the time vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function inSigFreq = correctFrequToSim(f, simTime)
    factor = round(f*simTime);
    inSigFreq = factor*1/simTime;
end

