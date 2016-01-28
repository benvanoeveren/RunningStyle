function smoothsignal=lowpassfilterzerolag(lpfreq,sampfreq,butorder,thesignal)
for i = 1: size(thesignal,2)
    [num,den]=butter(butorder,lpfreq./(sampfreq/2));%5 Hz lowpass filter 100 hz filter design, butterworth, 4th order
    smoothsignal(:,i) =(filtfilt(num,den,thesignal(:,i)));
end
end
