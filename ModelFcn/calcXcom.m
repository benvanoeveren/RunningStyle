function  [XCOM] = calcXcom(COM,fs)
l = mean(COM(:,3)); %pendulum length ==? height of COM?
w0 = sqrt(9.81/l);
Vcom = [nan(1,3);diff(COM).*fs];
XCOM = COM+ (Vcom./w0);
XCOM(:,3) = zeros(size(XCOM(:,3)));
end