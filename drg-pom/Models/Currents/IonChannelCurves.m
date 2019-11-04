% IonChannelCurves
% 24/09/2015

% Script to plot different ion channel gating IV curves

V = -100:0.1:80;

% INaV 1.7

a_m17 = 15.5./(1 + exp(-(V-5)/(12.08)));

b_m17 = 35.2./(1 + exp((V+72.7)/16.7));

m17_inf = a_m17./(a_m17 + b_m17);

a_h17 = 0.38685./(1 + exp((V+122.35)/15.29));

b_h17 = -0.00283 + 2.00283./(1 + exp(-(V+5.5266)/12.70195));

h17_inf = a_h17./(a_h17 + b_h17);

a_s17 = 0.00003 + 0.00092./(1 + exp((V+93.9)/16.6));

b_s17 = 132.05 - 132.05./(1 + exp((V-384.9)/28.5));

s17_inf = a_s17./(a_s17 + b_s17);


% INaV 1.8

a_m18 = 2.85 - 2.839./(1 + exp((V-1.159)/13.95));
b_m18 = 7.6205./(1 + exp((V+46.463)/8.8289));
m18_inf = a_m18./(a_m18 + b_m18);

tau_h18 = 1.218 + 42.043*exp(-((V+38.1).^2)./(2*15.19^2));

h18_inf = 1./(1 + exp((V+32.2)/4));

% INaV 1.8 Han
a_m18hw = 7.35 - 7.35./(1 + exp((V+1.38)/10.9));
b_m18hw = 5.97./(1 + exp((V+56.43)/18.26));

m18hw_inf = a_m18hw./(a_m18hw + b_m18hw);
m18hw_tau = 1./(a_m18hw + b_m18hw);

a_h18hw = 0.011 + 1.39./(1 + exp((V+78.04)/11.32));
b_h18hw = 0.56 - 0.56./(1 + exp((V-21.82)/20.03));

h18hw_inf = a_h18hw./(a_h18hw + b_h18hw);
h18hw_tau = 1./(a_h18hw + b_h18hw);

% INaV 1.9 Huang

a_m19 = 0.751./(1 + exp(-(V+32.26)/13.71));
b_m19 = 5.68./(1 + exp((V+123.71)/13.94));
m19_inf = a_m19./(a_m19 + b_m19);

a_h19 = 0.082./(1 + exp((V+113.69)/17.4));
b_h19 = 0.24./(1 + exp(-(V-10.1)/17.2));
h19_inf = a_h19./(a_h19 + b_h19);

a_s19 = 0.019./(1 + exp((V+154.51)/11.46));
b_s19 = 0.000376./(1 + exp(-(V+60.92)/15.79));
s19_inf = a_s19./(a_s19 + b_s19);

% INaV 1.9 Herzog (19Hz) (from Tigerholm)

a_m19hz = 1.032./(1 + exp(-(V+6.99)/14.87115));
b_m19hz = 5.79./(1 + exp((V+130.4)/22.9));
m19hz_inf = a_m19hz./(a_m19hz + b_m19hz);

a_h19hz = 0.06435./(1 + exp((V+73.26415)/3.71928));
b_h19hz = 0.13496./(1 + exp((V-10.27853)/29.09334));
h19hz_inf = a_h19hz./(a_h19hz + b_h19hz);

a_s19hz = 0.00000016*exp(-V/12);
b_s19hz = 0.0005./(1 + exp(-(V+32)/23));
s19hz_inf = a_s19hz./(a_s19hz + b_s19hz);

% --- Tigerholm Potassium Currents ---

% IKDr from Sheets (delayed rectifier)
nKdrSh_inf = 1./(1 + exp(-(V+45)/15.4));
IKdrSh = nKdrSh_inf.^4;

tau_n = 1000*(0.000688 + 1./(exp((V+75.2)/6.5) + exp(-(V-131.5)./(34.8)))); % Hump function, peaks around -50
tau_n2 = 0.16+0.8*exp(-0.0267*(V+11)); % Exp decay towards +Ve Vm

% IKM from Maigret (non-inactiVating)
nKM_inf = 1./(1 + exp(-(V+30)/6));

% IKA from Sheets with 15 mV shift in steady state Voltage dependence (transient fast-inactiVating)
nKASh_inf = 1./(1 + exp(-(V+5.4+15)/16.4));
hKASh_inf = 1./(1 + exp(V+49.9+15)/4.6);
IKASh = (nKASh_inf.^4).*hKASh_inf;

% fig; hold on;
% plot(V,IKdrSh,V,nKM_inf,'r',V,IKASh,'k')
% IKA

% Ih from KouranoVa
nHC = 1./(1+exp((V + 87.2)/9.7));

% IKNa
Nai = 5:0.1:15; %mM
wKNa = 1./(1+(38.7./Nai).^3.5);

% Leaks
% A bit of a hack? - Makes sure at RMP there's no net current by opposing
% currents at membrane potential?

