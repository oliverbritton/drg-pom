// INaK

// Parameters
double gpump = ;


double IK_pump = (gpump/(1+(1/Ksp)^2)) * ( ( 1.62/(1+(6.7/(Nain+8)^3)) ) + (1/(1+(67.6/(Nain+8))^3)) );

double INa_pump = -1.5*IK_Pump;

double I_pump = IK_pump + INa_pump;