function T = dynamic_g1_tt(T, y, x, params, steady_state, it_)
% function T = dynamic_g1_tt(T, y, x, params, steady_state, it_)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T             [#temp variables by 1]     double  vector of temporary terms to be filled by function
%   y             [#dynamic variables by 1]  double  vector of endogenous variables in the order stored
%                                                    in M_.lead_lag_incidence; see the Manual
%   x             [nperiods by M_.exo_nbr]   double  matrix of exogenous variables (in declaration order)
%                                                    for all simulation periods
%   steady_state  [M_.endo_nbr by 1]         double  vector of steady state values
%   params        [M_.param_nbr by 1]        double  vector of parameter values in declaration order
%   it_           scalar                     double  time period for exogenous variables for which
%                                                    to evaluate the model
%
% Output:
%   T           [#temp variables by 1]       double  vector of temporary terms
%

assert(length(T) >= 57);

T = KimModTheBuilder.dynamic_resid_tt(T, y, x, params, steady_state, it_);

T(25) = params(9)/2;
T(26) = (1-T(1))*(y(5)/(1-T(1)))^(1+params(8))+T(1)*(y(6)/T(1))^(1+params(8));
T(27) = T(7)^params(1);
T(28) = y(12)^(1-params(1));
T(29) = T(7)/y(12);
T(30) = (1-T(1))*1/(1-T(1))*getPowerDeriv(y(5)/(1-T(1)),1+params(8),1);
T(31) = getPowerDeriv(T(26),1/(1+params(8)),1);
T(32) = T(30)*T(31);
T(33) = ((1-T(1))*T(2)-y(5)*(1-T(1))*T(32))/((1-T(1))*T(2)*(1-T(1))*T(2));
T(34) = getPowerDeriv(y(5)/((1-T(1))*T(2)),params(8),1);
T(35) = (-(y(6)*T(1)*T(32)))/(T(1)*T(2)*T(1)*T(2));
T(36) = getPowerDeriv(y(6)/(T(1)*T(2)),params(8),1);
T(37) = (-y(6))/(y(1)*y(1));
T(38) = (-(T(25)*T(37)*2*(T(22)-1)));
T(39) = T(1)*1/T(1)*getPowerDeriv(y(6)/T(1),1+params(8),1);
T(40) = T(31)*T(39);
T(41) = (-(y(5)*(1-T(1))*T(40)))/((1-T(1))*T(2)*(1-T(1))*T(2));
T(42) = (T(1)*T(2)-y(6)*T(1)*T(40))/(T(1)*T(2)*T(1)*T(2));
T(43) = 1/y(1);
T(44) = (-(T(25)*2*(T(22)-1)*T(43)));
T(45) = (-y(15))/(y(6)*y(6));
T(46) = 2*y(15)/y(6);
T(47) = T(10)*y(17)*y(18)*T(16)*T(45)*T(46);
T(48) = 1/y(6);
T(49) = T(10)*y(17)*y(18)*T(16)*T(46)*T(48);
T(50) = getPowerDeriv(T(7),params(1),1);
T(51) = T(28)*T(50);
T(52) = 1/y(12);
T(53) = getPowerDeriv(T(29),params(1)-1,1);
T(54) = getPowerDeriv(T(29),params(1),1);
T(55) = getPowerDeriv(y(12),1-params(1),1);
T(56) = T(27)*T(55);
T(57) = (-T(7))/(y(12)*y(12));

end