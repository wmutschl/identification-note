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

assert(length(T) >= 68);

T = AnSchoModTheBuilder.dynamic_resid_tt(T, y, x, params, steady_state, it_);

T(28) = (1-params(5))^(1/params(4));
T(29) = (steady_state(6))*(y(13)/T(1))^params(6);
T(30) = (y(10)/T(2))^params(7);
T(31) = T(29)*T(30);
T(32) = y(10)/y(1);
T(33) = T(2)/T(3);
T(34) = T(32)/T(33);
T(35) = T(34)^params(8);
T(36) = (steady_state(7))^(1-params(17));
T(37) = (-y(10))/(y(1)*y(1))/T(33);
T(38) = getPowerDeriv(T(34),params(8),1);
T(39) = getPowerDeriv(T(13),1-params(9),1);
T(40) = 1/(steady_state(4));
T(41) = (-(T(9)/T(7)*y(17)))/(y(10)*y(10));
T(42) = 1/T(2);
T(43) = getPowerDeriv(y(10)/T(2),params(7),1);
T(44) = T(29)*T(42)*T(43);
T(45) = 1/y(1)/T(33);
T(46) = T(35)*T(44)+T(31)*T(38)*T(45);
T(47) = T(9)/T(7)/y(10);
T(48) = getPowerDeriv(y(11),(-params(4)),1);
T(49) = y(16)*T(48);
T(50) = getPowerDeriv(T(7),(-1),1);
T(51) = y(17)*(-(T(9)*T(49)))/(T(7)*T(7))/y(10);
T(52) = getPowerDeriv(y(18),(-params(4)),1);
T(53) = T(6)*T(52);
T(54) = y(17)*T(53)/T(7)/y(10);
T(55) = getPowerDeriv(y(2),params(9),1);
T(56) = (-(T(36)*getPowerDeriv(y(3),params(17),1)));
T(57) = T(36)*getPowerDeriv(y(13),params(17),1);
T(58) = (steady_state(6))*1/T(1)*getPowerDeriv(y(13)/T(1),params(6),1);
T(59) = T(30)*T(58);
T(60) = (-(T(28)*T(2)))/(T(3)*T(3));
T(61) = (-(T(32)*T(60)))/(T(33)*T(33));
T(62) = (-(T(28)*y(10)))/(T(2)*T(2));
T(63) = T(29)*T(43)*T(62);
T(64) = T(28)/T(3);
T(65) = (-(T(32)*T(64)))/(T(33)*T(33));
T(66) = y(17)*(-(T(4)*T(9)))/(T(7)*T(7))/y(10);
T(67) = y(17)*T(5)/T(7)/y(10);
T(68) = params(12)/100*T(22);

end