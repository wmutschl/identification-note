function T = dynamic_g2_tt(T, y, x, params, steady_state, it_)
% function T = dynamic_g2_tt(T, y, x, params, steady_state, it_)
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

assert(length(T) >= 78);

T = AnSchoModTheBuilder.dynamic_g1_tt(T, y, x, params, steady_state, it_);

T(71) = y(16)*getPowerDeriv(y(11),(-params(4)),2);
T(72) = T(6)*getPowerDeriv(y(18),(-params(4)),2);
T(73) = getPowerDeriv(T(7)/y(16),(-1),2);
T(74) = (-(T(37)*getPowerDeriv(y(3),params(17),2)));
T(75) = getPowerDeriv(T(35),params(8),2);
T(76) = getPowerDeriv(T(13),1-params(9),2);
T(77) = (T(36)*T(65)+T(32)*T(39)*T(67))*T(76);
T(78) = getPowerDeriv(y(10)/T(2),params(7),2);

end