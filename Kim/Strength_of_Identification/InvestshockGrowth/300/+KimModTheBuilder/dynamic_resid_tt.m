function T = dynamic_resid_tt(T, y, x, params, steady_state, it_)
% function T = dynamic_resid_tt(T, y, x, params, steady_state, it_)
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

assert(length(T) >= 21);

T(1) = (steady_state(3))/(steady_state(1));
T(2) = ((1-T(1))*(y(6)/(1-T(1)))^(1+params(8))+T(1)*(y(7)/T(1))^(1+params(8)))^(1/(1+params(8)));
T(3) = 0;
T(4) = 0;
T(5) = 0;
T(6) = 0;
T(7) = y(2);
T(8) = (y(6)/((1-T(1))*T(2)))^params(8);
T(9) = (y(6)-T(3))^(-1)*(1-T(4));
T(10) = 1/(1+params(2)/400);
T(11) = (y(14)-T(5))^(-1)*(-T(6));
T(12) = (y(7)/(T(1)*T(2)))^params(8);
T(13) = params(9)/2*(y(7)/y(1)-1)^2;
T(14) = params(9)*(y(7)/y(1)-1);
T(15) = y(19);
T(16) = params(9)*(y(15)/y(7)-1);
T(17) = T(7)^params(1);
T(18) = params(1)*T(7)^(params(1)-1);
T(19) = y(7)/y(1);
T(20) = (y(15)/y(7))^2;
T(21) = T(10)*y(17)*y(18)*T(15)*T(20);

end
