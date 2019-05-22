function T = static_g2_tt(T, y, x, params)
% function T = static_g2_tt(T, y, x, params)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T         [#temp variables by 1]  double   vector of temporary terms to be filled by function
%   y         [M_.endo_nbr by 1]      double   vector of endogenous variables in declaration order
%   x         [M_.exo_nbr by 1]       double   vector of exogenous variables in declaration order
%   params    [M_.param_nbr by 1]     double   vector of parameter values in declaration order
%
% Output:
%   T         [#temp variables by 1]  double   vector of temporary terms
%

assert(length(T) >= 42);

T = AnSchoModTheBuilder.static_g1_tt(T, y, x, params);

T(39) = getPowerDeriv(y(5),(-params(4)),2);
T(40) = getPowerDeriv(y(4)/T(2),params(7),2);
T(41) = getPowerDeriv(T(13),1-params(8),2);
T(42) = (y(6))*T(24)*T(38)*T(41);

end