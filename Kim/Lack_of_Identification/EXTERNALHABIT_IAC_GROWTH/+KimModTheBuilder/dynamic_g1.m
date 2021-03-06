function g1 = dynamic_g1(T, y, x, params, steady_state, it_, T_flag)
% function g1 = dynamic_g1(T, y, x, params, steady_state, it_, T_flag)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T             [#temp variables by 1]     double   vector of temporary terms to be filled by function
%   y             [#dynamic variables by 1]  double   vector of endogenous variables in the order stored
%                                                     in M_.lead_lag_incidence; see the Manual
%   x             [nperiods by M_.exo_nbr]   double   matrix of exogenous variables (in declaration order)
%                                                     for all simulation periods
%   steady_state  [M_.endo_nbr by 1]         double   vector of steady state values
%   params        [M_.param_nbr by 1]        double   vector of parameter values in declaration order
%   it_           scalar                     double   time period for exogenous variables for which
%                                                     to evaluate the model
%   T_flag        boolean                    boolean  flag saying whether or not to calculate temporary terms
%
% Output:
%   g1
%

if T_flag
    T = KimModTheBuilder.dynamic_g1_tt(T, y, x, params, steady_state, it_);
end
g1 = zeros(8, 17);
g1(1,5)=y(9)*T(28)*T(29)-(1-T(4))*getPowerDeriv(y(5)-T(3),(-1),1);
g1(1,12)=(-(T(10)*(-T(6))*getPowerDeriv(y(12)-T(5),(-1),1)));
g1(1,6)=y(9)*T(29)*T(36);
g1(1,9)=T(8);
g1(2,5)=y(9)*T(30)*T(31);
g1(2,1)=(-(y(9)*y(10)*T(13)*(T(33)-(T(15)*T(32)+T(20)*params(9)*T(32)))));
g1(2,6)=y(9)*T(31)*T(37)-(y(9)*y(10)*T(13)*(T(39)-(T(15)*T(38)+T(20)*params(9)*T(38)))+T(17)*T(42)+T(22)*params(9)*T(40));
g1(2,13)=(-(T(17)*T(44)+T(22)*params(9)*T(43)));
g1(2,9)=T(12)-(1-T(14)-T(20)*T(15))*y(10)*T(13);
g1(2,15)=(-(T(17)*T(21)*T(16)*T(10)*y(16)));
g1(2,10)=(-((1-T(14)-T(20)*T(15))*y(9)*T(13)));
g1(2,16)=(-(T(17)*T(21)*T(10)*y(15)*T(16)));
g1(3,14)=(-(T(10)*y(15)));
g1(3,9)=y(10);
g1(3,15)=(-(T(10)*(y(14)+y(16)*(1-params(3)))));
g1(3,10)=y(9);
g1(3,16)=(-(T(10)*y(15)*(1-params(3))));
g1(4,1)=(-(y(6)*T(13)*T(33)));
g1(4,6)=(-(y(6)*T(13)*T(39)+T(13)*(1-T(14))));
g1(4,2)=(-(1-params(3)));
g1(4,8)=1;
g1(5,4)=1;
g1(5,2)=(-(y(11)*T(45)));
g1(5,11)=(-T(18));
g1(6,7)=1;
g1(6,2)=(-(y(11)*T(46)));
g1(6,11)=(-T(19));
g1(7,4)=1;
g1(7,5)=(-T(27));
g1(7,6)=(-T(35));
g1(8,3)=(-(params(6)*1/y(3)));
g1(8,11)=1/y(11);
g1(8,17)=(-params(7));

end
