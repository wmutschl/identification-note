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
    T = AnSchoModTheBuilder.dynamic_g1_tt(T, y, x, params, steady_state, it_);
end
g1 = zeros(9, 23);
g1(1,9)=T(48);
g1(1,15)=(-(y(10)*T(9)*T(52)/y(16)/(T(11)*y(17))));
g1(1,10)=(-(T(9)*T(10)/y(16)/(T(11)*y(17))));
g1(1,16)=(-((-(y(10)*T(9)*T(10)))/(y(16)*y(16))/(T(11)*y(17))));
g1(1,17)=(-((-(y(10)*T(9)*T(10)/y(16)*T(11)))/(T(11)*y(17)*T(11)*y(17))));
g1(2,8)=T(9)*params(15)*y(16)*T(40);
g1(2,14)=T(9)*params(15)*y(16)*T(47);
g1(2,9)=(-(1/params(5)*(-(T(48)/T(5)*T(49)))-T(9)*params(15)*y(16)*T(51)));
g1(2,15)=T(9)*params(15)*y(16)*T(54);
g1(2,11)=(-(params(15)*(y(11)-T(12))+y(11)*params(15)-params(15)/(2*params(5))*2*(y(11)-T(12))));
g1(2,16)=T(9)*params(15)*(T(18)+y(16)*T(17));
g1(3,8)=1-T(16)*params(15)/2-(1-1/y(12));
g1(3,9)=(-1);
g1(3,11)=(-(y(8)*params(15)/2*2*(y(11)-T(12))));
g1(3,12)=(-(y(8)*(-((-1)/(y(12)*y(12))))));
g1(4,1)=(-(T(21)*T(20)*T(30)*T(35)*T(36)*T(37)));
g1(4,8)=(-(T(21)*T(20)*T(37)*T(45)));
g1(4,2)=(-(T(21)*T(19)*T(55)));
g1(4,10)=1;
g1(4,11)=(-(T(21)*T(20)*T(37)*T(34)*T(57)));
g1(4,3)=(-(T(21)*T(20)*T(37)*T(30)*T(36)*T(59)));
g1(4,12)=(-(T(21)*T(20)*T(37)*(T(34)*T(61)+T(30)*T(36)*T(63))));
g1(4,18)=(-(T(19)*T(20)*params(12)*T(21)));
g1(5,3)=(-(params(10)*1/y(3)));
g1(5,12)=1/y(12);
g1(5,19)=(-params(13));
g1(6,4)=(-(params(11)*1/y(4)));
g1(6,13)=1/y(13);
g1(6,20)=(-params(14));
g1(7,5)=1;
g1(7,1)=(-(100*(-(T(38)/T(23)))));
g1(7,8)=(-(100*T(38)/T(22)));
g1(7,13)=(-(100*1/(steady_state(9))/T(24)));
g1(7,21)=(-params(17));
g1(8,6)=1;
g1(8,11)=(-(400*1/(steady_state(7))/T(25)));
g1(8,22)=(-params(18));
g1(9,7)=1;
g1(9,10)=(-(400*1/(steady_state(6))/T(26)));
g1(9,23)=(-params(19));

end