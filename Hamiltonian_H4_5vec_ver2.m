% -------------------------------------------------------------------------
% Hks = Hks_12+ Hks_13 + Hks_14+ Hks_23 + Hks_34;
% Number of terms: 331
% -------------------------------------------------------------------------
%   78*s11 + 52*s12 + 78*s13 + 52*s14 - 52*s15 - 52*s16 - 52*s17 - 52*s18 
% + 78*s21 + 52*s22 + 78*s23 + 52*s24 - 52*s25 - 52*s26 - 52*s27 - 52*s28 
% + 78*s31 + 52*s32 + 78*s33 + 52*s34 - 52*s35 - 52*s36 - 52*s37 - 52*s38 
% + 78*s41 + 52*s42 + 78*s43 + 52*s44 - 52*s45 - 52*s46 - 52*s47 - 52*s48 
% - 52*s110 - 52*s210 - 52*s310 - 52*s410 + 20*s11*s12 + 20*s11*s13 
% + 20*s11*s14 + 20*s12*s13 - 40*s11*s15 - 40*s11*s16 - 40*s12*s15 
% + 20*s13*s14 - 40*s11*s17 - 40*s13*s16 - 40*s12*s18 - 40*s13*s18 
% - 40*s14*s17 + 6*s11*s21 + 2*s11*s22 + 2*s12*s21 + 2*s11*s23 + 4*s12*s22
% + 2*s13*s21 + 2*s11*s24 + 2*s12*s23 + 2*s13*s22 + 2*s14*s21 - 4*s11*s25 
% + 6*s13*s23 - 4*s15*s21 - 4*s11*s26 - 4*s12*s25 + 2*s13*s24 + 2*s14*s23 
% - 4*s15*s22 - 4*s16*s21 - 4*s11*s27 + 4*s14*s24 - 4*s17*s21 - 4*s13*s26 
% - 4*s16*s23 - 4*s12*s28 + 8*s15*s25 - 4*s18*s22 - 4*s13*s28 - 4*s14*s27 
% - 4*s17*s24 - 4*s18*s23 + 6*s11*s31 + 8*s16*s26 + 2*s11*s32 + 2*s12*s31
% + 20*s21*s22 + 2*s11*s33 + 4*s12*s32 + 2*s13*s31 + 8*s17*s27 + 20*s21*s23
% + 2*s11*s34 + 2*s12*s33 + 2*s13*s32 + 2*s14*s31 + 20*s21*s24 + 20*s22*s23
% - 4*s11*s35 + 6*s13*s33 - 4*s15*s31 + 8*s18*s28 - 40*s21*s25 - 4*s11*s36
% - 4*s12*s35 + 2*s13*s34 + 2*s14*s33 - 4*s15*s32 - 4*s16*s31 - 40*s21*s26
% - 40*s22*s25 + 20*s23*s24 - 4*s11*s37 + 4*s14*s34 - 4*s17*s31 - 40*s21*s27
% - 4*s13*s36 - 4*s16*s33 - 40*s23*s26 - 4*s12*s38 + 8*s15*s35 - 4*s18*s32
% - 40*s22*s28 - 4*s13*s38 - 4*s14*s37 - 4*s17*s34 - 4*s18*s33 - 40*s23*s28
% - 40*s24*s27 + 6*s11*s41 + 8*s16*s36 + 6*s21*s31 + 2*s11*s42 + 2*s12*s41 
% + 2*s21*s32 + 2*s22*s31 + 2*s11*s43 + 4*s12*s42 + 2*s13*s41 + 8*s17*s37 
% + 2*s21*s33 + 4*s22*s32 + 2*s23*s31 + 2*s11*s44 + 2*s12*s43 + 2*s13*s42
% + 2*s14*s41 + 2*s21*s34 + 2*s22*s33 + 2*s23*s32 + 2*s24*s31 - 4*s11*s45 
% + 6*s13*s43 - 4*s15*s41 + 8*s18*s38 - 4*s21*s35 + 6*s23*s33 - 4*s25*s31 
% - 4*s11*s46 - 4*s12*s45 + 2*s13*s44 + 2*s14*s43 - 4*s15*s42 - 4*s16*s41 
% - 4*s21*s36 - 4*s22*s35 + 2*s23*s34 + 2*s24*s33 - 4*s25*s32 - 4*s26*s31 
% - 4*s11*s47 + 4*s14*s44 - 4*s17*s41 - 4*s21*s37 + 4*s24*s34 - 4*s27*s31 
% - 4*s13*s46 - 4*s16*s43 - 4*s23*s36 - 4*s26*s33 - 4*s12*s48 + 8*s15*s45 
% - 4*s18*s42 - 4*s22*s38 + 8*s25*s35 - 4*s28*s32 - 4*s13*s48 - 4*s14*s47 
% - 4*s17*s44 - 4*s18*s43 - 4*s23*s38 - 4*s24*s37 - 4*s27*s34 - 4*s28*s33 
% + 8*s16*s46 + 6*s21*s41 + 8*s26*s36 + 2*s21*s42 + 2*s22*s41 + 20*s31*s32
% + 8*s17*s47 + 2*s21*s43 + 4*s22*s42 + 2*s23*s41 + 8*s27*s37 + 20*s31*s33
% + 2*s21*s44 + 2*s22*s43 + 2*s23*s42 + 2*s24*s41 + 20*s31*s34 + 20*s32*s33
% + 8*s18*s48 - 4*s21*s45 + 6*s23*s43 - 4*s25*s41 + 8*s28*s38 - 40*s31*s35 
% - 4*s21*s46 - 4*s22*s45 + 2*s23*s44 + 2*s24*s43 - 4*s25*s42 - 4*s26*s41 
% - 40*s31*s36 - 40*s32*s35 + 20*s33*s34 - 4*s21*s47 + 4*s24*s44 - 4*s27*s41 
% - 40*s31*s37 - 4*s23*s46 - 4*s26*s43 - 40*s33*s36 - 4*s22*s48 + 8*s25*s45
% - 4*s28*s42 - 40*s32*s38 - 4*s23*s48 - 4*s24*s47 - 4*s27*s44 - 4*s28*s43 
% - 40*s33*s38 - 40*s34*s37 + 8*s26*s46 + 6*s31*s41 + 2*s31*s42 + 2*s32*s41
% + 8*s27*s47 + 2*s31*s43 + 4*s32*s42 + 2*s33*s41 + 2*s31*s44 + 2*s32*s43
% + 2*s33*s42 + 2*s34*s41 + 8*s28*s48 - 4*s31*s45 + 6*s33*s43 - 4*s35*s41 
% - 4*s31*s46 - 4*s32*s45 + 2*s33*s44 + 2*s34*s43 - 4*s35*s42 - 4*s36*s41 
% - 4*s31*s47 + 4*s34*s44 - 4*s37*s41 - 4*s33*s46 - 4*s36*s43 - 4*s32*s48 
% + 8*s35*s45 - 4*s38*s42 - 4*s33*s48 - 4*s34*s47 - 4*s37*s44 - 4*s38*s43 
% + 8*s36*s46 + 20*s41*s42 + 8*s37*s47 + 20*s41*s43 + 20*s41*s44 + 20*s42*s43 
% + 8*s38*s48 - 40*s41*s45 - 40*s41*s46 - 40*s42*s45 + 20*s43*s44 - 40*s41*s47 
% - 40*s43*s46 - 40*s42*s48 - 40*s43*s48 - 40*s44*s47 - 40*s13*s110 
% - 40*s14*s110 - 4*s23*s110 - 4*s24*s110 - 4*s33*s110 - 4*s34*s110 
% - 4*s43*s110 - 4*s44*s110 - 4*s13*s210 - 4*s14*s210 - 40*s23*s210 
% - 40*s24*s210 - 4*s33*s210 - 4*s34*s210 - 4*s43*s210 - 4*s44*s210 
% + 8*s110*s210 - 4*s13*s310 - 4*s14*s310 - 4*s23*s310 - 4*s24*s310 
% - 40*s33*s310 - 40*s34*s310 - 4*s43*s310 - 4*s44*s310 + 8*s110*s310 
% - 4*s13*s410 - 4*s14*s410 - 4*s23*s410 - 4*s24*s410 - 4*s33*s410 
% - 4*s34*s410 - 40*s43*s410 - 40*s44*s410 + 8*s110*s410 + 8*s210*s310 
% + 8*s210*s410 + 8*s310*s410 + 1280
 