

%%%% This code is used to load the resting Ktrans values for the Dual-Bolus
%%%% method and the signal corrected T1 modeling method in several subjects
%%%% in order to determine the correlation between the 2 methods.


close all,clear all,clc

figure(456)
xlabel('Dual-Bolus Ktrans Estimates (ml/min/g)')
ylabel('T1 Ktrans Estimates (ml/min/g)')
axis([0 2 0 2]),axis square

%%% To initialize the vectors of T1 and Dual-Bolus (DB) Ktrans estimates
T1_REST=[];
DB_REST=[];

%%% Below is a subset of the rest data above, which has a corresponding Lexiscan and Adenosine stress scan
REST2=[];
LEX=[];
ADEN=[];


%%% Below is the per-slice and per-patient average results for the rest and stress data
T1_REST_avgSlice=[];T1_REST_avgPt=[];
DB_REST_avgSlice=[];DB_REST_avgPt=[];
REST2_avgSlice=[];REST2_avgPt=[];
LEX_avgSlice=[];LEX_avgPt=[];
ADEN_avgSlice=[];ADEN_avgPt=[];


%%%%% NOTE: to see results for 1 region per slice, change all the ".6."
%%%%% (ie. "point six point") to ".6." (ie. "point one point") using the
%%%%% Ctl-h 'Find & Replace' command.
%%%%% IF "1" REGION IS USED INSTEAD OF "6", 3 STRESS DATASETS MUST BE
%%%%% ADJUSTED TO MATCH SLICES:
%%%%% P043009, Line ~521, REST2=[REST2 T1_Ktrans_rest1(1:12)];--> change to (1:2)
%%%%% P061109, Line ~777, REST2=[REST2 T1_Ktrans_rest1(7:24)];--> change to (2:4)
%%%%% P081009, Line ~1973, REST2=[REST2 T1_Ktrans_rest1(13:18) T1_Ktrans_rest1(1:6)];--> change to (3) (1)





%%%---------P042809, scaled by 1/5 dose---------%%%
%%%% To load the Dual-Bolus Ktrans estimates
f='/v/raid1/npack/Processing/P042809_nufft_New/Output/flowvalues.study44.slice1.6.1_fixedDelay0.AIF_38_4_5.txt.full.txt';
Ktrans = getKtrans(f);
r1=Ktrans;
f='/v/raid1/npack/Processing/P042809_nufft_New/Output/flowvalues.study44.slice2.6.1_fixedDelay0.AIF_38_4_5.txt.full.txt';
Ktrans = getKtrans(f);
r2=Ktrans;
f='/v/raid1/npack/Processing/P042809_nufft_New/Output/flowvalues.study44.slice3.6.1_fixedDelay0.AIF_38_4_5.txt.full.txt';
Ktrans = getKtrans(f);
r3=Ktrans;
f='/v/raid1/npack/Processing/P042809_nufft_New/Output/flowvalues.study44.slice4.6.1_fixedDelay0.AIF_38_4_5.txt.full.txt';
Ktrans = getKtrans(f);
r4=Ktrans;
DB_Ktrans_rest1=[r1 r2 r3 r4];DB_REST_avgSlice=[DB_REST_avgSlice mean(r1) mean(r2) mean(r3) mean(r4)];DB_REST_avgPt=[DB_REST_avgPt mean(DB_Ktrans_rest1)];
DB_REST=[DB_REST DB_Ktrans_rest1];


%%%% To load the T1 Ktrans estimates
f='/v/raid1/npack/Processing/P042809_nufft_24rays/Output/flowvalues.study44.slice1.6.1_fixedDelay0.AIF_44_14_1.txt.full.txt';
Ktrans = getKtrans(f);
r1=Ktrans;
f='/v/raid1/npack/Processing/P042809_nufft_24rays/Output/flowvalues.study44.slice2.6.1_fixedDelay0.AIF_44_14_1.txt.full.txt';
Ktrans = getKtrans(f);
r2=Ktrans;
f='/v/raid1/npack/Processing/P042809_nufft_24rays/Output/flowvalues.study44.slice3.6.1_fixedDelay0.AIF_44_14_1.txt.full.txt';
Ktrans = getKtrans(f);
r3=Ktrans;
f='/v/raid1/npack/Processing/P042809_nufft_24rays/Output/flowvalues.study44.slice4.6.1_fixedDelay0.AIF_44_14_1.txt.full.txt';
Ktrans = getKtrans(f);
r4=Ktrans;
T1_Ktrans_rest1=[r1 r2 r3 r4];T1_REST_avgSlice=[T1_REST_avgSlice mean(r1) mean(r2) mean(r3) mean(r4)];T1_REST_avgPt=[T1_REST_avgPt mean(T1_Ktrans_rest1)];
T1_REST=[T1_REST T1_Ktrans_rest1];

figure(1),hold on,plot(DB_Ktrans_rest1,T1_Ktrans_rest1,'bo','LineWidth',3,'MarkerSize',8)
disp('Plot held for current patient; press key to continue.'),pause(.1)











%%%---------P042809, scaled by 1/10 dose---------%%%
%%%% To load the Dual-Bolus Ktrans estimates
f='/v/raid1/npack/Processing/P042809_nufft_New/Output/flowvalues.study44.slice1.6.1_fixedDelay0.AIF_34_1_10.txt.full.txt';
Ktrans = getKtrans(f);
r1=Ktrans;
f='/v/raid1/npack/Processing/P042809_nufft_New/Output/flowvalues.study44.slice2.6.1_fixedDelay0.AIF_34_1_10.txt.full.txt';
Ktrans = getKtrans(f);
r2=Ktrans;
f='/v/raid1/npack/Processing/P042809_nufft_New/Output/flowvalues.study44.slice3.6.1_fixedDelay0.AIF_34_1_10.txt.full.txt';
Ktrans = getKtrans(f);
r3=Ktrans;
f='/v/raid1/npack/Processing/P042809_nufft_New/Output/flowvalues.study44.slice4.6.1_fixedDelay0.AIF_34_1_10.txt.full.txt';
Ktrans = getKtrans(f);
r4=Ktrans;
DB_Ktrans_rest1=[r1 r2 r3 r4];
DB_REST_avgSlice=[DB_REST_avgSlice mean(r1) mean(r2) mean(r3) mean(r4)];
DB_REST_avgPt=[DB_REST_avgPt mean(DB_Ktrans_rest1)];
DB_REST=[DB_REST DB_Ktrans_rest1];


%%%% To load the T1 Ktrans estimates
f='/v/raid1/npack/Processing/P042809_nufft_24rays/Output/flowvalues.study44.slice1.6.1_fixedDelay0.AIF_44_11_1.txt.full.txt';
Ktrans = getKtrans(f);
r1=Ktrans;
f='/v/raid1/npack/Processing/P042809_nufft_24rays/Output/flowvalues.study44.slice2.6.1_fixedDelay0.AIF_44_11_1.txt.full.txt';
Ktrans = getKtrans(f);
r2=Ktrans;
f='/v/raid1/npack/Processing/P042809_nufft_24rays/Output/flowvalues.study44.slice3.6.1_fixedDelay0.AIF_44_11_1.txt.full.txt';
Ktrans = getKtrans(f);
r3=Ktrans;
f='/v/raid1/npack/Processing/P042809_nufft_24rays/Output/flowvalues.study44.slice4.6.1_fixedDelay0.AIF_44_11_1.txt.full.txt';
Ktrans = getKtrans(f);
r4=Ktrans;
T1_Ktrans_rest1=[r1 r2 r3 r4];
T1_REST_avgSlice=[T1_REST_avgSlice mean(r1) mean(r2) mean(r3) mean(r4)];
T1_REST_avgPt=[T1_REST_avgPt mean(T1_Ktrans_rest1)];
T1_REST=[T1_REST T1_Ktrans_rest1];

figure(1),hold on,plot(DB_Ktrans_rest1,T1_Ktrans_rest1,'bo','LineWidth',3,'MarkerSize',8)
disp('Plot held for current patient; press key to continue.'),pause(.1)












%%%---------P043009---------%%%
%%%% To load the Dual-Bolus Ktrans estimates
f='/v/raid1/npack/Processing/P043009_nufft_New/Output/flowvalues.study27.slice1.6.1_fixedDelay0.AIF_23_1_10.txt.full.txt';
Ktrans = getKtrans(f);
r1=Ktrans;
f='/v/raid1/npack/Processing/P043009_nufft_New/Output/flowvalues.study27.slice3.6.1_fixedDelay0.AIF_23_1_10.txt.full.txt';
Ktrans = getKtrans(f);
r3=Ktrans;
f='/v/raid1/npack/Processing/P043009_nufft_New/Output/flowvalues.study27.slice4.6.1_fixedDelay0.AIF_23_1_10.txt.full.txt';
Ktrans = getKtrans(f);
r4=Ktrans;
DB_Ktrans_rest1=[r1 r3 r4];DB_REST_avgSlice=[DB_REST_avgSlice mean(r1) mean(r3) mean(r4)];DB_REST_avgPt=[DB_REST_avgPt mean(DB_Ktrans_rest1)];
DB_REST=[DB_REST DB_Ktrans_rest1];


%%%% To load the T1 Ktrans estimates
f='/v/raid1/npack/Processing/P043009_nufft_24rays/Output/flowvalues.study27.slice1.6.1_fixedDelay0.AIF_27_11_1.txt.full.txt';
Ktrans = getKtrans(f);
r1=Ktrans;
f='/v/raid1/npack/Processing/P043009_nufft_24rays/Output/flowvalues.study27.slice3.6.1_fixedDelay0.AIF_27_11_1.txt.full.txt';
Ktrans = getKtrans(f);
r3=Ktrans;
f='/v/raid1/npack/Processing/P043009_nufft_24rays/Output/flowvalues.study27.slice4.6.1_fixedDelay0.AIF_27_11_1.txt.full.txt';
Ktrans = getKtrans(f);
r4=Ktrans;
T1_Ktrans_rest1=[r1 r3 r4];T1_REST_avgSlice=[T1_REST_avgSlice mean(r1) mean(r3) mean(r4)];T1_REST_avgPt=[T1_REST_avgPt mean(T1_Ktrans_rest1)];
T1_REST=[T1_REST T1_Ktrans_rest1];

%%%% To load the Adenosine Ktrans estimates
f='/v/raid1/npack/Processing/P043009_nufft_24rays/Output/flowvalues.study35.slice1.6.1_fixedDelay0.AIF_35_13_1.txt.full.txt';
Ktrans = getKtrans(f);
s1=Ktrans;
f='/v/raid1/npack/Processing/P043009_nufft_24rays/Output/flowvalues.study35.slice3.6.1_fixedDelay0.AIF_35_13_1.txt.full.txt';
Ktrans = getKtrans(f);
s3=Ktrans;
Aden_1=[s1 s3];ADEN_avgSlice=[ADEN_avgSlice mean(s1) mean(s3)];ADEN_avgPt=[ADEN_avgPt mean(Aden_1)];
ADEN=[ADEN Aden_1];

%%%% To load the Lexiscan Ktrans estimates
f='/v/raid1/npack/Processing/P043009_nufft_24rays/Output/flowvalues.study56.slice1.6.1_fixedDelay0.AIF_56_13_1.txt.full.txt';
Ktrans = getKtrans(f);
s1=Ktrans;
f='/v/raid1/npack/Processing/P043009_nufft_24rays/Output/flowvalues.study56.slice3.6.1_fixedDelay0.AIF_56_13_1.txt.full.txt';
Ktrans = getKtrans(f);
s3=Ktrans;
Lex_1=[s1 s3];LEX_avgSlice=[LEX_avgSlice mean(s1) mean(s3)];LEX_avgPt=[LEX_avgPt mean(Lex_1)];
LEX=[LEX Lex_1];

REST2=[REST2 T1_Ktrans_rest1(1:12)];REST2_avgSlice=[REST2_avgSlice mean(T1_Ktrans_rest1(1:6)) mean(T1_Ktrans_rest1(7:12))];REST2_avgPt=[REST2_avgPt mean(T1_Ktrans_rest1(1:12))];
figure(1),hold on,plot(DB_Ktrans_rest1,T1_Ktrans_rest1,'bo','LineWidth',3,'MarkerSize',8)
disp('Plot held for current patient; press key to continue.'),pause(.1)











%%%---------P061109--------%%%
%%%% To load the Dual-Bolus Ktrans estimates
f='/v/raid1/npack/Processing/P061109_nufft/Output/flowvalues.study29.slice1.6.1_fixedDelay0.AIF_17_1_5.txt.full.txt';
Ktrans = getKtrans(f);
r1=Ktrans;
f='/v/raid1/npack/Processing/P061109_nufft/Output/flowvalues.study29.slice2.6.1_fixedDelay0.AIF_17_1_5.txt.full.txt';
Ktrans = getKtrans(f);
r2=Ktrans;
f='/v/raid1/npack/Processing/P061109_nufft/Output/flowvalues.study29.slice4.6.1_fixedDelay0.AIF_17_1_5.txt.full.txt';
Ktrans = getKtrans(f);
r4=Ktrans;
f='/v/raid1/npack/Processing/P061109_nufft/Output/flowvalues.study29.slice5.6.1_fixedDelay0.AIF_17_1_5.txt.full.txt';
Ktrans = getKtrans(f);
r5=Ktrans;
DB_Ktrans_rest1=[r1 r2 r4 r5];DB_REST_avgSlice=[DB_REST_avgSlice mean(r1) mean(r2) mean(r4) mean(r5)];DB_REST_avgPt=[DB_REST_avgPt mean(DB_Ktrans_rest1)];
DB_REST=[DB_REST DB_Ktrans_rest1];


%%%% To load the T1 Ktrans estimates
f='/v/raid1/npack/Processing/P061109_nufft_24rays/Output/flowvalues.study29.slice1.6.1_fixedDelay0.AIF_29_11_1.txt.full.txt';
Ktrans = getKtrans(f);
r1=Ktrans;
f='/v/raid1/npack/Processing/P061109_nufft_24rays/Output/flowvalues.study29.slice2.6.1_fixedDelay0.AIF_29_11_1.txt.full.txt';
Ktrans = getKtrans(f);
r2=Ktrans;
f='/v/raid1/npack/Processing/P061109_nufft_24rays/Output/flowvalues.study29.slice4.6.1_fixedDelay0.AIF_29_11_1.txt.full.txt';
Ktrans = getKtrans(f);
r4=Ktrans;
f='/v/raid1/npack/Processing/P061109_nufft_24rays/Output/flowvalues.study29.slice5.6.1_fixedDelay0.AIF_29_11_1.txt.full.txt';
Ktrans = getKtrans(f);
r5=Ktrans;
T1_Ktrans_rest1=[r1 r2 r4 r5];T1_REST_avgSlice=[T1_REST_avgSlice mean(r1) mean(r2) mean(r4) mean(r5)];T1_REST_avgPt=[T1_REST_avgPt mean(T1_Ktrans_rest1)];
T1_REST=[T1_REST T1_Ktrans_rest1];

%%%% To load the Adenosine Ktrans estimates
f='/v/raid1/npack/Processing/P061109_nufft_24rays/Output/flowvalues.study35.slice1.6.1_fixedDelay0.AIF_35_11_1.txt.full.txt';
Ktrans = getKtrans(f);
s1=Ktrans;
f='/v/raid1/npack/Processing/P061109_nufft_24rays/Output/flowvalues.study35.slice2.6.1_fixedDelay0.AIF_35_11_1.txt.full.txt';
Ktrans = getKtrans(f);
s2=Ktrans;
f='/v/raid1/npack/Processing/P061109_nufft_24rays/Output/flowvalues.study35.slice3.6.1_fixedDelay0.AIF_35_11_1.txt.full.txt';
Ktrans = getKtrans(f);
s3=Ktrans;
Aden_1=[s1 s2 s3];ADEN_avgSlice=[ADEN_avgSlice mean(s1) mean(s2) mean(s3)];ADEN_avgPt=[ADEN_avgPt mean(Aden_1)];
ADEN=[ADEN Aden_1];

%%%% To load the Lexiscan Ktrans estimates
f='/v/raid1/npack/Processing/P061109_nufft_24rays/Output/flowvalues.study57.slice1.6.1_fixedDelay0.AIF_57_11_1.txt.full.txt';
Ktrans = getKtrans(f);
s1=Ktrans;
f='/v/raid1/npack/Processing/P061109_nufft_24rays/Output/flowvalues.study57.slice2.6.1_fixedDelay0.AIF_57_11_1.txt.full.txt';
Ktrans = getKtrans(f);
s2=Ktrans;
f='/v/raid1/npack/Processing/P061109_nufft_24rays/Output/flowvalues.study57.slice3.6.1_fixedDelay0.AIF_57_11_1.txt.full.txt';
Ktrans = getKtrans(f);
s3=Ktrans;
Lex_1=[s1 s2 s3];LEX_avgSlice=[LEX_avgSlice mean(s1) mean(s2) mean(s3)];LEX_avgPt=[LEX_avgPt mean(Lex_1)];
LEX=[LEX Lex_1];

REST2=[REST2 T1_Ktrans_rest1(7:24)]; REST2_avgSlice=[REST2_avgSlice mean(T1_Ktrans_rest1(7:12)) mean(T1_Ktrans_rest1(13:18)) mean(T1_Ktrans_rest1(19:24))];REST2_avgPt=[REST2_avgPt mean(T1_Ktrans_rest1(7:24))];
figure(1),hold on,plot(DB_Ktrans_rest1,T1_Ktrans_rest1,'bo','LineWidth',3,'MarkerSize',8)
disp('Plot held for current patient; press key to continue.'),pause(.1)











%%%---------P062309, scaled by 1/20 dose---------%%%
%%%% To load the Dual-Bolus Ktrans estimates
f='/v/raid1/npack/Processing/P062309_nufft/Output/flowvalues.study56.slice1.6.1_fixedDelay0.AIF_16_1_20.txt.full.txt';
Ktrans = getKtrans(f);
r1=Ktrans;
f='/v/raid1/npack/Processing/P062309_nufft/Output/flowvalues.study56.slice2.6.1_fixedDelay0.AIF_16_1_20.txt.full.txt';
Ktrans = getKtrans(f);
r2=Ktrans;
f='/v/raid1/npack/Processing/P062309_nufft/Output/flowvalues.study56.slice3.6.1_fixedDelay0.AIF_16_1_20.txt.full.txt';
Ktrans = getKtrans(f);
r3=Ktrans;
DB_Ktrans_rest1=[r1 r2 r3];DB_REST_avgSlice=[DB_REST_avgSlice mean(r1) mean(r2) mean(r3)];DB_REST_avgPt=[DB_REST_avgPt mean(DB_Ktrans_rest1)];
DB_REST=[DB_REST DB_Ktrans_rest1];


%%%% To load the T1 Ktrans estimates
f='/v/raid1/npack/Processing/P062309_nufft_24rays/Output/flowvalues.study56.slice1.6.1_fixedDelay0.AIF_56_11_1.txt.full.txt';
Ktrans = getKtrans(f);
r1=Ktrans;
f='/v/raid1/npack/Processing/P062309_nufft_24rays/Output/flowvalues.study56.slice2.6.1_fixedDelay0.AIF_56_11_1.txt.full.txt';
Ktrans = getKtrans(f);
r2=Ktrans;
f='/v/raid1/npack/Processing/P062309_nufft_24rays/Output/flowvalues.study56.slice3.6.1_fixedDelay0.AIF_56_11_1.txt.full.txt';
Ktrans = getKtrans(f);
r3=Ktrans;
T1_Ktrans_rest1=[r1 r2 r3];T1_REST_avgSlice=[T1_REST_avgSlice mean(r1) mean(r2) mean(r3)];T1_REST_avgPt=[T1_REST_avgPt mean(T1_Ktrans_rest1)];
T1_REST=[T1_REST T1_Ktrans_rest1];

figure(1),hold on,plot(DB_Ktrans_rest1,T1_Ktrans_rest1,'bo','LineWidth',3,'MarkerSize',8)
disp('Plot held for current patient; press key to continue.'),pause(.1)











%%%---------P062309, scaled by 1/5 dose---------%%%
%%%% To load the Dual-Bolus Ktrans estimates
f='/v/raid1/npack/Processing/P062309_nufft/Output/flowvalues.study56.slice1.6.1_fixedDelay0.AIF_36_1_5.txt.full.txt';
Ktrans = getKtrans(f);
r1=Ktrans;
f='/v/raid1/npack/Processing/P062309_nufft/Output/flowvalues.study56.slice2.6.1_fixedDelay0.AIF_36_1_5.txt.full.txt';
Ktrans = getKtrans(f);
r2=Ktrans;
f='/v/raid1/npack/Processing/P062309_nufft/Output/flowvalues.study56.slice3.6.1_fixedDelay0.AIF_36_1_5.txt.full.txt';
Ktrans = getKtrans(f);
r3=Ktrans;
DB_Ktrans_rest1=[r1 r2 r3];DB_REST_avgSlice=[DB_REST_avgSlice mean(r1) mean(r2) mean(r3)];DB_REST_avgPt=[DB_REST_avgPt mean(DB_Ktrans_rest1)];
DB_REST=[DB_REST DB_Ktrans_rest1];


%%%% To load the T1 Ktrans estimates
f='/v/raid1/npack/Processing/P062309_nufft_24rays/Output/flowvalues.study56.slice1.6.1_fixedDelay0.AIF_56_11_1.txt.full.txt';
Ktrans = getKtrans(f);
r1=Ktrans;
f='/v/raid1/npack/Processing/P062309_nufft_24rays/Output/flowvalues.study56.slice2.6.1_fixedDelay0.AIF_56_11_1.txt.full.txt';
Ktrans = getKtrans(f);
r2=Ktrans;
f='/v/raid1/npack/Processing/P062309_nufft_24rays/Output/flowvalues.study56.slice3.6.1_fixedDelay0.AIF_56_11_1.txt.full.txt';
Ktrans = getKtrans(f);
r3=Ktrans;
T1_Ktrans_rest1=[r1 r2 r3];T1_REST_avgSlice=[T1_REST_avgSlice mean(r1) mean(r2) mean(r3)];T1_REST_avgPt=[T1_REST_avgPt mean(T1_Ktrans_rest1)];
T1_REST=[T1_REST T1_Ktrans_rest1];

figure(1),hold on,plot(DB_Ktrans_rest1,T1_Ktrans_rest1,'bo','LineWidth',3,'MarkerSize',8)
disp('Plot held for current patient; press key to continue.'),pause(.1)












%%%---------P070909---------%%%
%%%% To load the Dual-Bolus Ktrans estimates
f='/v/raid1/npack/Processing/P070909_nufft/Output/flowvalues.study21.slice1.6.1_fixedDelay0.AIF_16_3_5.txt.full.txt';
Ktrans = getKtrans(f);
r1=Ktrans;
f='/v/raid1/npack/Processing/P070909_nufft/Output/flowvalues.study21.slice2.6.1_fixedDelay0.AIF_16_3_5.txt.full.txt';
Ktrans = getKtrans(f);
r2=Ktrans;
f='/v/raid1/npack/Processing/P070909_nufft/Output/flowvalues.study21.slice3.6.1_fixedDelay0.AIF_16_3_5.txt.full.txt';
Ktrans = getKtrans(f);
r3=Ktrans;
DB_Ktrans_rest1=[r1 r2 r3];DB_REST_avgSlice=[DB_REST_avgSlice mean(r1) mean(r2) mean(r3)];DB_REST_avgPt=[DB_REST_avgPt mean(DB_Ktrans_rest1)];
DB_REST=[DB_REST DB_Ktrans_rest1];


%%%% To load the T1 Ktrans estimates
f='/v/raid1/npack/Processing/P070909_nufft_24rays/Output/flowvalues.study21.slice1.6.1_fixedDelay0.AIF_21_13_1.txt.full.txt';
Ktrans = getKtrans(f);
r1=Ktrans;
f='/v/raid1/npack/Processing/P070909_nufft_24rays/Output/flowvalues.study21.slice2.6.1_fixedDelay0.AIF_21_13_1.txt.full.txt';
Ktrans = getKtrans(f);
r2=Ktrans;
f='/v/raid1/npack/Processing/P070909_nufft_24rays/Output/flowvalues.study21.slice3.6.1_fixedDelay0.AIF_21_13_1.txt.full.txt';
Ktrans = getKtrans(f);
r3=Ktrans;
T1_Ktrans_rest1=[r1 r2 r3];T1_REST_avgSlice=[T1_REST_avgSlice mean(r1) mean(r2) mean(r3)];T1_REST_avgPt=[T1_REST_avgPt mean(T1_Ktrans_rest1)];
T1_REST=[T1_REST T1_Ktrans_rest1];

%%%% To load the Adenosine Ktrans estimates
f='/v/raid1/npack/Processing/P070909_nufft_24rays/Output/flowvalues.study26.slice1.6.1_fixedDelay0.AIF_26_13_1.txt.full.txt';
Ktrans = getKtrans(f);
s1=Ktrans;
f='/v/raid1/npack/Processing/P070909_nufft_24rays/Output/flowvalues.study26.slice2.6.1_fixedDelay0.AIF_26_13_1.txt.full.txt';
Ktrans = getKtrans(f);
s2=Ktrans;
f='/v/raid1/npack/Processing/P070909_nufft_24rays/Output/flowvalues.study26.slice3.6.1_fixedDelay0.AIF_26_13_1.txt.full.txt';
Ktrans = getKtrans(f);
s3=Ktrans;
Aden_1=[s1 s2 s3];ADEN_avgSlice=[ADEN_avgSlice mean(s1) mean(s2) mean(s3)];ADEN_avgPt=[ADEN_avgPt mean(Aden_1)];
ADEN=[ADEN Aden_1];

%%%% To load the Lexiscan Ktrans estimates
f='/v/raid1/npack/Processing/P070909_nufft_24rays/Output/flowvalues.study60.slice1.6.1_fixedDelay0.AIF_60_13_1.txt.full.txt';
Ktrans = getKtrans(f);
s1=Ktrans;
f='/v/raid1/npack/Processing/P070909_nufft_24rays/Output/flowvalues.study60.slice2.6.1_fixedDelay0.AIF_60_13_1.txt.full.txt';
Ktrans = getKtrans(f);
s2=Ktrans;
f='/v/raid1/npack/Processing/P070909_nufft_24rays/Output/flowvalues.study60.slice3.6.1_fixedDelay0.AIF_60_13_1.txt.full.txt';
Ktrans = getKtrans(f);
s3=Ktrans;
Lex_1=[s1 s2 s3];LEX_avgSlice=[LEX_avgSlice mean(s1) mean(s2) mean(s3)];LEX_avgPt=[LEX_avgPt mean(Lex_1)];
LEX=[LEX Lex_1];

REST2=[REST2 T1_Ktrans_rest1];REST2_avgSlice=[REST2_avgSlice mean(T1_Ktrans_rest1(1:6)) mean(T1_Ktrans_rest1(7:12)) mean(T1_Ktrans_rest1(13:18))];REST2_avgPt=[REST2_avgPt mean(T1_Ktrans_rest1)];
figure(1),hold on,plot(DB_Ktrans_rest1,T1_Ktrans_rest1,'bo','LineWidth',3,'MarkerSize',8)
disp('Plot held for current patient; press key to continue.'),pause(.1)











%%%---------P071709---------%%%
%%%% To load the Dual-Bolus Ktrans estimates
f='/v/raid1/npack/Processing/P071709_nufft/Output/flowvalues.study37.slice1.6.1_fixedDelay0.AIF_32_3_5.txt.full.txt';
Ktrans = getKtrans(f);
r1=Ktrans;
f='/v/raid1/npack/Processing/P071709_nufft/Output/flowvalues.study37.slice2.6.1_fixedDelay0.AIF_32_3_5.txt.full.txt';
Ktrans = getKtrans(f);
r2=Ktrans;
f='/v/raid1/npack/Processing/P071709_nufft/Output/flowvalues.study37.slice3.6.1_fixedDelay0.AIF_32_3_5.txt.full.txt';
Ktrans = getKtrans(f);
r3=Ktrans;
DB_Ktrans_rest1=[r1 r2 r3];DB_REST_avgSlice=[DB_REST_avgSlice mean(r1) mean(r2) mean(r3)];DB_REST_avgPt=[DB_REST_avgPt mean(DB_Ktrans_rest1)];
DB_REST=[DB_REST DB_Ktrans_rest1];


%%%% To load the T1 Ktrans estimates
f='/v/raid1/npack/Processing/P071709_nufft_24rays/Output/flowvalues.study37.slice1.6.1_fixedDelay0.AIF_37_13_1.txt.full.txt';
Ktrans = getKtrans(f);
r1=Ktrans;
f='/v/raid1/npack/Processing/P071709_nufft_24rays/Output/flowvalues.study37.slice2.6.1_fixedDelay0.AIF_37_13_1.txt.full.txt';
Ktrans = getKtrans(f);
r2=Ktrans;
f='/v/raid1/npack/Processing/P071709_nufft_24rays/Output/flowvalues.study37.slice3.6.1_fixedDelay0.AIF_37_13_1.txt.full.txt';
Ktrans = getKtrans(f);
r3=Ktrans;
T1_Ktrans_rest1=[r1 r2 r3];T1_REST_avgSlice=[T1_REST_avgSlice mean(r1) mean(r2) mean(r3)];T1_REST_avgPt=[T1_REST_avgPt mean(T1_Ktrans_rest1)];
T1_REST=[T1_REST T1_Ktrans_rest1];

%%%% To load the Adenosine Ktrans estimates
f='/v/raid1/npack/Processing/P071709_nufft_24rays/Output/flowvalues.study41.slice1.6.1_fixedDelay0.AIF_41_12_1.txt.full.txt';
Ktrans = getKtrans(f);
s1=Ktrans;
f='/v/raid1/npack/Processing/P071709_nufft_24rays/Output/flowvalues.study41.slice2.6.1_fixedDelay0.AIF_41_12_1.txt.full.txt';
Ktrans = getKtrans(f);
s2=Ktrans;
f='/v/raid1/npack/Processing/P071709_nufft_24rays/Output/flowvalues.study41.slice3.6.1_fixedDelay0.AIF_41_12_1.txt.full.txt';
Ktrans = getKtrans(f);
s3=Ktrans;
Aden_1=[s1 s2 s3];ADEN_avgSlice=[ADEN_avgSlice mean(s1) mean(s2) mean(s3)];ADEN_avgPt=[ADEN_avgPt mean(Aden_1)];
ADEN=[ADEN Aden_1];

%%%% To load the Lexiscan Ktrans estimates
f='/v/raid1/npack/Processing/P071709_nufft_24rays/Output/flowvalues.study72.slice1.6.1_fixedDelay0.AIF_72_12_1.txt.full.txt';
Ktrans = getKtrans(f);
s1=Ktrans;
f='/v/raid1/npack/Processing/P071709_nufft_24rays/Output/flowvalues.study72.slice2.6.1_fixedDelay0.AIF_72_12_1.txt.full.txt';
Ktrans = getKtrans(f);
s2=Ktrans;
f='/v/raid1/npack/Processing/P071709_nufft_24rays/Output/flowvalues.study72.slice3.6.1_fixedDelay0.AIF_72_12_1.txt.full.txt';
Ktrans = getKtrans(f);
s3=Ktrans;
Lex_1=[s1 s2 s3];LEX_avgSlice=[LEX_avgSlice mean(s1) mean(s2) mean(s3)];LEX_avgPt=[LEX_avgPt mean(Lex_1)];
LEX=[LEX Lex_1];

REST2=[REST2 T1_Ktrans_rest1];REST2_avgSlice=[REST2_avgSlice mean(T1_Ktrans_rest1(1:6)) mean(T1_Ktrans_rest1(7:12)) mean(T1_Ktrans_rest1(13:18))];REST2_avgPt=[REST2_avgPt mean(T1_Ktrans_rest1)];
figure(1),hold on,plot(DB_Ktrans_rest1,T1_Ktrans_rest1,'bo','LineWidth',3,'MarkerSize',8)
disp('Plot held for current patient; press key to continue.'),pause(.1)











%%% BAD AIFs..., no back-check valve in place, possibly along with other injection errors
% %%%---------P072309--->AIFs from each slice separately--------%%%
% %%%% To load the Dual-Bolus Ktrans estimates
% f='/v/raid1/npack/Processing/P072309_nufft/Output/flowvalues.study21.slice1.6.1_fixedDelay0.AIF_18_1_7.5.txt.full.txt';
% fid=fopen(f,'r');
% a=fscanf(fid,'%s',1);
% indexflow=1;
% while ~isempty(a)
%     a=fscanf(fid,'%s',1);
%     switch(a)
%         case {'Ktrans'}
%             a =fscanf(fid,'%s',1);
%             a =fscanf(fid,'%s',1);   % argh, must be a cleaner way!
%             aa=str2num(a);
%             Ktrans(indexflow)=aa; indexflow=indexflow+1;
%         otherwise
%     end
% end
% r1=Ktrans;
% f='/v/raid1/npack/Processing/P072309_nufft/Output/flowvalues.study21.slice2.6.1_fixedDelay0.AIF_18_2_7.5.txt.full.txt';
% fid=fopen(f,'r');
% a=fscanf(fid,'%s',1);
% indexflow=1;
% while ~isempty(a)
%     a=fscanf(fid,'%s',1);
%     switch(a)
%         case {'Ktrans'}
%             a =fscanf(fid,'%s',1);
%             a =fscanf(fid,'%s',1);   % argh, must be a cleaner way!
%             aa=str2num(a);
%             Ktrans(indexflow)=aa; indexflow=indexflow+1;
%         otherwise
%     end
% end
% r2=Ktrans;
% f='/v/raid1/npack/Processing/P072309_nufft/Output/flowvalues.study21.slice3.6.1_fixedDelay0.AIF_18_3_7.5.txt.full.txt';
% fid=fopen(f,'r');
% a=fscanf(fid,'%s',1);
% indexflow=1;
% while ~isempty(a)
%     a=fscanf(fid,'%s',1);
%     switch(a)
%         case {'Ktrans'}
%             a =fscanf(fid,'%s',1);
%             a =fscanf(fid,'%s',1);   % argh, must be a cleaner way!
%             aa=str2num(a);
%             Ktrans(indexflow)=aa; indexflow=indexflow+1;
%         otherwise
%     end
% end
% r3=Ktrans;
% f='/v/raid1/npack/Processing/P072309_nufft/Output/flowvalues.study21.slice4.6.1_fixedDelay0.AIF_18_4_7.5.txt.full.txt';
% fid=fopen(f,'r');
% a=fscanf(fid,'%s',1);
% indexflow=1;
% while ~isempty(a)
%     a=fscanf(fid,'%s',1);
%     switch(a)
%         case {'Ktrans'}
%             a =fscanf(fid,'%s',1);
%             a =fscanf(fid,'%s',1);   % argh, must be a cleaner way!
%             aa=str2num(a);
%             Ktrans(indexflow)=aa; indexflow=indexflow+1;
%         otherwise
%     end
% end
% r4=Ktrans;
% DB_Ktrans_rest1=[r1 r2 r3 r4];DB_REST_avgSlice=[DB_REST_avgSlice mean(r1) mean(r2) mean(r3) mean(r4)];DB_REST_avgPt=[DB_REST_avgPt mean(DB_Ktrans_rest1)];
% DB_REST=[DB_REST DB_Ktrans_rest1];
% 
% 
% %%%% To load the T1 Ktrans estimates
% f='/v/raid1/npack/Processing/P072309_nufft_24rays/Output/flowvalues.study21.slice1.6.1_fixedDelay0.AIF_21_11_1.txt.full.txt';
% fid=fopen(f,'r');
% a=fscanf(fid,'%s',1);
% indexflow=1;
% while ~isempty(a)
%     a=fscanf(fid,'%s',1);
%     switch(a)
%         case {'Ktrans'}
%             a =fscanf(fid,'%s',1);
%             a =fscanf(fid,'%s',1);   % argh, must be a cleaner way!
%             aa=str2num(a);
%             Ktrans(indexflow)=aa; indexflow=indexflow+1;
%         otherwise
%     end
% end
% r1=Ktrans;
% f='/v/raid1/npack/Processing/P072309_nufft_24rays/Output/flowvalues.study21.slice2.6.1_fixedDelay0.AIF_21_12_1.txt.full.txt';
% fid=fopen(f,'r');
% a=fscanf(fid,'%s',1);
% indexflow=1;
% while ~isempty(a)
%     a=fscanf(fid,'%s',1);
%     switch(a)
%         case {'Ktrans'}
%             a =fscanf(fid,'%s',1);
%             a =fscanf(fid,'%s',1);   % argh, must be a cleaner way!
%             aa=str2num(a);
%             Ktrans(indexflow)=aa; indexflow=indexflow+1;
%         otherwise
%     end
% end
% r2=Ktrans;
% f='/v/raid1/npack/Processing/P072309_nufft_24rays/Output/flowvalues.study21.slice3.6.1_fixedDelay0.AIF_21_13_1.txt.full.txt';
% fid=fopen(f,'r');
% a=fscanf(fid,'%s',1);
% indexflow=1;
% while ~isempty(a)
%     a=fscanf(fid,'%s',1);
%     switch(a)
%         case {'Ktrans'}
%             a =fscanf(fid,'%s',1);
%             a =fscanf(fid,'%s',1);   % argh, must be a cleaner way!
%             aa=str2num(a);
%             Ktrans(indexflow)=aa; indexflow=indexflow+1;
%         otherwise
%     end
% end
% r3=Ktrans;
% f='/v/raid1/npack/Processing/P072309_nufft_24rays/Output/flowvalues.study21.slice4.6.1_fixedDelay0.AIF_21_14_1.txt.full.txt';
% fid=fopen(f,'r');
% a=fscanf(fid,'%s',1);
% indexflow=1;
% while ~isempty(a)
%     a=fscanf(fid,'%s',1);
%     switch(a)
%         case {'Ktrans'}
%             a =fscanf(fid,'%s',1);
%             a =fscanf(fid,'%s',1);   % argh, must be a cleaner way!
%             aa=str2num(a);
%             Ktrans(indexflow)=aa; indexflow=indexflow+1;
%         otherwise
%     end
% end
% r4=Ktrans;
% T1_Ktrans_rest1=[r1 r2 r3 r4];T1_REST_avgSlice=[T1_REST_avgSlice mean(r1) mean(r2) mean(r3) mean(r4)];T1_REST_avgPt=[T1_REST_avgPt mean(T1_Ktrans_rest1)];
% T1_REST=[T1_REST T1_Ktrans_rest1];
% 
% figure(1),hold on,plot(DB_Ktrans_rest1,T1_Ktrans_rest1,'bo','LineWidth',3,'MarkerSize',8)
% disp('Plot held for current patient; press key to continue.'),pause













% %%%---------P080609---------%%% %%%% This subject is omitted, since GD was accidentally injected onto the floor during adenosine stress (rest values are okay)
% %%%% To load the Dual-Bolus Ktrans estimates
% f='/v/raid1/npack/Processing/P080609_nufft/Output/flowvalues.study21.slice1.6.1_fixedDelay0.AIF_16_3_5.txt.full.txt';
% fid=fopen(f,'r');
% a=fscanf(fid,'%s',1);
% indexflow=1;
% while ~isempty(a)
%     a=fscanf(fid,'%s',1);
%     switch(a)
%         case {'Ktrans'}
%             a =fscanf(fid,'%s',1);
%             a =fscanf(fid,'%s',1);   % argh, must be a cleaner way!
%             aa=str2num(a);
%             Ktrans(indexflow)=aa; indexflow=indexflow+1;
%         otherwise
%     end
% end
% r1=Ktrans;
% f='/v/raid1/npack/Processing/P080609_nufft/Output/flowvalues.study21.slice3.6.1_fixedDelay0.AIF_16_3_5.txt.full.txt';
% fid=fopen(f,'r');
% a=fscanf(fid,'%s',1);
% indexflow=1;
% while ~isempty(a)
%     a=fscanf(fid,'%s',1);
%     switch(a)
%         case {'Ktrans'}
%             a =fscanf(fid,'%s',1);
%             a =fscanf(fid,'%s',1);   % argh, must be a cleaner way!
%             aa=str2num(a);
%             Ktrans(indexflow)=aa; indexflow=indexflow+1;
%         otherwise
%     end
% end
% r3=Ktrans;
% DB_Ktrans_rest1=[r1 r3];DB_REST_avgSlice=[DB_REST_avgSlice mean(r1) mean(r3)];DB_REST_avgPt=[DB_REST_avgPt mean(DB_Ktrans_rest1)];
% DB_REST=[DB_REST DB_Ktrans_rest1];
% 
% 
% %%%% To load the T1 Ktrans estimates
% f='/v/raid1/npack/Processing/P080609_nufft_24rays/Output/flowvalues.study21.slice1.6.1_fixedDelay0.AIF_21_13_1.txt.full.txt';
% fid=fopen(f,'r');
% a=fscanf(fid,'%s',1);
% indexflow=1;
% while ~isempty(a)
%     a=fscanf(fid,'%s',1);
%     switch(a)
%         case {'Ktrans'}
%             a =fscanf(fid,'%s',1);
%             a =fscanf(fid,'%s',1);   % argh, must be a cleaner way!
%             aa=str2num(a);
%             Ktrans(indexflow)=aa; indexflow=indexflow+1;
%         otherwise
%     end
% end
% r1=Ktrans;
% f='/v/raid1/npack/Processing/P080609_nufft_24rays/Output/flowvalues.study21.slice3.6.1_fixedDelay0.AIF_21_13_1.txt.full.txt';
% fid=fopen(f,'r');
% a=fscanf(fid,'%s',1);
% indexflow=1;
% while ~isempty(a)
%     a=fscanf(fid,'%s',1);
%     switch(a)
%         case {'Ktrans'}
%             a =fscanf(fid,'%s',1);
%             a =fscanf(fid,'%s',1);   % argh, must be a cleaner way!
%             aa=str2num(a);
%             Ktrans(indexflow)=aa; indexflow=indexflow+1;
%         otherwise
%     end
% end
% r3=Ktrans;
% T1_Ktrans_rest1=[r1 r3];T1_REST_avgSlice=[T1_REST_avgSlice mean(r1) mean(r3)];T1_REST_avgPt=[T1_REST_avgPt mean(T1_Ktrans_rest1)];
% T1_REST=[T1_REST T1_Ktrans_rest1];
% 
% %%%% To load the Adenosine Ktrans estimates
% f='/v/raid1/npack/Processing/P080609_nufft_24rays/Output/flowvalues.study26.slice1.6.1_fixedDelay0.AIF_26_12_1.txt.full.txt';
% fid=fopen(f,'r');
% a=fscanf(fid,'%s',1);
% indexflow=1;
% while ~isempty(a)
%     a=fscanf(fid,'%s',1);
%     switch(a)
%         case {'Ktrans'}
%             a =fscanf(fid,'%s',1);
%             a =fscanf(fid,'%s',1);   % argh, must be a cleaner way!
%             aa=str2num(a);
%             Ktrans(indexflow)=aa; indexflow=indexflow+1;
%         otherwise
%     end
% end
% s1=Ktrans;
% f='/v/raid1/npack/Processing/P080609_nufft_24rays/Output/flowvalues.study26.slice2.6.1_fixedDelay0.AIF_26_12_1.txt.full.txt';
% fid=fopen(f,'r');
% a=fscanf(fid,'%s',1);
% indexflow=1;
% while ~isempty(a)
%     a=fscanf(fid,'%s',1);
%     switch(a)
%         case {'Ktrans'}
%             a =fscanf(fid,'%s',1);
%             a =fscanf(fid,'%s',1);   % argh, must be a cleaner way!
%             aa=str2num(a);
%             Ktrans(indexflow)=aa; indexflow=indexflow+1;
%         otherwise
%     end
% end
% s2=Ktrans;
% Aden_1=[s1 s2];ADEN_avgSlice=[ADEN_avgSlice mean(s1) mean(s2)];ADEN_avgPt=[ADEN_avgPt mean(Aden_1)];
% ADEN=[ADEN Aden_1];
% 
% %%%% To load the Lexiscan Ktrans estimates
% f='/v/raid1/npack/Processing/P080609_nufft_24rays/Output/flowvalues.study37.slice1.6.1_fixedDelay0.AIF_37_11_1.txt.full.txt';
% fid=fopen(f,'r');
% a=fscanf(fid,'%s',1);
% indexflow=1;
% while ~isempty(a)
%     a=fscanf(fid,'%s',1);
%     switch(a)
%         case {'Ktrans'}
%             a =fscanf(fid,'%s',1);
%             a =fscanf(fid,'%s',1);   % argh, must be a cleaner way!
%             aa=str2num(a);
%             Ktrans(indexflow)=aa; indexflow=indexflow+1;
%         otherwise
%     end
% end
% s1=Ktrans;
% f='/v/raid1/npack/Processing/P080609_nufft_24rays/Output/flowvalues.study37.slice2.6.1_fixedDelay0.AIF_37_11_1.txt.full.txt';
% fid=fopen(f,'r');
% a=fscanf(fid,'%s',1);
% indexflow=1;
% while ~isempty(a)
%     a=fscanf(fid,'%s',1);
%     switch(a)
%         case {'Ktrans'}
%             a =fscanf(fid,'%s',1);
%             a =fscanf(fid,'%s',1);   % argh, must be a cleaner way!
%             aa=str2num(a);
%             Ktrans(indexflow)=aa; indexflow=indexflow+1;
%         otherwise
%     end
% end
% s2=Ktrans;
% Lex_1=[s1 s2];LEX_avgSlice=[LEX_avgSlice mean(s1) mean(s2)];LEX_avgPt=[LEX_avgPt mean(Lex_1)];
% LEX=[LEX Lex_1];
% 
% REST2=[REST2 fliplr(T1_Ktrans_rest1)]; REST2_avgSlice=[REST2_avgSlice mean(T1_Ktrans_rest1(7:12)) mean(T1_Ktrans_rest1(1:6))];REST2_avgPt=[REST2_avgPt mean(T1_Ktrans_rest1)];% switched order to properly match stress slices
% figure(1),hold on,plot(DB_Ktrans_rest1,T1_Ktrans_rest1,'bo','LineWidth',3,'MarkerSize',8)
% disp('Plot held for current patient; press key to continue.'),pause












%%%---------P081009---------%%%
%%%% To load the Dual-Bolus Ktrans estimates
f='/v/raid1/npack/Processing/P081009_nufft/Output/flowvalues.study22.slice1.6.1_fixedDelay0.AIF_17_1_5.txt.full.txt';
Ktrans = getKtrans(f);
r1=Ktrans;
f='/v/raid1/npack/Processing/P081009_nufft/Output/flowvalues.study22.slice2.6.1_fixedDelay0.AIF_17_1_5.txt.full.txt';
Ktrans = getKtrans(f);
r2=Ktrans;
f='/v/raid1/npack/Processing/P081009_nufft/Output/flowvalues.study22.slice3.6.1_fixedDelay0.AIF_17_1_5.txt.full.txt';
Ktrans = getKtrans(f);
r3=Ktrans;
DB_Ktrans_rest1=[r1 r2 r3];DB_REST_avgSlice=[DB_REST_avgSlice mean(r1) mean(r2) mean(r3)];DB_REST_avgPt=[DB_REST_avgPt mean(DB_Ktrans_rest1)];
DB_REST=[DB_REST DB_Ktrans_rest1];


%%%% To load the T1 Ktrans estimates
f='/v/raid1/npack/Processing/P081009_nufft_24rays/Output/flowvalues.study22.slice1.6.1_fixedDelay0.AIF_22_11_1.txt.full.txt';
Ktrans = getKtrans(f);
r1=Ktrans;
f='/v/raid1/npack/Processing/P081009_nufft_24rays/Output/flowvalues.study22.slice2.6.1_fixedDelay0.AIF_22_11_1.txt.full.txt';
Ktrans = getKtrans(f);
r2=Ktrans;
f='/v/raid1/npack/Processing/P081009_nufft_24rays/Output/flowvalues.study22.slice3.6.1_fixedDelay0.AIF_22_11_1.txt.full.txt';
Ktrans = getKtrans(f);
r3=Ktrans;
T1_Ktrans_rest1=[r1 r2 r3];T1_REST_avgSlice=[T1_REST_avgSlice mean(r1) mean(r2) mean(r3)];T1_REST_avgPt=[T1_REST_avgPt mean(T1_Ktrans_rest1)];
T1_REST=[T1_REST T1_Ktrans_rest1];

%%%% To load the Adenosine Ktrans estimates
f='/v/raid1/npack/Processing/P081009_nufft_24rays/Output/flowvalues.study26.slice1.6.1_fixedDelay0.AIF_26_12_1.txt.full.txt';
Ktrans = getKtrans(f);
s1=Ktrans;
f='/v/raid1/npack/Processing/P081009_nufft_24rays/Output/flowvalues.study26.slice2.6.1_fixedDelay0.AIF_26_12_1.txt.full.txt';
Ktrans = getKtrans(f);
s2=Ktrans;
Aden_1=[s1 s2];ADEN_avgSlice=[ADEN_avgSlice mean(s1) mean(s2)];ADEN_avgPt=[ADEN_avgPt mean(Aden_1)];
ADEN=[ADEN Aden_1];

%%%% To load the Lexiscan Ktrans estimates
f='/v/raid1/npack/Processing/P081009_nufft_24rays/Output/flowvalues.study50.slice1.6.1_fixedDelay0.AIF_50_12_1.txt.full.txt';
Ktrans = getKtrans(f);
s1=Ktrans;
f='/v/raid1/npack/Processing/P081009_nufft_24rays/Output/flowvalues.study50.slice2.6.1_fixedDelay0.AIF_50_12_1.txt.full.txt';
Ktrans = getKtrans(f);
s2=Ktrans;
Lex_1=[s1 s2];LEX_avgSlice=[LEX_avgSlice mean(s1) mean(s2)];LEX_avgPt=[LEX_avgPt mean(Lex_1)];
LEX=[LEX Lex_1];

REST2=[REST2 T1_Ktrans_rest1(13:18) T1_Ktrans_rest1(1:6)];  REST2_avgSlice=[REST2_avgSlice mean(T1_Ktrans_rest1(13:18)) mean(T1_Ktrans_rest1(1:6))];REST2_avgPt=[REST2_avgPt (mean(T1_Ktrans_rest1(13:18)) + mean(T1_Ktrans_rest1(1:6)))/2];
figure(1),hold on,plot(DB_Ktrans_rest1,T1_Ktrans_rest1,'bo','LineWidth',3,'MarkerSize',8)
disp('Plot held for current patient; press key to continue.'),pause(.1)













%%%---------P081309---------%%%
%%%% To load the Dual-Bolus Ktrans estimates
f='/v/raid1/npack/Processing/P081309_nufft/Output/flowvalues.study23.slice1.6.1_fixedDelay0.AIF_18_1_5.txt.full.txt';
Ktrans = getKtrans(f);
r1=Ktrans;
f='/v/raid1/npack/Processing/P081309_nufft/Output/flowvalues.study23.slice2.6.1_fixedDelay0.AIF_18_1_5.txt.full.txt';
Ktrans = getKtrans(f);
r2=Ktrans;
f='/v/raid1/npack/Processing/P081309_nufft/Output/flowvalues.study23.slice3.6.1_fixedDelay0.AIF_18_1_5.txt.full.txt';
Ktrans = getKtrans(f);
r3=Ktrans;
DB_Ktrans_rest1=[r1 r2 r3];DB_REST_avgSlice=[DB_REST_avgSlice mean(r1) mean(r2) mean(r3)];DB_REST_avgPt=[DB_REST_avgPt mean(DB_Ktrans_rest1)];
DB_REST=[DB_REST DB_Ktrans_rest1];


%%%% To load the T1 Ktrans estimates
f='/v/raid1/npack/Processing/P081309_nufft_24rays/Output/flowvalues.study23.slice1.6.1_fixedDelay0.AIF_23_11_1.txt.full.txt';
Ktrans = getKtrans(f);
r1=Ktrans;
f='/v/raid1/npack/Processing/P081309_nufft_24rays/Output/flowvalues.study23.slice2.6.1_fixedDelay0.AIF_23_11_1.txt.full.txt';
Ktrans = getKtrans(f);
r2=Ktrans;
f='/v/raid1/npack/Processing/P081309_nufft_24rays/Output/flowvalues.study23.slice3.6.1_fixedDelay0.AIF_23_11_1.txt.full.txt';
Ktrans = getKtrans(f);
r3=Ktrans;
T1_Ktrans_rest1=[r1 r2 r3];T1_REST_avgSlice=[T1_REST_avgSlice mean(r1) mean(r2) mean(r3)];T1_REST_avgPt=[T1_REST_avgPt mean(T1_Ktrans_rest1)];
T1_REST=[T1_REST T1_Ktrans_rest1];

%%%% To load the Adenosine Ktrans estimates
f='/v/raid1/npack/Processing/P081309_nufft_24rays/Output/flowvalues.study28.slice1.6.1_fixedDelay0.AIF_28_11_1.txt.full.txt';
Ktrans = getKtrans(f);
s1=Ktrans;
f='/v/raid1/npack/Processing/P081309_nufft_24rays/Output/flowvalues.study28.slice2.6.1_fixedDelay0.AIF_28_11_1.txt.full.txt';
Ktrans = getKtrans(f);
s2=Ktrans;
f='/v/raid1/npack/Processing/P081309_nufft_24rays/Output/flowvalues.study28.slice3.6.1_fixedDelay0.AIF_28_11_1.txt.full.txt';
Ktrans = getKtrans(f);
s3=Ktrans;
Aden_1=[s1 s2 s3];ADEN_avgSlice=[ADEN_avgSlice mean(s1) mean(s2) mean(s3)];ADEN_avgPt=[ADEN_avgPt mean(Aden_1)];
ADEN=[ADEN Aden_1];

%%%% To load the Lexiscan Ktrans estimates
f='/v/raid1/npack/Processing/P081309_nufft_24rays/Output/flowvalues.study55.slice1.6.1_fixedDelay0.AIF_55_11_1.txt.full.txt';
Ktrans = getKtrans(f);
s1=Ktrans;
f='/v/raid1/npack/Processing/P081309_nufft_24rays/Output/flowvalues.study55.slice2.6.1_fixedDelay0.AIF_55_11_1.txt.full.txt';
Ktrans = getKtrans(f);
s2=Ktrans;
f='/v/raid1/npack/Processing/P081309_nufft_24rays/Output/flowvalues.study55.slice3.6.1_fixedDelay0.AIF_55_11_1.txt.full.txt';
Ktrans = getKtrans(f);
s3=Ktrans;
Lex_1=[s1 s2 s3];LEX_avgSlice=[LEX_avgSlice mean(s1) mean(s2) mean(s3)];LEX_avgPt=[LEX_avgPt mean(Lex_1)];
LEX=[LEX Lex_1];

REST2=[REST2 T1_Ktrans_rest1];REST2_avgSlice=[REST2_avgSlice mean(T1_Ktrans_rest1(1:6)) mean(T1_Ktrans_rest1(7:12)) mean(T1_Ktrans_rest1(13:18))];REST2_avgPt=[REST2_avgPt mean(T1_Ktrans_rest1)];
figure(1),hold on,plot(DB_Ktrans_rest1,T1_Ktrans_rest1,'bo','LineWidth',3,'MarkerSize',8)
disp('Plot held for current patient; press key to continue.'),pause(.1)












%%%---------P081809---------%%%
%%%% To load the Dual-Bolus Ktrans estimates
f='/v/raid1/npack/Processing/P081809_nufft/Output/flowvalues.study25.slice1.6.1_fixedDelay0.AIF_20_1_5.txt.full.txt';
Ktrans = getKtrans(f);
r1=Ktrans;
f='/v/raid1/npack/Processing/P081809_nufft/Output/flowvalues.study25.slice2.6.1_fixedDelay0.AIF_20_1_5.txt.full.txt';
Ktrans = getKtrans(f);
r2=Ktrans;
f='/v/raid1/npack/Processing/P081809_nufft/Output/flowvalues.study25.slice3.6.1_fixedDelay0.AIF_20_1_5.txt.full.txt';
Ktrans = getKtrans(f);
r3=Ktrans;
DB_Ktrans_rest1=[r1 r2 r3];DB_REST_avgSlice=[DB_REST_avgSlice mean(r1) mean(r2) mean(r3)];DB_REST_avgPt=[DB_REST_avgPt mean(DB_Ktrans_rest1)];
DB_REST=[DB_REST DB_Ktrans_rest1];


%%%% To load the T1 Ktrans estimates
f='/v/raid1/npack/Processing/P081809_nufft_24rays/Output/flowvalues.study25.slice1.6.1_fixedDelay0.AIF_25_11_1.txt.full.txt';
Ktrans = getKtrans(f);
r1=Ktrans;
f='/v/raid1/npack/Processing/P081809_nufft_24rays/Output/flowvalues.study25.slice2.6.1_fixedDelay0.AIF_25_11_1.txt.full.txt';
Ktrans = getKtrans(f);
r2=Ktrans;
f='/v/raid1/npack/Processing/P081809_nufft_24rays/Output/flowvalues.study25.slice3.6.1_fixedDelay0.AIF_25_11_1.txt.full.txt';
Ktrans = getKtrans(f);
r3=Ktrans;
T1_Ktrans_rest1=[r1 r2 r3];T1_REST_avgSlice=[T1_REST_avgSlice mean(r1) mean(r2) mean(r3)];T1_REST_avgPt=[T1_REST_avgPt mean(T1_Ktrans_rest1)];
T1_REST=[T1_REST T1_Ktrans_rest1];

%%%% To load the Adenosine Ktrans estimates
f='/v/raid1/npack/Processing/P081809_nufft_24rays/Output/flowvalues.study29.slice1.6.1_fixedDelay0.AIF_29_11_1.txt.full.txt';
Ktrans = getKtrans(f);
s1=Ktrans;
f='/v/raid1/npack/Processing/P081809_nufft_24rays/Output/flowvalues.study29.slice2.6.1_fixedDelay0.AIF_29_11_1.txt.full.txt';
Ktrans = getKtrans(f);
s2=Ktrans;
f='/v/raid1/npack/Processing/P081809_nufft_24rays/Output/flowvalues.study29.slice3.6.1_fixedDelay0.AIF_29_11_1.txt.full.txt';
Ktrans = getKtrans(f);
s3=Ktrans;
Aden_1=[s1 s2 s3];ADEN_avgSlice=[ADEN_avgSlice mean(s1) mean(s2) mean(s3)];ADEN_avgPt=[ADEN_avgPt mean(Aden_1)];
ADEN=[ADEN Aden_1];

%%%% To load the Lexiscan Ktrans estimates
f='/v/raid1/npack/Processing/P081809_nufft_24rays/Output/flowvalues.study56.slice1.6.1_fixedDelay0.AIF_56_11_1.txt.full.txt';
Ktrans = getKtrans(f);
s1=Ktrans;
f='/v/raid1/npack/Processing/P081809_nufft_24rays/Output/flowvalues.study56.slice2.6.1_fixedDelay0.AIF_56_11_1.txt.full.txt';
Ktrans = getKtrans(f);
s2=Ktrans;
f='/v/raid1/npack/Processing/P081809_nufft_24rays/Output/flowvalues.study56.slice3.6.1_fixedDelay0.AIF_56_11_1.txt.full.txt';
Ktrans = getKtrans(f);
s3=Ktrans;
Lex_1=[s1 s2 s3];LEX_avgSlice=[LEX_avgSlice mean(s1) mean(s2) mean(s3)];LEX_avgPt=[LEX_avgPt mean(Lex_1)];
LEX=[LEX Lex_1];

REST2=[REST2 T1_Ktrans_rest1];REST2_avgSlice=[REST2_avgSlice mean(T1_Ktrans_rest1(1:6)) mean(T1_Ktrans_rest1(7:12)) mean(T1_Ktrans_rest1(13:18))];REST2_avgPt=[REST2_avgPt mean(T1_Ktrans_rest1)];
figure(1),hold on,plot(DB_Ktrans_rest1,T1_Ktrans_rest1,'bo','LineWidth',3,'MarkerSize',8)
disp('Plot held for current patient; press key to continue.'),pause(.1)












%%%---------P082609, 1st low-dose scan---------%%%
%%%% To load the Dual-Bolus Ktrans estimates
f='/v/raid1/npack/Processing/P082609_nufft/Output/flowvalues.study33.slice1.6.1_fixedDelay0.AIF_28_3_5.txt.full.txt';
Ktrans = getKtrans(f);
r1=Ktrans;
f='/v/raid1/npack/Processing/P082609_nufft/Output/flowvalues.study33.slice2.6.1_fixedDelay0.AIF_28_3_5.txt.full.txt';
Ktrans = getKtrans(f);
r2=Ktrans;
f='/v/raid1/npack/Processing/P082609_nufft/Output/flowvalues.study33.slice3.6.1_fixedDelay0.AIF_28_3_5.txt.full.txt';
Ktrans = getKtrans(f);
r3=Ktrans;
DB_Ktrans_rest1=[r1 r2 r3];DB_REST_avgSlice=[DB_REST_avgSlice mean(r1) mean(r2) mean(r3)];DB_REST_avgPt=[DB_REST_avgPt mean(DB_Ktrans_rest1)];
DB_REST=[DB_REST DB_Ktrans_rest1];


%%%% To load the T1 Ktrans estimates
f='/v/raid1/npack/Processing/P082609_nufft_24rays/Output/flowvalues.study33.slice1.6.1_fixedDelay0.AIF_33_13_1.txt.full.txt';
Ktrans = getKtrans(f);
r1=Ktrans;
f='/v/raid1/npack/Processing/P082609_nufft_24rays/Output/flowvalues.study33.slice2.6.1_fixedDelay0.AIF_33_13_1.txt.full.txt';
Ktrans = getKtrans(f);
r2=Ktrans;
f='/v/raid1/npack/Processing/P082609_nufft_24rays/Output/flowvalues.study33.slice3.6.1_fixedDelay0.AIF_33_13_1.txt.full.txt';
Ktrans = getKtrans(f);
r3=Ktrans;
T1_Ktrans_rest1=[r1 r2 r3];T1_REST_avgSlice=[T1_REST_avgSlice mean(r1) mean(r2) mean(r3)];T1_REST_avgPt=[T1_REST_avgPt mean(T1_Ktrans_rest1)];
T1_REST=[T1_REST T1_Ktrans_rest1];

figure(1),hold on,plot(DB_Ktrans_rest1,T1_Ktrans_rest1,'bo','LineWidth',3,'MarkerSize',8)
disp('Plot held for current patient; press key to continue.'),pause(.1)










%%%---------P082609, 2nd low-dose scan---------%%%
%%%% To load the Dual-Bolus Ktrans estimates
f='/v/raid1/npack/Processing/P082609_nufft/Output/flowvalues.study33.slice1.6.1_fixedDelay0.AIF_38_1_5.txt.full.txt';
Ktrans = getKtrans(f);
r1=Ktrans;
f='/v/raid1/npack/Processing/P082609_nufft/Output/flowvalues.study33.slice2.6.1_fixedDelay0.AIF_38_1_5.txt.full.txt';
Ktrans = getKtrans(f);
r2=Ktrans;
f='/v/raid1/npack/Processing/P082609_nufft/Output/flowvalues.study33.slice3.6.1_fixedDelay0.AIF_38_1_5.txt.full.txt';
Ktrans = getKtrans(f);
r3=Ktrans;
DB_Ktrans_rest1=[r1 r2 r3];DB_REST_avgSlice=[DB_REST_avgSlice mean(r1) mean(r2) mean(r3)];DB_REST_avgPt=[DB_REST_avgPt mean(DB_Ktrans_rest1)];
DB_REST=[DB_REST DB_Ktrans_rest1];


%%%% To load the T1 Ktrans estimates
f='/v/raid1/npack/Processing/P082609_nufft_24rays/Output/flowvalues.study33.slice1.6.1_fixedDelay0.AIF_33_11_1.txt.full.txt';
Ktrans = getKtrans(f);
r1=Ktrans;
f='/v/raid1/npack/Processing/P082609_nufft_24rays/Output/flowvalues.study33.slice2.6.1_fixedDelay0.AIF_33_11_1.txt.full.txt';
Ktrans = getKtrans(f);
r2=Ktrans;
f='/v/raid1/npack/Processing/P082609_nufft_24rays/Output/flowvalues.study33.slice3.6.1_fixedDelay0.AIF_33_11_1.txt.full.txt';
Ktrans = getKtrans(f);
r3=Ktrans;
T1_Ktrans_rest1=[r1 r2 r3];T1_REST_avgSlice=[T1_REST_avgSlice mean(r1) mean(r2) mean(r3)];T1_REST_avgPt=[T1_REST_avgPt mean(T1_Ktrans_rest1)];
T1_REST=[T1_REST T1_Ktrans_rest1];

figure(1),hold on,plot(DB_Ktrans_rest1,T1_Ktrans_rest1,'bo','LineWidth',3,'MarkerSize',8)
disp('Plot held for current patient; press key to continue.'),pause(.1)










% %%%---------P082609, 0.06 high-dose scan---------%%%
% %%%% To load the Dual-Bolus Ktrans estimates
% f='/v/raid1/npack/Processing/P082609_nufft/Output/flowvalues.study43.slice1.6.1_fixedDelay0.AIF_28_1_15.txt.full.txt';
% fid=fopen(f,'r');
% a=fscanf(fid,'%s',1);
% indexflow=1;
% while ~isempty(a)
%     a=fscanf(fid,'%s',1);
%     switch(a)
%         case {'Ktrans'}
%             a =fscanf(fid,'%s',1);
%             a =fscanf(fid,'%s',1);   % argh, must be a cleaner way!
%             aa=str2num(a);
%             Ktrans(indexflow)=aa; indexflow=indexflow+1;
%         otherwise
%     end
% end
% r1=Ktrans;
% f='/v/raid1/npack/Processing/P082609_nufft/Output/flowvalues.study43.slice2.6.1_fixedDelay0.AIF_28_1_15.txt.full.txt';
% fid=fopen(f,'r');
% a=fscanf(fid,'%s',1);
% indexflow=1;
% while ~isempty(a)
%     a=fscanf(fid,'%s',1);
%     switch(a)
%         case {'Ktrans'}
%             a =fscanf(fid,'%s',1);
%             a =fscanf(fid,'%s',1);   % argh, must be a cleaner way!
%             aa=str2num(a);
%             Ktrans(indexflow)=aa; indexflow=indexflow+1;
%         otherwise
%     end
% end
% r2=Ktrans;
% f='/v/raid1/npack/Processing/P082609_nufft/Output/flowvalues.study43.slice3.6.1_fixedDelay0.AIF_28_1_15.txt.full.txt';
% fid=fopen(f,'r');
% a=fscanf(fid,'%s',1);
% indexflow=1;
% while ~isempty(a)
%     a=fscanf(fid,'%s',1);
%     switch(a)
%         case {'Ktrans'}
%             a =fscanf(fid,'%s',1);
%             a =fscanf(fid,'%s',1);   % argh, must be a cleaner way!
%             aa=str2num(a);
%             Ktrans(indexflow)=aa; indexflow=indexflow+1;
%         otherwise
%     end
% end
% r3=Ktrans;
% DB_Ktrans_rest1=[r1 r2 r3];DB_REST_avgSlice=[DB_REST_avgSlice mean(r1) mean(r2) mean(r3)];DB_REST_avgPt=[DB_REST_avgPt mean(DB_Ktrans_rest1)];
% DB_REST=[DB_REST DB_Ktrans_rest1];
% 
% 
% %%%% To load the T1 Ktrans estimates
% f='/v/raid1/npack/Processing/P082609_nufft_24rays/Output/flowvalues.study43.slice1.6.1_fixedDelay0.AIF_43_11_1.txt.full.txt';
% fid=fopen(f,'r');
% a=fscanf(fid,'%s',1);
% indexflow=1;
% while ~isempty(a)
%     a=fscanf(fid,'%s',1);
%     switch(a)
%         case {'Ktrans'}
%             a =fscanf(fid,'%s',1);
%             a =fscanf(fid,'%s',1);   % argh, must be a cleaner way!
%             aa=str2num(a);
%             Ktrans(indexflow)=aa; indexflow=indexflow+1;
%         otherwise
%     end
% end
% r1=Ktrans;
% f='/v/raid1/npack/Processing/P082609_nufft_24rays/Output/flowvalues.study43.slice2.6.1_fixedDelay0.AIF_43_11_1.txt.full.txt';
% fid=fopen(f,'r');
% a=fscanf(fid,'%s',1);
% indexflow=1;
% while ~isempty(a)
%     a=fscanf(fid,'%s',1);
%     switch(a)
%         case {'Ktrans'}
%             a =fscanf(fid,'%s',1);
%             a =fscanf(fid,'%s',1);   % argh, must be a cleaner way!
%             aa=str2num(a);
%             Ktrans(indexflow)=aa; indexflow=indexflow+1;
%         otherwise
%     end
% end
% r2=Ktrans;
% f='/v/raid1/npack/Processing/P082609_nufft_24rays/Output/flowvalues.study43.slice3.6.1_fixedDelay0.AIF_43_11_1.txt.full.txt';
% fid=fopen(f,'r');
% a=fscanf(fid,'%s',1);
% indexflow=1;
% while ~isempty(a)
%     a=fscanf(fid,'%s',1);
%     switch(a)
%         case {'Ktrans'}
%             a =fscanf(fid,'%s',1);
%             a =fscanf(fid,'%s',1);   % argh, must be a cleaner way!
%             aa=str2num(a);
%             Ktrans(indexflow)=aa; indexflow=indexflow+1;
%         otherwise
%     end
% end
% r3=Ktrans;
% T1_Ktrans_rest1=[r1 r2 r3];T1_REST_avgSlice=[T1_REST_avgSlice mean(r1) mean(r2) mean(r3)];T1_REST_avgPt=[T1_REST_avgPt mean(T1_Ktrans_rest1)];
% T1_REST=[T1_REST T1_Ktrans_rest1];
% 
% figure(1),hold on,plot(DB_Ktrans_rest1,T1_Ktrans_rest1,'rs','LineWidth',3,'MarkerSize',8)
% disp('Plot held for current patient; press key to continue.'),pause










%%%---------P082809---------%%%
%%%% To load the Dual-Bolus Ktrans estimates
f='/v/raid1/npack/Processing/P082809_nufft/Output/flowvalues.study35.slice1.6.1_fixedDelay0.AIF_25_1_5.txt.full.txt';
Ktrans = getKtrans(f);
r1=Ktrans;
f='/v/raid1/npack/Processing/P082809_nufft/Output/flowvalues.study35.slice2.6.1_fixedDelay0.AIF_25_1_5.txt.full.txt';
Ktrans = getKtrans(f);
r2=Ktrans;
f='/v/raid1/npack/Processing/P082809_nufft/Output/flowvalues.study35.slice3.6.1_fixedDelay0.AIF_25_1_5.txt.full.txt';
Ktrans = getKtrans(f);
r3=Ktrans;
DB_Ktrans_rest1=[r1 r2 r3];DB_REST_avgSlice=[DB_REST_avgSlice mean(r1) mean(r2) mean(r3)];DB_REST_avgPt=[DB_REST_avgPt mean(DB_Ktrans_rest1)];
DB_REST=[DB_REST DB_Ktrans_rest1];


%%%% To load the T1 Ktrans estimates
f='/v/raid1/npack/Processing/P082809_nufft_24rays/Output/flowvalues.study35.slice1.6.1_fixedDelay0.AIF_35_11_1.txt.full.txt';
Ktrans = getKtrans(f);
r1=Ktrans;
f='/v/raid1/npack/Processing/P082809_nufft_24rays/Output/flowvalues.study35.slice2.6.1_fixedDelay0.AIF_35_11_1.txt.full.txt';
Ktrans = getKtrans(f);
r2=Ktrans;
f='/v/raid1/npack/Processing/P082809_nufft_24rays/Output/flowvalues.study35.slice3.6.1_fixedDelay0.AIF_35_11_1.txt.full.txt';
Ktrans = getKtrans(f);
r3=Ktrans;
T1_Ktrans_rest1=[r1 r2 r3];T1_REST_avgSlice=[T1_REST_avgSlice mean(r1) mean(r2) mean(r3)];T1_REST_avgPt=[T1_REST_avgPt mean(T1_Ktrans_rest1)];
T1_REST=[T1_REST T1_Ktrans_rest1];

figure(1),hold on,plot(DB_Ktrans_rest1,T1_Ktrans_rest1,'bo','LineWidth',3,'MarkerSize',8)
disp('Plot held for current patient; press key to continue.'),pause(.1)










% %%%---------P082809,  0.06 high-dose injection---------%%%
% %%%% To load the Dual-Bolus Ktrans estimates
% f='/v/raid1/npack/Processing/P082809_nufft/Output/flowvalues.study39.slice1.6.1_fixedDelay0.AIF_25_1_15.txt.full.txt';
% fid=fopen(f,'r');
% a=fscanf(fid,'%s',1);
% indexflow=1;
% while ~isempty(a)
%     a=fscanf(fid,'%s',1);
%     switch(a)
%         case {'Ktrans'}
%             a =fscanf(fid,'%s',1);
%             a =fscanf(fid,'%s',1);   % argh, must be a cleaner way!
%             aa=str2num(a);
%             Ktrans(indexflow)=aa; indexflow=indexflow+1;
%         otherwise
%     end
% end
% r1=Ktrans;
% f='/v/raid1/npack/Processing/P082809_nufft/Output/flowvalues.study39.slice2.6.1_fixedDelay0.AIF_25_1_15.txt.full.txt';
% fid=fopen(f,'r');
% a=fscanf(fid,'%s',1);
% indexflow=1;
% while ~isempty(a)
%     a=fscanf(fid,'%s',1);
%     switch(a)
%         case {'Ktrans'}
%             a =fscanf(fid,'%s',1);
%             a =fscanf(fid,'%s',1);   % argh, must be a cleaner way!
%             aa=str2num(a);
%             Ktrans(indexflow)=aa; indexflow=indexflow+1;
%         otherwise
%     end
% end
% r2=Ktrans;
% f='/v/raid1/npack/Processing/P082809_nufft/Output/flowvalues.study39.slice3.6.1_fixedDelay0.AIF_25_1_15.txt.full.txt';
% fid=fopen(f,'r');
% a=fscanf(fid,'%s',1);
% indexflow=1;
% while ~isempty(a)
%     a=fscanf(fid,'%s',1);
%     switch(a)
%         case {'Ktrans'}
%             a =fscanf(fid,'%s',1);
%             a =fscanf(fid,'%s',1);   % argh, must be a cleaner way!
%             aa=str2num(a);
%             Ktrans(indexflow)=aa; indexflow=indexflow+1;
%         otherwise
%     end
% end
% r3=Ktrans;
% DB_Ktrans_rest1=[r1 r2 r3];DB_REST_avgSlice=[DB_REST_avgSlice mean(r1) mean(r2) mean(r3)];DB_REST_avgPt=[DB_REST_avgPt mean(DB_Ktrans_rest1)];
% DB_REST=[DB_REST DB_Ktrans_rest1];
% 
% 
% %%%% To load the T1 Ktrans estimates
% f='/v/raid1/npack/Processing/P082809_nufft_24rays/Output/flowvalues.study39.slice1.6.1_fixedDelay0.AIF_39_11_1.txt.full.txt';
% fid=fopen(f,'r');
% a=fscanf(fid,'%s',1);
% indexflow=1;
% while ~isempty(a)
%     a=fscanf(fid,'%s',1);
%     switch(a)
%         case {'Ktrans'}
%             a =fscanf(fid,'%s',1);
%             a =fscanf(fid,'%s',1);   % argh, must be a cleaner way!
%             aa=str2num(a);
%             Ktrans(indexflow)=aa; indexflow=indexflow+1;
%         otherwise
%     end
% end
% r1=Ktrans;
% f='/v/raid1/npack/Processing/P082809_nufft_24rays/Output/flowvalues.study39.slice2.6.1_fixedDelay0.AIF_39_11_1.txt.full.txt';
% fid=fopen(f,'r');
% a=fscanf(fid,'%s',1);
% indexflow=1;
% while ~isempty(a)
%     a=fscanf(fid,'%s',1);
%     switch(a)
%         case {'Ktrans'}
%             a =fscanf(fid,'%s',1);
%             a =fscanf(fid,'%s',1);   % argh, must be a cleaner way!
%             aa=str2num(a);
%             Ktrans(indexflow)=aa; indexflow=indexflow+1;
%         otherwise
%     end
% end
% r2=Ktrans;
% f='/v/raid1/npack/Processing/P082809_nufft_24rays/Output/flowvalues.study39.slice3.6.1_fixedDelay0.AIF_39_11_1.txt.full.txt';
% fid=fopen(f,'r');
% a=fscanf(fid,'%s',1);
% indexflow=1;
% while ~isempty(a)
%     a=fscanf(fid,'%s',1);
%     switch(a)
%         case {'Ktrans'}
%             a =fscanf(fid,'%s',1);
%             a =fscanf(fid,'%s',1);   % argh, must be a cleaner way!
%             aa=str2num(a);
%             Ktrans(indexflow)=aa; indexflow=indexflow+1;
%         otherwise
%     end
% end
% r3=Ktrans;
% T1_Ktrans_rest1=[r1 r2 r3];T1_REST_avgSlice=[T1_REST_avgSlice mean(r1) mean(r2) mean(r3)];T1_REST_avgPt=[T1_REST_avgPt mean(T1_Ktrans_rest1)];
% T1_REST=[T1_REST T1_Ktrans_rest1];
% 
% figure(1),hold on,plot(DB_Ktrans_rest1,T1_Ktrans_rest1,'rs','LineWidth',3,'MarkerSize',8)
% disp('Plot held for current patient; press key to continue.'),pause














%%%---------P090309---------%%%
%%%% To load the Dual-Bolus Ktrans estimates
f='/v/raid1/npack/Processing/P090309_nufft/Output/flowvalues.study27.slice1.6.1_fixedDelay0.AIF_22_1_5.txt.full.txt';
Ktrans = getKtrans(f);
r1=Ktrans;
f='/v/raid1/npack/Processing/P090309_nufft/Output/flowvalues.study27.slice2.6.1_fixedDelay0.AIF_22_1_5.txt.full.txt';
Ktrans = getKtrans(f);
r2=Ktrans;
f='/v/raid1/npack/Processing/P090309_nufft/Output/flowvalues.study27.slice3.6.1_fixedDelay0.AIF_22_1_5.txt.full.txt';
Ktrans = getKtrans(f);
r3=Ktrans;
DB_Ktrans_rest1=[r1 r2 r3];DB_REST_avgSlice=[DB_REST_avgSlice mean(r1) mean(r2) mean(r3)];DB_REST_avgPt=[DB_REST_avgPt mean(DB_Ktrans_rest1)];
DB_REST=[DB_REST DB_Ktrans_rest1];


%%%% To load the T1 Ktrans estimates
f='/v/raid1/npack/Processing/P090309_nufft_24rays/Output/flowvalues.study27.slice1.6.1_fixedDelay0.AIF_27_11_1.txt.full.txt';
Ktrans = getKtrans(f);
r1=Ktrans;
f='/v/raid1/npack/Processing/P090309_nufft_24rays/Output/flowvalues.study27.slice2.6.1_fixedDelay0.AIF_27_11_1.txt.full.txt';
Ktrans = getKtrans(f);
r2=Ktrans;
f='/v/raid1/npack/Processing/P090309_nufft_24rays/Output/flowvalues.study27.slice3.6.1_fixedDelay0.AIF_27_11_1.txt.full.txt';
Ktrans = getKtrans(f);
r3=Ktrans;
T1_Ktrans_rest1=[r1 r2 r3];T1_REST_avgSlice=[T1_REST_avgSlice mean(r1) mean(r2) mean(r3)];T1_REST_avgPt=[T1_REST_avgPt mean(T1_Ktrans_rest1)];
T1_REST=[T1_REST T1_Ktrans_rest1];

%%%% To load the Adenosine Ktrans estimates
f='/v/raid1/npack/Processing/P090309_nufft_24rays/Output/flowvalues.study31.slice1.6.1_fixedDelay0.AIF_31_12_1.txt.full.txt';
Ktrans = getKtrans(f);
s1=Ktrans;
f='/v/raid1/npack/Processing/P090309_nufft_24rays/Output/flowvalues.study31.slice2.6.1_fixedDelay0.AIF_31_12_1.txt.full.txt';
Ktrans = getKtrans(f);
s2=Ktrans;
Aden_1=[s1 s2];ADEN_avgSlice=[ADEN_avgSlice mean(s1) mean(s2)];ADEN_avgPt=[ADEN_avgPt mean(Aden_1)]; %sl 1 removed to PVCs w/ bad gating
ADEN=[ADEN Aden_1];

%%%% To load the Lexiscan Ktrans estimates
f='/v/raid1/npack/Processing/P090309_nufft_24rays/Output/flowvalues.study64.slice1.6.1_fixedDelay0.AIF_64_12_1.txt.full.txt';
Ktrans = getKtrans(f);
s1=Ktrans;
f='/v/raid1/npack/Processing/P090309_nufft_24rays/Output/flowvalues.study64.slice2.6.1_fixedDelay0.AIF_64_12_1.txt.full.txt';
Ktrans = getKtrans(f);
s2=Ktrans;
Lex_1=[s1 s2];LEX_avgSlice=[LEX_avgSlice mean(s1) mean(s2)];LEX_avgPt=[LEX_avgPt mean(Lex_1)]; %sl 1 removed to PVCs w/ bad gating
LEX=[LEX Lex_1];

REST2=[REST2 (T1_Ktrans_rest1(13:18)) (T1_Ktrans_rest1(7:12)) ];REST2_avgSlice=[REST2_avgSlice mean(T1_Ktrans_rest1(13:18)) mean(T1_Ktrans_rest1(7:12)) ];REST2_avgPt=[REST2_avgPt mean(T1_Ktrans_rest1(7:18))]; %sl 2 removed to PVCs w/ bad gating
figure(1),hold on,plot(DB_Ktrans_rest1,T1_Ktrans_rest1,'bo','LineWidth',3,'MarkerSize',8)
disp('Plot held for current patient; press key to continue.'),pause(.1)










%%%---------P092909,  low-dose 1/5*0.02 injection---------%%%
%%%% To load the Dual-Bolus Ktrans estimates
f='/v/raid1/npack/Processing/P092909_nufft/Output/flowvalues.study55.slice1.6.1_fixedDelay0.AIF_51_1_5.txt.full.txt';
Ktrans = getKtrans(f);
r1=Ktrans;
f='/v/raid1/npack/Processing/P092909_nufft/Output/flowvalues.study55.slice2.6.1_fixedDelay0.AIF_51_1_5.txt.full.txt';
Ktrans = getKtrans(f);
r2=Ktrans;
f='/v/raid1/npack/Processing/P092909_nufft/Output/flowvalues.study55.slice3.6.1_fixedDelay0.AIF_51_1_5.txt.full.txt';
Ktrans = getKtrans(f);
r3=Ktrans;
DB_Ktrans_rest1=[r1 r2 r3];DB_REST_avgSlice=[DB_REST_avgSlice mean(r1) mean(r2) mean(r3)];DB_REST_avgPt=[DB_REST_avgPt mean(DB_Ktrans_rest1)];
DB_REST=[DB_REST DB_Ktrans_rest1];


%%%% To load the T1 Ktrans estimates
f='/v/raid1/npack/Processing/P092909_nufft_24rays/Output/flowvalues.study55.slice1.6.1_fixedDelay0.AIF_55_11_1.txt.full.txt';
Ktrans = getKtrans(f);
r1=Ktrans;
f='/v/raid1/npack/Processing/P092909_nufft_24rays/Output/flowvalues.study55.slice2.6.1_fixedDelay0.AIF_55_11_1.txt.full.txt';
Ktrans = getKtrans(f);
r2=Ktrans;
f='/v/raid1/npack/Processing/P092909_nufft_24rays/Output/flowvalues.study55.slice3.6.1_fixedDelay0.AIF_55_11_1.txt.full.txt';
Ktrans = getKtrans(f);
r3=Ktrans;
T1_Ktrans_rest1=[r1 r2 r3];T1_REST_avgSlice=[T1_REST_avgSlice mean(r1) mean(r2) mean(r3)];T1_REST_avgPt=[T1_REST_avgPt mean(T1_Ktrans_rest1)];
T1_REST=[T1_REST T1_Ktrans_rest1];

figure(1),hold on,plot(DB_Ktrans_rest1,T1_Ktrans_rest1,'bo','LineWidth',3,'MarkerSize',8)
disp('Plot held for current patient; press key to continue.'),pause(.1)










% %%%---------P092909,  0.06 high-dose injection---------%%%
% %%%% To load the Dual-Bolus Ktrans estimates
% f='/v/raid1/npack/Processing/P092909_nufft/Output/flowvalues.study59.slice1.6.1_fixedDelay0.AIF_51_1_15.txt.full.txt';
% fid=fopen(f,'r');
% a=fscanf(fid,'%s',1);
% indexflow=1;
% while ~isempty(a)
%     a=fscanf(fid,'%s',1);
%     switch(a)
%         case {'Ktrans'}
%             a =fscanf(fid,'%s',1);
%             a =fscanf(fid,'%s',1);   % argh, must be a cleaner way!
%             aa=str2num(a);
%             Ktrans(indexflow)=aa; indexflow=indexflow+1;
%         otherwise
%     end
% end
% r1=Ktrans;
% f='/v/raid1/npack/Processing/P092909_nufft/Output/flowvalues.study59.slice2.6.1_fixedDelay0.AIF_51_1_15.txt.full.txt';
% fid=fopen(f,'r');
% a=fscanf(fid,'%s',1);
% indexflow=1;
% while ~isempty(a)
%     a=fscanf(fid,'%s',1);
%     switch(a)
%         case {'Ktrans'}
%             a =fscanf(fid,'%s',1);
%             a =fscanf(fid,'%s',1);   % argh, must be a cleaner way!
%             aa=str2num(a);
%             Ktrans(indexflow)=aa; indexflow=indexflow+1;
%         otherwise
%     end
% end
% r2=Ktrans;
% f='/v/raid1/npack/Processing/P092909_nufft/Output/flowvalues.study59.slice3.6.1_fixedDelay0.AIF_51_1_15.txt.full.txt';
% fid=fopen(f,'r');
% a=fscanf(fid,'%s',1);
% indexflow=1;
% while ~isempty(a)
%     a=fscanf(fid,'%s',1);
%     switch(a)
%         case {'Ktrans'}
%             a =fscanf(fid,'%s',1);
%             a =fscanf(fid,'%s',1);   % argh, must be a cleaner way!
%             aa=str2num(a);
%             Ktrans(indexflow)=aa; indexflow=indexflow+1;
%         otherwise
%     end
% end
% r3=Ktrans;
% DB_Ktrans_rest1=[r1 r2 r3];DB_REST_avgSlice=[DB_REST_avgSlice mean(r1) mean(r2) mean(r3)];DB_REST_avgPt=[DB_REST_avgPt mean(DB_Ktrans_rest1)];
% DB_REST=[DB_REST DB_Ktrans_rest1];
% 
% 
% %%%% To load the T1 Ktrans estimates
% f='/v/raid1/npack/Processing/P092909_nufft_24rays/Output/flowvalues.study59.slice1.6.1_fixedDelay0.AIF_59_11_1.txt.full.txt';
% fid=fopen(f,'r');
% a=fscanf(fid,'%s',1);
% indexflow=1;
% while ~isempty(a)
%     a=fscanf(fid,'%s',1);
%     switch(a)
%         case {'Ktrans'}
%             a =fscanf(fid,'%s',1);
%             a =fscanf(fid,'%s',1);   % argh, must be a cleaner way!
%             aa=str2num(a);
%             Ktrans(indexflow)=aa; indexflow=indexflow+1;
%         otherwise
%     end
% end
% r1=Ktrans;
% f='/v/raid1/npack/Processing/P092909_nufft_24rays/Output/flowvalues.study59.slice2.6.1_fixedDelay0.AIF_59_11_1.txt.full.txt';
% fid=fopen(f,'r');
% a=fscanf(fid,'%s',1);
% indexflow=1;
% while ~isempty(a)
%     a=fscanf(fid,'%s',1);
%     switch(a)
%         case {'Ktrans'}
%             a =fscanf(fid,'%s',1);
%             a =fscanf(fid,'%s',1);   % argh, must be a cleaner way!
%             aa=str2num(a);
%             Ktrans(indexflow)=aa; indexflow=indexflow+1;
%         otherwise
%     end
% end
% r2=Ktrans;
% f='/v/raid1/npack/Processing/P092909_nufft_24rays/Output/flowvalues.study59.slice3.6.1_fixedDelay0.AIF_59_11_1.txt.full.txt';
% fid=fopen(f,'r');
% a=fscanf(fid,'%s',1);
% indexflow=1;
% while ~isempty(a)
%     a=fscanf(fid,'%s',1);
%     switch(a)
%         case {'Ktrans'}
%             a =fscanf(fid,'%s',1);
%             a =fscanf(fid,'%s',1);   % argh, must be a cleaner way!
%             aa=str2num(a);
%             Ktrans(indexflow)=aa; indexflow=indexflow+1;
%         otherwise
%     end
% end
% r3=Ktrans;
% T1_Ktrans_rest1=[r1 r2 r3];T1_REST_avgSlice=[T1_REST_avgSlice mean(r1) mean(r2) mean(r3)];T1_REST_avgPt=[T1_REST_avgPt mean(T1_Ktrans_rest1)];
% T1_REST=[T1_REST T1_Ktrans_rest1];
% 
% figure(1),hold on,plot(DB_Ktrans_rest1,T1_Ktrans_rest1,'rs','LineWidth',3,'MarkerSize',8)
% disp('Plot held for current patient; press key to continue.'),pause










%%%---------P100609,  low-dose 1/5*0.02 injection---------%%%
%%%% To load the Dual-Bolus Ktrans estimates
f='/v/raid1/npack/Processing/P100609_nufft/Output/flowvalues.study47.slice1.6.1_fixedDelay0.AIF_43_1_5.txt.full.txt';
Ktrans = getKtrans(f);
r1=Ktrans;
f='/v/raid1/npack/Processing/P100609_nufft/Output/flowvalues.study47.slice2.6.1_fixedDelay0.AIF_43_1_5.txt.full.txt';
Ktrans = getKtrans(f);
r2=Ktrans;
f='/v/raid1/npack/Processing/P100609_nufft/Output/flowvalues.study47.slice3.6.1_fixedDelay0.AIF_43_1_5.txt.full.txt';
Ktrans = getKtrans(f);
r3=Ktrans;
DB_Ktrans_rest1=[r1 r2 r3];DB_REST_avgSlice=[DB_REST_avgSlice mean(r1) mean(r2) mean(r3)];DB_REST_avgPt=[DB_REST_avgPt mean(DB_Ktrans_rest1)];
DB_REST=[DB_REST DB_Ktrans_rest1];


%%%% To load the T1 Ktrans estimates
f='/v/raid1/npack/Processing/P100609_nufft_24rays/Output/flowvalues.study47.slice1.6.1_fixedDelay0.AIF_47_11_1.txt.full.txt';
Ktrans = getKtrans(f);
r1=Ktrans;
f='/v/raid1/npack/Processing/P100609_nufft_24rays/Output/flowvalues.study47.slice2.6.1_fixedDelay0.AIF_47_11_1.txt.full.txt';
Ktrans = getKtrans(f);
r2=Ktrans;
f='/v/raid1/npack/Processing/P100609_nufft_24rays/Output/flowvalues.study47.slice3.6.1_fixedDelay0.AIF_47_11_1.txt.full.txt';
Ktrans = getKtrans(f);
r3=Ktrans;
T1_Ktrans_rest1=[r1 r2 r3];T1_REST_avgSlice=[T1_REST_avgSlice mean(r1) mean(r2) mean(r3)];T1_REST_avgPt=[T1_REST_avgPt mean(T1_Ktrans_rest1)];
T1_REST=[T1_REST T1_Ktrans_rest1];

figure(1),hold on,plot(DB_Ktrans_rest1,T1_Ktrans_rest1,'bo','LineWidth',3,'MarkerSize',8)
disp('Plot held for current patient; press key to continue.'),pause(.1)










% %%%---------P100609,  0.06 high-dose injection---------%%%
% %%%% To load the Dual-Bolus Ktrans estimates
% f='/v/raid1/npack/Processing/P100609_nufft/Output/flowvalues.study51.slice1.6.1_fixedDelay0.AIF_43_1_15.txt.full.txt';
% fid=fopen(f,'r');
% a=fscanf(fid,'%s',1);
% indexflow=1;
% while ~isempty(a)
%     a=fscanf(fid,'%s',1);
%     switch(a)
%         case {'Ktrans'}
%             a =fscanf(fid,'%s',1);
%             a =fscanf(fid,'%s',1);   % argh, must be a cleaner way!
%             aa=str2num(a);
%             Ktrans(indexflow)=aa; indexflow=indexflow+1;
%         otherwise
%     end
% end
% r1=Ktrans;
% f='/v/raid1/npack/Processing/P100609_nufft/Output/flowvalues.study51.slice2.6.1_fixedDelay0.AIF_43_1_15.txt.full.txt';
% fid=fopen(f,'r');
% a=fscanf(fid,'%s',1);
% indexflow=1;
% while ~isempty(a)
%     a=fscanf(fid,'%s',1);
%     switch(a)
%         case {'Ktrans'}
%             a =fscanf(fid,'%s',1);
%             a =fscanf(fid,'%s',1);   % argh, must be a cleaner way!
%             aa=str2num(a);
%             Ktrans(indexflow)=aa; indexflow=indexflow+1;
%         otherwise
%     end
% end
% r2=Ktrans;
% f='/v/raid1/npack/Processing/P100609_nufft/Output/flowvalues.study51.slice3.6.1_fixedDelay0.AIF_43_1_15.txt.full.txt';
% fid=fopen(f,'r');
% a=fscanf(fid,'%s',1);
% indexflow=1;
% while ~isempty(a)
%     a=fscanf(fid,'%s',1);
%     switch(a)
%         case {'Ktrans'}
%             a =fscanf(fid,'%s',1);
%             a =fscanf(fid,'%s',1);   % argh, must be a cleaner way!
%             aa=str2num(a);
%             Ktrans(indexflow)=aa; indexflow=indexflow+1;
%         otherwise
%     end
% end
% r3=Ktrans;
% DB_Ktrans_rest1=[r1 r2 r3];DB_REST_avgSlice=[DB_REST_avgSlice mean(r1) mean(r2) mean(r3)];DB_REST_avgPt=[DB_REST_avgPt mean(DB_Ktrans_rest1)];
% DB_REST=[DB_REST DB_Ktrans_rest1];
% 
% 
% %%%% To load the T1 Ktrans estimates
% f='/v/raid1/npack/Processing/P100609_nufft_24rays/Output/flowvalues.study51.slice1.6.1_fixedDelay0.AIF_51_11_1.txt.full.txt';
% fid=fopen(f,'r');
% a=fscanf(fid,'%s',1);
% indexflow=1;
% while ~isempty(a)
%     a=fscanf(fid,'%s',1);
%     switch(a)
%         case {'Ktrans'}
%             a =fscanf(fid,'%s',1);
%             a =fscanf(fid,'%s',1);   % argh, must be a cleaner way!
%             aa=str2num(a);
%             Ktrans(indexflow)=aa; indexflow=indexflow+1;
%         otherwise
%     end
% end
% r1=Ktrans;
% f='/v/raid1/npack/Processing/P100609_nufft_24rays/Output/flowvalues.study51.slice2.6.1_fixedDelay0.AIF_51_11_1.txt.full.txt';
% fid=fopen(f,'r');
% a=fscanf(fid,'%s',1);
% indexflow=1;
% while ~isempty(a)
%     a=fscanf(fid,'%s',1);
%     switch(a)
%         case {'Ktrans'}
%             a =fscanf(fid,'%s',1);
%             a =fscanf(fid,'%s',1);   % argh, must be a cleaner way!
%             aa=str2num(a);
%             Ktrans(indexflow)=aa; indexflow=indexflow+1;
%         otherwise
%     end
% end
% r2=Ktrans;
% f='/v/raid1/npack/Processing/P100609_nufft_24rays/Output/flowvalues.study51.slice3.6.1_fixedDelay0.AIF_51_11_1.txt.full.txt';
% fid=fopen(f,'r');
% a=fscanf(fid,'%s',1);
% indexflow=1;
% while ~isempty(a)
%     a=fscanf(fid,'%s',1);
%     switch(a)
%         case {'Ktrans'}
%             a =fscanf(fid,'%s',1);
%             a =fscanf(fid,'%s',1);   % argh, must be a cleaner way!
%             aa=str2num(a);
%             Ktrans(indexflow)=aa; indexflow=indexflow+1;
%         otherwise
%     end
% end
% r3=Ktrans;
% T1_Ktrans_rest1=[r1 r2 r3];T1_REST_avgSlice=[T1_REST_avgSlice mean(r1) mean(r2) mean(r3)];T1_REST_avgPt=[T1_REST_avgPt mean(T1_Ktrans_rest1)];
% T1_REST=[T1_REST T1_Ktrans_rest1];
% 
% figure(1),hold on,plot(DB_Ktrans_rest1,T1_Ktrans_rest1,'rs','LineWidth',3,'MarkerSize',8)
% disp('Plot held for current patient; press key to continue.'),pause


















%%%---------P111309---------%%%
%%%% To load the Dual-Bolus Ktrans estimates
f='/v/raid1/npack/Processing_72ray_subset123/P111309_nufft_24rays/Output/flowvalues.study32.slice1.6.1_fixedDelay0.AIF_23_1_5.txt.full.txt';
Ktrans = getKtrans(f);
r1=Ktrans;
f='/v/raid1/npack/Processing_72ray_subset123/P111309_nufft_24rays/Output/flowvalues.study32.slice2.6.1_fixedDelay0.AIF_23_2_5.txt.full.txt';
Ktrans = getKtrans(f);
r2=Ktrans;
f='/v/raid1/npack/Processing_72ray_subset123/P111309_nufft_24rays/Output/flowvalues.study32.slice3.6.1_fixedDelay0.AIF_23_3_5.txt.full.txt';
Ktrans = getKtrans(f);
r3=Ktrans;
f='/v/raid1/npack/Processing_72ray_subset123/P111309_nufft_24rays/Output/flowvalues.study32.slice4.6.1_fixedDelay0.AIF_23_4_5.txt.full.txt';
Ktrans = getKtrans(f);
r4=Ktrans;
DB_Ktrans_rest1=[r1 r2 r3 r4];DB_REST_avgSlice=[DB_REST_avgSlice mean(r1) mean(r2) mean(r3) mean(r4)];DB_REST_avgPt=[DB_REST_avgPt mean(DB_Ktrans_rest1)];
DB_REST=[DB_REST DB_Ktrans_rest1];


%%%% To load the T1 Ktrans estimates
f='/v/raid1/npack/Processing_72ray_subset123/P111309_nufft_24rays/Output/flowvalues.study32.slice1.6.1_fixedDelay0.AIF_32_11_1.txt.full.txt';
Ktrans = getKtrans(f);
r1=Ktrans;
f='/v/raid1/npack/Processing_72ray_subset123/P111309_nufft_24rays/Output/flowvalues.study32.slice2.6.1_fixedDelay0.AIF_32_12_1.txt.full.txt';
Ktrans = getKtrans(f);
r2=Ktrans;
f='/v/raid1/npack/Processing_72ray_subset123/P111309_nufft_24rays/Output/flowvalues.study32.slice3.6.1_fixedDelay0.AIF_32_13_1.txt.full.txt';
Ktrans = getKtrans(f);
r3=Ktrans;
f='/v/raid1/npack/Processing_72ray_subset123/P111309_nufft_24rays/Output/flowvalues.study32.slice4.6.1_fixedDelay0.AIF_32_14_1.txt.full.txt';
Ktrans = getKtrans(f);
r4=Ktrans;
T1_Ktrans_rest1=[r1 r2 r3 r4];T1_REST_avgSlice=[T1_REST_avgSlice mean(r1) mean(r2) mean(r3) mean(r4)];T1_REST_avgPt=[T1_REST_avgPt mean(T1_Ktrans_rest1)];
T1_REST=[T1_REST T1_Ktrans_rest1];

%%%% To load the Adenosine Ktrans estimates
f='/v/raid1/npack/Processing_72ray_subset123/P111309_nufft_24rays/Output/flowvalues.study39.slice1.6.1_fixedDelay0.AIF_39_11_1.txt.full.txt';
Ktrans = getKtrans(f);
s1=Ktrans;
f='/v/raid1/npack/Processing_72ray_subset123/P111309_nufft_24rays/Output/flowvalues.study39.slice2.6.1_fixedDelay0.AIF_39_12_1.txt.full.txt';
Ktrans = getKtrans(f);
s2=Ktrans;
f='/v/raid1/npack/Processing_72ray_subset123/P111309_nufft_24rays/Output/flowvalues.study39.slice3.6.1_fixedDelay0.AIF_39_13_1.txt.full.txt';
Ktrans = getKtrans(f);
s3=Ktrans;
Aden_1=[s1 s2 s3];ADEN_avgSlice=[ADEN_avgSlice mean(s1) mean(s2) mean(s3)];ADEN_avgPt=[ADEN_avgPt mean(Aden_1)];
ADEN=[ADEN Aden_1];

%%%% To load the Lexiscan Ktrans estimates
f='/v/raid1/npack/Processing_72ray_subset123/P111309_nufft_24rays/Output/flowvalues.study55.slice1.6.1_fixedDelay0.AIF_55_11_1.txt.full.txt';
Ktrans = getKtrans(f);
s1=Ktrans;
f='/v/raid1/npack/Processing_72ray_subset123/P111309_nufft_24rays/Output/flowvalues.study55.slice2.6.1_fixedDelay0.AIF_55_12_1.txt.full.txt';
Ktrans = getKtrans(f);
s2=Ktrans;
f='/v/raid1/npack/Processing_72ray_subset123/P111309_nufft_24rays/Output/flowvalues.study55.slice3.6.1_fixedDelay0.AIF_55_13_1.txt.full.txt';
Ktrans = getKtrans(f);
s3=Ktrans;
Lex_1=[s1 s2 s3];LEX_avgSlice=[LEX_avgSlice mean(s1) mean(s2) mean(s3)];LEX_avgPt=[LEX_avgPt mean(Lex_1)];
LEX=[LEX Lex_1];

REST2=[REST2 T1_Ktrans_rest1(7:24)];REST2_avgSlice=[REST2_avgSlice mean(T1_Ktrans_rest1(7:12)) mean(T1_Ktrans_rest1(13:18)) mean(T1_Ktrans_rest1(19:24))];REST2_avgPt=[REST2_avgPt mean(T1_Ktrans_rest1)];
figure(1),hold on,plot(DB_Ktrans_rest1,T1_Ktrans_rest1,'bo','LineWidth',3,'MarkerSize',8)
disp('Plot held for current patient; press key to continue.'),pause(.1)



















%%%---------P111609---------%%%
%%%% To load the Dual-Bolus Ktrans estimates
f='/v/raid1/npack/Processing_72ray_subset123/P111609_nufft_24rays/Output/flowvalues.study44.slice1.6.1_fixedDelay0.AIF_25_1_5.txt.full.txt';
Ktrans = getKtrans(f);
r1=Ktrans;
f='/v/raid1/npack/Processing_72ray_subset123/P111609_nufft_24rays/Output/flowvalues.study44.slice2.6.1_fixedDelay0.AIF_25_2_5.txt.full.txt';
Ktrans = getKtrans(f);
r2=Ktrans;
f='/v/raid1/npack/Processing_72ray_subset123/P111609_nufft_24rays/Output/flowvalues.study44.slice3.6.1_fixedDelay0.AIF_25_3_5.txt.full.txt';
Ktrans = getKtrans(f);
r3=Ktrans;
DB_Ktrans_rest1=[r1 r2 r3];DB_REST_avgSlice=[DB_REST_avgSlice mean(r1) mean(r2) mean(r3)];DB_REST_avgPt=[DB_REST_avgPt mean(DB_Ktrans_rest1)];
DB_REST=[DB_REST DB_Ktrans_rest1];


%%%% To load the T1 Ktrans estimates
f='/v/raid1/npack/Processing_72ray_subset123/P111609_nufft_24rays/Output/flowvalues.study44.slice1.6.1_fixedDelay0.AIF_44_11_1.txt.full.txt';
Ktrans = getKtrans(f);
r1=Ktrans;
f='/v/raid1/npack/Processing_72ray_subset123/P111609_nufft_24rays/Output/flowvalues.study44.slice2.6.1_fixedDelay0.AIF_44_12_1.txt.full.txt';
Ktrans = getKtrans(f);
r2=Ktrans;
f='/v/raid1/npack/Processing_72ray_subset123/P111609_nufft_24rays/Output/flowvalues.study44.slice3.6.1_fixedDelay0.AIF_44_13_1.txt.full.txt';
Ktrans = getKtrans(f);
r3=Ktrans;
T1_Ktrans_rest1=[r1 r2 r3];T1_REST_avgSlice=[T1_REST_avgSlice mean(r1) mean(r2) mean(r3)];T1_REST_avgPt=[T1_REST_avgPt mean(T1_Ktrans_rest1)];
T1_REST=[T1_REST T1_Ktrans_rest1];

%%%% To load the Adenosine Ktrans estimates
f='/v/raid1/npack/Processing_72ray_subset123/P111609_nufft_24rays/Output/flowvalues.study47.slice1.6.1_fixedDelay0.AIF_47_11_1.txt.full.txt';
Ktrans = getKtrans(f);
s1=Ktrans;
f='/v/raid1/npack/Processing_72ray_subset123/P111609_nufft_24rays/Output/flowvalues.study47.slice2.6.1_fixedDelay0.AIF_47_12_1.txt.full.txt';
Ktrans = getKtrans(f);
s2=Ktrans;
Aden_1=[s1 s2];ADEN_avgSlice=[ADEN_avgSlice mean(s1) mean(s2)];ADEN_avgPt=[ADEN_avgPt mean(Aden_1)];
ADEN=[ADEN Aden_1];

%%%% To load the Lexiscan Ktrans estimates
f='/v/raid1/npack/Processing_72ray_subset123/P111609_nufft_24rays/Output/flowvalues.study68.slice1.6.1_fixedDelay0.AIF_68_11_1.txt.full.txt';
Ktrans = getKtrans(f);
s1=Ktrans;
f='/v/raid1/npack/Processing_72ray_subset123/P111609_nufft_24rays/Output/flowvalues.study68.slice2.6.1_fixedDelay0.AIF_68_12_1.txt.full.txt';
Ktrans = getKtrans(f);
s2=Ktrans;
Lex_1=[s1 s2];LEX_avgSlice=[LEX_avgSlice mean(s1) mean(s2)];LEX_avgPt=[LEX_avgPt mean(Lex_1)];
LEX=[LEX Lex_1];

REST2=[REST2 T1_Ktrans_rest1(1:12)];REST2_avgSlice=[REST2_avgSlice mean(T1_Ktrans_rest1(1:6)) mean(T1_Ktrans_rest1(7:12))];REST2_avgPt=[REST2_avgPt mean(T1_Ktrans_rest1(1:12))];
figure(1),hold on,plot(DB_Ktrans_rest1,T1_Ktrans_rest1,'bo','LineWidth',3,'MarkerSize',8)
disp('Plot held for current patient; press key to continue.'),pause(.1)




















%%%---------P120109---------%%%
%%%% To load the Dual-Bolus Ktrans estimates
f='/v/raid1/npack/Processing/P120109_nufft_24rays/Output/flowvalues.study24.slice1.6.1_fixedDelay0.AIF_17_1_5.txt.full.txt';
Ktrans = getKtrans(f);
r1=Ktrans;
f='/v/raid1/npack/Processing/P120109_nufft_24rays/Output/flowvalues.study24.slice2.6.1_fixedDelay0.AIF_17_1_5.txt.full.txt';
Ktrans = getKtrans(f);
r2=Ktrans;
f='/v/raid1/npack/Processing/P120109_nufft_24rays/Output/flowvalues.study24.slice3.6.1_fixedDelay0.AIF_17_1_5.txt.full.txt';
Ktrans = getKtrans(f);
r3=Ktrans;
DB_Ktrans_rest1=[r1 r2 r3];DB_REST_avgSlice=[DB_REST_avgSlice mean(r1) mean(r2) mean(r3)];DB_REST_avgPt=[DB_REST_avgPt mean(DB_Ktrans_rest1)];
DB_REST=[DB_REST DB_Ktrans_rest1];


%%%% To load the T1 Ktrans estimates
f='/v/raid1/npack/Processing/P120109_nufft_24rays/Output/flowvalues.study24.slice1.6.1_fixedDelay0.AIF_24_11_1.txt.full.txt';
Ktrans = getKtrans(f);
r1=Ktrans;
f='/v/raid1/npack/Processing/P120109_nufft_24rays/Output/flowvalues.study24.slice2.6.1_fixedDelay0.AIF_24_11_1.txt.full.txt';
Ktrans = getKtrans(f);
r2=Ktrans;
f='/v/raid1/npack/Processing/P120109_nufft_24rays/Output/flowvalues.study24.slice3.6.1_fixedDelay0.AIF_24_11_1.txt.full.txt';
Ktrans = getKtrans(f);
r3=Ktrans;
T1_Ktrans_rest1=[r1 r2 r3];T1_REST_avgSlice=[T1_REST_avgSlice mean(r1) mean(r2) mean(r3)];T1_REST_avgPt=[T1_REST_avgPt mean(T1_Ktrans_rest1)];
T1_REST=[T1_REST T1_Ktrans_rest1];

%%%% To load the Adenosine Ktrans estimates
f='/v/raid1/npack/Processing/P120109_nufft_24rays/Output/flowvalues.study28.slice1.6.1_fixedDelay0.AIF_28_11_1.txt.full.txt';
Ktrans = getKtrans(f);
s1=Ktrans;
f='/v/raid1/npack/Processing/P120109_nufft_24rays/Output/flowvalues.study28.slice2.6.1_fixedDelay0.AIF_28_11_1.txt.full.txt';
Ktrans = getKtrans(f);
s2=Ktrans;
f='/v/raid1/npack/Processing/P120109_nufft_24rays/Output/flowvalues.study28.slice3.6.1_fixedDelay0.AIF_28_11_1.txt.full.txt';
Ktrans = getKtrans(f);
s3=Ktrans;
Aden_1=[s1 s2 s3];ADEN_avgSlice=[ADEN_avgSlice mean(s1) mean(s2) mean(s3)];ADEN_avgPt=[ADEN_avgPt mean(Aden_1)];
ADEN=[ADEN Aden_1];

%%%% To load the Lexiscan Ktrans estimates
f='/v/raid1/npack/Processing/P120109_nufft_24rays/Output/flowvalues.study49.slice1.6.1_fixedDelay0.AIF_49_13_1.txt.full.txt';
Ktrans = getKtrans(f);
s1=Ktrans;
f='/v/raid1/npack/Processing/P120109_nufft_24rays/Output/flowvalues.study49.slice2.6.1_fixedDelay0.AIF_49_13_1.txt.full.txt';
Ktrans = getKtrans(f);
s2=Ktrans;
f='/v/raid1/npack/Processing/P120109_nufft_24rays/Output/flowvalues.study49.slice3.6.1_fixedDelay0.AIF_49_13_1.txt.full.txt';
Ktrans = getKtrans(f);
s3=Ktrans;
Lex_1=[s1 s2 s3];LEX_avgSlice=[LEX_avgSlice mean(s1) mean(s2) mean(s3)];LEX_avgPt=[LEX_avgPt mean(Lex_1)];
LEX=[LEX Lex_1];

REST2=[REST2 T1_Ktrans_rest1];REST2_avgSlice=[REST2_avgSlice mean(T1_Ktrans_rest1(1:6)) mean(T1_Ktrans_rest1(7:12)) mean(T1_Ktrans_rest1(13:18))];REST2_avgPt=[REST2_avgPt mean(T1_Ktrans_rest1)];
figure(1),hold on,plot(DB_Ktrans_rest1,T1_Ktrans_rest1,'bo','LineWidth',3,'MarkerSize',8)
disp('Plot held for current patient; press key to continue.'),pause(.1)












%%%---------P120409---------%%%
%%%% To load the Dual-Bolus Ktrans estimates
f='/v/raid1/npack/Processing/P120409_nufft_24rays/Output/flowvalues.study19.slice1.6.1_fixedDelay0.AIF_16_3_5.txt.full.txt';
Ktrans = getKtrans(f);
r1=Ktrans;
f='/v/raid1/npack/Processing/P120409_nufft_24rays/Output/flowvalues.study19.slice2.6.1_fixedDelay0.AIF_16_3_5.txt.full.txt';
Ktrans = getKtrans(f);
r2=Ktrans;
f='/v/raid1/npack/Processing/P120409_nufft_24rays/Output/flowvalues.study19.slice3.6.1_fixedDelay0.AIF_16_3_5.txt.full.txt';
Ktrans = getKtrans(f);
r3=Ktrans;
DB_Ktrans_rest1=[r1 r2 r3];DB_REST_avgSlice=[DB_REST_avgSlice mean(r1) mean(r2) mean(r3)];DB_REST_avgPt=[DB_REST_avgPt mean(DB_Ktrans_rest1)];
DB_REST=[DB_REST DB_Ktrans_rest1];


%%%% To load the T1 Ktrans estimates
f='/v/raid1/npack/Processing/P120409_nufft_24rays/Output/flowvalues.study19.slice1.6.1_fixedDelay0.AIF_19_13_1.txt.full.txt';
Ktrans = getKtrans(f);
r1=Ktrans;
f='/v/raid1/npack/Processing/P120409_nufft_24rays/Output/flowvalues.study19.slice2.6.1_fixedDelay0.AIF_19_13_1.txt.full.txt';
Ktrans = getKtrans(f);
r2=Ktrans;
f='/v/raid1/npack/Processing/P120409_nufft_24rays/Output/flowvalues.study19.slice3.6.1_fixedDelay0.AIF_19_13_1.txt.full.txt';
Ktrans = getKtrans(f);
r3=Ktrans;
T1_Ktrans_rest1=[r1 r2 r3];T1_REST_avgSlice=[T1_REST_avgSlice mean(r1) mean(r2) mean(r3)];T1_REST_avgPt=[T1_REST_avgPt mean(T1_Ktrans_rest1)];
T1_REST=[T1_REST T1_Ktrans_rest1];

%%%% To load the Adenosine Ktrans estimates
f='/v/raid1/npack/Processing/P120409_nufft_24rays/Output/flowvalues.study56.slice1.6.1_fixedDelay0.AIF_56_13_1.txt.full.txt';
Ktrans = getKtrans(f);
s1=Ktrans;
f='/v/raid1/npack/Processing/P120409_nufft_24rays/Output/flowvalues.study56.slice2.6.1_fixedDelay0.AIF_56_13_1.txt.full.txt';
Ktrans = getKtrans(f);
s2=Ktrans;
f='/v/raid1/npack/Processing/P120409_nufft_24rays/Output/flowvalues.study56.slice3.6.1_fixedDelay0.AIF_56_13_1.txt.full.txt';
Ktrans = getKtrans(f);
s3=Ktrans;
Aden_1=[s1 s2 s3];ADEN_avgSlice=[ADEN_avgSlice mean(s1) mean(s2) mean(s3)];ADEN_avgPt=[ADEN_avgPt mean(Aden_1)];
ADEN=[ADEN Aden_1];

%%%% To load the Lexiscan Ktrans estimates
f='/v/raid1/npack/Processing/P120409_nufft_24rays/Output/flowvalues.study26.slice1.6.1_fixedDelay0.AIF_26_13_1.txt.full.txt';
Ktrans = getKtrans(f);
s1=Ktrans;
f='/v/raid1/npack/Processing/P120409_nufft_24rays/Output/flowvalues.study26.slice2.6.1_fixedDelay0.AIF_26_13_1.txt.full.txt';
Ktrans = getKtrans(f);
s2=Ktrans;
f='/v/raid1/npack/Processing/P120409_nufft_24rays/Output/flowvalues.study26.slice3.6.1_fixedDelay0.AIF_26_13_1.txt.full.txt';
Ktrans = getKtrans(f);
s3=Ktrans;
Lex_1=[s1 s2 s3];LEX_avgSlice=[LEX_avgSlice mean(s1) mean(s2) mean(s3)];LEX_avgPt=[LEX_avgPt mean(Lex_1)];
LEX=[LEX Lex_1];

REST2=[REST2 T1_Ktrans_rest1];REST2_avgSlice=[REST2_avgSlice mean(T1_Ktrans_rest1(1:6)) mean(T1_Ktrans_rest1(7:12)) mean(T1_Ktrans_rest1(13:18))];REST2_avgPt=[REST2_avgPt mean(T1_Ktrans_rest1)];
figure(1),hold on,plot(DB_Ktrans_rest1,T1_Ktrans_rest1,'bo','LineWidth',3,'MarkerSize',8)
disp('Plot held for current patient; press key to continue.'),pause(.1)


















%%%---------P120709---------%%%
%%%% To load the Dual-Bolus Ktrans estimates
f='/v/raid1/npack/Processing/P120709_nufft_24rays/Output/flowvalues.study27.slice1.6.1_fixedDelay0.AIF_20_1_5.txt.full.txt';
Ktrans = getKtrans(f);
r1=Ktrans;
f='/v/raid1/npack/Processing/P120709_nufft_24rays/Output/flowvalues.study27.slice2.6.1_fixedDelay0.AIF_20_1_5.txt.full.txt';
Ktrans = getKtrans(f);
r2=Ktrans;
f='/v/raid1/npack/Processing/P120709_nufft_24rays/Output/flowvalues.study27.slice3.6.1_fixedDelay0.AIF_20_1_5.txt.full.txt';
Ktrans = getKtrans(f);
r3=Ktrans;
DB_Ktrans_rest1=[r1 r2 r3];DB_REST_avgSlice=[DB_REST_avgSlice mean(r1) mean(r2) mean(r3)];DB_REST_avgPt=[DB_REST_avgPt mean(DB_Ktrans_rest1)];
DB_REST=[DB_REST DB_Ktrans_rest1];


%%%% To load the T1 Ktrans estimates
f='/v/raid1/npack/Processing/P120709_nufft_24rays/Output/flowvalues.study27.slice1.6.1_fixedDelay0.AIF_27_11_1.txt.full.txt';
Ktrans = getKtrans(f);
r1=Ktrans;
f='/v/raid1/npack/Processing/P120709_nufft_24rays/Output/flowvalues.study27.slice2.6.1_fixedDelay0.AIF_27_11_1.txt.full.txt';
Ktrans = getKtrans(f);
r2=Ktrans;
f='/v/raid1/npack/Processing/P120709_nufft_24rays/Output/flowvalues.study27.slice3.6.1_fixedDelay0.AIF_27_11_1.txt.full.txt';
Ktrans = getKtrans(f);
r3=Ktrans;
T1_Ktrans_rest1=[r1 r2 r3];T1_REST_avgSlice=[T1_REST_avgSlice mean(r1) mean(r2) mean(r3)];T1_REST_avgPt=[T1_REST_avgPt mean(T1_Ktrans_rest1)];
T1_REST=[T1_REST T1_Ktrans_rest1];

%%%% To load the Adenosine Ktrans estimates
f='/v/raid1/npack/Processing/P120709_nufft_24rays/Output/flowvalues.study33.slice1.6.1_fixedDelay0.AIF_33_11_1.txt.full.txt';
Ktrans = getKtrans(f);
s1=Ktrans;
f='/v/raid1/npack/Processing/P120709_nufft_24rays/Output/flowvalues.study33.slice2.6.1_fixedDelay0.AIF_33_11_1.txt.full.txt';
Ktrans = getKtrans(f);
s2=Ktrans;
Aden_1=[s1 s2];ADEN_avgSlice=[ADEN_avgSlice mean(s1) mean(s2)];ADEN_avgPt=[ADEN_avgPt mean(Aden_1)];
ADEN=[ADEN Aden_1];

%%%% To load the Lexiscan Ktrans estimates
f='/v/raid1/npack/Processing/P120709_nufft_24rays/Output/flowvalues.study60.slice1.6.1_fixedDelay0.AIF_60_11_1.txt.full.txt';
Ktrans = getKtrans(f);
s1=Ktrans;
f='/v/raid1/npack/Processing/P120709_nufft_24rays/Output/flowvalues.study60.slice2.6.1_fixedDelay0.AIF_60_11_1.txt.full.txt';
Ktrans = getKtrans(f);
s2=Ktrans;
Lex_1=[s1 s2];LEX_avgSlice=[LEX_avgSlice mean(s1) mean(s2)];LEX_avgPt=[LEX_avgPt mean(Lex_1)];
LEX=[LEX Lex_1];

REST2=[REST2 T1_Ktrans_rest1(1:6) T1_Ktrans_rest1(7:12)];REST2_avgSlice=[REST2_avgSlice mean(T1_Ktrans_rest1(1:6)) mean(T1_Ktrans_rest1(7:12))];REST2_avgPt=[REST2_avgPt mean(T1_Ktrans_rest1(1:12))];
figure(1),hold on,plot(DB_Ktrans_rest1,T1_Ktrans_rest1,'bo','LineWidth',3,'MarkerSize',8)
disp('Plot held for current patient; press key to continue.'),pause(.1)


















%%%---------P121509---------%%%
%%%% To load the Dual-Bolus Ktrans estimates
f='/v/raid1/npack/Processing/P121509_nufft_24rays/Output/flowvalues.study24.slice1.6.1_fixedDelay0.AIF_19_1_5.txt.full.txt';
Ktrans = getKtrans(f);
r1=Ktrans;
f='/v/raid1/npack/Processing/P121509_nufft_24rays/Output/flowvalues.study24.slice2.6.1_fixedDelay0.AIF_19_1_5.txt.full.txt';
Ktrans = getKtrans(f);
r2=Ktrans;
f='/v/raid1/npack/Processing/P121509_nufft_24rays/Output/flowvalues.study24.slice3.6.1_fixedDelay0.AIF_19_1_5.txt.full.txt';
Ktrans = getKtrans(f);
r3=Ktrans;
f='/v/raid1/npack/Processing/P121509_nufft_24rays/Output/flowvalues.study24.slice4.6.1_fixedDelay0.AIF_19_1_5.txt.full.txt';
Ktrans = getKtrans(f);
r4=Ktrans;
DB_Ktrans_rest1=[r1 r2 r3 r4];DB_REST_avgSlice=[DB_REST_avgSlice mean(r1) mean(r2) mean(r3) mean(r4)];DB_REST_avgPt=[DB_REST_avgPt mean(DB_Ktrans_rest1)];
DB_REST=[DB_REST DB_Ktrans_rest1];


%%%% To load the T1 Ktrans estimates
f='/v/raid1/npack/Processing/P121509_nufft_24rays/Output/flowvalues.study24.slice1.6.1_fixedDelay0.AIF_24_11_1.txt.full.txt';
Ktrans = getKtrans(f);
r1=Ktrans;
f='/v/raid1/npack/Processing/P121509_nufft_24rays/Output/flowvalues.study24.slice2.6.1_fixedDelay0.AIF_24_11_1.txt.full.txt';
Ktrans = getKtrans(f);
r2=Ktrans;
f='/v/raid1/npack/Processing/P121509_nufft_24rays/Output/flowvalues.study24.slice3.6.1_fixedDelay0.AIF_24_11_1.txt.full.txt';
Ktrans = getKtrans(f);
r3=Ktrans;
f='/v/raid1/npack/Processing/P121509_nufft_24rays/Output/flowvalues.study24.slice4.6.1_fixedDelay0.AIF_24_11_1.txt.full.txt';
Ktrans = getKtrans(f);
r4=Ktrans;
T1_Ktrans_rest1=[r1 r2 r3 r4];T1_REST_avgSlice=[T1_REST_avgSlice mean(r1) mean(r2) mean(r3) mean(r4)];T1_REST_avgPt=[T1_REST_avgPt mean(T1_Ktrans_rest1)];
T1_REST=[T1_REST T1_Ktrans_rest1];

%%%% To load the Adenosine Ktrans estimates
f='/v/raid1/npack/Processing/P121509_nufft_24rays/Output/flowvalues.study32.slice1.6.1_fixedDelay0.AIF_32_11_1.txt.full.txt';
Ktrans = getKtrans(f);
s1=Ktrans;
f='/v/raid1/npack/Processing/P121509_nufft_24rays/Output/flowvalues.study32.slice2.6.1_fixedDelay0.AIF_32_11_1.txt.full.txt';
Ktrans = getKtrans(f);
s2=Ktrans;
f='/v/raid1/npack/Processing/P121509_nufft_24rays/Output/flowvalues.study32.slice3.6.1_fixedDelay0.AIF_32_11_1.txt.full.txt';
Ktrans = getKtrans(f);
s3=Ktrans;
Aden_1=[s1 s2 s3];ADEN_avgSlice=[ADEN_avgSlice mean(s1) mean(s2) mean(s3)];ADEN_avgPt=[ADEN_avgPt mean(Aden_1)];
ADEN=[ADEN Aden_1];

%%%% To load the Lexiscan Ktrans estimates
f='/v/raid1/npack/Processing/P121509_nufft_24rays/Output/flowvalues.study55.slice1.6.1_fixedDelay0.AIF_55_11_1.txt.full.txt';
Ktrans = getKtrans(f);
s1=Ktrans;
f='/v/raid1/npack/Processing/P121509_nufft_24rays/Output/flowvalues.study55.slice2.6.1_fixedDelay0.AIF_55_11_1.txt.full.txt';
Ktrans = getKtrans(f);
s2=Ktrans;
f='/v/raid1/npack/Processing/P121509_nufft_24rays/Output/flowvalues.study55.slice3.6.1_fixedDelay0.AIF_55_11_1.txt.full.txt';
Ktrans = getKtrans(f);
s3=Ktrans;
Lex_1=[s1 s2 s3];LEX_avgSlice=[LEX_avgSlice mean(s1) mean(s2) mean(s3)];LEX_avgPt=[LEX_avgPt mean(Lex_1)];
LEX=[LEX Lex_1];

REST2=[REST2 T1_Ktrans_rest1(7:12) T1_Ktrans_rest1(19:24) T1_Ktrans_rest1(13:18)];REST2_avgSlice=[REST2_avgSlice mean(T1_Ktrans_rest1(7:12)) mean(T1_Ktrans_rest1(19:24)) mean(T1_Ktrans_rest1(13:18))];REST2_avgPt=[REST2_avgPt mean(T1_Ktrans_rest1(7:24))];
figure(1),hold on,plot(DB_Ktrans_rest1,T1_Ktrans_rest1,'bo','LineWidth',3,'MarkerSize',8)
disp('Plot held for current patient; press key to continue.'),pause(.1)















%%%---------P121609---------%%%
%%%% To load the Dual-Bolus Ktrans estimates
f='/v/raid1/npack/Processing/P121609_nufft_24rays/Output/flowvalues.study23.slice1.6.1_fixedDelay0.AIF_19_1_5.txt.full.txt';
Ktrans = getKtrans(f);
r1=Ktrans;
f='/v/raid1/npack/Processing/P121609_nufft_24rays/Output/flowvalues.study23.slice2.6.1_fixedDelay0.AIF_19_1_5.txt.full.txt';
Ktrans = getKtrans(f);
r2=Ktrans;
f='/v/raid1/npack/Processing/P121609_nufft_24rays/Output/flowvalues.study23.slice3.6.1_fixedDelay0.AIF_19_1_5.txt.full.txt';
Ktrans = getKtrans(f);
r3=Ktrans;
DB_Ktrans_rest1=[r1 r2 r3];DB_REST_avgSlice=[DB_REST_avgSlice mean(r1) mean(r2) mean(r3)];DB_REST_avgPt=[DB_REST_avgPt mean(DB_Ktrans_rest1)];
DB_REST=[DB_REST DB_Ktrans_rest1];


%%%% To load the T1 Ktrans estimates
f='/v/raid1/npack/Processing/P121609_nufft_24rays/Output/flowvalues.study23.slice1.6.1_fixedDelay0.AIF_23_11_1.txt.full.txt';
Ktrans = getKtrans(f);
r1=Ktrans;
f='/v/raid1/npack/Processing/P121609_nufft_24rays/Output/flowvalues.study23.slice2.6.1_fixedDelay0.AIF_23_11_1.txt.full.txt';
Ktrans = getKtrans(f);
r2=Ktrans;
f='/v/raid1/npack/Processing/P121609_nufft_24rays/Output/flowvalues.study23.slice3.6.1_fixedDelay0.AIF_23_11_1.txt.full.txt';
Ktrans = getKtrans(f);
r3=Ktrans;
T1_Ktrans_rest1=[r1 r2 r3];T1_REST_avgSlice=[T1_REST_avgSlice mean(r1) mean(r2) mean(r3)];T1_REST_avgPt=[T1_REST_avgPt mean(T1_Ktrans_rest1)];
T1_REST=[T1_REST T1_Ktrans_rest1];

%%%% To load the Adenosine Ktrans estimates
f='/v/raid1/npack/Processing/P121609_nufft_24rays/Output/flowvalues.study26.slice1.6.1_fixedDelay0.AIF_26_11_1.txt.full.txt';
Ktrans = getKtrans(f);
s1=Ktrans;
f='/v/raid1/npack/Processing/P121609_nufft_24rays/Output/flowvalues.study26.slice2.6.1_fixedDelay0.AIF_26_11_1.txt.full.txt';
Ktrans = getKtrans(f);
s2=Ktrans;
f='/v/raid1/npack/Processing/P121609_nufft_24rays/Output/flowvalues.study26.slice3.6.1_fixedDelay0.AIF_26_11_1.txt.full.txt';
Ktrans = getKtrans(f);
s3=Ktrans;
Aden_1=[s1 s2 s3];ADEN_avgSlice=[ADEN_avgSlice mean(s1) mean(s2) mean(s3)];ADEN_avgPt=[ADEN_avgPt mean(Aden_1)];
ADEN=[ADEN Aden_1];

%%%% To load the Lexiscan Ktrans estimates
f='/v/raid1/npack/Processing/P121609_nufft_24rays/Output/flowvalues.study49.slice1.6.1_fixedDelay0.AIF_49_11_1.txt.full.txt';
Ktrans = getKtrans(f);
s1=Ktrans;
f='/v/raid1/npack/Processing/P121609_nufft_24rays/Output/flowvalues.study49.slice2.6.1_fixedDelay0.AIF_49_11_1.txt.full.txt';
Ktrans = getKtrans(f);
s2=Ktrans;
f='/v/raid1/npack/Processing/P121609_nufft_24rays/Output/flowvalues.study49.slice3.6.1_fixedDelay0.AIF_49_11_1.txt.full.txt';
Ktrans = getKtrans(f);
s3=Ktrans;
Lex_1=[s1 s2 s3];LEX_avgSlice=[LEX_avgSlice mean(s1) mean(s2) mean(s3)];LEX_avgPt=[LEX_avgPt mean(Lex_1)];
LEX=[LEX Lex_1];

REST2=[REST2 T1_Ktrans_rest1];REST2_avgSlice=[REST2_avgSlice mean(T1_Ktrans_rest1(1:6)) mean(T1_Ktrans_rest1(7:12)) mean(T1_Ktrans_rest1(13:18))];REST2_avgPt=[REST2_avgPt mean(T1_Ktrans_rest1)];
figure(1),hold on,plot(DB_Ktrans_rest1,T1_Ktrans_rest1,'bo','LineWidth',3,'MarkerSize',8)
disp('Plot held for current patient; press key to continue.'),pause(.1)
















%%%---------P010710---------%%%
%%%% To load the Dual-Bolus Ktrans estimates
f='/v/raid1/npack/Processing/P010710_nufft_24rays/Output/flowvalues.study15.slice1.6.1_fixedDelay0.AIF_13_1_5.txt.full.txt';
Ktrans = getKtrans(f);
r1=Ktrans;
f='/v/raid1/npack/Processing/P010710_nufft_24rays/Output/flowvalues.study15.slice2.6.1_fixedDelay0.AIF_13_1_5.txt.full.txt';
Ktrans = getKtrans(f);
r2=Ktrans;
DB_Ktrans_rest1=[r1 r2];DB_REST_avgSlice=[DB_REST_avgSlice mean(r1) mean(r2)];DB_REST_avgPt=[DB_REST_avgPt mean(DB_Ktrans_rest1)];
DB_REST=[DB_REST DB_Ktrans_rest1];


%%%% To load the T1 Ktrans estimates
f='/v/raid1/npack/Processing/P010710_nufft_24rays/Output/flowvalues.study15.slice1.6.1_fixedDelay0.AIF_15_11_1.txt.full.txt';
Ktrans = getKtrans(f);
r1=Ktrans;
f='/v/raid1/npack/Processing/P010710_nufft_24rays/Output/flowvalues.study15.slice2.6.1_fixedDelay0.AIF_15_11_1.txt.full.txt';
Ktrans = getKtrans(f);
r2=Ktrans;
T1_Ktrans_rest1=[r1 r2];T1_REST_avgSlice=[T1_REST_avgSlice mean(r1) mean(r2)];T1_REST_avgPt=[T1_REST_avgPt mean(T1_Ktrans_rest1)];
T1_REST=[T1_REST T1_Ktrans_rest1];

%%%% To load the Adenosine Ktrans estimates
f='/v/raid1/npack/Processing/P010710_nufft_24rays/Output/flowvalues.study20.slice1.6.1_fixedDelay0.AIF_20_11_1.txt.full.txt';
Ktrans = getKtrans(f);
s1=Ktrans;
f='/v/raid1/npack/Processing/P010710_nufft_24rays/Output/flowvalues.study20.slice2.6.1_fixedDelay0.AIF_20_11_1.txt.full.txt';
Ktrans = getKtrans(f);
s2=Ktrans;
Aden_1=[s1 s2];ADEN_avgSlice=[ADEN_avgSlice mean(s1) mean(s2)];ADEN_avgPt=[ADEN_avgPt mean(Aden_1)];
ADEN=[ADEN Aden_1];

%%%% To load the Lexiscan Ktrans estimates
f='/v/raid1/npack/Processing/P010710_nufft_24rays/Output/flowvalues.study57.slice1.6.1_fixedDelay0.AIF_57_11_1.txt.full.txt';
Ktrans = getKtrans(f);
s1=Ktrans;
f='/v/raid1/npack/Processing/P010710_nufft_24rays/Output/flowvalues.study57.slice2.6.1_fixedDelay0.AIF_57_11_1.txt.full.txt';
Ktrans = getKtrans(f);
s2=Ktrans;
Lex_1=[s1 s2];LEX_avgSlice=[LEX_avgSlice mean(s1) mean(s2)];LEX_avgPt=[LEX_avgPt mean(Lex_1)];
LEX=[LEX Lex_1];

REST2=[REST2 T1_Ktrans_rest1];REST2_avgSlice=[REST2_avgSlice mean(T1_Ktrans_rest1(1:6)) mean(T1_Ktrans_rest1(7:12))];REST2_avgPt=[REST2_avgPt mean(T1_Ktrans_rest1)];
figure(1),hold on,plot(DB_Ktrans_rest1,T1_Ktrans_rest1,'bo','LineWidth',3,'MarkerSize',8)
disp('Plot held for current patient; press key to continue.'),pause(.1)











%figure(1),hold on,mpi_twoplot(T1_REST,DB_REST,1);

figure(2),hold on
plot(DB_REST,'b*')
plot(T1_REST,'g sq')
axis([0 500 0 5])
legend('Ktrans (Dual Bolus)','Ktrans (T1)')
xlabel('Number of Regions (all subjects)')
ylabel('Ktrans (ml/min/g)')

disp('Aggregate Dual-Bolus rest Ktrans (ml/min/g)')
mean(DB_REST),std(DB_REST)
disp('Aggregate T1 Equation rest Ktrans (ml/min/g)')
mean(T1_REST),std(T1_REST)



% Ktrans plots
figure(3)
subplot(1,3,1),plot(DB_REST,T1_REST,'bo','LineWidth',3,'MarkerSize',8)
title('Ktrans--6 regions per slice')
xlabel('Dual-Bolus Ktrans Estimates (ml/min/g)')
ylabel('T1 Ktrans Estimates (ml/min/g)')
axis([0 2 0 2]),axis square
subplot(1,3,2),plot(DB_REST_avgSlice,T1_REST_avgSlice,'bo','LineWidth',3,'MarkerSize',8)
title('Ktrans--1 region per slice')
xlabel('Dual-Bolus Ktrans Estimates (ml/min/g)')
ylabel('T1 Ktrans Estimates (ml/min/g)')
axis([0 2 0 2]),axis square
subplot(1,3,3),plot(DB_REST_avgPt,T1_REST_avgPt,'bo','LineWidth',3,'MarkerSize',8)
title('Ktrans--1 region per subject')
xlabel('Dual-Bolus Ktrans Estimates (ml/min/g)')
ylabel('T1 Ktrans Estimates (ml/min/g)')
axis([0 2 0 2]),axis square

MPR_Aden=ADEN./REST2;
MPR_Lex=LEX./REST2;
disp('Aggregate Adenosine MPR')
mean(MPR_Aden),std(MPR_Aden)
disp('Aggregate Lexiscan MPR')
mean(MPR_Lex),std(MPR_Lex)

 
MPR_Aden_avgSlice=ADEN_avgSlice./REST2_avgSlice;
MPR_Lex_avgSlice=LEX_avgSlice./REST2_avgSlice;
MPR_Aden_avgPt=ADEN_avgPt./REST2_avgPt;
MPR_Lex_avgPt=LEX_avgPt./REST2_avgPt;

% MPR plots
figure(4)
subplot(1,3,1),plot(MPR_Aden,MPR_Lex,'bo','LineWidth',3,'MarkerSize',8)
title('MPR--6 regions per slice')
xlabel('Dual-Bolus MPR Estimates')
ylabel('T1 MPR Estimates')
axis([0 10 0 10]),axis square
subplot(1,3,2),plot(MPR_Aden_avgSlice,MPR_Lex_avgSlice,'bo','LineWidth',3,'MarkerSize',8)
title('MPR--1 region per slice')
xlabel('Dual-Bolus MPR Estimates')
ylabel('T1 MPR Estimates')
axis([0 10 0 10]),axis square
subplot(1,3,3),plot(MPR_Aden_avgPt,MPR_Lex_avgPt,'bo','LineWidth',3,'MarkerSize',8)
title('MPR--1 region per subject')
xlabel('Dual-Bolus MPR Estimates')
ylabel('T1 MPR Estimates')
axis([0 10 0 10]),axis square




figure(5),plot(LEX,ADEN,'bo','MarkerSize',6,'LineWidth',3),axis([0 8 0 8]),axis square
xlabel('Lexiscan Ktrans (ml/min/g)'),ylabel('Adenoscan Ktrans (ml/min/g)')

%%% To separate the stress/MPR results by gender
ADEN_male=[ADEN(13:30) ADEN(67:96) ADEN(115:174) ADEN(205:240)];
ADEN_female=[ADEN(1:12) ADEN(31:66) ADEN(97:114) ADEN(175:204) ADEN(241:252)];
LEX_male=[LEX(13:30) LEX(67:96) LEX(115:174) LEX(205:240)];
LEX_female=[LEX(1:12) LEX(31:66) LEX(97:114) LEX(175:204) LEX(241:252)];

figure(6),plot(LEX_male,ADEN_male,'bo','MarkerSize',6,'LineWidth',3)
hold on,plot(LEX_female,ADEN_female,'gs','MarkerSize',6,'LineWidth',3)
legend('Male','Female'),axis([0 5 0 5]),axis square
xlabel('Lexiscan Ktrans (ml/min/g)'),ylabel('Adenoscan Ktrans (ml/min/g)')

MPR_Aden_male=[MPR_Aden(13:30) MPR_Aden(67:96) MPR_Aden(115:174) MPR_Aden(205:240)];
MPR_Aden_female=[MPR_Aden(1:12) MPR_Aden(31:66) MPR_Aden(97:114) MPR_Aden(175:204) MPR_Aden(241:252)];
MPR_Lex_male=[MPR_Lex(13:30) MPR_Lex(67:96) MPR_Lex(115:174) MPR_Lex(205:240)];
MPR_Lex_female=[MPR_Lex(1:12) MPR_Lex(31:66) MPR_Lex(97:114) MPR_Lex(175:204) MPR_Lex(241:252)];

figure(7),plot(MPR_Lex_male,MPR_Aden_male,'bo','MarkerSize',6,'LineWidth',3)
hold on,plot(MPR_Lex_female,MPR_Aden_female,'gs','MarkerSize',6,'LineWidth',3)
legend('Male','Female'),axis([0 8 0 8]),axis square
xlabel('Lexiscan MPR'),ylabel('Adenoscan MPR')

% %%% To separate the stress/MPR results by weight
Weight=zeros(size(LEX));
Weight(1:12)=132;
Weight(13:30)=135;
Weight(31:48)=145;
Weight(49:66)=240;
Weight(67:78)=210;
Weight(79:96)=275;
Weight(97:114)=220;
Weight(115:126)=230;
Weight(127:144)=175;
Weight(145:156)=252;
Weight(157:174)=260;
Weight(175:192)=135;
Weight(193:204)=180;
Weight(205:222)=205;
Weight(223:240)=226;
Weight(241:252)=160;

%%% For now, I've arbitrarily separated the weight categories into 4 ranges as follows:
heavy1=find(Weight>=240);
heavy2=find(Weight>=200 & Weight<240);
light1=find(Weight>145 & Weight<200);
light2=find(Weight<=145);

ADEN_heavy1=ADEN(heavy1);
ADEN_heavy2=ADEN(heavy2);
ADEN_light1=ADEN(light1);
ADEN_light2=ADEN(light2);
MPR_Aden_heavy1=MPR_Aden(heavy1);
MPR_Aden_heavy2=MPR_Aden(heavy2);
MPR_Aden_light1=MPR_Aden(light1);
MPR_Aden_light2=MPR_Aden(light2);

LEX_heavy1=LEX(heavy1);
LEX_heavy2=LEX(heavy2);
LEX_light1=LEX(light1);
LEX_light2=LEX(light2);
MPR_Lex_heavy1=MPR_Lex(heavy1);
MPR_Lex_heavy2=MPR_Lex(heavy2);
MPR_Lex_light1=MPR_Lex(light1);
MPR_Lex_light2=MPR_Lex(light2);

figure(8),plot(LEX_heavy1,ADEN_heavy1,'bo','MarkerSize',8,'LineWidth',3)
hold on,plot(LEX_heavy2,ADEN_heavy2,'gs','MarkerSize',8,'LineWidth',3)
hold on,plot(LEX_light1,ADEN_light1,'r*','MarkerSize',8,'LineWidth',3)
hold on,plot(LEX_light2,ADEN_light2,'c>','MarkerSize',8,'LineWidth',3)
legend('wt: > 240#','wt: 200# -> 240#','wt: 145# -> 200#','wt: < 145#'),axis([0 5 0 5]),axis square
xlabel('Lexiscan Ktrans (ml/min/g)','FontSize',24),ylabel('Adenoscan Ktrans (ml/min/g)','FontSize',24)

figure(100),hold on
myweights = Weight/max(Weight);
mycolors = gray;
for i=1:length(Weight)
    plot(LEX(i),ADEN(i),'d','Color',[mycolors(ceil(myweights(i)*64),1) mycolors(ceil(myweights(i)*64),2) mycolors(ceil(myweights(i)*64),3)]);
end
colormap(mycolors);
colorbar;

figure(101),hold on
for i=1:length(Weight)
    plot(MPR_Lex(i),MPR_Aden(i),'d','Color',[mycolors(ceil(myweights(i)*64),1) mycolors(ceil(myweights(i)*64),2) mycolors(ceil(myweights(i)*64),3)]);
end
colormap(mycolors);
colorbar;
pause();
figure(9),plot(MPR_Lex_heavy1,MPR_Aden_heavy1,'bo','MarkerSize',8,'LineWidth',3)
hold on,plot(MPR_Lex_heavy2,MPR_Aden_heavy2,'gs','MarkerSize',8,'LineWidth',3)
hold on,plot(MPR_Lex_light1,MPR_Aden_light1,'r*','MarkerSize',8,'LineWidth',3)
hold on,plot(MPR_Lex_light2,MPR_Aden_light2,'c>','MarkerSize',8,'LineWidth',3)
legend('wt: > 240#','wt: 200# -> 240#','wt: 145# -> 200#','wt: < 145#'),axis([0 8 0 8]),axis square
xlabel('Lexiscan MPR','FontSize',24),ylabel('Adenoscan MPR','FontSize',24)

figure(10)
subplot(2,2,1),plot(Weight,ADEN,'bo'),xlabel('Subject Weight (lbs)','FontSize',24),ylabel('Adenoscan Flow (ml/min/g)','FontSize',24),axis([120 280 0 8])
subplot(2,2,2),plot(Weight,MPR_Aden,'go'),xlabel('Subject Weight (lbs)','FontSize',24),ylabel('Adenoscan MPR','FontSize',24),axis([120 280 0 10])
subplot(2,2,3),plot(Weight,LEX,'bo'),xlabel('Subject Weight (lbs)','FontSize',24),ylabel('Lexiscan Flow (ml/min/g)','FontSize',24),axis([120 280 0 8])
subplot(2,2,4),plot(Weight,MPR_Lex,'go'),xlabel('Subject Weight (lbs)','FontSize',24),ylabel('Lexiscan MPR)','FontSize',24),axis([120 280 0 10])



%%% Also to separate the stress/MPR results by patient BMI
BMI=zeros(size(LEX));
BMI(1:12)=21.3;
BMI(13:30)=20.8;
BMI(31:48)=25.7;
BMI(49:66)=46.6;
BMI(67:78)=29.3;
BMI(79:96)=39.5;
BMI(97:114)=37.8;
BMI(115:126)=30.5;
BMI(127:144)=25.2;
BMI(145:156)=35.2;
BMI(157:174)=30.9;
BMI(175:192)=19.7;
BMI(193:204)=35.2;
BMI(205:222)=27.1;
BMI(223:240)=27.6;
BMI(241:252)=29.3;

%%% Also for now, I've arbitrarily separated the BMI categories into 4 ranges as follows:
BMIheavy1=find(BMI>=35);
BMIheavy2=find(BMI>=30 & BMI<35);
BMIlight1=find(BMI>25 & BMI<30);
BMIlight2=find(BMI<=25);

ADEN_BMIheavy1=ADEN(BMIheavy1);
ADEN_BMIheavy2=ADEN(BMIheavy2);
ADEN_BMIlight1=ADEN(BMIlight1);
ADEN_BMIlight2=ADEN(BMIlight2);
MPR_Aden_BMIheavy1=MPR_Aden(BMIheavy1);
MPR_Aden_BMIheavy2=MPR_Aden(BMIheavy2);
MPR_Aden_BMIlight1=MPR_Aden(BMIlight1);
MPR_Aden_BMIlight2=MPR_Aden(BMIlight2);

LEX_BMIheavy1=LEX(BMIheavy1);
LEX_BMIheavy2=LEX(BMIheavy2);
LEX_BMIlight1=LEX(BMIlight1);
LEX_BMIlight2=LEX(BMIlight2);
MPR_Lex_BMIheavy1=MPR_Lex(BMIheavy1);
MPR_Lex_BMIheavy2=MPR_Lex(BMIheavy2);
MPR_Lex_BMIlight1=MPR_Lex(BMIlight1);
MPR_Lex_BMIlight2=MPR_Lex(BMIlight2);

figure(11),plot(LEX_BMIheavy1,ADEN_BMIheavy1,'bo','MarkerSize',8,'LineWidth',3)
hold on,plot(LEX_BMIheavy2,ADEN_BMIheavy2,'gs','MarkerSize',8,'LineWidth',3)
hold on,plot(LEX_BMIlight1,ADEN_BMIlight1,'r*','MarkerSize',8,'LineWidth',3)
hold on,plot(LEX_BMIlight2,ADEN_BMIlight2,'c>','MarkerSize',8,'LineWidth',3)
legend('wt: > 240#','wt: 200# -> 240#','wt: 145# -> 200#','wt: < 145#'),axis([0 8 0 8]),axis square
xlabel('Lexiscan Ktrans (ml/min/g)','FontSize',24),ylabel('Adenoscan Ktrans (ml/min/g)','FontSize',24)

figure(12),plot(MPR_Lex_BMIheavy1,MPR_Aden_BMIheavy1,'bo','MarkerSize',8,'LineWidth',3)
hold on,plot(MPR_Lex_BMIheavy2,MPR_Aden_BMIheavy2,'gs','MarkerSize',8,'LineWidth',3)
hold on,plot(MPR_Lex_BMIlight1,MPR_Aden_BMIlight1,'r*','MarkerSize',8,'LineWidth',3)
hold on,plot(MPR_Lex_BMIlight2,MPR_Aden_BMIlight2,'c>','MarkerSize',8,'LineWidth',3)
legend('BMI: > 35','BMI: 30 -> 35#','BMI: 25 -> 30','BMI: < 25'),axis([0 8 0 8]),axis square
xlabel('Lexiscan MPR','FontSize',24),ylabel('Adenoscan MPR','FontSize',24)

figure(13)
subplot(2,2,1),plot(BMI,ADEN,'bo'),xlabel('Subject BMI','FontSize',24),ylabel('Adenoscan Flow (ml/min/g)','FontSize',24),axis([15 50 0 8])
subplot(2,2,2),plot(BMI,MPR_Aden,'go'),xlabel('Subject BMI','FontSize',24),ylabel('Adenoscan MPR','FontSize',24),axis([15 50 0 10])
subplot(2,2,3),plot(BMI,LEX,'bo'),xlabel('Subject BMI','FontSize',24),ylabel('Lexiscan Flow (ml/min/g)','FontSize',24),axis([15 50 0 8])
subplot(2,2,4),plot(BMI,MPR_Lex,'go'),xlabel('Subject BMI','FontSize',24),ylabel('Lexiscan MPR)','FontSize',24),axis([15 50 0 10])
% 






%%%%% Below is only used if ALL subjects are included.  Since P080609 had a
%%%%% technical problem with the Gd injection, it's excluded here.
% %%% To separate the stress/MPR results by gender
% ADEN_male=[ADEN(13:30) ADEN(67:108) ADEN(127:186) ADEN(217:252)];
% ADEN_female=[ADEN(1:12) ADEN(31:66) ADEN(109:126) ADEN(187:216) ADEN(253:264)];
% LEX_male=[LEX(13:30) LEX(67:108) LEX(127:186) LEX(217:252)];
% LEX_female=[LEX(1:12) LEX(31:66) LEX(109:126) LEX(187:216) LEX(253:264)];
% 
% figure,plot(LEX_male,ADEN_male,'bo','MarkerSize',6,'LineWidth',3)
% hold on,plot(LEX_female,ADEN_female,'gs','MarkerSize',6,'LineWidth',3)
% legend('Male','Female'),axis([0 8 0 8]),axis square
% xlabel('Lexiscan Ktrans (ml/min/g)'),ylabel('Adenoscan Ktrans (ml/min/g)')
% 
% MPR_Aden_male=[MPR_Aden(13:30) MPR_Aden(67:108) MPR_Aden(127:186) MPR_Aden(217:252)];
% MPR_Aden_female=[MPR_Aden(1:12) MPR_Aden(31:66) MPR_Aden(109:126) MPR_Aden(187:216) MPR_Aden(253:264)];
% MPR_Lex_male=[MPR_Lex(13:30) MPR_Lex(67:108) MPR_Lex(127:186) MPR_Lex(217:252)];
% MPR_Lex_female=[MPR_Lex(1:12) MPR_Lex(31:66) MPR_Lex(109:126) MPR_Lex(187:216) MPR_Lex(253:264)];
% 
% figure,plot(MPR_Lex_male,MPR_Aden_male,'bo','MarkerSize',6,'LineWidth',3)
% hold on,plot(MPR_Lex_female,MPR_Aden_female,'gs','MarkerSize',6,'LineWidth',3)
% legend('Male','Female'),axis([0 8 0 8]),axis square
% xlabel('Lexiscan MPR'),ylabel('Adenoscan MPR')
% 
% % %%% To separate the stress/MPR results by weight
% Weight=zeros(size(LEX));
% Weight(1:12)=132;
% Weight(13:30)=135;
% Weight(31:48)=145;
% Weight(49:66)=240;
% Weight(67:78)=210;
% Weight(79:90)=210;
% Weight(91:108)=275;
% Weight(109:126)=220;
% Weight(127:138)=230;
% Weight(139:156)=175;
% Weight(157:168)=252;
% Weight(169:186)=260;
% Weight(187:204)=135;
% Weight(205:216)=180;
% Weight(217:234)=205;
% Weight(235:252)=226;
% Weight(253:264)=160;
% 
% %%% For now, I've arbitrarily separated the weight categories into 4 ranges as follows:
% heavy1=find(Weight>=240);
% heavy2=find(Weight>=200 & Weight<240);
% light1=find(Weight>145 & Weight<200);
% light2=find(Weight<=145);
% 
% ADEN_heavy1=ADEN(heavy1);
% ADEN_heavy2=ADEN(heavy2);
% ADEN_light1=ADEN(light1);
% ADEN_light2=ADEN(light2);
% MPR_Aden_heavy1=MPR_Aden(heavy1);
% MPR_Aden_heavy2=MPR_Aden(heavy2);
% MPR_Aden_light1=MPR_Aden(light1);
% MPR_Aden_light2=MPR_Aden(light2);
% 
% LEX_heavy1=LEX(heavy1);
% LEX_heavy2=LEX(heavy2);
% LEX_light1=LEX(light1);
% LEX_light2=LEX(light2);
% MPR_Lex_heavy1=MPR_Lex(heavy1);
% MPR_Lex_heavy2=MPR_Lex(heavy2);
% MPR_Lex_light1=MPR_Lex(light1);
% MPR_Lex_light2=MPR_Lex(light2);
% 
% figure,plot(LEX_heavy1,ADEN_heavy1,'bo','MarkerSize',8,'LineWidth',3)
% hold on,plot(LEX_heavy2,ADEN_heavy2,'gs','MarkerSize',8,'LineWidth',3)
% hold on,plot(LEX_light1,ADEN_light1,'r*','MarkerSize',8,'LineWidth',3)
% hold on,plot(LEX_light2,ADEN_light2,'c>','MarkerSize',8,'LineWidth',3)
% legend('wt: > 240#','wt: 200# -> 240#','wt: 145# -> 200#','wt: < 145#'),axis([0 8 0 8]),axis square
% xlabel('Lexiscan Ktrans (ml/min/g)','FontSize',24),ylabel('Adenoscan Ktrans (ml/min/g)','FontSize',24)
% 
% figure,plot(MPR_Lex_heavy1,MPR_Aden_heavy1,'bo','MarkerSize',8,'LineWidth',3)
% hold on,plot(MPR_Lex_heavy2,MPR_Aden_heavy2,'gs','MarkerSize',8,'LineWidth',3)
% hold on,plot(MPR_Lex_light1,MPR_Aden_light1,'r*','MarkerSize',8,'LineWidth',3)
% hold on,plot(MPR_Lex_light2,MPR_Aden_light2,'c>','MarkerSize',8,'LineWidth',3)
% legend('wt: > 240#','wt: 200# -> 240#','wt: 145# -> 200#','wt: < 145#'),axis([0 8 0 8]),axis square
% xlabel('Lexiscan MPR','FontSize',24),ylabel('Adenoscan MPR','FontSize',24)
% 
% figure
% subplot(2,2,1),plot(Weight,ADEN,'bo'),xlabel('Subject Weight (lbs)','FontSize',24),ylabel('Adenoscan Flow (ml/min/g)','FontSize',24),axis([120 280 0 8])
% subplot(2,2,2),plot(Weight,MPR_Aden,'go'),xlabel('Subject Weight (lbs)','FontSize',24),ylabel('Adenoscan MPR','FontSize',24),axis([120 280 0 10])
% subplot(2,2,3),plot(Weight,LEX,'bo'),xlabel('Subject Weight (lbs)','FontSize',24),ylabel('Lexiscan Flow (ml/min/g)','FontSize',24),axis([120 280 0 8])
% subplot(2,2,4),plot(Weight,MPR_Lex,'go'),xlabel('Subject Weight (lbs)','FontSize',24),ylabel('Lexiscan MPR)','FontSize',24),axis([120 280 0 10])




%%% Bland-Altman plots
figure(14)
subplot(1,3,1),plot((MPR_Aden + MPR_Lex)/2,(MPR_Aden - MPR_Lex),'bo')
% a=mean(MPR_Aden-MPR_Lex);
% a2=1.96*std(a);
% a3=0:0.1:10;
% a4=
title('MPR--6 regions per slice')
xlabel('(MPR_L_e_x + MPR_A_d_e_n)/2')
ylabel('MPR_A_d_e_n - MPR_L_e_x')
axis([0 10 -6 6]),axis square




subplot(1,3,2),plot((MPR_Aden_avgSlice + MPR_Lex_avgSlice)/2,(MPR_Aden_avgSlice - MPR_Lex_avgSlice),'bo')
title('MPR--1 region per slice')
xlabel('(MPR_L_e_x + MPR_A_d_e_n)/2')
ylabel('MPR_A_d_e_n - MPR_L_e_x')
axis([0 10 -6 6]),axis square
subplot(1,3,3),plot((MPR_Aden_avgPt + MPR_Lex_avgPt)/2,(MPR_Aden_avgPt - MPR_Lex_avgPt),'bo')
title('MPR--1 region per subject')
xlabel('(MPR_L_e_x + MPR_A_d_e_n)/2')
ylabel('MPR_A_d_e_n - MPR_L_e_x')
axis([0 10 -6 6]),axis square

disp('All high dose studies (ie. 0.06mmol/kg) omitted,')
disp('P080609 is omitted since Gd spilled and late Aden response')
disp('Removed 1 slice from P090309 since pt. had PVCs and bad response due to gating')


%%%% Note the equations of best fit below and r-values were taken from the
%%%% Excel spreadsheet "Astellas_Ktrans-Values_011910.xls"
figure(15),hold on
plot(MPR_Lex,MPR_Aden,'bo','MarkerSize',8,'LineWidth',3)
xval=0:0.1:10;    yval=0.8321*xval+0.5534;
plot(xval,yval,'k','LineWidth',2);
rval=sqrt(.4868);
text(5,1.5,'y=0.83x + 0.55')
text(5,1,'r=0.70')
axis([0 11 0 11]),axis square
title('MPR: All regions')
xlabel('MPR (Regad)')
ylabel('MPR (Adeno)')


figure(16),hold on
plot(MPR_Lex_avgPt,MPR_Aden_avgPt,'bo','MarkerSize',8,'LineWidth',3)
xval=0:0.1:10;    yval=0.9267*xval+0.5012;
plot(xval,yval,'k','LineWidth',2)
rval=sqrt(.5012);
text(2.5,1,'y=0.93x + 0.42')
text(2.5,0.7,'r=0.71','FontSize',20,'FontWeight','bold')
axis([0 7 0 7]),axis square
title('MPR: 1 Mean value per subject')
xlabel('MPR (Regad)')
ylabel('MPR (Adeno)')




figure(17),axis([0.8 2.2 0 8]),ylabel('MPR')
x=1:2;
for i=1:length(MPR_Lex_avgPt),tmp1(i,:)=[MPR_Lex_avgPt(i) MPR_Aden_avgPt(i)];end
hold on,plot(x,tmp1,'-*','LineWidth',4,'MarkerSize',10)


RPP_Ad=[7992 8308 8580 15484 19140 9792 16720 13965 11431 21625 15224 11931 15000 12432 13310]; % no value for P121609, also P080609 is omitted due to spilled Gd
RPP_Lx=[7770 8848 14580 19760 15865 7980 16102 15984 13760 20406 12096 13362 16256 12584 15400]; % no value for P121609, also P080609 is omitted due to spilled Gd
figure(18),axis([0.8 2.2 5000 25000]),ylabel('RPP')
x=1:2;
for i=1:length(RPP_Ad),tmp2(i,:)=[RPP_Ad(i) RPP_Lx(i)];end
hold on,plot(x,tmp2,'-*','LineWidth',4,'MarkerSize',10)
