%% Data Import

fname = '07092024_Fn at 50 ug per ml at pH 4 + dSF pH4 then pH 7 PBS wash on func. Au at 25C_n=4_KL.csv';
%file name listed above is known to work with the script
dataImport = importdata(fname);

%% Data Organization
idx=[1,14,15,29,30,44,45,59,60,74,75,89,90,104,105]; 
selected_data=dataImport.data(:,idx); % data for important parameters (time, changes in dissipation and frequency)

% Determining if data is missing
threshold_dataselected=(selected_data(:,:));
thresdata=true(size(threshold_dataselected));
thresdata(isnan(threshold_dataselected))=false; %if data is there, value should not equal zero
sumCols=sum(thresdata,1);

[numRows,numCols]=size(threshold_dataselected);

for i=1:numCols
    if sumCols(1,i)==numRows
        dataFinal(:,i)=selected_data(:,i);
    else 
        disp('There is missing data for column:')
        disp(i)
    end
end


%dataFinal=selected_data(:,nonzeroLocation);
dataFinal(:,1)=dataFinal(:,1)./60; % converting seconds to minutes for time

%p = baseline
%q= ECM protein or precursor protein (Fn or Col-II) 
%r = Precursor film postwash 
%w = optional buffer pH change (pH 4 --> pH 7)
%s = BSA
%t = BSA postwash
%u = DSF
%v = DSF postwash

figure
plot (dataFinal(:,1),dataFinal(:,5),'r-')

title ('Frequency data for 3rd overtone against time')
xlabel('Time (s)')
ylabel('Frequency change')


%colHeadersfreq=['delta f1';'delta f3';'delta f5';'delta f7';'delta f9';'delta f11';'delta f13']';
%frequency changes in these columns
freqCols=[3,5,7,9,11,13,15];
freq_df=[dataFinal(:,freqCols)];

%time in this column
time=dataFinal(:,1);


%dissipation changes in these columns
dissCols=[2,4,6,8,10,12,14];
diss_df=[dataFinal(:,dissCols)];
%% Change in Frequency to find where steps change
%finding the changes in freq drops

freqFind=[freq_df(:,2)];
figure
plot(time,freqFind, 'b-')
hold on
title('Frequency Against Time with Overlay of Time Stamps for Steps')
xlabel('Time(s)')


baseline=1:125;% ranges corrpeond to row numbers
ECM=418:447;
precursor=607:647;
bufferChange=765:831;
BSA=995:1019;
BSAwash=1152:1189;
dSF=1484:1546;
dSFwash=1679:1727;
% baseline=1:64; 
% ECM=490:520;
% precursor=672:712;
% bufferChange=832:903;
% BSA=1076:1093;
% BSAwash=1190:1249;
% dSF=1540:1606;
% dSFwash=1769:1817;

plot(time(baseline),freqFind(baseline),'r-')
plot(time(ECM),freqFind(ECM),'g-')
plot(time(precursor),freqFind(precursor),'k-')
plot(time(bufferChange),freqFind(bufferChange),'c-')
plot(time(BSA),freqFind(BSA),'y-')
plot(time(BSAwash),freqFind(BSAwash),'r-')
plot(time(dSF),freqFind(dSF),'m-')
plot(time(dSFwash),freqFind(dSFwash),'bx')

baselineFreq=freq_df(baseline,:);
ECMfreq=freq_df(ECM,:);
precursorFreq=freq_df(precursor,:);
bufferChangeFreq=freq_df(bufferChange,:);
BSAfreq=freq_df(BSA,:);
BSAwashFreq=freq_df(BSAwash,:);
dSFfreq=freq_df(dSF,:);
dSFwashFreq=freq_df(dSFwash,:);

baselineDis=diss_df(baseline,:);
ECMdis=diss_df(ECM,:);
precursorDis=diss_df(precursor,:);
bufferChangeDis=diss_df(bufferChange,:);
BSAdis=diss_df(BSA,:);
BSAwashDis=diss_df(BSAwash,:);
dSFdis=diss_df(dSF,:);
dSFwashDis=diss_df(dSFwash,:);

%% Frequency Calculations

%baseline 
p=mean(baselineFreq,1);
p_std=std(baselineFreq,1);

%Fn
q=mean(ECMfreq,1);
q_std=std(ECMfreq,1);

%wash
r=mean(precursorFreq,1);
r_std=std(precursorFreq,1);

%buffer change
w=mean(bufferChangeFreq,1);
w_std=std(bufferChangeFreq,1);

%BSA 
s=mean(BSAfreq,1);
s_std=std(BSAfreq,1);

%BSA post wash
t=mean(BSAwashFreq,1);
t_std=std(BSAwashFreq,1);

%dSF
u=mean(dSFfreq,1);
u_std=std(dSFfreq,1);

%dSF wash
v=mean(dSFwashFreq,1);
v_std=std(dSFwashFreq,1);

%amount of Fn adsorption prewash and postwash
baselineAdsorptionAvg=p;
baselineAdsorptionSTD=p_std;

precursor_film_adsorption_prewash=q-p;
precursor_film_adsorption_postwash=r-p;

BSA_adsorption_prewash=s-w;
BSA_adsorption_postwash=t-w;

dSF_adsorption_prewash=u-t;
dSF_adsorption_postwash=v-t;


overtones= [1,3,5,7,9,11,13];

df=[overtones;baselineAdsorptionAvg;precursor_film_adsorption_prewash;precursor_film_adsorption_postwash;...
    BSA_adsorption_prewash;BSA_adsorption_postwash;dSF_adsorption_prewash;dSF_adsorption_postwash];
%% Dissipation Calculations

%baseline 
pp=mean(baselineDis,1);
pp_std=std(baselineDis,1);

%Fn
qq=mean(ECMdis,1);
qq_std=std(ECMdis,1);

%wash
rr=mean(precursorDis,1);
rr_std=std(precursorDis,1);

%buffer change
ww=mean(bufferChangeDis,1);
ww_std=std(bufferChangeDis,1);

%BSA 
ss=mean(BSAdis,1);
ss_std=std(BSAdis,1);

%BSA post wash
tt=mean(BSAwashDis,1);
tt_std=std(BSAwashDis,1);

%dSF
uu=mean(dSFdis,1);
uu_std=std(dSFdis,1);

%dSF wash
vv=mean(dSFwashDis,1);
vv_std=std(dSFwashDis,1);

%amount of Fn adsorption prewash and postwash
baselineDiss=pp;

precursor_film_diss_prewash=qq-pp;
precursor_film_diss_postwash=rr-pp;

BSA_diss_prewash=ss-ww;
BSA_diss_postwash=tt-ww;

dSF_diss_prewash=uu-tt;
dSF_diss_postwash=vv-tt;


dd=[baselineDiss;precursor_film_diss_prewash;precursor_film_diss_postwash;...
    BSA_diss_prewash;BSA_diss_postwash;dSF_diss_prewash;dSF_diss_postwash];
%% Sauerbrey Mass
%fundamental freq
fundFreq=freq_df(:,1);
mf_fundamental=-17.7.*(fundFreq/1);

%3rd overtone
thirdFreq=freq_df(:,2);
mf_third=-17.7.*(thirdFreq/3);

%5th overtone
fifthFreq=freq_df(:,3);
mf_fifth=-17.7.*(fifthFreq/5);

%7th overtone
seventhFreq=freq_df(:,4);
mf_seventh=-17.7.*(seventhFreq/7);

%9th overtone
ninthFreq=freq_df(:,5);
mf_ninth=-17.7.*(ninthFreq/9);

%11th overtone
eleventhFreq=freq_df(:,6);
mf_eleventh=-17.7.*(eleventhFreq/11);

%13th overtone
thirteenthFreq=freq_df(:,7);
mf_thirteenth=-17.7.*(thirteenthFreq/13);


%fundamental avg std mf
baseline_mf_fund_mean=mean(mf_fundamental(baseline,:));
ECM_mf_fund_mean=mean(mf_fundamental(ECM,:));
precursor_mf_fund_mean=mean(mf_fundamental(precursor,:));
bufferChange_mf_fund_mean=mean(mf_fundamental(bufferChange,:));
BSA_mf_fund_mean=mean(mf_fundamental(BSA,:));
BSAwash_mf_fund_mean=mean(mf_fundamental(BSAwash,:));
dSF_mf_fund_mean=mean(mf_fundamental(dSF,:));
dSFwash_mf_fund_mean=mean(mf_fundamental(dSFwash,:));

baseline_mf_fund_std=std(mf_fundamental(baseline,:));
ECM_mf_fund_std=std(mf_fundamental(ECM,:));
precursor_mf_fund_std=std(mf_fundamental(precursor,:));
bufferChange_mf_fund_std=std(mf_fundamental(bufferChange,:));
BSA_mf_fund_std=std(mf_fundamental(BSA,:));
BSAwash_mf_fund_std=std(mf_fundamental(BSAwash,:));
dSF_mf_fund_std=std(mf_fundamental(dSF,:));
dSFwash_mf_fund_std=std(mf_fundamental(dSFwash,:));

%3rd overtone
baseline_mf_third_mean=mean(mf_third(baseline,:));
ECM_mf_third_mean=mean(mf_third(ECM,:));
precursor_mf_third_mean=mean(mf_third(precursor,:));
bufferChange_mf_third_mean=mean(mf_third(bufferChange,:));
BSA_mf_third_mean=mean(mf_third(BSA,:));
BSAwash_mf_third_mean=mean(mf_third(BSAwash,:));
dSF_mf_third_mean=mean(mf_third(dSF,:));
dSFwash_mf_third_mean=mean(mf_third(dSFwash,:));

baseline_mf_third_std=std(mf_third(baseline,:));
ECM_mf_third_std=std(mf_third(ECM,:));
precursor_mf_third_std=std(mf_third(precursor,:));
bufferChange_mf_third_std=std(mf_third(bufferChange,:));
BSA_mf_third_std=std(mf_third(BSA,:));
BSAwash_mf_third_std=std(mf_third(BSAwash,:));
dSF_mf_third_std=std(mf_third(dSF,:));
dSFwash_mf_third_std=std(mf_third(dSFwash,:));

%5th overtone
baseline_mf_fifth_mean=mean(mf_fifth(baseline,:));
ECM_mf_fifth_mean=mean(mf_fifth(ECM,:));
precursor_mf_fifth_mean=mean(mf_fifth(precursor,:));
bufferChange_mf_fifth_mean=mean(mf_fifth(bufferChange,:));
BSA_mf_fifth_mean=mean(mf_fifth(BSA,:));
BSAwash_mf_fifth_mean=mean(mf_fifth(BSAwash,:));
dSF_mf_fifth_mean=mean(mf_fifth(dSF,:));
dSFwash_mf_fifth_mean=mean(mf_fifth(dSFwash,:));

baseline_mf_fifth_std=std(mf_fifth(baseline,:));
ECM_mf_fifth_std=std(mf_fifth(ECM,:));
precursor_mf_fifth_std=std(mf_fifth(precursor,:));
bufferChange_mf_fifth_std=std(mf_fifth(bufferChange,:));
BSA_mf_fifth_std=std(mf_fifth(BSA,:));
BSAwash_mf_fifth_std=std(mf_fifth(BSAwash,:));
dSF_mf_fifth_std=std(mf_fifth(dSF,:));
dSFwash_mf_fifth_std=std(mf_fifth(dSFwash,:));


%7th overtone
baseline_mf_seventh_mean=mean(mf_seventh(baseline,:));
ECM_mf_seventh_mean=mean(mf_seventh(ECM,:));
precursor_mf_seventh_mean=mean(mf_seventh(precursor,:));
bufferChange_mf_seventh_mean=mean(mf_seventh(bufferChange,:));
BSA_mf_seventh_mean=mean(mf_seventh(BSA,:));
BSAwash_mf_seventh_mean=mean(mf_seventh(BSAwash,:));
dSF_mf_seventh_mean=mean(mf_seventh(dSF,:));
dSFwash_mf_seventh_mean=mean(mf_fifth(dSFwash,:));

baseline_mf_seventh_std=std(mf_seventh(baseline,:));
ECM_mf_seventh_std=std(mf_seventh(ECM,:));
precursor_mf_seventh_std=std(mf_seventh(precursor,:));
bufferChange_mf_seventh_std=std(mf_seventh(bufferChange,:));
BSA_mf_seventh_std=std(mf_seventh(BSA,:));
BSAwash_mf_seventh_std=std(mf_seventh(BSAwash,:));
dSF_mf_seventh_std=std(mf_seventh(dSF,:));
dSFwash_mf_seventh_std=std(mf_seventh(dSFwash,:));


%9th overtone
baseline_mf_ninth_mean=mean(mf_ninth(baseline,:));
ECM_mf_ninth_mean=mean(mf_ninth(ECM,:));
precursor_mf_ninth_mean=mean(mf_ninth(precursor,:));
bufferChange_mf_ninth_mean=mean(mf_ninth(bufferChange,:));
BSA_mf_ninth_mean=mean(mf_ninth(BSA,:));
BSAwash_mf_ninth_mean=mean(mf_ninth(BSAwash,:));
dSF_mf_ninth_mean=mean(mf_ninth(dSF,:));
dSFwash_mf_ninth_mean=mean(mf_ninth(dSFwash,:));

baseline_mf_ninth_std=std(mf_ninth(baseline,:));
ECM_mf_ninth_std=std(mf_ninth(ECM,:));
precursor_mf_ninth_std=std(mf_ninth(precursor,:));
bufferChange_mf_ninth_std=std(mf_ninth(bufferChange,:));
BSA_mf_ninth_std=std(mf_ninth(BSA,:));
BSAwash_mf_ninth_std=std(mf_ninth(BSAwash,:));
dSF_mf_ninth_std=std(mf_ninth(dSF,:));
dSFwash_mf_ninth_std=std(mf_ninth(dSFwash,:));


%11th overtone
baseline_mf_eleventh_mean=mean(mf_eleventh(baseline,:));
ECM_mf_eleventh_mean=mean(mf_eleventh(ECM,:));
precursor_mf_eleventh_mean=mean(mf_eleventh(precursor,:));
bufferChange_mf_eleventh_mean=mean(mf_eleventh(bufferChange,:));
BSA_mf_eleventh_mean=mean(mf_eleventh(BSA,:));
BSAwash_mf_eleventh_mean=mean(mf_eleventh(BSAwash,:));
dSF_mf_eleventh_mean=mean(mf_eleventh(dSF,:));
dSFwash_mf_eleventh_mean=mean(mf_eleventh(dSFwash,:));

baseline_mf_eleventh_std=std(mf_eleventh(baseline,:));
ECM_mf_eleventh_std=std(mf_eleventh(ECM,:));
precursor_mf_eleventh_std=std(mf_eleventh(precursor,:));
bufferChange_mf_eleventh_std=std(mf_eleventh(bufferChange,:));
BSA_mf_eleventh_std=std(mf_eleventh(BSA,:));
BSAwash_mf_eleventh_std=std(mf_eleventh(BSAwash,:));
dSF_mf_eleventh_std=std(mf_eleventh(dSF,:));
dSFwash_mf_eleventh_std=std(mf_eleventh(dSFwash,:));


%13th overtone 
baseline_mf_thirteenth_mean=mean(mf_thirteenth(baseline,:));
ECM_mf_thirteenth_mean=mean(mf_thirteenth(ECM,:));
precursor_mf_thirteenth_mean=mean(mf_thirteenth(precursor,:));
bufferChange_mf_thirteenth_mean=mean(mf_thirteenth(bufferChange,:));
BSA_mf_thirteenth_mean=mean(mf_thirteenth(BSA,:));
BSAwash_mf_thirteenth_mean=mean(mf_thirteenth(BSAwash,:));
dSF_mf_thirteenth_mean=mean(mf_thirteenth(dSF,:));
dSFwash_mf_thirteenth_mean=mean(mf_thirteenth(dSFwash,:));

baseline_mf_thirteenth_std=std(mf_thirteenth(baseline,:));
ECM_mf_thirteenth_std=std(mf_thirteenth(ECM,:));
precursor_mf_thirteenth_std=std(mf_thirteenth(precursor,:));
bufferChange_mf_thirteenth_std=std(mf_thirteenth(bufferChange,:));
BSA_mf_thirteenth_std=std(mf_thirteenth(BSA,:));
BSAwash_mf_thirteenth_std=std(mf_thirteenth(BSAwash,:));
dSF_mf_thirteenth_std=std(mf_thirteenth(dSF,:));
dSFwash_mf_thirteenth_std=std(mf_thirteenth(dSFwash,:));


%Sauerbrey mass=(-SLOPE(Y8:Y14,X8:X14)*17.7)
% baseline df
% ECM
% precursor
% bufferChange
% BSA
% BSAwash
% dSF
% dSFwash

overtone_select=[1,2,3,4,6,7]; %because ninth freq lost, have to find way to not hard code here

%Fn prewash mass
X=df(1,overtone_select); %overtones
Y=df(3,overtone_select);
z = polyfit(X,Y,1);
mf_slope_Fnprewash=z(1);

Sauerbrey_mass_adsorption_precursor_film_prewash = mf_slope_Fnprewash*(-17.7); 

%Fn postwash mass
X=df(1,overtone_select); %overtones
Y=df(4,overtone_select); 
z2 = polyfit(X,Y,1);
mf_slope_Fnpostwash=z2(1);

Sauerbrey_mass_adsorption_precursor_film_postwash = mf_slope_Fnpostwash*(-17.7); 


%BSA prewash mass
X=df(1,overtone_select); %overtones
Y=df(5,overtone_select); 
z3 = polyfit(X,Y,1);
mf_slope_BSAprewash=z3(1);

Sauerbrey_mass_adsorption_BSA_film_prewash = mf_slope_BSAprewash*(-17.7); 

%BSA postwash mass
X=df(1,overtone_select); %overtones
Y=df(6,overtone_select); 
z4 = polyfit(X,Y,1);
mf_slope_BSApostwash=z4(1);

Sauerbrey_mass_adsorption_BSA_film_postwash = mf_slope_BSApostwash*(-17.7); 

%dSF prewash mass
X=df(1,overtone_select); %overtones
Y=df(7,overtone_select); 
z5 = polyfit(X,Y,1);
mf_slope_dSFprewash=z5(1);

Sauerbrey_mass_adsorption_dSF_film_prewash = mf_slope_dSFprewash*(-17.7); 

%dSF postwash mass
X=df(1,overtone_select); %overtones
Y=df(8,overtone_select); 
z6 = polyfit(X,Y,1);
mf_slope_dSFpostwash=z6(1);

Sauerbrey_mass_adsorption_dSF_film_postwash = mf_slope_dSFpostwash*(-17.7); 

Mf_titles={'Sauerbrey Mass Fn prewash (ng/cm^2)' 'Sauerbrey Mass Fn postwash (ng/cm^2)'...
    'Sauerbrey Mass dSF prewash (ng/cm^2)' 'Sauebrey Mass dSF postwash (ng/cm^2)'};
Mf_data=[Sauerbrey_mass_adsorption_precursor_film_prewash, Sauerbrey_mass_adsorption_precursor_film_postwash,...
    Sauerbrey_mass_adsorption_dSF_film_prewash,Sauerbrey_mass_adsorption_dSF_film_postwash];
Mf= [Mf_titles; num2cell(Mf_data)];


writecell(Mf,'SauerbreyMass.csv'); %converts Mf cell to a csv file in your current folder

disp(['Sauerbrey Mass Fn prewash',' ',num2str(Sauerbrey_mass_adsorption_precursor_film_prewash),' ','(ng/cm^2)'])
disp(['Sauerbrey Mass Fn postwash',' ',num2str(Sauerbrey_mass_adsorption_precursor_film_postwash),' ','(ng/cm^2)'])
disp(['Sauerbrey Mass dSF prewash',' ',num2str(Sauerbrey_mass_adsorption_dSF_film_prewash),' ','(ng/cm^2)'])
disp(['Sauerbrey Mass dSF postwash',' ',num2str(Sauerbrey_mass_adsorption_dSF_film_postwash),' ','(ng/cm^2)'])
%% Jf prime Calculations
%Jf= slope/(2*pi*n*1)*10^-3
%slope of del bandwidth & ndelF

%ndelF=n*df
ndel_fund=1*(df(2:end,1));
ndel_third=3*(df(2:end,2));
ndel_fifth=5*(df(2:end,3));
ndel_seventh=7*(df(2:end,4));
ndel_ninth=9*(df(2:end,5));
ndel_eleventh=11*(df(2:end,6));
ndel_thirteenth=13*(df(2:end,7));


del_bandwidth_fund=1*(dd(:,1))*(5*10^6/2);
del_bandwidth_third=3*(dd(:,2))*(5*10^6/2);
del_bandwidth_fifth=5*(dd(:,3))*(5*10^6/2);
del_bandwidth_seventh=7*(dd(:,4))*(5*10^6/2);
del_bandwidth_ninth=9*(dd(:,5))*(5*10^6/2);
del_bandwidth_eleventh=11*(dd(:,6))*(5*10^6/2);
del_bandwidth_thirteenth=13*(dd(:,7))*(5*10^6/2);

del_freqn_fund=(df(2:end,1))/1;
del_freqn_third=(df(2:end,2))/3;
del_freqn_fifth=(df(2:end,3))/5;
del_freqn_seventh=(df(2:end,4))/7;
del_freqn_ninth=(df(2:end,5))/9;
del_freqn_eleventh=(df(2:end,6))/11;
del_freqn_thirteenth=(df(2:end,7))/13;

del_bandn_fund=(del_bandwidth_fund)/1;
del_bandn_third=(del_bandwidth_third)/3;
del_bandn_fifth=(del_bandwidth_fifth)/5;
del_bandn_seventh=(del_bandwidth_seventh)/7;
del_bandn_ninth=(del_bandwidth_ninth)/9;
del_bandn_eleventh=(del_bandwidth_eleventh)/11;
del_bandn_thirteenth=(del_bandwidth_thirteenth)/13;

%if any dissipation values are missing (i.e. missing data from data
%organization section) it is important that it not be included in the
%del_bandcondition or ndel_condition as it will skew the Jf prime results

%Fn prewash
del_bandFnprewash=[del_bandwidth_third(2,1),del_bandwidth_fifth(2,1),...
    del_bandwidth_seventh(2,1),del_bandwidth_eleventh(2,1),del_bandwidth_thirteenth(2,1)];
ndel_Fnprewash=[ndel_third(2,1),ndel_fifth(2,1),ndel_seventh(2,1),ndel_eleventh(2,1),ndel_thirteenth(2,1)];
Fnprewash= polyfit(ndel_Fnprewash,del_bandFnprewash,1);
slope_Fnprewash=Fnprewash(1);
Jf_Fnprewash=-(slope_Fnprewash/((2*pi*5))*(10^-3));

%Fn postwash
del_bandFnpostwash=[del_bandwidth_third(3,1),del_bandwidth_fifth(3,1),...
    del_bandwidth_seventh(3,1),del_bandwidth_eleventh(3,1),del_bandwidth_thirteenth(3,1)];
ndel_Fnpostwash=[ndel_third(3,1),ndel_fifth(3,1),ndel_seventh(3,1),ndel_eleventh(3,1),ndel_thirteenth(3,1)];
Fnpostwash= polyfit(ndel_Fnpostwash,del_bandFnpostwash,1);
slope_Fnpostwash=Fnpostwash(1);
Jf_Fnpostwash=-(slope_Fnpostwash/((2*pi*5))*(10^-3));

%dSF prewash
del_bandDSFprewash=[del_bandwidth_third(6,1),del_bandwidth_fifth(6,1),...
    del_bandwidth_seventh(6,1),del_bandwidth_eleventh(6,1),del_bandwidth_thirteenth(6,1)];
ndel_DSFprewash=[ndel_third(6,1),ndel_fifth(6,1),ndel_seventh(6,1),ndel_eleventh(6,1),ndel_thirteenth(6,1)];
DSFprewash= polyfit(ndel_DSFprewash,del_bandDSFprewash,1);
slope_DSFprewash=DSFprewash(1);
Jf_DSFprewash=-(slope_DSFprewash/((2*pi*5))*(10^-3));


%dSF postwash
del_bandDSFpostwash=[del_bandwidth_third(7,1),del_bandwidth_fifth(7,1),...
    del_bandwidth_seventh(7,1),del_bandwidth_eleventh(7,1),del_bandwidth_thirteenth(7,1)];
ndel_DSFpostwash=[ndel_third(7,1),ndel_fifth(7,1),ndel_seventh(7,1),ndel_eleventh(7,1),ndel_thirteenth(7,1)];
DSFpostwash= polyfit(ndel_DSFpostwash,del_bandDSFpostwash,1);
slope_DSFpostwash=DSFpostwash(1);
Jf_DSFpostwash=-(slope_DSFpostwash/((2*pi*5))*(10^-3));


Jf_titles={'Jf prime Fn prewash (1/Pa)' 'Jf prime Fn postwash (1/Pa)' 'Jf prime dSF prewash (1/Pa)' 'Jf prime dSF postwash (1/Pa)'};
Jf_data=[Jf_Fnprewash, Jf_Fnpostwash,Jf_DSFprewash,Jf_DSFpostwash];
Jf= [Jf_titles; num2cell(Jf_data)];


writecell(Jf,'ThinFilmCompliance.csv'); %converts Jf cell to a csv file in your current folder

disp(['Jf prime Fn prewash ',' ',num2str(Jf_Fnprewash),' ','(1/Pa)'])
disp(['Jf prime Fn postwash',' ',num2str(Jf_Fnpostwash),' ','(1/Pa)'])
disp(['Jf prime dSF prewash',' ',num2str(Jf_DSFprewash),' ','(1/Pa)'])
disp(['Jf prime dSF postwash',' ',num2str(Jf_DSFpostwash),' ','(1/Pa)'])



