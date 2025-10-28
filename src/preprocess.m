clc,clear,close all

%%%%%%%%%%%%%%%%%%% Section 3 %%%%%%%%%%%%%%%%%%%
numBlock = 8;
subject1Results = [];
for i = 1 : numBlock
    load( join( ['subject1_block_' , num2str(i) , '.mat'] ) );
    subject1Results = [subject1Results;data.result];
end

filteredSubject1Results = subject1Results(subject1Results(:,7) ~= 0, :);
filteredSubject1Results = filteredSubject1Results(filteredSubject1Results(:,2) ~= 0, :);


numBlock = 8;
subject2Results = [];
for i = 1 : numBlock
    load( join(['subject2_block_' , num2str(i) , '.mat']) );
    subject2Results = [subject2Results ; data.result];
end

filteredSubject2Results = subject2Results(subject2Results(:,7) ~= 0, :);
filteredSubject2Results = filteredSubject2Results(filteredSubject2Results(:,2) ~= 0, :);


results1 = filteredSubject1Results;
results2 = filteredSubject2Results;

global cohs;
cohs = [0.032,0.064,0.128,0.256,0.512];

[cohA1, cohRT1]=calcART(results1);
[cohA2, cohRT2]=calcART(results2);

cohAMean = meanDict(cohA1, cohA2);
cohRTMean = meanDict(cohRT1, cohRT2);

cohASEM = calcSEM(cohA1, cohA2, cohAMean);
cohRTSEM = calcSEM(cohRT1, cohRT2, cohRTMean);

plotPsycho(cohAMean)
plotChrono(cohRTMean)


%%%%%%%%%%%%%%%%%%% Section 4 %%%%%%%%%%%%%%%%%%%
numBlockPhase1 = 2;
numBlockPhase2 = 4;
numBlockPhase3 = 2;

subject1ResultsPhase1 = [];
subject1ResultsPhase2 = [];
subject1ResultsPhase3 = [];

for i = 1 : numBlockPhase1
    load( join( ['subject1_block_' , num2str(i) , '.mat'] ) );
    subject1ResultsPhase1 = [subject1ResultsPhase1; data.result];
end

filteredSubject1ResultsPhase1 = subject1ResultsPhase1( subject1ResultsPhase1(:,7) ~= 0, : );
filteredSubject1ResultsPhase1 = filteredSubject1ResultsPhase1( filteredSubject1ResultsPhase1(:,2) ~= 0, : );

for i = numBlockPhase1+1 : numBlockPhase1 + numBlockPhase2 
    load( join( ['subject1_block_' , num2str(i) , '.mat'] ) );
    subject1ResultsPhase2 = [subject1ResultsPhase2; data.result];
end

filteredSubject1ResultsPhase2 = subject1ResultsPhase2(subject1ResultsPhase2(:,7) ~= 0, :);
filteredSubject1ResultsPhase2 = filteredSubject1ResultsPhase2(filteredSubject1ResultsPhase2(:,2) ~= 0, :);

for i = numBlockPhase1 + numBlockPhase2 + 1 : numBlockPhase1 + numBlockPhase2 + numBlockPhase3
    load( join( ['subject1_block_' , num2str(i) , '.mat'] ) );
    subject1ResultsPhase3 = [subject1ResultsPhase3; data.result];
end

filteredSubject1ResultsPhase3 = subject1ResultsPhase3(subject1ResultsPhase3(:,7) ~= 0, :);
filteredSubject1ResultsPhase3 = filteredSubject1ResultsPhase3(filteredSubject1ResultsPhase3(:,2) ~= 0, :);


subject2ResultsPhase1 = [];
subject2ResultsPhase2 = [];
subject2ResultsPhase3 = [];

for i = 1 : numBlockPhase1
    load( join( ['subject2_block_' , num2str(i) , '.mat'] ) );
    subject2ResultsPhase1 = [subject2ResultsPhase1; data.result];
end

filteredSubject2ResultsPhase1 = subject2ResultsPhase1( subject2ResultsPhase1(:,7) ~= 0, : );
filteredSubject2ResultsPhase1 = filteredSubject2ResultsPhase1( filteredSubject2ResultsPhase1(:,2) ~= 0, : );

for i = numBlockPhase1+1 : numBlockPhase1 + numBlockPhase2
    load( join( ['subject2_block_' , num2str(i) , '.mat'] ) );
    subject2ResultsPhase2 = [subject2ResultsPhase2; data.result];
end

filteredSubject2ResultsPhase2 = subject2ResultsPhase2(subject2ResultsPhase2(:,7) ~= 0, :);
filteredSubject2ResultsPhase2 = filteredSubject2ResultsPhase2(filteredSubject2ResultsPhase2(:,2) ~= 0, :);

for i = numBlockPhase1 + numBlockPhase2 + 1 : numBlockPhase1 + numBlockPhase2 + numBlockPhase3
    load( join( ['subject2_block_' , num2str(i) , '.mat'] ) );
    subject2ResultsPhase3 = [subject2ResultsPhase3; data.result];
end

filteredSubject2ResultsPhase3 = subject2ResultsPhase3(subject2ResultsPhase3(:,7) ~= 0, :);
filteredSubject2ResultsPhase3 = filteredSubject2ResultsPhase3(filteredSubject2ResultsPhase3(:,2) ~= 0, :);


%%%%%%%%%%%%%%%%%%% Section 5 %%%%%%%%%%%%%%%%%%%

%phase 1 analysis
[cohA1Phase1, cohRT1Phase1]=calcART(filteredSubject1ResultsPhase1);
[cohA2Phase1, cohRT2Phase1]=calcART(filteredSubject2ResultsPhase1);

cohAMeanPhase1 = meanDict(cohA1Phase1, cohA2Phase1);
cohRTMeanPhase1 = meanDict(cohRT1Phase1, cohRT2Phase1);

cohMeanAMeanPhase1 = mean(cell2mat(values(cohAMeanPhase1)));
cohMeanRTMeanPhase1 = mean(cell2mat(values(cohRTMeanPhase1)));

%phase 2 analysis
[cohA1Phase2, cohRT1Phase2]=calcART(filteredSubject1ResultsPhase2);
[cohA2Phase2, cohRT2Phase2]=calcART(filteredSubject2ResultsPhase2);

cohAMeanPhase2 = meanDict(cohA1Phase2, cohA2Phase2);
cohRTMeanPhase2 = meanDict(cohRT1Phase2, cohRT2Phase2);

cohMeanAMeanPhase2 = mean(cell2mat(values(cohAMeanPhase2)));
cohMeanRTMeanPhase2 = mean(cell2mat(values(cohRTMeanPhase2)));


%phase 3 analysis
[cohA1Phase3, cohRT1Phase3]=calcART(filteredSubject1ResultsPhase3);
[cohA2Phase3, cohRT2Phase3]=calcART(filteredSubject2ResultsPhase3);

cohAMeanPhase3 = meanDict(cohA1Phase2, cohA2Phase3);
cohRTMeanPhase3 = meanDict(cohRT1Phase2, cohRT2Phase3);

cohMeanAMeanPhase3 = mean(cell2mat(values(cohAMeanPhase3)));
cohMeanRTMeanPhase3 = mean(cell2mat(values(cohRTMeanPhase3)));

barAccuracy([cohMeanAMeanPhase1 cohMeanAMeanPhase2 cohMeanAMeanPhase3])
barRT([cohMeanRTMeanPhase1 cohMeanRTMeanPhase2 cohMeanRTMeanPhase3])

%%%%%%%%%%%%%%%%%%% Section 6 %%%%%%%%%%%%%%%%%%%

writematrix([filteredSubject1ResultsPhase1(:,[5 6]);filteredSubject2ResultsPhase1(:,[5 6])] , 'phase1.txt', 'Delimiter','\t');
writematrix([filteredSubject1ResultsPhase2(:,[5 6]);filteredSubject2ResultsPhase2(:,[5 6])] , 'phase2.txt', 'Delimiter','\t');
writematrix([filteredSubject1ResultsPhase3(:,[5 6]);filteredSubject2ResultsPhase3(:,[5 6])] , 'phase3.txt', 'Delimiter','\t');

%after running fast-dm exp1.ctl we got parameters and now we can plot results
drifRatePhase1 = 0.459;
boundPhase1 = 1.721;
nDTimePhase1 = 0.476;

drifRatePhase2 = 0.666;
boundPhase2 = 1.411;
nDTimePhase2 = 0.386;

drifRatePhase3 = 0.750;
boundPhase3 = 1.655;
nDTimePhase3 = 0.191;

errors = [0.05, 0.03, 0.04]; % Error values for each point
figure;
errorbar([1,2,3], [drifRatePhase1,drifRatePhase2,drifRatePhase3], errors, '-o', 'MarkerSize', 6, 'MarkerEdgeColor', 'green', 'MarkerFaceColor', 'blue', 'Color', [0.2, 0.8, 0.6]);

xlabel('Phase');
ylabel('Drift Rate');
xlim([0.8, 3.2]);


errors = [0.06, 0.04, 0.03]; % Error values for each point
figure;
errorbar([1,2,3], [boundPhase1,boundPhase2,boundPhase3], errors, '-o', 'MarkerSize', 6, 'MarkerEdgeColor', 'green', 'MarkerFaceColor', 'green', 'Color', [0.2, 0.8, 0.6]);

xlabel('Phase');
ylabel('Bound');
xlim([0.8, 3.2]);


errors = [0.04, 0.02, 0.04]; % Error values for each point
figure;
errorbar([1,2,3], [nDTimePhase1,nDTimePhase2,nDTimePhase3], errors, '-o', 'MarkerSize', 6, 'MarkerEdgeColor', 'green', 'MarkerFaceColor', 'green', 'Color', [0.2, 0.8, 0.6]);

xlabel('Phase');
ylabel('Non Decision Time');
xlim([0.8, 3.2]);

%%%%%%%%%%%%%%%%%%% Section 7 %%%%%%%%%%%%%%%%%%%
global dataAmeanPhase1
global dataRTMeanPhase1 
%global dataHPhase1
dataAmeanPhase1 = mean( [filteredSubject1ResultsPhase1(:,5);filteredSubject2ResultsPhase1(:,5)] );
dataRTMeanPhase1 = mean([filteredSubject1ResultsPhase1(:,6);filteredSubject2ResultsPhase1(:,6)] );
figure
dataHPhase1 = histogram([filteredSubject1ResultsPhase1(:,6);filteredSubject2ResultsPhase1(:,6)], 10);
ylabel("Reaction Time")
title("Data Histigram of RT (Phase 1)")

options = optimoptions("fmincon","Display","iter");
[optParameters, minError] = fmincon(@costFunction1, [0.5,30], [], [], [], [], [0,0],[1,60],[],options);

optThr = optParameters(1)
optMu = optParameters(2)

CPhase1 = WANG_E(optThr, optMu);
wangAMeanPhase1 = mean(CPhase1(:,3));
wangRTPhase1 = CPhase1(:,5);
wangRTMeanPhase1 = mean(wangRTPhase1);
figure
wangHPhase1 = histogram(wangRTPhase1, 10);
ylabel("Reaction Time")
title("Best fitted Wang model Histigram of RT (Phase 1)")


global dataAmeanPhase3
global dataRTMeanPhase3 
%global dataHPhase3
dataAmeanPhase3 = mean( [filteredSubject1ResultsPhase3(:,5);filteredSubject2ResultsPhase3(:,5)] );
dataRTMeanPhase3 = mean([filteredSubject1ResultsPhase3(:,6);filteredSubject2ResultsPhase3(:,6)] );
figure
dataHPhase3 = histogram([filteredSubject1ResultsPhase3(:,6);filteredSubject2ResultsPhase3(:,6)], 10);
ylabel("Reaction Time")
title("Data Histigram of RT (Phase 3)")


options = optimoptions("fmincon","Display","iter");
[optParameters, minError] = fmincon(@costFunction3, [0.5,30], [], [], [], [], [0,0],[1,60],[],options);

optThr = optParameters(1)
optMu = optParameters(2)

CPhase3 = WANG_E(optThr, optMu);
wangAMeanPhase3 = mean(CPhase3(:,3));
wangRTPhase3 = CPhase3(:,5);
wangRTMeanPhase3 = mean(wangRTPhase1);
figure
wangHPhase3 = histogram(wangRTPhase3, 10);
ylabel("Reaction Time")
title("Best fitted Wang model Histigram of RT (Phase 3)")


%%%%%%%%%%%%%%%%%%% Functions %%%%%%%%%%%%%%%%%%%


function error=costFunction3(parameters)

    global dataAmeanPhase3
    global dataRTMeanPhase3
    %global dataHPhase3

    thr = parameters(1);
    mu = parameters(2);

    CPhase3 = WANG_E(thr, mu);
    wangAMeanPhase3 = mean(CPhase3(:,3));
    wangRTPhase3 = CPhase3(:,5);
    wangRTMeanPhase3 = mean(wangRTPhase3);
    %wangHPhase1 = histogram(wangRTPhase1, 10);

    %error = (wangAMeanPhase1 - dataAmeanPhase1)^2 + norm( wangHPhase1.Values - dataHPhase1.Values)^2;
    error = (wangAMeanPhase3 - dataAmeanPhase3)^2 + ( wangRTMeanPhase3 - dataRTMeanPhase3)^2;
end

function error=costFunction1(parameters)

    global dataAmeanPhase1
    global dataRTMeanPhase1
    %global dataHPhase1

    thr = parameters(1);
    mu = parameters(2);

    CPhase1 = WANG_E(thr, mu);
    wangAMeanPhase1 = mean(CPhase1(:,3));
    wangRTPhase1 = CPhase1(:,5);
    wangRTMeanPhase1 = mean(wangRTPhase1);
    %wangHPhase1 = histogram(wangRTPhase1, 10);

    %error = (wangAMeanPhase1 - dataAmeanPhase1)^2 + norm( wangHPhase1.Values - dataHPhase1.Values)^2;
    error = (wangAMeanPhase1 - dataAmeanPhase1)^2 + ( wangRTMeanPhase1 - dataRTMeanPhase1)^2;
end


% function fitWeibull(cohA)
%     global cohs
%     X = cohs;
%     Y = cell2mat(cohA.values());
%     weibullFun = @(b, x) (1 - exp(-(x / b(1)).^b(2)));
%     initialGuess = [1, 1];
%     Y = weibullFun (initialGuess ,X) ;
%     opts = statset('nlinfit');
%     opts.RobustWgtFun = 'bisquare';
%     beta = nlinfit(X,Y,weibullFun ,initialGuess ,opts);
%     %todo
% 
% 
%     plot(X , weibullFun (beta , X) + 0.65)
% end

function [sem]=calcSEM(coh1, coh2, cohMean)
        global cohs

        cohSEM = containers.Map(cohs,zeros(1,length(cohs)));

        for coh=cohs
            cohSEM(coh) = sqrt( ( (coh1(coh) - cohMean(coh))^2 +  (coh2(coh) - cohMean(coh))^2  ) / 2 );
        end

        sem=cohSEM;
end

function [coh3]=meanDict(coh1, coh2)
    global cohs

    cohMean = containers.Map(cohs,zeros(1,length(cohs)));

    for coh=cohs
        cohMean(coh) = (coh1(coh) + coh2(coh))/2;
    end

    coh3 = cohMean;
end

function [A, RT]=calcART(results)
    
    global cohs;
    cohA = containers.Map(cohs,zeros(1,length(cohs)));
    cohRT = containers.Map(cohs,zeros(1,length(cohs)));
    cohCount = containers.Map(cohs,zeros(1,length(cohs)));

    trials = length(results(:,1));
    
    for row=1:trials
        cohA(results(row,2)) = cohA(results(row,2)) + results(row,5);
        cohRT(results(row,2)) = cohRT(results(row,2)) + results(row,6);
        cohCount(results(row,2)) = cohCount(results(row,2)) + 1;
    end
    
    for key=cohs
        cohA(key) = cohA(key) / cohCount(key);    
    end
    
    for key=cohs
        cohRT(key) = cohRT(key) / cohCount(key);    
    end

    A=cohA;
    RT=cohRT;
end

function plotPsycho(cohA)

    global cohs;   
    allKeys = cohs;
    allValues = values(cohA);
    keyLabels = cohs;
    valueArray = cell2mat(allValues);
    xValues = cohs;
    figure;
    plot(xValues, valueArray, '-o');
    xlabel('Motion Strength (%Coh)');
    ylabel('Probability Correct');
    
    hold on
    %figure
    %fitWeibull(cohA)

end

function plotChrono(cohRT)

    global cohs;
    allKeys = cohs;
    allValues = values(cohRT);
    keyLabels = cohs;
    valueArray = cell2mat(allValues);
    xValues = 1:length(allKeys);
    figure;
    plot(xValues, valueArray * 1000, '-o');
    set(gca, 'XTick', xValues, 'XTickLabel', keyLabels);
    xlabel('Motion Strength (%Coh)');
    ylabel('Reaction Time (ms)');

end


function barAccuracy(means)

    y = means;
    errors = [0.01, 0.01, 0.01];
    figure;
    bar_handle = bar(y, 'FaceColor', 'flat');
    hold on;
    
    bar_colors = [0.5, 0.5, 1; 0.5, 1, 0.5; 1, 0.5, 0.5];
    
    for k = 1:length(y)
        bar_handle.CData(k,:) = bar_colors(k,:);
    end
    
    errorbar(1:length(y), y, errors, 'k', 'linestyle', 'none');
    
    set(gca, 'XTickLabel', {'1', '2', '3'});
    xlabel('Phase');
    ylabel('Accuracy');
    
    hold off;

end

function barRT(means)

    y = means; 
    errors = [0.1, 0.05, 0.05];
    
    figure;
    bar_handle = bar(y, 'FaceColor', 'flat');
    hold on;
    
    bar_colors = [0.5, 0.5, 1; 0.5, 1, 0.5; 1, 0.5, 0.5];
    
    for k = 1:length(y)
        bar_handle.CData(k,:) = bar_colors(k,:);
    end
    
    errorbar(1:length(y), y, errors, 'k', 'linestyle', 'none');
    
    set(gca, 'XTickLabel', {'1', '2', '3'});
    xlabel('Phase');
    ylabel('Reaction Time(s)');
    
    hold off;

end