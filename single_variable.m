data = readtable('sample_file.txt');
for i = 1:size(data, 2)
     if isnumeric(data{:, i})
         data{:, i}(isnan(data{:, i})) = Inf;
     end
 end

sc = creditscorecard(data,'GoodLabel',0,'ResponseVar','overdue');
sc = autobinning(sc, 'Algorithm', 'Monotone', 'AlgorithmOptions', {'InitialNumBins',10});



data = readtable('sample_file.txt');

data.lending_time=datetime(data.lending_time);
data = sortrows(data,'lending_time','ascend');
data.lending_time = [];
data = data(randperm(size(data, 1)), :);
cut = round(0.7*size(data, 1));

%取训练样本、测试样本
train_data = data(1:cut, :);
test_data = data((cut + 1):size(data, 1), :);

%sc = creditscorecard(train_data,'GoodLabel',0,'ResponseVar','overdue');
sc = creditscorecard(data,'GoodLabel',0,'ResponseVar','overdue');

% test_sc = creditscorecard(test_data,'GoodLabel',0,'ResponseVar','overdue');

sc = autobinning(sc, 'Algorithm', 'Monotone', 'AlgorithmOptions', {'InitialNumBins',10});

%[bi, cp] = bininfo(sc, 'close_mobile_3days_times');
%cp(5) = [];
%sc = modifybins(sc,'close_mobile_3days_times','CutPoints',cp);

iv_table = table(zeros(size(sc.PredictorVars, 2), 1), cell(size(sc.PredictorVars, 2), 1), 'VariableNames',{'iv' 'PredictorVars'});
for i = 1:size(sc.PredictorVars, 2)
    iv = ivcalc(sc,i);
    iv_table.iv(i) = iv;
    iv_table.PredictorVars{i} = sc.PredictorVars{i};
end
iv_table = sortrows(iv_table,'iv', 'descend');
InModelVars = {};
for i = 1:size(iv_table, 1)
    if iv_table.iv(i) >= 0.02 && ~strcmp(iv_table.PredictorVars{i}, 'province')
        InModelVars = [InModelVars  iv_table.PredictorVars{i}];
    end
end

train_data_for_psi = train_data;
train_data_for_psi.overdue(1:round(size(train_data_for_psi, 1) / 2)) = 0;
train_data_for_psi.overdue(round(size(train_data_for_psi, 1)/2+1):end) = 1;

psi_sc = creditscorecard(train_data_for_psi,'GoodLabel',1,'ResponseVar','overdue', 'PredictorVars', InModelVars);
%psi_sc = creditscorecard(train_data_for_psi,'GoodLabel',1,'ResponseVar','overdue');
psi_sc = autobinning(psi_sc, 'Algorithm', 'EqualFrequency', 'AlgorithmOptions', {'NumBins',10});
psi_table = table(zeros(size(psi_sc.PredictorVars, 2), 1), cell(size(psi_sc.PredictorVars, 2), 1), 'VariableNames',{'psi' 'PredictorVars'});
for i = 1:size(psi_sc.PredictorVars, 2)
    bi = bininfo(psi_sc, psi_sc.PredictorVars{i});
    psi = sum(bi.InfoValue(bi.InfoValue(1:end-1) < Inf)); %bi.InfoValue(end);
    psi_table.psi(i) = psi;
    psi_table.PredictorVars{i} = psi_sc.PredictorVars{i};
end
psi_table = sortrows(psi_table,'psi', 'ascend');
InModelVars = {};
for i = 1:size(psi_table, 1)
    if psi_table.psi(i) < 0.1 
        InModelVars = [InModelVars  psi_table.PredictorVars{i}];
    end
end

%模型拟合
model_sc = fitmodel(sc, 'PredictorVars', InModelVars);

%sc = formatpoints(sc,'Missing','MinPoints');
TargetPoints = 500;
TargetOdds = 2;
PDO = 50; % Points to double the odds
model_sc = formatpoints(model_sc,'PointsOddsAndPDO',[TargetPoints TargetOdds PDO]);
[Stats,T]=validatemodel(model_sc, test_data,'Plot',{'ROC','KS'});

p2 = displaypoints(sc);
disp(p2);

[Scores,Points] = score(sc, data);

fid=fopen('validation_sample.txt', 'wt'); %打开文件
for i = 1:size(Scores, 1)
    fprintf(fid, '%d,%f\n', data.overdue(i), Scores(i));
end

fclose(fid);



