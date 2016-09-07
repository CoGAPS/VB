% load snippet of Handel's Messiah included with Matlab.
% (replace temp.y with your favorite 1-d audio signal.)
temp = load('handel');
% compute its power spectrogram.
handelsgram = abs(specgram(temp.y)).^2;

% create gapnmfclass object with parameters alpha=1, a=0.1, b=0.1, K=50.
% feel free to play around with these parameters, particularly K.
obj = gapnmfclass(handelsgram, 1, 0.1, 0.1, 50);

% iterate updates until convergence.
score = -inf;
criterion = 0.0005;
for iteration = 1:1000,
    % does one update of the variational parameters.
    obj = obj.update();

    % displays the current expectations of the model parameters and
    % the reconstruction obtained from them.
    obj.figures();

    % compute the current variational bound, and see if we should quit.
    lastscore = score;
    score = obj.bound();
    improvement = (score - lastscore) / abs(lastscore);
    fprintf('iteration %d:  bound=%f  (%f improvement)\n', iteration, ...
        score, improvement);
    if (improvement < criterion)
        break;
    end
end