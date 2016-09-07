classdef gapnmfclass
    properties
        % data
        x;
        
        % hyperparams
        alpha;
        a;
        b;

        % dimensions
        K;
        M;
        N;
        
        % variational parameters
        rhow;
        tauw;
        rhoh;
        tauh;
        rhot;
        taut;

        % parameter expectations
        Ew;
        Ewinv;
        Ewinvinv;
        Eh;
        Ehinv;
        Ehinvinv;
        Et;
        Etinv;
        Etinvinv;
    end

    methods
        % runs one iteration of the optimization
        function obj = update(obj, varargin)
            verbose = 0;
            verbose2 = 0;
            score = 0;
            lastscore = 0;    
            if (~isempty(varargin) && strcmp(varargin{1}, 'verbose'))
                verbose = 1;
            end
            if (~isempty(varargin) && strcmp(varargin{1}, 'verbose2'))
                verbose = 1;
                verbose2 = 'verbose';
            end

            if verbose, lastscore = score; score = obj.bound(verbose2); fprintf('update:  %f  initial\n', score-lastscore); end
            obj = obj.updateh();
            if verbose, lastscore = score; score = obj.bound(verbose2); fprintf('update:  %f  H\n', score-lastscore); end
            if (score < lastscore), error('oops.'); end;
            obj = obj.updatew();
            if verbose, lastscore = score; score = obj.bound(verbose2); fprintf('update:  %f  W\n', score-lastscore); end
            if (score < lastscore), error('oops.'); end;
            obj = obj.updatetheta();
            if verbose, lastscore = score; score = obj.bound(verbose2); fprintf('update:  %f  Theta\n', score-lastscore); end
            if (score < lastscore), error('oops.'); end;

            obj = obj.clearbadk();
%             if verbose, lastscore = score; score = obj.bound(verbose2); fprintf('update:  %f  clearbadk\n', score-lastscore); end
%             if (score < lastscore), error('oops.'); end;
        end

        % These functions (updatewslow, updatehslow, and updatethetaslow)
        % are equivalent to their counterparts updatew, updateh, and 
        % updatetheta, but are slower (in matlab) and use much more memory.
        % They represent the bound parameters phi and omega explicitly, and
        % may therefore be easier to follow.
        function obj = updatewslow(obj)
            goodk = obj.goodk();
            phi = zeros(obj.M, obj.N, obj.K);
            for k = goodk,
                phi(:, :, k) = obj.Ewinvinv(:, k) * obj.Etinvinv(k) * obj.Ehinvinv(k, :);
            end
            phinorm = 1./sum(phi, 3);
            for k = goodk,
                phi(:, :, k) = phi(:, :, k) .* phinorm;
            end
            omega = obj.Ew(:, goodk) * diag(obj.Et(goodk)) * obj.Eh(goodk, :);

            obj.rhow(:, goodk) = obj.a + (1./omega) * obj.Eh(goodk, :)' * ...
                diag(obj.Et(goodk));
            for k = goodk,
                obj.tauw(:, k) = phi(:, :, k).^2 .* obj.x * ...
                    obj.Ehinv(k, :)' * obj.Etinv(k);
            end
            obj.tauw(obj.tauw < 1e-100) = 0;
            [obj.Ew(:, goodk) obj.Ewinv(:, goodk)] = computegigexpectations(obj.a, obj.rhow(:, goodk), obj.tauw(:, goodk));
            obj.Ewinvinv(:, goodk) = obj.Ewinv(:, goodk).^-1;
        end

        function obj = updatehslow(obj)
            goodk = obj.goodk();
            phi = zeros(obj.M, obj.N, obj.K);
            for k = goodk,
                phi(:, :, k) = obj.Ewinvinv(:, k) * obj.Etinvinv(k) * obj.Ehinvinv(k, :);
            end
            phinorm = 1./sum(phi, 3);
            for k = goodk,
                phi(:, :, k) = phi(:, :, k) .* phinorm;
            end
            omega = obj.Ew(:, goodk) * diag(obj.Et(goodk)) * obj.Eh(goodk, :);

            obj.rhoh(goodk, :) = obj.b + diag(obj.Et(goodk)) * ...
                obj.Ew(:, goodk)' * (1./omega);
            for k = goodk,
                obj.tauh(k, :) = obj.Ewinv(:, k)' * obj.Etinv(k) * ...
                    (phi(:, :, k).^2 .* obj.x);
            end
            obj.tauh(obj.tauh < 1e-100) = 0;
            [obj.Eh(goodk, :) obj.Ehinv(goodk, :)] = computegigexpectations(obj.b, obj.rhoh(goodk, :), obj.tauh(goodk, :));
            obj.Ehinvinv(goodk, :) = obj.Ehinv(goodk, :).^-1;
        end
        
        function obj = updatethetaslow(obj)
            goodk = obj.goodk();
            phi = zeros(obj.M, obj.N, obj.K);
            for k = goodk,
                phi(:, :, k) = obj.Ewinvinv(:, k) * obj.Etinvinv(k) * obj.Ehinvinv(k, :);
            end
            phinorm = 1./sum(phi, 3);
            for k = goodk,
                phi(:, :, k) = phi(:, :, k) .* phinorm;
            end
            omega = obj.Ew(:, goodk) * diag(obj.Et(goodk)) * obj.Eh(goodk, :);

            obj.rhot(goodk) = obj.alpha + sum(obj.Ew(:, goodk)' * ...
                (1./omega) .* obj.Eh(goodk, :), 2);
            for k = goodk,
                obj.taut(k) = sum(obj.Ewinv(:, k)' * ...
                    (phi(:, :, k).^2 .* obj.x) .* obj.Ehinv(k, :), 2);
            end
            obj.taut(obj.taut < 1e-100) = 0;
            [obj.Et(goodk) obj.Etinv(goodk)] = computegigexpectations(obj.alpha/obj.K, obj.rhot(goodk), obj.taut(goodk));
            obj.Etinvinv(goodk) = obj.Etinv(goodk).^-1;
        end

        % These versions of updatew, updateh, and updatetheta are faster
        % than (at least in matlab) and equivalent to their counterparts
        % updatewslow, updatehslow, and updatethetaslow, and use much less
        % memory. The reason is that they don't explicitly compute the
        % bound parameters phi and omega.
        function obj = updatew(obj)
            goodk = obj.goodk();
            xxtwidinvsq = obj.x .* obj.xtwid(goodk).^-2;
            xbarinv = obj.xbar(goodk).^-1;
            dEt = diag(obj.Et(goodk));
            dEtinvinv = diag(obj.Etinvinv(goodk));
            obj.rhow(:, goodk) = obj.a + xbarinv * obj.Eh(goodk, :)' * dEt;
            obj.tauw(:, goodk) = obj.Ewinvinv(:, goodk).^2 .* ...
                (xxtwidinvsq * obj.Ehinvinv(goodk, :)' * dEtinvinv);
            obj.tauw(obj.tauw < 1e-100) = 0;
            [obj.Ew(:, goodk) obj.Ewinv(:, goodk)] = computegigexpectations(obj.a, obj.rhow(:, goodk), obj.tauw(:, goodk));
            obj.Ewinvinv(:, goodk) = obj.Ewinv(:, goodk).^-1;
        end

        function obj = updateh(obj)
            goodk = obj.goodk();
            xxtwidinvsq = obj.x .* obj.xtwid(goodk).^-2;
            xbarinv = obj.xbar(goodk).^-1;
            dEt = diag(obj.Et(goodk));
            dEtinvinv = diag(obj.Etinvinv(goodk));
            obj.rhoh(goodk, :) = obj.b + dEt * (obj.Ew(:, goodk)' * xbarinv);
            obj.tauh(goodk, :) = obj.Ehinvinv(goodk, :).^2 .* ...
                (dEtinvinv * (obj.Ewinvinv(:, goodk)' * xxtwidinvsq));
            obj.tauh(obj.tauh < 1e-100) = 0;
            [obj.Eh(goodk, :) obj.Ehinv(goodk, :)] = computegigexpectations(obj.b, obj.rhoh(goodk, :), obj.tauh(goodk, :));
            obj.Ehinvinv(goodk, :) = obj.Ehinv(goodk, :).^-1;
        end
        
        function obj = updatetheta(obj)
            goodk = obj.goodk();
            xxtwidinvsq = obj.x .* obj.xtwid(goodk).^-2;
            xbarinv = obj.xbar(goodk).^-1;
            obj.rhot(goodk) = obj.alpha + sum(obj.Ew(:, goodk)' * xbarinv .* obj.Eh(goodk, :), 2);
            obj.taut(goodk) = obj.Etinvinv(goodk).^2 .* ...
                sum(obj.Ewinvinv(:, goodk)' * xxtwidinvsq .* obj.Ehinvinv(goodk, :), 2);
            obj.taut(obj.taut < 1e-100) = 0;
            [obj.Et(goodk) obj.Etinv(goodk)] = computegigexpectations(obj.alpha/obj.K, obj.rhot(goodk), obj.taut(goodk));
            obj.Etinvinv(goodk) = obj.Etinv(goodk).^-1;
        end
        
        function goodk = goodk(obj, cutoff)
            if (nargin < 2),
                cutoff = 1e-10 * max(obj.x(:));
            end
            powers = obj.Et .* max(obj.Ew, [], 1)' .* max(obj.Eh, [], 2);
            [temp sorted] = sort(powers, 'descend');
            goodk = sorted(1:find(temp/max(temp) > cutoff, 1, 'last'))';
            if (powers(goodk(end)) < cutoff),
                goodk(end) = [];
            end
        end
        
        % This function sets unused components' posteriors equal to their
        % priors.
        function obj = clearbadk(obj)
            goodk = obj.goodk();
            badk = setdiff(1:obj.K, goodk);
            obj.rhow(:, badk) = obj.a;
            obj.tauw(:, badk) = 0;
            obj.rhoh(badk, :) = obj.b;
            obj.tauh(badk, :) = 0;
            obj = obj.recomputeexpectations();
        end
        
        function obj = recomputeexpectations(obj)
            [obj.Ew obj.Ewinv] = computegigexpectations(obj.a, obj.rhow, obj.tauw);
            obj.Ewinvinv = obj.Ewinv.^-1;
            [obj.Eh obj.Ehinv] = computegigexpectations(obj.b, obj.rhoh, obj.tauh);
            obj.Ehinvinv = obj.Ehinv.^-1;
            [obj.Et obj.Etinv] = computegigexpectations(obj.alpha/obj.K, obj.rhot, obj.taut);
            obj.Etinvinv = obj.Etinv.^-1;
        end
        
        function figures(obj)
            subplot(3, 2, 1);
            imagesc(log(obj.Ew));
            title('E[W]')
            xlabel('component index')
            ylabel('frequency')
            subplot(3, 2, 2);
            imagesc(log(obj.Eh));
            title('E[H]')
            xlabel('time')
            ylabel('component index')
            subplot(3, 2, 3);
            bar(1:obj.K, obj.Et);
%             bar(1:obj.K, sort(obj.Et, 'descend'));
            title('E[\theta]')
            xlabel('component index')
            ylabel('E[\theta]')
            subplot(3, 2, 5);
            imagesc(log(obj.x));
            title('Original Spectrogram')
            xlabel('time')
            ylabel('frequency')
            subplot(3, 2, 6);
            imagesc(log(obj.xbar()));
            title('Reconstructed Spectrogram')
            xlabel('time')
            ylabel('frequency')
            pause(0.000001);
        end
        
        function score = bound(obj, varargin)
            verbose = 0;
            if (~isempty(varargin) && strcmp(varargin{1}, 'verbose'))
                verbose = 1;
                lastscore = 0;
            end
            
            score = 0;
            
            goodk = obj.goodk();
            
            xbar = obj.xbar(goodk);
            xtwid = obj.xtwid(goodk);
            score = score - sum(obj.x(:) ./ xtwid(:) + log(xbar(:)));
            if (verbose)
                fprintf('%f:  X\n', score - lastscore);
                lastscore = score;
            end
            
            score = score + giggammaterm(obj.Ew, obj.Ewinv, obj.rhow, obj.tauw, obj.a, obj.a);
            if (verbose)
                fprintf('%f:  W\n', score - lastscore);
                lastscore = score;
            end
            score = score + giggammaterm(obj.Eh, obj.Ehinv, obj.rhoh, obj.tauh, obj.b, obj.b);
            if (verbose)
                fprintf('%f:  H\n', score - lastscore);
                lastscore = score;
            end
            score = score + giggammaterm(obj.Et, obj.Etinv, obj.rhot, obj.taut, obj.alpha/obj.K, obj.alpha);
            if (verbose)
                fprintf('%f:  Theta\n', score - lastscore);
                lastscore = score;
            end
        end
        
        function xbar = xbar(obj, goodk)
            if (nargin < 2)
                goodk = 1:obj.K;
            end
            xbar = obj.Ew(:, goodk) * diag(obj.Et(goodk)) * obj.Eh(goodk, :);
        end

        function xtwid = xtwid(obj, goodk)
            if (nargin < 2)
                goodk = 1:obj.K;
            end
            xtwid = obj.Ewinvinv(:, goodk) * diag(obj.Etinvinv(goodk)) * obj.Ehinvinv(goodk, :);
        end

        function obj = gapnmfclass(x, alpha, a, b, K, smoothness)
            if (nargin < 6)
                smoothness = 100;
            end
            
            obj.x = x / mean(x(:));
            obj.alpha = alpha;
            obj.a = a;
            obj.b = b;
            obj.K = K;
            M = size(x, 1);
            N = size(x, 2);
            obj.M = size(x, 1);
            obj.N = size(x, 2);

            obj.rhow = 10000*gamrnd(smoothness, 1/smoothness, M, K);
            obj.tauw = 10000*gamrnd(smoothness, 1/smoothness, M, K);
            obj.rhoh = 10000*gamrnd(smoothness, 1/smoothness, K, N);
            obj.tauh = 10000*gamrnd(smoothness, 1/smoothness, K, N);
            obj.rhot = K*10000*gamrnd(smoothness, 1/smoothness, K, 1);
            obj.taut = 1/K*10000*gamrnd(smoothness, 1/smoothness, K, 1);
            
            [obj.Ew obj.Ewinv] = computegigexpectations(a, obj.rhow, obj.tauw);
            obj.Ewinvinv = obj.Ewinv.^-1;
            [obj.Eh obj.Ehinv] = computegigexpectations(b, obj.rhoh, obj.tauh);
            obj.Ehinvinv = obj.Ehinv.^-1;
            [obj.Et obj.Etinv] = computegigexpectations(alpha/K, obj.rhot, obj.taut);
            obj.Etinvinv = obj.Etinv.^-1;
        end
    end
end
