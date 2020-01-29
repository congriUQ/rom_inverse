classdef ROM_SPDE < handle
    %Class for ROM SPDE model
    
    properties  %public
        %% Finescale data specifications
        %Finescale system size
        nElFX = 256
        nElFY = 256
        %Finescale conductivities, binary material
        lowerConductivity = 1
        upperConductivity = 10
        %Conductivity field distribution type
        nBochnerBasis = 3  % Set to false if no fixed basis is desired
        conductivityDistributionType = "continuous"     %binary (default) or continuous
        conductivityDistribution = "squaredExponential"
        randFieldGenMatrix      %matrix for quick generation of discretized conductivity fields
        %Boundary condition functions
        %evaluate those on boundaries to get boundary conditions
        boundaryTemperature
        boundaryHeatFlux
        bcArrayTrain   %boundary condition cell array for training data
        bcArrayTest    %boundary condition cell array for test data
        
        naturalNodes
        %Directory where finescale data is stored; specify basename here
        fineScaleDataPath
        %matfile handle
        trainingDataMatfile
        testDataMatfile
        %Finescale Mesh object
        fineMesh
        %Array holding fine scale data output; possibly large
        fineScaleDataOutput
        %number of samples per generated matfile
        nSets = [4096 4096]
        %Output data characteristics
        outputVariance
        outputMean
        meanOutputVariance
        coarseCellMask          %Mask to get fine scale pixels of coarse cell
        neighborDictionary %Gives neighbors of macrocells
        %% Model training parameters
        nStart = 1              %first training data sample in file
        nTrain = 16             %number of samples used for training
        mode = 'none'           %useNeighbor, useLocalNeighbor, useDiagNeighbor,
                                %useLocalDiagNeighbor, useLocal, global
                                %global: take whole microstructure as feature 
                                %function input, not only local window (only 
                                %recommended for pooling)
        %E-step inference method. variationalInference or monteCarlo
        inferenceMethod = 'variationalInference'
        
        %Principal components of the whole microstructure used as features 
        globalPcaComponents = 3
        %Principal components of single macro-cell used as features
        localPcaComponents = 7
        pcaSamples = 4096
        mix_S = 0             %To slow down convergence of S
        mix_theta = 0
        mix_sigma = 0
        fixSigmaInit = 0       %Fix sigma in initial iterations
        
        %% Prior specifications
        sigmaPriorType = 'none'    %none, delta, expSigSq
        sigmaPriorHyperparam = 1
        thetaPriorType = 'RVM'
        thetaPriorHyperparam = [0 1e-30]
        
        %% Model parameters
        theta_c
        theta_cf
        free_W = false
        %Cell array containing local feature function handles
        featureFunctions
        %cell array with handles to global feature functions
        globalFeatureFunctions
        %transformation of finescale conductivity to real axis
        conductivityTransformation
        latentDim = 0              %If autoencoder is used
        sumPhiTPhi             %Design matrix precomputation
        
        %% Feature function rescaling parameters
        %'standardize' for zero mean and unit variance of features, 'rescale' 
        %to have all between 0 and 1, 'normalize' to have unit variance only
        %(and same mean as before)
%         featureScaling = 'normalize'
        featureFunctionMean
        featureFunctionSqMean
        featureFunctionMin
        featureFunctionMax
        loadDesignMatrix = false

        
        %% Prediction parameters
        nSamples_p_c = 1000
        %Use Laplace approx around MAP for prediction?
        useLaplaceApproximation = false
        testSamples = 17:22       %pick out specific test samples here
        trainingSamples   %pick out specific training samples here
        
        %% Prediction outputs
        predMeanArray
        predVarArray
        meanPredMeanOutput                %mean predicted mean of output field
        %mean squared distance of predicted mean to true solution
        squaredDistance
        normError
        meanSquaredDistance
        meanSquaredDistanceField
        meanSquaredDistanceError          %Monte Carlo error
        meanLogLikelihood
        meanLogPerplexity
        meanPerplexity
        meanMahalanobisError
        meanEffCond
        
        %% Finescale data
        lambdak
        xk
        conductivityBasisFunctions      %handle to conductivity basis functions
        conductivityCoefficients
        
        %% Computational quantities
        varExpect_p_cf_exp_mean
        XMean                  %Expected value of X under q
        XSqMean                %<X^2>_q
        mean_ufucT
        mean_ucucT
        thetaArray
        thetaHyperparamArray
        sigmaArray
        
        EM_iterations = 1
        epoch = 0      %how often every data point has been seen
        epoch_old = 0  %To check if epoch has changed
        maxEpochs      %Maximum number of epochs
    end
    
    
    properties(SetAccess = private)
        %% finescale data specifications
        %delta for fixed length scale, lognormal for rand
        conductivityLengthScaleDist = "delta"
        
        %{volumeFraction, correlationLength, sigma_f2}
        %for log normal length scale, the
        %length scale parameters are log normal mu and
        %sigma
        conductivityDistributionParams = {-1 [.2 .2] 9}
        
        %Coefficients giving boundary conditions, specify as string
        boundaryConditions = '[0 1 1 0]'
        boundaryConditionVariance = [0 0 0 0]
        
        %% Coarse model specifications
        coarseMesh
        coarseGridVectorX = (1/4)*ones(1, 4)
        coarseGridVectorY = (1/4)*ones(1, 4)
        
        %Design matrices. Cell index gives data point, row index coarse cell, 
        %and column index feature function
        designMatrix
        originalDesignMatrix    %design matrix without any locality mode
        testDesignMatrix    %design matrices on independent test set
    end
    
    
    methods
        function self = ROM_SPDE(mode)
            %Constructor
            addpath('./FEMgradient')
            %Create data directory
            if ~exist('./data/', 'dir')
                mkdir('./data/');
            end
            
            %add needed paths
            addpath('./aux');
            %set up path
            self.genFineScaleDataPath;
            %Set handles to boundary condition functions
            self.genBoundaryConditionFunctions;
            self.naturalNodes = [2:(2*self.nElFX + 2*self.nElFY)];
            %Set up coarseMesh; must be done after b.c.'s are set up
            self.genCoarseMesh;
            %prealloc
            self.XMean = zeros(self.coarseMesh.nEl, self.nTrain);
            self.XSqMean = ones(self.coarseMesh.nEl, self.nTrain);
            %Set up default value for test samples
            self.nStart = randi(self.nSets(1) - self.nTrain, 1);
            self.trainingSamples = self.nStart:(self.nStart + self.nTrain - 1);
            %Set conductivity transformation
            self.conductivityTransformation.anisotropy = false;
            self.conductivityTransformation.type = 'log';
            if strcmp(self.conductivityTransformation.type, 'log')
                self.conductivityTransformation.limits = [1e-6 1e6];
            elseif strcmp(self.conductivityTransformation.type, 'logit')
                self.conductivityTransformation.limits = [(1 - 1e-2)*self.lowerConductivity,...
                                                            (1 + 1e-2)*self.upperConductivity];
            else
                self.conductivityTransformation.limits = [1e-8 1e8];
            end
            conductivityTransformation = self.conductivityTransformation;
            save('./data/conductivityTransformation', 'conductivityTransformation');
            
            if ~strcmp(mode, 'genData')
                %Load fine scale domain and set boundary conditions
                self.loadTrainingData;
                %Set up feature function handles
                %Prealloc
                self.varExpect_p_cf_exp_mean = zeros(self.fineMesh.nNodes, 1);
            end
            
            %Check
            if(strcmp(self.thetaPriorType, 'sharedRVM') && ~strcmp(self.mode, 'useLocal'))
                error('sharedRVM prior only allowed for useLocal mode')
            end
        end
        
        function genFineScaleDomain(self)
            %Generate finescale domain
            tic
            disp('Generate fine scale domain...')
            addpath('./heatFEM')    %to find Domain class
            self.fineMesh = Domain(self.nElFX, self.nElFY);
            %Only fix lower left corner as essential node
            self.fineMesh = setBoundaries(self.fineMesh, self.naturalNodes, self.boundaryTemperature,...
                                          self.boundaryHeatFlux);
            disp('done')
            domain_generating_time = toc
        end
        
        function genFineScaleData(self, boundaryConditions, condDistParams)
            %Function to generate and save finescale data
            disp('Generate fine scale data...')
            
            if(nargin > 1)
                self.setBoundaryConditions(boundaryConditions);
            end
            
            if(nargin > 2)
                self = self.setConductivityDistributionParams(condDistParams);
            end
            
            %for boundary condition functions
            if(isempty(self.boundaryTemperature) || isempty(self.boundaryHeatFlux))
                self.genBoundaryConditionFunctions;
            end
            
            %Generate fine scale domain
            self.genFineScaleDomain();
            
            if ~exist(self.fineScaleDataPath, 'dir')
                mkdir(self.fineScaleDataPath);
            end
            
            %Generate finescale conductivity samples and solve FEM
            for n = 1:numel(self.nSets)
                filename = sprintf(self.fineScaleDataPath + "set%u-samples=%u.mat", n, self.nSets(n));
                cond = self.generateConductivityField(n, filename);
                self.solveFEM(n, cond, filename);
            end
            
            %save params
            fineMesh = self.fineMesh;
            fineMesh = fineMesh.shrink();
            save(self.fineScaleDataPath + "fineMesh.mat", "fineMesh");
            disp('done')
        end
        
        function genBoundaryConditionFunctions(self)
            %Set up boundary condition functions
            if isempty(self.boundaryConditions)
                error('No string specified for boundary conditions')
            end
            bc = str2num(self.boundaryConditions);
            self.boundaryTemperature = @(x) bc(1) + bc(2)*x(1) + bc(3)*x(2) + bc(4)*x(1)*x(2);
            self.boundaryHeatFlux{1} = @(x) -(bc(3) + bc(4)*x);    %lower bound
            self.boundaryHeatFlux{2} = @(y) (bc(2) + bc(4)*y);     %right bound
            self.boundaryHeatFlux{3} = @(x) (bc(3) + bc(4)*x);     %upper bound
            self.boundaryHeatFlux{4} = @(y) -(bc(2) + bc(4)*y);    %left bound
        end
        
        function [cond, xi] = generateConductivityField(self, nSet, filename)
            %nSet is the number of the data set
            %nSet is the set (file number) index
            addpath('./genConductivity')    %to find genBochnerSamples
            addpath('./rom');
            
            % Draw conductivity/log conductivity
            fprintf('Generating finescale conductivity field...\n')
            tic
            
            l = self.conductivityDistributionParams{2};
            if self.nBochnerBasis
                %Fixed random field basis functions; the only random parameter is the coefficient xi
                fprintf('Using fixed basis of %u Bochner samples. Randomized length scale not possible.\n',...
                    self.nBochnerBasis)
                self.conductivityBasisFunctions = genBochnerSamplesFixedBasis(...
                    l, self.conductivityDistributionParams{3}, self.nBochnerBasis, self.conductivityDistribution);
                xi = normrnd(0, 1, self.nSets(nSet), self.nBochnerBasis);
            else
                error('Only fixed random field basis functions is implemented.')
            end
            basisFuns = self.conductivityBasisFunctions;
            save(filename, 'xi', 'basisFuns', '-v7.3')
            
            cond = self.getDiscretizedConductivityField(xi);
            save(filename, 'cond', '-append')
            fprintf('...done. Time: %.2f\n', toc)
        end
        
        function cond = getDiscretizedConductivityField(self, xi)
            %ATTENTION: only isotrop. distributions (length scales) possible. Compute coordinates of element centers
            if isempty(self.randFieldGenMatrix)
                self.randFieldGenMatrix = self.conductivityBasisFunctions(self.fineMesh.getCenterCoordinates())';
            end
            cond = self.lowerConductivity + exp(self.randFieldGenMatrix*xi');
        end
        
        function solveFEM(self, nSet, cond, filename)
            %Solve finite element model
            addpath('./computation')        %parpool
            fprintf('Solving finescale problem...\n')
            tic
            %cell array for parfor
            uf = mat2cell(zeros(self.fineMesh.nNodes, self.nSets(nSet)), ...
                                self.fineMesh.nNodes, ones(1, self.nSets(nSet)));

            %To avoid broadcasting overhead
            fineMesh = self.fineMesh;
            bcMean = str2num(self.boundaryConditions);
            bcVariance = self.boundaryConditionVariance;
            naturalNodes = self.naturalNodes;
            if(any(bcVariance))
                for i = 1:self.nSets(nSet)
                    bc{i} = mvnrnd(bcMean, diag(bcVariance));
                end
            else
                bc = [];
            end
            
            parPoolInit(self.nSets(nSet));
            for i = 1:self.nSets(nSet)
                if(any(bcVariance))
                    [bcPressure, bcFlux] = bcCoeffs2Fun(bc{i});
                    fineMeshTemp = fineMesh.setBoundaries(naturalNodes, bcPressure, bcFlux);
                else
                    fineMeshTemp = fineMesh;
                end
                FEMout = heat2d(fineMeshTemp, cond(:, i));
                %Store fine temperatures as a vector uf, use reshape(uf(:, i), domain.nElX + 1, domain.nElY + 1)
                %and then transpose result to reconvert it to original temperature field
                uf{i} = FEMout.u;
            end
            uf = cell2mat(uf);
            fprintf('FEM systems solved. Time: %.2fs\n', toc)
            
            if(nargin > 2)
                fprintf("saving finescale data to " + filename + "\n")
                save(filename,'uf', '-append')
                if(any(bcVariance))
                    save(filename, 'bc', '-append')
                end
                fprintf('done.')
            end
        end
        
        function genFineScaleDataPath(self)
            volFrac = self.conductivityDistributionParams{1};
            sigma_f2 = self.conductivityDistributionParams{3};
            if string(java.net.InetAddress.getLocalHost.getHostName) == "workstation1-room0436"
                %workstation
                self.fineScaleDataPath = "~/cluster/matlab/data/fineData/";
            else
                %cluster
                self.fineScaleDataPath = "~/matlab/data/fineData/";
            end
            self.fineScaleDataPath = sprintf(self.fineScaleDataPath + "systemSize=%ux%u/", self.nElFX, self.nElFY);
            %Type of conductivity distribution
            if(self.conductivityDistribution == "squaredExponential" || ...
               self.conductivityDistribution == "ornsteinUhlenbeck" || ...
               self.conductivityDistribution == "sincCov" || ...
               self.conductivityDistribution == "sincSqCov" || ...
               self.conductivityDistribution == "matern")
                if self.conductivityLengthScaleDist == "delta"
                    if(self.conductivityDistributionParams{2}(1) == self.conductivityDistributionParams{2}(2))
                        corrLength = self.conductivityDistributionParams{2}(1);
                    else
                        corrLength = self.conductivityDistributionParams{2};
                    end
                    self.fineScaleDataPath = sprintf(self.fineScaleDataPath + self.conductivityDistribution +...
                        '/l=%.2f_sigmafSq=%.2f/volumeFraction=%.2f/locond=%.2f',...
                        corrLength, sigma_f2, volFrac, self.lowerConductivity)
                    if self.conductivityDistributionType == "continuous"
                        self.fineScaleDataPath = self.fineScaleDataPath + "_upcond=continuous";
                    else
                        self.fineScaleDataPath =...
                            sprintf(self.fineScaleDataPath + "_upcond=%.2f", self.upperConductivity);
                    end
                    self.fineScaleDataPath = strcat(self.fineScaleDataPath, '/BCcoeffs=', self.boundaryConditions, '/');
                elseif strcmp(self.conductivityLengthScaleDist, 'lognormal')
                    corrLength1 = self.conductivityDistributionParams{2}(1);
                    corrLength2 = self.conductivityDistributionParams{2}(2);
                    self.fineScaleDataPath = sprintf(self.fineScaleDataPath + self.conductivityDistribution + ...
                        "/l=lognormal_mu=%.2f_sigma=%.2f_sigmafSq=.2f/volumeFraction=%.2f/locond=%.2f_upcond=%.2f" + ...
                        "/BCcoeffs=", self.boundaryConditions, '/', corrLength1, corrLength2, sigma_f2, volFrac,...
                        self.lowerConductivity, self.upperConductivity);
                else
                    error('Unknown length scale distribution')
                end
            elseif strcmp(cond_distribution, 'binary')
                self.fineScaleDataPath = sprintf(self.fineScaleDataPath + self.conductivityDistribution +...
                    '/volumeFraction=%.2f/locond=%.2f_upcond=%.2f/BCcoeffs='+ self.boundaryConditions + '/',...
                    volFrac, self.lowerConductivity, self.upperConductivity);
            else
                error('Unknown conductivity distribution')
            end
            
            if(any(self.boundaryConditionVariance))
                self.fineScaleDataPath = strcat(self.fineScaleDataPath, 'BCvar=[',...
                                                num2str(self.boundaryConditionVariance), ']/');
            end
            
            %Name of training data file
            trainFileName = sprintf("set1-samples=%u.mat", self.nSets(1));
            self.trainingDataMatfile = matfile(strcat(self.fineScaleDataPath, trainFileName));
            testFileName = sprintf("set2-samples=%u.mat", self.nSets(2));
            self.testDataMatfile = matfile(strcat(self.fineScaleDataPath, testFileName));
        end
        
        function loadTrainingData(self)
            %load data params; warning for variable FD can be ignored
            fprintf("\nloading data from \n" + self.fineScaleDataPath + "\n\n")
            load(self.fineScaleDataPath + "fineMesh.mat", "fineMesh");
            self.fineMesh = fineMesh;
            %for finescale domain class
            addpath('./heatFEM')
            %for boundary condition functions
            if(isempty(self.boundaryTemperature) || isempty(self.boundaryHeatFlux))
                self = self.genBoundaryConditionFunctions;
            end
            
            %there is no cum_lEl (cumulated finite element length) in old files
            if(~numel(self.fineMesh.cum_lElX) || ~numel(self.fineMesh.cum_lElX))
                self.fineMesh.cum_lElX = linspace(0, 1, self.fineMesh.nElX + 1);
                self.fineMesh.cum_lElY = linspace(0, 1, self.fineMesh.nElY + 1);
            end
            
            self.conductivityBasisFunctions = self.trainingDataMatfile.basisFuns;
%             self.conductivityCoefficients = self.trainingDataMatfile.xi(self.trainingSamples, :);
            
            %load finescale temperatures partially
            self.fineScaleDataOutput = self.trainingDataMatfile.uf(:, self.nStart:(self.nStart + self.nTrain - 1));
        end
        
        function genCoarseMesh(self)
            %Generate coarse domain object
            nX = length(self.coarseGridVectorX);
            nY = length(self.coarseGridVectorY);
            addpath('./heatFEM')        %to find Domain class
            self.coarseMesh = Domain(nX, nY, self.coarseGridVectorX, self.coarseGridVectorY);
            %ATTENTION: natural nodes have to be set manually
            %and consistently in coarse and fine scale domain!!
            self.coarseMesh.compute_grad = true;
            self.coarseMesh = setBoundaries(self.coarseMesh, [2:(2*nX + 2*nY)],...
                self.boundaryTemperature, self.boundaryHeatFlux);
            
            %Legacy, for predictions
            if ~exist('./data/', 'dir')
                mkdir('./data/');
            end
            filename = './data/coarseMesh.mat';
            coarseMesh = self.coarseMesh;
            save(filename, 'coarseMesh');
        end
        
        function estimateDataVariance(self)
            
            uftemp = self.trainingDataMatfile.uf(:, 1);
            uf_true_mean = zeros(size(uftemp));
            uf_true_sq_mean = zeros(size(uftemp));
            nSamples = self.nSets(1);
            window = nSamples;
            nWindows = ceil(nSamples/window);
            tic
            for i = 1:nWindows
                initial = 1 + (i - 1)*window;
                final = i*window;
                if final > nSamples
                    final = nSamples;
                end
                uftemp = self.trainingDataMatfile.uf(:, initial:final);
                uf_mean_temp = mean(uftemp, 2);
                uf_sq_mean_temp = mean(uftemp.^2, 2);
                clear uftemp;
                uf_true_mean = ((i - 1)/i)*uf_true_mean + (1/i)*uf_mean_temp;
                uf_true_sq_mean = ...
                    ((i - 1)/i)*uf_true_sq_mean + (1/i)*uf_sq_mean_temp;
            end
            
            uf_true_var = uf_true_sq_mean - uf_true_mean.^2;
            self.outputMean = uf_true_mean;
            self.outputVariance = uf_true_var;
            self.meanOutputVariance = mean(uf_true_var);
            toc
            
            %Mean log likelihood
            uf_true_var(self.fineMesh.essentialNodes) = NaN;
            nNatNodes = self.fineMesh.nNodes -...
                numel(self.fineMesh.essentialNodes);
            Lm = -.5*log(2*pi) - .5 - .5*(1/nNatNodes)*...
                sum(log(uf_true_var), 'omitnan');
            sv = true;
            if sv
                %     savedir = '~/matlab/data/trueMC/';
                savedir = self.fineScaleDataPath;
                if ~exist(savedir, 'dir')
                    mkdir(savedir);
                end
                save(strcat(savedir, 'trueMC', '_nSamples=',...
                    num2str(nSamples), '.mat'), 'uf_true_mean',...
                    'uf_true_var', 'Lm')
            end
        end
        
        function [meanLogLikelihood, err]= estimateLogL(self, nTrainingData, uf)
            
            natNodes = true(self.fineMesh.nNodes, 1);
            natNodes(self.fineMesh.essentialNodes) = false;
            nNatNodes = sum(natNodes);
            uf = uf(natNodes, :);
            meanLogLikelihood = 0;
            meanLogLikelihoodSq = 0;
            converged = false;
            i = 1;
            while(~converged)
                randSamples = randperm(self.nSets(1));
                randSamples_params = randSamples(1:nTrainingData);
                randSamples_samples =...
                    randSamples((nTrainingData + 1):(2*nTrainingData));
                uftemp = uf(:, randSamples_params);
                mu_data = mean(uftemp, 2);
                var_data = var(uftemp')';
                
                term1 = .5*log(geomean(var_data));
                term2 = .5*mean(mean(...
                    (uf(:, randSamples_samples) - mu_data).^2, 2)./var_data);
                
                meanLogLikelihood =...
                    ((i - 1)/i)*meanLogLikelihood + (1/i)*(term1 + term2);
                meanLogLikelihoodSq =...
                    ((i - 1)/i)*meanLogLikelihoodSq + (1/i)*(term1 + term2)^2;
                if(mod(i, 1000) == 0)
                    i
                    meanLogLikelihood
                    err = sqrt((meanLogLikelihoodSq - meanLogLikelihood^2)/i);
                    relErr = abs(err/meanLogLikelihood)
                    if((relErr < 1e-2 && i > 2e3) || i > 1e5)
                        converged = true;
                    end
                end
                i = i + 1;
            end
            meanLogLikelihood = meanLogLikelihood + .5*log(2*pi);
            meanLogLikelihood = meanLogLikelihood + .87; %remove this!
            err = sqrt((meanLogLikelihoodSq - meanLogLikelihood^2)/i);
            
        end
        
        function [Xopt, LambdaOpt, s2] = detOpt_p_cf(self, nStart, nTrain)
            %Det. optimization of log(p_cf) to check capabilities of model
            
            %don't change these!
            theta_cfOptim.S = 1;
            theta_cfOptim.sumLogS = 0;
            theta_cfOptim.Sinv = 1;
            theta_cfOptim.Sinv_vec = ones(self.fineMesh.nNodes, 1);
            theta_cfOptim.W = self.theta_cf.W;
            theta_cfOptim.WTSinv = self.theta_cf.WTSinv;
            
            options = optimoptions(@fminunc,'Display','iter',...
                'Algorithm', 'trust-region', 'SpecifyObjectiveGradient', true);
            Xinit = 0*ones(self.coarseMesh.nEl, 1);
            Xopt = zeros(self.coarseMesh.nEl, nTrain);
            LambdaOpt = Xopt;
            s2 = zeros(1, nTrain);
            j = 1;
            addpath('./tests/detOptP_cf')
            for i = nStart:(nStart + nTrain -1)
                uf = self.trainingDataMatfile.uf(:, i);
                objFun = @(X) objective(X, uf, self.coarseMesh,...
                    self.conductivityTransformation, theta_cfOptim);
                [XoptTemp, fvalTemp] = fminunc(objFun, Xinit, options);
                LambdaOptTemp = conductivityBackTransform(XoptTemp,...
                    self.conductivityTransformation);
                Xopt(:, j) = XoptTemp;
                LambdaOpt(:, j) = LambdaOptTemp;
                
                %s2 is the squared distance of truth to optimal 
                %coarse averaged over all nodes
                s2(j) = fvalTemp/self.fineMesh.nNodes
                j = j + 1;
            end
        end
        
        function loadTrainedParams(self)
            %Load trained model parameters from disk to workspace
            
            
            if exist(strcat('./data/coarseMesh.mat'), 'file')
                load(strcat('./data/coarseMesh.mat'));
                self.coarseMesh = coarseMesh;
            else
                warning(strcat('No coarse domain file found.',...
                    'Take b.c.s from finescale data and regenerate.',...
                    'Please make sure everything is correct!'))
                self.genCoarseDomain;
            end
            %Load trained params from disk
            disp('Loading optimal parameters from disk...')
            self.thetaPriorHyperparam = dlmread('./data/thetaPriorHyperparam');
            self.thetaPriorHyperparam = self.thetaPriorHyperparam(end, :);
            self.theta_c.theta = dlmread('./data/theta');
            self.theta_c.theta = self.theta_c.theta(end, :)';
            self.theta_c.Sigma = dlmread('./data/sigma');
            if(numel(self.theta_c.Sigma) == self.coarseMesh.nEl)
                self.theta_c.Sigma = diag(self.theta_c.Sigma(end, :));
            else
                self.theta_c.Sigma = diag(self.theta_c.Sigma(end, :));
            end
            self.theta_cf.S = dlmread('./data/S')';
            W = dlmread('./data/Wmat');
            W = reshape(W, length(W)/3, 3)';
            self.theta_cf.W = sparse(W(1, :), W(2, :), W(3, :));
            self.theta_cf.mu = dlmread('./data/mu')';
            disp('done')
            
            disp('Loading data normalization data...')
            try
                self.featureFunctionMean =...
                    dlmread('./data/featureFunctionMean');
                self.featureFunctionSqMean = ...
                    dlmread('./data/featureFunctionSqMean');
            catch
                warning('featureFunction(Sq)Mean not found, setting it to 0.')
                self.featureFunctionMean = 0;
                self.featureFunctionSqMean = 0;
            end
            
            try
                self.featureFunctionMin = dlmread('./data/featureFunctionMin');
                self.featureFunctionMax = dlmread('./data/featureFunctionMax');
            catch
                warning('featureFunctionMin/Max not found, setting it to 0.')
                self.featureFunctionMin = 0;
                self.featureFunctionMax = 0;
            end
            disp('done')
            
            if(isempty(self.coarseMesh) || isempty(self.fineMesh))
                disp('Loading fine and coarse domain objects...')
                addpath('./heatFEM')        %to find Domain class
                try
                    load(strcat(self.fineScaleDataPath, 'fineMesh.mat'));
                    self.fineMesh = fineMesh;
                catch
                    temp = load(strcat(self.fineScaleDataPath, 'romObj.mat'));
                    self.fineMesh = temp.obj.fineMesh;
                end
                disp('done')
            end
        end
        
        function M_step(self)
            tic
            fprintf("M-step: find optimal params...")
            %Optimal S (decelerated convergence)
            lowerBoundS = eps;
            self.theta_cf.S = (1 - self.mix_S)*self.varExpect_p_cf_exp_mean...
                + self.mix_S*self.theta_cf.S + lowerBoundS*ones(self.fineMesh.nNodes, 1);

            if self.free_W
                self.mean_ucucT(self.coarseMesh.essentialNodes, :) = [];
                self.mean_ucucT(:, self.coarseMesh.essentialNodes) = [];
                if isempty(self.fineMesh.essentialNodes)
                    warning("No essential nodes stored in fineMesh. Setting first node to be essential.")
                    self.fineMesh.essentialNodes = 1;
                end
                self.mean_ufucT(self.fineMesh.essentialNodes, :) = [];
                self.mean_ufucT(:, self.coarseMesh.essentialNodes) = [];
                W_temp = self.mean_ufucT/self.mean_ucucT;
                
                natNodesFine = 1:self.fineMesh.nNodes;
                natNodesFine(self.fineMesh.essentialNodes) = [];
                natNodesCoarse = 1:self.coarseMesh.nNodes;
                natNodesCoarse(self.coarseMesh.essentialNodes) = [];
                self.theta_cf.W(natNodesFine, natNodesCoarse) = W_temp;
            end
            
            self.theta_cf.Sinv = sparse(1:self.fineMesh.nNodes,...
                1:self.fineMesh.nNodes, 1./self.theta_cf.S);
            self.theta_cf.Sinv_vec = 1./self.theta_cf.S;
            %Precomputation for efficiency
            self.theta_cf.WTSinv = self.theta_cf.W'*self.theta_cf.Sinv;
            
            %optimal theta_c and sigma
            Sigma_old = self.theta_c.Sigma;
            theta_old = self.theta_c.theta;
            
            self.updateTheta_c;

            self.theta_c.Sigma =...
                (1 - self.mix_sigma)*self.theta_c.Sigma +...
                self.mix_sigma*Sigma_old;
            self.theta_c.theta =...
                (1 - self.mix_theta)*self.theta_c.theta +...
                self.mix_theta*theta_old;
            
            fprintf("   M-step done. Time: %.2fs\n", toc)
        end
        
        function updateTheta_c(self)
            %Find optimal theta_c and sigma
            dim_theta = numel(self.theta_c.theta);
            
            %% Solve self-consistently: compute optimal sigma2, then theta, 
            %then sigma2 again and so on
            %Start from previous best estimate
            theta = self.theta_c.theta;
            I = speye(dim_theta);
            Sigma = self.theta_c.Sigma;
            
            %sum_i Phi_i^T Sigma^-1 <X^i>_qi
            sumPhiTSigmaInvXmean = 0;
            SigmaInv = self.theta_c.SigmaInv;
            SigmaInvXMean = SigmaInv*self.XMean;
            sumPhiTSigmaInvPhi = 0;
            PhiThetaMat = zeros(self.coarseMesh.nEl, self.nTrain);
            
            for n = 1:self.nTrain
                sumPhiTSigmaInvXmean = sumPhiTSigmaInvXmean +...
                    self.designMatrix{n}'*SigmaInvXMean(:, n);
                sumPhiTSigmaInvPhi = sumPhiTSigmaInvPhi +...
                    self.designMatrix{n}'*SigmaInv*self.designMatrix{n};
                PhiThetaMat(:, n) = self.designMatrix{n}*self.theta_c.theta;
            end
            
            stabilityParam = 1e-2;    %for stability in matrix inversion
            if(strcmp(self.thetaPriorType, 'adaptiveGaussian') ||...
                    strcmp(self.thetaPriorType, 'RVM') || ...
                    strcmp(self.thetaPriorType, 'sharedRVM'))
                %Find prior hyperparameter by max marginal likelihood
                converged = false;
                iter = 0;
                if strcmp(self.thetaPriorType, 'adaptiveGaussian')
                    self.thetaPriorHyperparam = 1;
                elseif(strcmp(self.thetaPriorType, 'RVM') ||...
                        strcmp(self.thetaPriorType, 'sharedRVM'))
                    if(numel(self.thetaPriorHyperparam) ~= dim_theta)
                        warning('resizing theta hyperparam')
                        self.thetaPriorHyperparam = 1e-4*ones(dim_theta, 1);
                    end
					nElc = self.coarseMesh.nEl;
					nFeatures = dim_theta/nElc; %for shared RVM
                end
                while(~converged)
                    if strcmp(self.thetaPriorType, 'adaptiveGaussian')
                        SigmaTilde = inv(sumPhiTSigmaInvPhi +...
                            (self.thetaPriorHyperparam + stabilityParam)*I);
                        muTilde = SigmaTilde*sumPhiTSigmaInvXmean;
                        theta_prior_hyperparam_old = self.thetaPriorHyperparam;
                        self.thetaPriorHyperparam = dim_theta/...
                            (muTilde'*muTilde + trace(SigmaTilde));
                    elseif(strcmp(self.thetaPriorType, 'RVM') ||...
                            strcmp(self.thetaPriorType, 'sharedRVM'))
                        SigmaTilde = inv(sumPhiTSigmaInvPhi + diag(self.thetaPriorHyperparam));
                        muTilde = (sumPhiTSigmaInvPhi + diag(self.thetaPriorHyperparam))\sumPhiTSigmaInvXmean;
                        theta_prior_hyperparam_old = self.thetaPriorHyperparam;
                        if strcmp(self.thetaPriorType, 'RVM')
                            self.thetaPriorHyperparam = 1./(muTilde.^2 + diag(SigmaTilde));
                        elseif strcmp(self.thetaPriorType, 'sharedRVM')
                            muTildeSq = muTilde.^2;
                            varTilde = diag(SigmaTilde);
                            lambdaInv = (1/nElc)*(sum(reshape(muTildeSq, nFeatures, nElc), 2) +...
                                sum(reshape(varTilde, nFeatures, nElc), 2));
                            self.thetaPriorHyperparam = repmat(1./lambdaInv, nElc, 1);
                        end
                        self.thetaPriorHyperparam = self.thetaPriorHyperparam + stabilityParam;
                    end
                    crit = norm(1./self.thetaPriorHyperparam -...
                        1./theta_prior_hyperparam_old)/norm(1./self.thetaPriorHyperparam);
                    if(crit < 1e-5 || iter >= 5)
                        converged = true;
                    elseif(any(~isfinite(self.thetaPriorHyperparam)) || any(self.thetaPriorHyperparam <= 0))
                        converged = true;
                        muTilde.^2
                        self.thetaPriorHyperparam
                        self.thetaPriorHyperparam = ones(dim_theta, 1);
                        warning(strcat('Gaussian hyperparameter precision',...
                            ' is negative or not a number. Setting it to 1.'))
                    end
                    iter = iter + 1;
                end
%                 lambda_end = [obj.thetaPriorHyperparam (1:100)']
            end
            
            linsolveOpts.SYM = true;
            linsolveOpts.POSDEF = true;
            iter = 0;
            converged = false;
            while(~converged)
                theta_old = theta;  %to check for iterative convergence
                
                %Matrix M is pos. def., invertible even if badly conditioned
                %warning('off', 'MATLAB:nearlySingularMatrix');
                if strcmp(self.thetaPriorType, 'hierarchical_laplace')
                    offset = 1e-30;
                    U = diag(sqrt((abs(theta) + offset)/...
                        self.thetaPriorHyperparam(1)));
                elseif strcmp(self.thetaPriorType, 'hierarchical_gamma')
                    U = diag(sqrt((.5*abs(theta).^2 +...
                        self.thetaPriorHyperparam(2))./...
                        (self.thetaPriorHyperparam(1) + .5)));
                elseif(strcmp(self.thetaPriorType, 'gaussian') ||...
                        strcmp(self.thetaPriorType, 'adaptiveGaussian'))
                    sumPhiTSigmaInvPhi = sumPhiTSigmaInvPhi +...
                        (self.thetaPriorHyperparam(1) + stabilityParam)*I;
                elseif(strcmp(self.thetaPriorType, 'RVM') ||...
                        strcmp(self.thetaPriorType, 'sharedRVM'))
                    sumPhiTSigmaInvPhi =...
                        sumPhiTSigmaInvPhi + diag(self.thetaPriorHyperparam);
                elseif strcmp(self.thetaPriorType, 'none')
                else
                    error('Unknown prior on theta_c')
                end
                
                if (strcmp(self.thetaPriorType, 'gaussian') ||...
                        strcmp(self.thetaPriorType, 'RVM') ||...
                        strcmp(self.thetaPriorType, 'none') ||...
                        strcmp(self.thetaPriorType, 'adaptiveGaussian') || ...
                        strcmp(self.thetaPriorType, 'sharedRVM'))
                    theta_temp = sumPhiTSigmaInvPhi\sumPhiTSigmaInvXmean;
                    %is this true? we do not need to iteratively maximize theta
                    converged = true;
                else
                    %theta_temp =...
                    %   U*((U*sumPhiTSigmaInvPhi*U + I)\U)*sumPhiTSigmaInvXmean;
                    %A= linsolve((U*sumPhiTSigmaInvPhi*U + I), U, linsolveOpts);
                    theta_temp = U*linsolve((U*sumPhiTSigmaInvPhi*U + I),...
                        U, linsolveOpts)*sumPhiTSigmaInvXmean;
                end
                [~, msgid] = lastwarn;     %to catch nearly singular matrix
                
                if(strcmp(msgid, 'MATLAB:singularMatrix') ||...
                        strcmp(msgid, 'MATLAB:nearlySingularMatrix')...
                        || strcmp(msgid, 'MATLAB:illConditionedMatrix') ||...
                        norm(theta_temp)/length(theta) > 1e8)
                    warning(strcat('theta_c is assuming unusually large' , ...
                        ' values. Only go small step.'))
                    theta = .5*(theta + .1*(norm(theta)/...
                        norm(theta_temp))*theta_temp)
                    if any(~isfinite(theta))
                        %restart from 0
                        warning(strcat('Some components of theta are not',...
                            ' finite. Restarting from theta = 0...'))
                        theta = 0*theta;
                    end
                else
                    theta = theta_temp;
                end
                theta= theta_temp;
                
                PhiThetaMat = zeros(self.coarseMesh.nEl, self.nTrain);
                for n = 1:self.nTrain
                    PhiThetaMat(:, n) = self.designMatrix{n}*theta;
                end
                if(self.EM_iterations < self.fixSigmaInit)
                    sigma_prior_type = 'delta';
                else
                    sigma_prior_type = self.sigmaPriorType;
                end
                if strcmp(sigma_prior_type, 'none')
                    if self.theta_c.full_Sigma
                        sumPhiThetaPhiTThetaT = 0;
                        sumPhiThetaXT = 0;
                        sumXXT = 0;
                        for n = 1:self.nTrain
                            sumPhiThetaPhiTThetaT = sumPhiThetaPhiTThetaT +...
                                PhiThetaMat(:, n)*PhiThetaMat(:, n)';
                            sumPhiThetaXT = sumPhiThetaXT +...
                                PhiThetaMat(:, n)*self.XMean(:, n)';
                            sumXXT= sumXXT + self.XMean(:, n)*self.XMean(:, n)';
                        end
                        Sigma = diag(mean(self.XSqMean, 2)) +...
                            (sumXXT - diag(diag(sumXXT)))/self.nTrain +...
                            (sumPhiThetaPhiTThetaT/self.nTrain) -...
                            (sumPhiThetaXT + sumPhiThetaXT')/self.nTrain;
                    else
                        Sigma = sparse(1:self.coarseMesh.nEl,...
                            1:self.coarseMesh.nEl,...
                            mean(self.XSqMean -...
                            2*(PhiThetaMat.*self.XMean)+PhiThetaMat.^2, 2));
                    end
                    %Sigma(Sigma < 0) = eps; %for numerical stability
                    %Variances must be positive
                    Sigma(logical(eye(size(Sigma)))) =...
                        abs(Sigma(logical(eye(size(Sigma))))) + eps;
                    
                    
                    %sum_i Phi_i^T Sigma^-1 <X^i>_qi
                    sumPhiTSigmaInvXmean = 0;
                    %Only valid for diagonal Sigma
%                     s = diag(Sigma);
%                     SigmaInv = sparse(diag(1./s));
                    SigmaInv = inv(Sigma);
%                     SigmaInvXMean = SigmaInv*obj.XMean;
                    SigmaInvXMean = Sigma\self.XMean;
                    sumPhiTSigmaInvPhi = 0;
                    
                    for n = 1:self.nTrain
                        sumPhiTSigmaInvXmean = sumPhiTSigmaInvXmean +...
                            self.designMatrix{n}'*SigmaInvXMean(:, n);
                        sumPhiTSigmaInvPhi = sumPhiTSigmaInvPhi +...
                            self.designMatrix{n}'*(Sigma\self.designMatrix{n});
                    end
                elseif strcmp(sigma_prior_type, 'expSigSq')
                    %Vector of variances
                    sigmaSqVec = .5*(1./self.sigmaPriorHyperparam).*...
                        (-.5*self.nTrain +...
                        sqrt(2*self.nTrain*self.sigmaPriorHyperparam.*...
                        mean(self.XSqMean - 2*(PhiThetaMat.*self.XMean) +...
                        PhiThetaMat.^2, 2) + .25*self.nTrain^2));
                    
                    Sigma = sparse(1:self.coarseMesh.nEl,...
                        1:self.coarseMesh.nEl, sigmaSqVec);
                    Sigma(Sigma < 0) = eps; %for numerical stability
                    sumPhiTSigmaInvXmean = 0;
                    %Only valid for diagonal Sigma
                    s = diag(Sigma);
                    SigmaInv = sparse(diag(1./s));
                    SigmaInvXMean = SigmaInv*self.XMean;
                    sumPhiTSigmaInvPhi = 0;
                    
                    for n = 1:self.nTrain
                        sumPhiTSigmaInvXmean = sumPhiTSigmaInvXmean +...
                            self.designMatrix{n}'*SigmaInvXMean(:, n);
                        sumPhiTSigmaInvPhi = sumPhiTSigmaInvPhi +...
                            self.designMatrix{n}'*SigmaInv*self.designMatrix{n};
                    end
                elseif strcmp(sigma_prior_type, 'delta')
                    %Don't change sigma
                    %sum_i Phi_i^T Sigma^-1 <X^i>_qi
                    sumPhiTSigmaInvXmean = 0;
                    %Only valid for diagonal Sigma
                    s = diag(Sigma);
                    SigmaInv = sparse(diag(1./s));
                    SigmaInvXMean = SigmaInv*self.XMean;
                    sumPhiTSigmaInvPhi = 0;
                    
                    for n = 1:self.nTrain
                        sumPhiTSigmaInvXmean = sumPhiTSigmaInvXmean +...
                            self.designMatrix{n}'*SigmaInvXMean(:, n);
                        sumPhiTSigmaInvPhi = sumPhiTSigmaInvPhi +...
                            self.designMatrix{n}'*SigmaInv*self.designMatrix{n};
                    end
                end
                
                iter = iter + 1;
                thetaDiffRel=norm(theta_old - theta)/(norm(theta)*numel(theta));
                if((iter > 5 && thetaDiffRel < 1e-8) || iter > 200)
                    converged = true;
                end
            end
            
            self.theta_c.theta = theta;
            self.theta_c.Sigma = Sigma;
            self.theta_c.SigmaInv = SigmaInv;
            self.thetaPriorHyperparam = self.thetaPriorHyperparam;
%             noPriorSigma = mean(obj.XSqMean - 2*(PhiThetaMat.*obj.XMean) +...
%                 PhiThetaMat.^2, 2);
%             save('./data/noPriorSigma.mat', 'noPriorSigma');
            
        end
        
        function dispCurrentParams(self)
            disp('Current params:')
            [~, index] = sort(abs(self.theta_c.theta));
            if strcmp(self.mode, 'useNeighbor')
                feature = mod((index - 1), numel(self.featureFunctions)) + 1;
                %counted counterclockwise from right to lower neighbor
                neighborElement=floor((index - 1)/numel(self.featureFunctions));
                curr_theta = [self.theta_c.theta(index) feature neighborElement]
            elseif strcmp(self.mode, 'useDiagNeighbor')
                feature = mod((index - 1), numel(self.featureFunctions)) + 1;
                %counted counterclockwise from right to lower right neighbor
                neighborElement=floor((index - 1)/numel(self.featureFunctions));
                curr_theta = [self.theta_c.theta(index) feature neighborElement]
            elseif strcmp(self.mode, 'useLocal')
                feature = mod((index - 1), size(self.featureFunctions, 2) +...
                    size(self.globalFeatureFunctions, 2)) + 1;
                Element= floor((index - 1)/(size(self.featureFunctions, 2) +...
                    size(self.globalFeatureFunctions, 2))) + 1;
                curr_theta = [self.theta_c.theta(index) feature Element]
            elseif(strcmp(self.mode, 'useLocalNeighbor') ||...
                    strcmp(self.mode, 'useLocalDiagNeighbor'))
                disp('theta feature coarseElement neighbor')
                curr_theta = [self.theta_c.theta(index), self.neighborDictionary(index, 1),...
                    self.neighborDictionary(index, 2), self.neighborDictionary(index, 3)]
            else
                curr_theta = [self.theta_c.theta(index) index]
            end
            
            curr_sigma = self.theta_c.Sigma
            if self.theta_c.full_Sigma
                diag_sigma = diag(self.theta_c.Sigma)
            end
            fprintf("Mean reconstruction variance <S> = %f", mean(self.theta_cf.S))
            %curr_theta_hyperparam = obj.thetaPriorHyperparam
        end
        
        function [lambda_theta_c, lambda_log_s2, lambda_log_sigma2] = laplaceApproximation(self)
            %Computes parameter precisions based on second derivatives 
            %of posterior lower bound F
            
            %p_cf
            %precision of variances s_n^2
            lambda_log_s2 = .5*self.nTrain*ones(1, self.fineMesh.nNodes);
            
            %p_c
            %precision of variances sigma_k^2
            lambda_log_sigma2 = .5*ones(1, self.coarseMesh.nEl);
            
            %precision on the theta_c's
            lambda_theta_c = zeros(length(self.theta_c.theta),...
                length(self.theta_c.theta));
            if isempty(self.designMatrix)
                load('./persistentData/trainDesignMatrix.mat');
                self.designMatrix = designMatrix;
            end
            for n = 1:self.nTrain
                lambda_theta_c = lambda_theta_c + ...
                    self.designMatrix{n}'*...
                    (self.theta_c.Sigma\self.designMatrix{n});
            end
            %To ensure symmetry on machine precision; 
            %this shouldn't change anything
            lambda_theta_c = .5*(lambda_theta_c + lambda_theta_c');
            
            %Contribution from prior;
            %the contribution from a Laplacian prior is 0
            if strcmp(self.thetaPriorType, 'RVM')
                lambda_theta_c =...
                    lambda_theta_c + diag(self.thetaPriorHyperparam);
            elseif(strcmp(self.thetaPriorType, 'gaussian') ||...
                    strcmp(self.thetaPriorType, 'adaptiveGaussian'))
                lambda_theta_c = lambda_theta_c +...
                    self.thetaPriorHyperparam*eye(length(self.theta_c.theta));
            elseif strcmp(self.thetaPriorType, 'hierarchical_laplace')
                warning(strcat('Hessian of Laplace distribution',...
                    ' is ill-defined. Ignoring contribution from prior.'))
            elseif strcmp(self.thetaPriorType, 'hierarchical_gamma')
                lambda_theta_c = lambda_theta_c + ...
                    (self.thetaPriorHyperparam(1) + .5)*...
                    diag(1./(self.thetaPriorHyperparam(2) +...
                .5*self.theta_c.theta.^2) - (self.theta_c.theta.^2)./...
                    ((self.thetaPriorHyperparam(2) +...
                    .5*self.theta_c.theta.^2).^2));
            end
        end
        
        function p_ROM(xi, uf)
            %Takes fine scale basis function coefficients xi and computes
            %the probability density p(uf|xi) for a fixed uf
            
        end
        
        function predict(self, mode, boundaryConditions)
            %Function to predict finescale output from generative model
            
            if(nargin < 2)
                mode = "test";
            end
            
            if(nargin > 2)
                %Predict on different boundary conditions
                self.setBoundaryConditions(boundaryConditions);
                
                %Set up coarseMesh;
                %must be done after boundary conditions are set up
                self.genCoarseDomain;
            end
                
            %Load test file
            if strcmp(mode, 'self')
                %Self-prediction on training set
                uf = self.trainingDataMatfile.uf(:, self.nStart:(self.nStart + self.nTrain - 1));
                tempVars = whos(self.trainingDataMatfile);
                bcVar = ismember('bc', {tempVars.name});
            else
                uf = self.testDataMatfile.uf(:, self.testSamples);
                tempVars = whos(self.testDataMatfile);
                bcVar = ismember('bc', {tempVars.name});
            end
            self.loadTrainedParams;
            if strcmp(mode, 'self')
                self.computeDesignMatrix('train');
                designMatrixPred = self.designMatrix;
            else
                self.computeDesignMatrix('test');
                designMatrixPred = self.testDesignMatrix;
            end
            
            %% Sample from p_c
            disp('Sampling from p_c...')
            if strcmp(mode, 'self')
                nTest = self.nTrain;
            else
                nTest = numel(self.testSamples);
            end
            
            %short hand notation/ avoiding broadcast overhead
            nElc = self.coarseMesh.nEl;
            nSamples = self.nSamples_p_c;
            %bcVar = any(obj.boundaryConditionVariance);
            nFineNodes = self.fineMesh.nNodes;
            
            Xsamples = zeros(nElc, nSamples, nTest);
            LambdaSamples{1} = zeros(nElc, nSamples);
            LambdaSamples = repmat(LambdaSamples, nTest, 1);
            self.meanEffCond = zeros(nElc, nTest);
            
            if self.useLaplaceApproximation
                [precisionTheta, precisionLogS, precisionLogSigma] =...
                    self.laplaceApproximation;
                SigmaTheta = inv(precisionTheta) +...
                    eps*eye(numel(self.theta_c.theta));
%                 stdLogS = sqrt(1./precisionLogS)';
%                 stdLogSigma = sqrt(1./precisionLogSigma);
                stdLogS = 0;
                stdLogSigma = 0;
            else
                stdLogS = [];   %for parfor
            end
            
            for i = 1:nTest
                %Samples from p_c
                if self.useLaplaceApproximation
                    %First sample theta from Laplace approx, then sample X
                    theta = mvnrnd(self.theta_c.theta', SigmaTheta, nSamples)';
                    for j = 1:nSamples
                        Sigma2 = exp(normrnd(log(diag(self.theta_c.Sigma))',...
                            stdLogSigma));
%                         Sigma2 = diag(obj.theta_c.Sigma)
                        Xsamples(:, j, i) =...
                            mvnrnd((designMatrixPred{i}*theta(:, j))', Sigma2)';
                    end
                else
                    Xsamples(:, :, i) =...
                        mvnrnd((designMatrixPred{i}*self.theta_c.theta)',...
                        self.theta_c.Sigma, nSamples)';
                end
                %Conductivities
                LambdaSamples{i} =...
                    conductivityBackTransform(Xsamples(1:nElc, :, i),...
                    self.conductivityTransformation);
                if(strcmp(self.conductivityTransformation.type, 'log'))
                    self.meanEffCond(:, i) = exp(self.testDesignMatrix{i}*...
                        self.theta_c.theta + .5*diag(self.theta_c.Sigma));
                else
                    self.meanEffCond(:, i) = mean(LambdaSamples{i}, 2);
                end
            end
            disp('done')
            
            %% Run coarse model and sample from p_cf
            disp('Solving coarse model and sample from p_cf...')
            ufMeanArray{1} = zeros(self.fineMesh.nNodes, 1);
            ufMeanArray = repmat(ufMeanArray, nTest, 1);
            ufVarArray = ufMeanArray;
            uf_sq_mean = ufMeanArray;
            
            if(bcVar)
                %Set coarse domain for data with different boundary conditions
                nX = self.coarseMesh.nElX;
                nY = self.coarseMesh.nElY;
                if strcmp(mode, 'self')
                    bc = self.trainingDataMatfile.bc;
                else
                    bc = self.testDataMatfile.bc;
                end
                for j = 1:nTest
                    if strcmp(mode, 'self')
                        i = self.trainingSamples(j);
                    else
                        i = self.testSamples(j);
                    end
                    bcT = @(x) bc{i}(1) + bc{i}(2)*x(1) +...
                        bc{i}(3)*x(2) + bc{i}(4)*x(1)*x(2);
                    bcQ{1} = @(x) -(bc{i}(3) + bc{i}(4)*x);      %lower bound
                    bcQ{2} = @(y) (bc{i}(2) + bc{i}(4)*y);       %right bound
                    bcQ{3} = @(x) (bc{i}(3) + bc{i}(4)*x);       %upper bound
                    bcQ{4} = @(y) -(bc{i}(2) + bc{i}(4)*y);      %left bound
                    cd(i) = self.coarseMesh;
                    cd(i) = cd(i).setBoundaries([2:(2*nX + 2*nY)], bcT, bcQ);
                end
            else
                cd = self.coarseMesh;
            end
            t_cf = self.theta_cf;
            lapAp = self.useLaplaceApproximation;
            %             t_c = obj.theta_c;
            natNodes = true(self.fineMesh.nNodes, 1);
            natNodes(self.fineMesh.essentialNodes) = false;
            nNatNodes = sum(natNodes);
            addpath('./heatFEM');
            parfor j = 1:nTest
                if bcVar
                    coarseDomain = cd(j);
                else
                    coarseDomain = cd;
                end
                for i = 1:nSamples
                    FEMout = heat2d(coarseDomain, LambdaSamples{j}(:, i));
                    uctemp = FEMout.u';
                    
                    %sample from p_cf
                    mu_cf = t_cf.mu + t_cf.W*uctemp(:);
                    %only for diagonal S!!
                    %Sequentially compute mean and <uf^2> to save memory
                    %U_f-integration can be done analyt.
                    ufMeanArray{j} = ((i - 1)/i)*ufMeanArray{j} + (1/i)*mu_cf;
                    uf_sq_mean{j} = ((i - 1)/i)*uf_sq_mean{j} + (1/i)*mu_cf.^2;
                end
                if lapAp
%                     S = exp(normrnd(log(t_cf.S), stdLogS));
                    S = t_cf.S;
                    uf_sq_mean{j} = uf_sq_mean{j} + S;
                else
                    uf_sq_mean{j} = uf_sq_mean{j} + t_cf.S;
                end
                %abs to avoid negative variance due to numerical error
                uf_var = abs(uf_sq_mean{j} - ufMeanArray{j}.^2);
                meanuf_meanMCErr = mean(sqrt(uf_var/nSamples))
                ufVarArray{j} = uf_var;
                
                meanMahaErrTemp{j} = mean(sqrt(abs((1./(uf_var)).*(uf(:, j) - ufMeanArray{j}).^2)));
                sqDist{j} = (uf(:, j) - ufMeanArray{j}).^2;
                meanSqDistTemp{j} = mean(sqDist{j});
                
                uf_var_nat = uf_var(natNodes);
                logLikelihood{j} = -.5*nNatNodes*log(2*pi) - .5*sum(log(uf_var_nat), 'omitnan') - ...
                    .5*sum(sqDist{j}(natNodes)./uf_var_nat, 'omitnan');
                logPerplexity{j} = -(1/(nNatNodes))*logLikelihood{j};
            end
            
            self.meanPredMeanOutput = mean(cell2mat(ufMeanArray'), 2);
            self.meanMahalanobisError = mean(cell2mat(meanMahaErrTemp));
            self.meanSquaredDistanceField = mean(cell2mat(sqDist), 2);
            self.meanSquaredDistance = mean(cell2mat(meanSqDistTemp));
            self.squaredDistance = meanSqDistTemp;
            norm_data = sqrt(sum(uf.^2));
            self.normError = cell2mat(meanSqDistTemp)./norm_data;
            meanSqDistSq = mean(cell2mat(meanSqDistTemp).^2);
            self.meanSquaredDistanceError = sqrt((meanSqDistSq -...
                self.meanSquaredDistance^2)/nTest);
            self.meanLogLikelihood = mean(cell2mat(logLikelihood))/nNatNodes;
            self.meanLogPerplexity = mean(cell2mat(logPerplexity));
            self.meanPerplexity = exp(self.meanLogPerplexity);
            storeArray = false;
            if storeArray
                self.predMeanArray = ufMeanArray;
                self.predVarArray = ufVarArray;
            end
            
            plotPrediction = true;
            if plotPrediction
                f = figure('units','normalized','outerposition',[0 0 1 1]);
                pstart = 1;
                j = 1;
                max_uf = max(max(uf(:, pstart:(pstart + 5))));
                min_uf = min(min(uf(:, pstart:(pstart + 5))));
                if strcmp(mode, 'self')
                    cond = self.trainingDataMatfile.cond(:,...
                        self.nStart:(self.nStart + self.nTrain - 1));
                    cond = cond(:, pstart:(pstart + 5));
                else
                    cond = self.testDataMatfile.cond(:,...
                        self.testSamples(1):(self.testSamples(1) + 5));
                end
                %to use same color scale
                cond = ((min_uf - max_uf)/(min(min(cond)) - max(max(cond))))*cond + max_uf - ...
                    ((min_uf - max_uf)/(min(min(cond)) - max(max(cond))))*max(max(cond));
                for i = pstart:(pstart + 5)
                    ax = subplot(2, 3, j);
                    uf_i_min = min(uf(:, i))
                    s(j, 1) = surf(reshape(uf(:, i) - uf_i_min, (self.nElFX + 1), (self.nElFY + 1)));
                    s(j, 1).LineStyle = 'none';
                    hold on;
                    s(j, 2) = surf(reshape(ufMeanArray{i} - uf_i_min, (self.nElFX + 1), (self.nElFY + 1)));
                    s(j, 2).LineStyle = 'none';
                    s(j, 2).FaceColor = 'b';
                    s(j, 3) = surf(reshape(ufMeanArray{i} - uf_i_min, (self.nElFX + 1), (self.nElFY + 1)) +...
                        sqrt(reshape(ufVarArray{i}, (self.nElFX + 1), (self.nElFY + 1))));
                    s(j, 3).LineStyle = 'none';
                    s(j, 3).FaceColor = [.85 .85 .85];
                    s(j, 4) = surf(reshape(ufMeanArray{i} - uf_i_min, (self.nElFX + 1), (self.nElFY + 1)) -...
                        sqrt(reshape(ufVarArray{i}, (self.nElFX + 1), (self.nElFY + 1))));
                    s(j, 4).LineStyle = 'none';
                    s(j, 4).FaceColor = [.85 .85 .85];
                    ax.FontSize = 30;
                    ax.View = [-125, -5];
                    im(j) = imagesc(reshape(cond(:, i), self.fineMesh.nElX, self.fineMesh.nElY));
                    xticks([0 64 128 192 256]);
                    yticks([0 64 128 192 256]);
%                     zticks(100:100:800)
                    xticklabels({});
                    yticklabels({});
%                     zticklabels({});
                    axis tight;
                    axis square;
                    box on;
                    zlim([min_uf max_uf]);
                    caxis([min_uf max_uf]);
                    j = j + 1;
                end
%                 print(f, './predictions', '-dpng', '-r300')
            end
        end

        function [predMean, predSqMean] = randThetaPredict(self, mode, nSamplesTheta, boundaryConditions)
            %Function to predict finescale output from generative model
            
            if(nargin < 2)
                mode = 'test';
            end
            
            if(nargin > 3)
                %Predict on different boundary conditions
                self.setBoundaryConditions(boundaryConditions);
                
                %Set up coarseMesh;
                %must be done after boundary conditions are set up
                self.genCoarseDomain;
            end
                
            %Load test file
            if strcmp(mode, 'self')
                %Self-prediction on training set
                uf = self.trainingDataMatfile.uf(:, self.nStart:(self.nStart + self.nTrain - 1));
                tempVars = whos(self.trainingDataMatfile);
                bcVar = ismember('bc', {tempVars.name});
            else
                uf = self.testDataMatfile.uf(:, self.testSamples);
                tempVars = whos(self.testDataMatfile);
                bcVar = ismember('bc', {tempVars.name});
            end
            self.loadTrainedParams;
            if strcmp(mode, 'self')
                self.computeDesignMatrix('train');
                designMatrixPred = self.designMatrix;
            else
                self.computeDesignMatrix('test');
                designMatrixPred = self.testDesignMatrix;
            end
            
            %% Sample from p_c
            disp('Sampling from p_c...')
            if strcmp(mode, 'self')
                nTest = self.nTrain;
            else
                nTest = numel(self.testSamples);
            end
            
            %short hand notation/ avoiding broadcast overhead
            nElc = self.coarseMesh.nEl;
            nSamplesLambda_c = self.nSamples_p_c;
            
            XsamplesTemp = zeros(nElc, nSamplesLambda_c);
            
            LambdaSamples{1} = zeros(nElc, nSamplesLambda_c);
            LambdaSamples = repmat(LambdaSamples, nTest, nSamplesTheta);
            self.meanEffCond = zeros(nElc, nTest);
            
            if self.useLaplaceApproximation
                precisionTheta = self.laplaceApproximation;
                SigmaTheta = inv(precisionTheta) + eps*eye(numel(...
                    self.theta_c.theta));
            else
                SigmaTheta = zeros(numel(self.theta_c.theta));
                nSamplesTheta = 1;
            end
            
            for n = 1:nTest
                %Samples from p_c
                %First sample theta from Laplace approx, then sample X
                theta = mvnrnd(self.theta_c.theta', SigmaTheta, nSamplesTheta)';
                for t = 1:nSamplesTheta
                    for k = 1:nSamplesLambda_c
                        Sigma2 = self.theta_c.Sigma;
                        XsamplesTemp(:, k) =...
                            mvnrnd((designMatrixPred{n}*theta(:, t))', Sigma2)';
                    end
                %Conductivities
                LambdaSamples{n, t} = conductivityBackTransform(XsamplesTemp,...
                    self.conductivityTransformation);
                end
            end
            clear Xsamples; %save memory
            disp('done')
            
            %% Run coarse model and sample from p_cf
            disp('Solving coarse model and sample from p_cf...')
            ufMean{1} = zeros(self.fineMesh.nNodes, nTest);
            ufMean = repmat(ufMean, nSamplesTheta, 1);
            ufVar = ufMean;
            uf_sq_mean = ufMean;
            
            if(bcVar)
                %Set coarse domain for data with different boundary conditions
                nX = self.coarseMesh.nElX;
                nY = self.coarseMesh.nElY;
                if strcmp(mode, 'self')
                    bc = self.trainingDataMatfile.bc;
                else
                    bc = self.testDataMatfile.bc;
                end
                for n = 1:nTest
                    if strcmp(mode, 'self')
                        i = self.trainingSamples(n);
                    else
                        i = self.testSamples(n);
                    end
                    bcT = @(x) bc{i}(1) + bc{i}(2)*x(1) +...
                        bc{i}(3)*x(2) + bc{i}(4)*x(1)*x(2);
                    bcQ{1} = @(x) -(bc{i}(3) + bc{i}(4)*x);      %lower bound
                    bcQ{2} = @(y) (bc{i}(2) + bc{i}(4)*y);       %right bound
                    bcQ{3} = @(x) (bc{i}(3) + bc{i}(4)*x);       %upper bound
                    bcQ{4} = @(y) -(bc{i}(2) + bc{i}(4)*y);      %left bound
                    cd(n) = self.coarseMesh;
                    cd(n) = cd(n).setBoundaries([2:(2*nX + 2*nY)], bcT, bcQ);
                end
            else
                cd = self.coarseMesh;
            end
            t_cf = self.theta_cf;
            natNodes = true(self.fineMesh.nNodes, 1);
            natNodes(self.fineMesh.essentialNodes) = false;
            addpath('./heatFEM');
            parfor t = 1:nSamplesTheta
                for n = 1:nTest
                    if bcVar
                        coarseDomain = cd(n);
                    else
                        coarseDomain = cd;
                    end
                    for i = 1:nSamplesLambda_c
                        D = zeros(2, 2, coarseDomain.nEl);
                        for e = 1:coarseDomain.nEl
                            D(:, :, e) = LambdaSamples{n, t}(e, i)*eye(2);
                        end
                        
                        FEMout= heat2d(coarseDomain, LambdaSamples{n, t}(:, i));
                        uctemp = FEMout.uff';
                        
                        %sample from p_cf
                        mu_cf = t_cf.mu + t_cf.W*uctemp(:);
                        %only for diagonal S!!
                        %Sequentially compute mean and <uf^2> to save memory
                        %U_f-integration can be done analyt.
                        ufMean{t}(:, n) =...
                            ((i - 1)/i)*ufMean{t}(:, n) + (1/i)*mu_cf;
                        uf_sq_mean{t}(:, n) =...
                            ((i - 1)/i)*uf_sq_mean{t}(:, n) + (1/i)*mu_cf.^2;
                    end
                    
                    uf_sq_mean{t}(:, n) = uf_sq_mean{t}(:, n) + t_cf.S;
                end
                predMean{t} = mean(ufMean{t}, 2);
                predSqMean{t} = mean(uf_sq_mean{t}, 2);
                ufMean{t} = []; %save memory
                uf_sq_mean{t} = [];
            end
            

        end
        
        
        
        %% Design matrix functions
        function getCoarseElement(self)
            debug = false;
            self.coarseCellMask = {zeros(self.fineMesh.nEl, 1, 'logical')};
            self.coarseCellMask = repmat(self.coarseCellMask, self.coarseMesh.nEl, 1);
            e = 1;  %element number
            for row_fine = 1:self.fineMesh.nElY
                %coordinate of lower boundary of fine element
                y_coord = self.fineMesh.cum_lElY(row_fine);
                row_coarse = sum(y_coord >= self.coarseMesh.cum_lElY);
                for col_fine = 1:self.fineMesh.nElX
                    %coordinate of left boundary of fine element
                    x_coord = self.fineMesh.cum_lElX(col_fine);
                    col_coarse = sum(x_coord >= self.coarseMesh.cum_lElX);
                    self.coarseCellMask{(row_coarse - 1)*self.coarseMesh.nElX + col_coarse}(e) = true;
                    e = e + 1;
                end
            end
            
            if debug
                figure
%                 imagesc(E)
                pause
            end
        end

        function [lambdak, xk] = get_coarseElementConductivities(self, conductivity)
            %Cuts out conductivity fields from macro-cells
			addpath('./rom');    
            
            nData = size(conductivity, 2);
            %Mapping from fine cell index to coarse cell index
            self.getCoarseElement;
                                    
            %prealloc
            lambdak = cell(nData, self.coarseMesh.nEl);
            if(nargout > 1)
                xk = lambdak;
            end
            
            for s = 1:nData
                %inputs belonging to same coarse element are in the same column of xk. They are ordered in x-direction.
                %Get conductivity fields in coarse cell windows. Might be wrong for non-square fine scale domains
                conductivityMat = reshape(conductivity(:, s), self.fineMesh.nElX, self.fineMesh.nElY);
                for e = 1:self.coarseMesh.nEl
%                     lambdakTemp = conductivityMat.*self.coarseCellMask{e};
                    lambdakTemp = conductivityMat(self.coarseCellMask{e});
                    %Cut elements from matrix that do not belong to coarse cell
                    lambdakTemp(~any(lambdakTemp, 2), :) = [];
                    lambdakTemp(:, ~any(lambdakTemp, 1)) = [];
                    lambdak{s, e} = lambdakTemp;
                    if(nargout > 1)
                        xk{s, e} = conductivityTransform(lambdak{s, e}, self.conductivityTransformation);
                    end
                end
            end
        end
        
        function [Phi, d_Phi] = computeDesignMatrix(self, mode, recompute, xi)
            %Actual computation of design matrix set recompute to true if design matrices have to be
            %recomputed during optimization (parametric features)
            if(self.loadDesignMatrix && ~recompute && mode == "train")
                load(strcat('./persistentData/', mode, 'DesignMatrix.mat'));
                self.designMatrix = designMatrix;
            else
                if mode ~= "inverse"
                    fprintf('Computing design matrices...\n')
                    tic
                end
                
                if mode == "train"
                    dataFile = self.trainingDataMatfile;
                    dataSamples = self.trainingSamples;
                    d_Phi = [];
                    nData = numel(dataSamples);
                    self.getCoarseElement;
                elseif mode == "test"
                    dataFile = self.testDataMatfile;
                    dataSamples = self.testSamples;
                    d_Phi = [];
                    nData = numel(dataSamples);
                    self.getCoarseElement;
                elseif mode == "inverse"
                    %do nothing. Microstructures are not given by data
                    nData = 1;
                else
                    error('Compute design matrices for train or test data?')
                end

                %Generate discretized conductivity field from basis functions
                if nargin < 4
                    xi = dataFile.xi(dataSamples, :);
                end
                conductivity = self.getDiscretizedConductivityField(xi);

                %set feature function handles
                if isempty(self.featureFunctions) && isempty(self.globalFeatureFunctions)
                    [phi, phiGlobal] = self.setFeatureFunctions;
                else
                    phi = self.featureFunctions;
                    phiGlobal = self.globalFeatureFunctions;
                    self.writeFeatureFunctionList(phi, phiGlobal);
                end

                nFeatureFunctions = size(self.featureFunctions, 2);
                nGlobalFeatureFunctions = size(self.globalFeatureFunctions, 2);
                nFeaturesTot = nFeatureFunctions + nGlobalFeatureFunctions;

                Phi{1} = zeros(self.coarseMesh.nEl, nFeaturesTot);
%                 [lambdak, xk] = self.get_coarseElementConductivities(conductivity);
                if mode == "inverse"
                    %ONLY WORKS IF MACRO CELLS LAMBDAK HAVE IDENTICAL SIZES
%                     d_Phi = zeros(self.coarseMesh.nEl*(nFeaturesTot), numel(lambdak{1, 1}));
                    d_Phi = spalloc(self.coarseMesh.nEl*nFeaturesTot, self.fineMesh.nEl, numel(self.coarseCellMask{1}));
                end
                Phi = repmat(Phi, nData, 1);

                %for cheap features, serial evaluation might be more efficient
                for s = 1:nData
                    %inputs belonging to same coarse element are in 
                    %the same column of xk. They are ordered in x-direction.
                    
                    %construct conductivity design matrix
                    for k = 1:self.coarseMesh.nEl
                        %local features
                        for j = 1:nFeatureFunctions
                            %only take pixels of corresponding macro-cell as input for features
                            if mode == "inverse"
                                %FIXME: needs to be generalized for global feature functions
%                                 [Phi{s}(k, j), d_Phi(((k - 1)*nFeaturesTot) + j, :)] = phi{k, j}(lambdak{s, k});
                                [Phi{s}(k, j), d_Phi(((k - 1)*nFeaturesTot) + j, :)] =...
                                    phi{j}(conductivity(:, s), self.coarseCellMask{k});
                            else
                                Phi{s}(k, j) = phi{j}(conductivity(:, s), self.coarseCellMask{k});
                            end
                        end
                        %global features
                        for j = 1:nGlobalFeatureFunctions
                            %Take whole microstructure as input for feature function can be wrong
                            %for non-sq. fine domains
                            conductivityMat = reshape(conductivity(:, s), self.fineMesh.nElX, self.fineMesh.nElY);
                            Phi{s}(k, nFeatureFunctions + j) = phiGlobal{k, j}(conductivityMat);
                        end
                    end
                end
                
                self.designMatrixTest(Phi);
                
                if mode == "train"
                    self.designMatrix = Phi;
                elseif mode == "test"
                    self.testDesignMatrix = Phi;
                elseif mode == "inverse"
                    %pass
                else
                    error('Wrong design matrix computation model')
                end
                                
                %Design matrix is always stored in its original form. Local modes are applied after loading.
                if strcmp(mode, 'train')
                    designMatrix = self.designMatrix;
                    self.originalDesignMatrix = self.designMatrix;
                    save(strcat('./persistentData/', mode, 'DesignMatrix.mat'), 'designMatrix')
                end
            end
            
            self.applyLocalityModeToDesignMatrix(Phi, mode);
            if mode ~= "inverse"
                self.computeSumPhiTPhi;
                fprintf('design matrices computed.\n')
                fprintf("Design matrix computation time: %.2f\n", toc)
            end
        end
        
        function applyLocalityModeToDesignMatrix(self, Phi, mode)
            %Use specific nonlocality mode
            if strcmp(self.mode, 'useNeighbor')
                %use feature function information from nearest neighbors
                self.includeNearestNeighborFeatures(Phi, mode);
            elseif strcmp(self.mode, 'useLocalNeighbor')
                self.includeLocalNearestNeighborFeatures(Phi, mode);
            elseif strcmp(self.mode, 'useLocalDiagNeighbor')
                self.includeLocalDiagNeighborFeatures(Phi, mode);
            elseif strcmp(self.mode, 'useDiagNeighbor')
                %use feature function info from nearest and diagonal neighbors
                self.includeDiagNeighborFeatures(Phi, mode);
            elseif strcmp(self.mode, 'useLocal')
                %Use separate parameters for every macro-cell
                self.localTheta_c(Phi, mode);
            else
                self.originalDesignMatrix = [];
            end
        end
        
        function designMatrixTest(self, Phi)
            %Check for real finite inputs
            for i = 1:numel(Phi)
                if(~all(all(all(isfinite(Phi{i})))))
                    [coarseElement, featureFunction] = ind2sub(size(Phi{i}), find(~isfinite(Phi{i})));
                    warning("Non-finite design matrix at data point %u, element %u, feature %u" +...
                        "Setting non-finite component to 0.", i, coarseElement, featureFunction)
                    Phi{i}(~isfinite(Phi{i})) = 0;
                elseif(~all(all(all(isreal(Phi{i})))))
                    [coarseElement, featureFunction] = ind2sub(size(Phi{i}), find(imag(Phi{i})));
                    warning("Complex design matrix at data point %u, element %u, feature %u" +...
                        "Setting imaginary part to 0.", i, coarseElement, featureFunction)
                    Phi{i} = real(Phi{i});
                end
            end
        end
        
        function applyDesignMatrixNormalization(self, mode)
            %Normalize design matrices
            if self.featureScaling == "standardize"
                self.standardizeDesignMatrix(mode, PhiCell);
                frintf("Feature functions standardized.\n")
                warning("Feature normalization not possible for inverse problems.")
            elseif self.featureScaling == "rescale"
                self.rescaleDesignMatrix(mode, PhiCell);
                frintf("Feature functions rescaled.\n")
                warning("Feature normalization not possible for inverse problems.")
            elseif self.featureScaling == "normalize"
                self.normalizeDesignMatrix(mode, PhiCell);
                frintf("Feature Functions normalized.\n")
                warning("Feature normalization not possible for inverse problems.")
            else
                frintf("No feature normalization used.\n")
            end
        end
        
        function writeFeatureFunctionList(self, phi, phiGlobal)
            for j = 1:size(phi, 2)
                if(j == 1)
                    dlmwrite('./data/features', func2str(phi{1, j}), 'delimiter', '');
                else
                    dlmwrite('./data/features', func2str(phi{1, j}), 'delimiter', '', '-append');
                end
            end
            for j = 1:size(phiGlobal, 2)
                dlmwrite('./data/features', func2str(phiGlobal{1, j}), 'delimiter', '', '-append');
            end
        end
        
        function computeFeatureFunctionMean(self)
            %Must be executed BEFORE useLocal etc.
            self.featureFunctionMean = 0;
            for n = 1:numel(self.designMatrix)
                self.featureFunctionMean = self.featureFunctionMean + mean(self.designMatrix{n}, 1);
            end
            self.featureFunctionMean = self.featureFunctionMean/numel(self.designMatrix);
        end

        function computeFeatureFunctionSqMean(self)
            featureFunctionSqSum = 0;
            for i = 1:numel(self.designMatrix)
                featureFunctionSqSum =...
                    featureFunctionSqSum + sum(self.designMatrix{i}.^2, 1);
            end
            self.featureFunctionSqMean = featureFunctionSqSum/...
                (numel(self.designMatrix)*size(self.designMatrix{1}, 1));
        end

        function standardizeDesignMatrix(self, mode, designMatrix)
            %Standardize covariates to have 0 mean and unit variance
            disp('Standardize design matrix...')
            %Compute std
            if strcmp(mode, 'test')
                featureFunctionStd = sqrt(self.featureFunctionSqMean -...
                    self.featureFunctionMean.^2);
            else
                self.computeFeatureFunctionMean;
                self.computeFeatureFunctionSqMean;
                featureFunctionStd = sqrt(self.featureFunctionSqMean -...
                    self.featureFunctionMean.^2);
                if(any(~isreal(featureFunctionStd)))
                    warning('Imaginary standard deviation. Setting it to 0.')
                    featureFunctionStd = real(featureFunctionStd);
                end
            end
            
            %Check if there is a constant feature function
            for i = 1:length(featureFunctionStd)
                if(featureFunctionStd == 0)
                    i
                    featureFunctionStd(i) = self.featureFunctionMean(i);
                    warning(strcat('At least one feature always has the', ...
                        'same output. It will be rescaled to one.'))
                    break;
                end
            end
            
            %centralize
            for i = 1:numel(designMatrix)
                designMatrix{i} = designMatrix{i} - self.featureFunctionMean;
            end
            
            %normalize
            for i = 1:numel(designMatrix)
                designMatrix{i} = designMatrix{i}./featureFunctionStd;
            end
            
            %Check for finiteness
            for i = 1:numel(designMatrix)
                if(~all(all(all(isfinite(designMatrix{i})))))
                    warning('Non-fin. des. mat. Setting non-fin. comp. to 0.')
                    designMatrix{i}(~isfinite(designMatrix{i})) = 0;
                elseif(~all(all(all(isreal(designMatrix{i})))))
                    warning('Complex feature function output:')
                    dataPoint = i
                    [coarseElement, featureFunction] = ...
                        ind2sub(size(designMatrix{i}),...
                        find(imag(designMatrix{i})))
                    disp('Ignoring imaginary part...')
                    designMatrix{i} = real(designMatrix{i});
                end
            end
            if strcmp(mode, 'train')
                self.designMatrix = designMatrix;
            elseif strcmp(mode, 'test')
                self.testDesignMatrix = designMatrix;
            else
                error('Unknown mode')
            end
            self.saveNormalization('standardization');
            disp('done')
        end
        
        function normalizeDesignMatrix(self, mode, designMatrix)
            %Standardize covariates to unit variance
            disp('Normalize design matrix...')
            %Compute std
            if strcmp(mode, 'test')
                featureFunctionStd = sqrt(self.featureFunctionSqMean -...
                    self.featureFunctionMean.^2);
            else
                self.computeFeatureFunctionMean;
                self.computeFeatureFunctionSqMean;
                featureFunctionStd = sqrt(self.featureFunctionSqMean -...
                    self.featureFunctionMean.^2);
                if(any(~isreal(featureFunctionStd)))
                    warning('Imaginary standard deviation. Setting it to 0.')
                    featureFunctionStd = real(featureFunctionStd);
                end
            end
            
            %Check if there is a constant feature function
            for i = 1:length(featureFunctionStd)
                if(featureFunctionStd(i) == 0)
                    i
                    featureFunctionStd(i) = self.featureFunctionMean(i);
                    warning(strcat('At least one feature always has the',...
                        'same output. It will be rescaled to one.'))
                    break;
                end
            end
            
            %normalize
            for i = 1:numel(designMatrix)
                designMatrix{i} = designMatrix{i}./featureFunctionStd;
            end
            
            %Check for finiteness
            for i = 1:numel(designMatrix)
                if(~all(all(all(isfinite(designMatrix{i})))))
                    warning(strcat('Non-finite design matrix.',...
                        'Setting non-finite component to 0.'))
                    designMatrix{i}(~isfinite(designMatrix{i})) = 0;
                elseif(~all(all(all(isreal(designMatrix{i})))))
                    warning('Complex feature function output:')
                    dataPoint = i
                    [coarseElement, featureFunction] =...
                        ind2sub(size(designMatrix{i}),...
                        find(imag(designMatrix{i})))
                    disp('Ignoring imaginary part...')
                    designMatrix{i} = real(designMatrix{i});
                end
            end
            if strcmp(mode, 'train')
                self.designMatrix = designMatrix;
            elseif strcmp(mode, 'test')
                self.testDesignMatrix = designMatrix;
            else
                error('wrong mode')
            end
            self.saveNormalization('standardization');
            disp('done')
        end
        
        function computeFeatureFunctionMinMax(self)
            %Computes min/max of feature function outputs over training data, 
            %separately for every macro cell
            self.featureFunctionMin = self.designMatrix{1};
            self.featureFunctionMax = self.designMatrix{1};
            for n = 1:numel(self.designMatrix)
                self.featureFunctionMin(self.featureFunctionMin > self.designMatrix{n}) =...
                    self.designMatrix{n}(self.featureFunctionMin > self.designMatrix{n});
                self.featureFunctionMax(self.featureFunctionMax < self.designMatrix{n}) =...
                    self.designMatrix{n}(self.featureFunctionMax < self.designMatrix{n});
            end
        end
        
        function rescaleDesignMatrix(self, mode, designMatrix)
            %Rescale design matrix s.t. outputs are between 0 and 1
            disp('Rescale design matrix...')
            if strcmp(mode, 'test')
                featFuncDiff= self.featureFunctionMax - self.featureFunctionMin;
                %to avoid irregularities due to rescaling 
                %(if every macro cell has the same feature function output)
                self.featureFunctionMin(featFuncDiff == 0) = 0;
                featFuncDiff(featFuncDiff == 0) = 1;
                for n = 1:numel(designMatrix)
                    designMatrix{n} = (designMatrix{n} - ...
                        self.featureFunctionMin)./(featFuncDiff);
                end
            else
                self.computeFeatureFunctionMinMax;
                featFuncDiff =self.featureFunctionMax - self.featureFunctionMin;
                %to avoid irregularities due to rescaling 
                %(if every macro cell has the same feature function output)
                self.featureFunctionMin(featFuncDiff == 0) = 0;
                featFuncDiff(featFuncDiff == 0) = 1;
                for n = 1:numel(designMatrix)
                    designMatrix{n} = (designMatrix{n} -...
                        self.featureFunctionMin)./(featFuncDiff);
                end
            end
            %Check for finiteness
            for n = 1:numel(designMatrix)
                if(~all(all(all(isfinite(designMatrix{n})))))
                    warning('Non-finite des. mat. Setting non-fin. comp. to 0.')
                    designMatrix{n}(~isfinite(designMatrix{n})) = 0;
                    dataPoint = n
                    [coarseElement, featureFunction] = ...
                        ind2sub(size(designMatrix{n}),...
                        find(~isfinite(designMatrix{n})))
                elseif(~all(all(all(isreal(designMatrix{n})))))
                    warning('Complex feature function output:')
                    dataPoint = n
                    [coarseElement, featureFunction] =...
                        ind2sub(size(designMatrix{n}),...
                        find(imag(designMatrix{n})))
                    disp('Ignoring imaginary part...')
                    designMatrix{n} = real(designMatrix{n});
                end
            end
            if strcmp(mode, 'train')
                self.designMatrix = designMatrix;
            elseif strcmp(mode, 'test')
                self.testDesignMatrix = designMatrix;
            else
                error('wrong mode')
            end
            self.saveNormalization('rescaling');
            disp('done')
        end
        
        function saveNormalization(self, type)
            disp('Saving design matrix normalization...')
            if(isempty(self.featureFunctionMean))
                self.computeFeatureFunctionMean;
            end
            if(isempty(self.featureFunctionSqMean))
                self.computeFeatureFunctionSqMean;
            end
            if ~exist('./data')
                mkdir('./data');
            end
            if strcmp(type, 'standardization')
                featureFunctionMean = self.featureFunctionMean;
                featureFunctionSqMean = self.featureFunctionSqMean;
                save('./data/featureFunctionMean', 'featureFunctionMean', '-ascii');
                save('./data/featureFunctionSqMean', 'featureFunctionSqMean', '-ascii');
            elseif strcmp(type, 'rescaling')
                featureFunctionMin = self.featureFunctionMin;
                featureFunctionMax = self.featureFunctionMax;
                save('./data/featureFunctionMin', 'featureFunctionMin', '-ascii');
                save('./data/featureFunctionMax', 'featureFunctionMax', '-ascii');
            else
                error('Which type of data normalization?')
            end
        end
        
        function includeNearestNeighborFeatures(self, designMatrix, mode)
            %Includes feature function information of neighboring cells
            %Can only be executed after standardization/rescaling!
            %nc/nf: coarse/fine elements in x/y direction
            disp('Including nearest neighbor feature function information...')
            nFeatureFunctionsTotal = size(designMatrix{1}, 2);
            PhiCell{1} = zeros(self.coarseMesh.nEl,...
                5*nFeatureFunctionsTotal);
            nData = numel(designMatrix);
            PhiCell = repmat(PhiCell, nData, 1);
            
            for n = 1:nData
                %The first columns contain feature function information of
                %the original cell
                PhiCell{n}(:, 1:nFeatureFunctionsTotal) = designMatrix{n};
                
                %Only assign nonzero values to design matrix for neighboring
                %elements if neighbor in respective direction exists
                for k = 1:self.coarseMesh.nEl
                    if(mod(k, self.coarseMesh.nElX) ~= 0)
                        %right neighbor of coarse element exists
                        PhiCell{n}(k, (nFeatureFunctionsTotal + 1):...
                            (2*nFeatureFunctionsTotal)) =...
                            designMatrix{n}(k + 1, :);
                    end
                    
                    if(k <= self.coarseMesh.nElX*...
                            (self.coarseMesh.nElY - 1))
                        %upper neighbor of coarse element exists
                        PhiCell{n}(k, (2*nFeatureFunctionsTotal + 1):...
                            (3*nFeatureFunctionsTotal)) =...
                            designMatrix{n}(k + self.coarseMesh.nElX, :);
                    end
                    
                    if(mod(k - 1, self.coarseMesh.nElX) ~= 0)
                        %left neighbor of coarse element exists
                        PhiCell{n}(k, (3*nFeatureFunctionsTotal + 1):...
                            (4*nFeatureFunctionsTotal)) =...
                            designMatrix{n}(k - 1, :);
                    end
                    
                    if(k > self.coarseMesh.nElX)
                        %lower neighbor of coarse element exists
                        PhiCell{n}(k, (4*nFeatureFunctionsTotal + 1):...
                            (5*nFeatureFunctionsTotal)) =...
                            designMatrix{n}(k - self.coarseMesh.nElX, :);
                    end
                end
            end
            if strcmp(mode, 'train')
                self.designMatrix = PhiCell;
            elseif strcmp(mode, 'test')
                self.testDesignMatrix = PhiCell;
            else
                error('wrong mode');
            end
            disp('done')
        end%includeNearestNeighborFeatures
        
        function includeLocalNearestNeighborFeatures(self, designMatrix, mode)
            %Includes feature function information of neighboring cells
            %Can only be executed after standardization/rescaling!
            %nc/nf: coarse/fine elements in x/y direction
            disp('Incl. n.n. feature func. information sep. for each cell...')
            nFeatureFunctionsTotal = size(designMatrix{1}, 2);
            PhiCell{1} =...
                zeros(self.coarseMesh.nEl, 5*nFeatureFunctionsTotal);
            nData = numel(designMatrix);
            PhiCell = repmat(PhiCell, nData, 1);
            
            for n = 1:nData
                %Only assign nonzero values to design matrix for neighboring
                %elements if neighbor in respective direction exists
                k = 0;
                for i = 1:self.coarseMesh.nEl
                    PhiCell{n}(i, (k*nFeatureFunctionsTotal + 1):((k + 1)*...
                        nFeatureFunctionsTotal)) = designMatrix{n}(i, :);
                    self.neighborDictionary((k*nFeatureFunctionsTotal + 1):...
                        ((k + 1)*nFeatureFunctionsTotal), 1) = ...
                        (1:nFeatureFunctionsTotal)'; %feature index
                    self.neighborDictionary((k*nFeatureFunctionsTotal + 1):...
                        ((k + 1)*nFeatureFunctionsTotal), 2) = ...
                        i; %coarse element index
                    self.neighborDictionary((k*nFeatureFunctionsTotal + 1):...
                        ((k + 1)*nFeatureFunctionsTotal), 3) = ...
                        0; %center element
                    k = k + 1;
                    if(mod(i, self.coarseMesh.nElX) ~= 0)
                        %right neighbor of coarse element exists
                        PhiCell{n}(i, (k*nFeatureFunctionsTotal + 1):...
                            ((k + 1)*nFeatureFunctionsTotal)) =...
                            designMatrix{n}(i + 1, :);
                        self.neighborDictionary((k*nFeatureFunctionsTotal+1):...
                            ((k + 1)*nFeatureFunctionsTotal), 1) = ...
                            (1:nFeatureFunctionsTotal)'; %feature index
                        self.neighborDictionary((k*nFeatureFunctionsTotal+1):...
                            ((k + 1)*nFeatureFunctionsTotal), 2) = ...
                            i; %coarse element index
                        self.neighborDictionary((k*nFeatureFunctionsTotal+1):...
                            ((k + 1)*nFeatureFunctionsTotal), 3) = ...
                            1; %right neighbor
                        k = k + 1;
                    end
                    
                    if(i <= self.coarseMesh.nElX*...
                            (self.coarseMesh.nElY - 1))
                        %upper neighbor of coarse element exists
                        PhiCell{n}(i, (k*nFeatureFunctionsTotal + 1):...
                            ((k + 1)*nFeatureFunctionsTotal)) =...
                            designMatrix{n}(i + self.coarseMesh.nElX, :);
                        self.neighborDictionary((k*nFeatureFunctionsTotal+1):...
                            ((k + 1)*nFeatureFunctionsTotal), 1) = ...
                            (1:nFeatureFunctionsTotal)'; %feature index
                        self.neighborDictionary((k*nFeatureFunctionsTotal+1):...
                            ((k + 1)*nFeatureFunctionsTotal), 2) = ...
                            i; %coarse element index
                        self.neighborDictionary((k*nFeatureFunctionsTotal+1):...
                            ((k + 1)*nFeatureFunctionsTotal), 3) = ...
                            2; %upper neighbor
                        k = k + 1;
                    end
                    
                    if(mod(i - 1, self.coarseMesh.nElX) ~= 0)
                        %left neighbor of coarse element exists
                        PhiCell{n}(i, (k*nFeatureFunctionsTotal + 1):...
                            ((k + 1)*nFeatureFunctionsTotal)) =...
                            designMatrix{n}(i - 1, :);
                        self.neighborDictionary((k*nFeatureFunctionsTotal+1):...
                            ((k + 1)*nFeatureFunctionsTotal), 1) = ...
                            (1:nFeatureFunctionsTotal)'; %feature index
                        self.neighborDictionary((k*nFeatureFunctionsTotal+1):...
                            ((k + 1)*nFeatureFunctionsTotal), 2) = ...
                            i; %coarse element index
                        self.neighborDictionary((k*nFeatureFunctionsTotal+1):...
                            ((k + 1)*nFeatureFunctionsTotal), 3) = ...
                            3; %left neighbor
                        k = k + 1;
                    end
                    
                    if(i > self.coarseMesh.nElX)
                        %lower neighbor of coarse element exists
                        PhiCell{n}(i, (k*nFeatureFunctionsTotal + 1):((k+ 1)*...
                            nFeatureFunctionsTotal)) =...
                            designMatrix{n}(i - self.coarseMesh.nElX, :);
                        self.neighborDictionary((k*nFeatureFunctionsTotal+1):...
                            ((k + 1)*nFeatureFunctionsTotal), 1) = ...
                            (1:nFeatureFunctionsTotal)'; %feature index
                        self.neighborDictionary((k*nFeatureFunctionsTotal+1):...
                            ((k + 1)*nFeatureFunctionsTotal), 2) = ...
                            i; %coarse element index
                        self.neighborDictionary((k*nFeatureFunctionsTotal+1):...
                            ((k + 1)*nFeatureFunctionsTotal), 3) = ...
                            4; %lower neighbor
                        k = k + 1;
                    end
                end
            end
            if strcmp(mode, 'train')
                self.designMatrix = PhiCell;
            elseif strcmp(mode, 'test')
                self.testDesignMatrix = PhiCell;
            else
                error('wrong mode');
            end
            disp('done')
        end%includeLocalNearestNeighborFeatures
        
        function includeDiagNeighborFeatures(self, designMatrix, mode)
            %includes feature function information of all other cells
            %Can only be executed after standardization/rescaling!
            %nc/nf: coarse/fine elements in x/y direction
            disp('Incl. nearest and diagonal neighbor feature function info...')
            nFeatureFunctionsTotal = size(designMatrix{1}, 2);
            PhiCell{1} = zeros(self.coarseMesh.nEl,...
                9*nFeatureFunctionsTotal);
            nData = numel(designMatrix);
            PhiCell = repmat(PhiCell, nData, 1);
            
            for n = 1:nData
                %The first columns contain feature function information 
                %of the original cell
                PhiCell{n}(:, 1:nFeatureFunctionsTotal) = designMatrix{n};
                
                %Only assign nonzero values to design matrix for neighboring 
                %elements if neighbor in respective direction exists
                for i = 1:self.coarseMesh.nEl
                    if(mod(i, self.coarseMesh.nElX) ~= 0)
                        %right neighbor of coarse element exists
                        PhiCell{n}(i, (nFeatureFunctionsTotal + 1):...
                            (2*nFeatureFunctionsTotal)) = ...
                            designMatrix{n}(i + 1, :);
                        if(i <= self.coarseMesh.nElX*...
                                (self.coarseMesh.nElY - 1))
                            %upper right neighbor of coarse element exists
                            PhiCell{n}(i, (2*nFeatureFunctionsTotal + 1):...
                                (3*nFeatureFunctionsTotal)) =...
                                designMatrix{n}(i +...
                                self.coarseMesh.nElX + 1, :);
                        end
                    end
                    
                    if(i <= self.coarseMesh.nElX*...
                            (self.coarseMesh.nElY - 1))
                        %upper neighbor of coarse element exists
                        PhiCell{n}(i, (3*nFeatureFunctionsTotal + 1):...
                            (4*nFeatureFunctionsTotal)) =...
                            designMatrix{n}(i + self.coarseMesh.nElX, :);
                        if(mod(i - 1, self.coarseMesh.nElX) ~= 0)
                            %upper left neighbor exists
                            PhiCell{n}(i, (4*nFeatureFunctionsTotal + 1):...
                                (5*nFeatureFunctionsTotal)) =...
                            designMatrix{n}(i + ...
                            self.coarseMesh.nElX - 1, :);
                        end
                    end
                    
                    if(mod(i - 1, self.coarseMesh.nElX) ~= 0)
                        %left neighbor of coarse element exists
                        PhiCell{n}(i, (5*nFeatureFunctionsTotal + 1):...
                            (6*nFeatureFunctionsTotal)) =...
                            designMatrix{n}(i - 1, :);
                        if(i > self.coarseMesh.nElX)
                            %lower left neighbor exists
                            PhiCell{n}(i, (6*nFeatureFunctionsTotal + 1):...
                                (7*nFeatureFunctionsTotal)) =...
                            designMatrix{n}(i - ...
                            self.coarseMesh.nElX - 1, :);
                        end
                    end
                    
                    if(i > self.coarseMesh.nElX)
                        %lower neighbor of coarse element exists
                        PhiCell{n}(i, (7*nFeatureFunctionsTotal + 1):...
                            (8*nFeatureFunctionsTotal)) =...
                            designMatrix{n}(i - self.coarseMesh.nElX, :);
                        if(mod(i, self.coarseMesh.nElX) ~= 0)
                            %lower right neighbor exists
                            PhiCell{n}(i, (8*nFeatureFunctionsTotal + 1):...
                                (9*nFeatureFunctionsTotal)) =...
                            designMatrix{n}...
                            (i - self.coarseMesh.nElX + 1, :);
                        end
                    end
                end
            end
            if strcmp(mode, 'train')
                self.designMatrix = PhiCell;
            elseif strcmp(mode, 'test')
                self.testDesignMatrix = PhiCell;
            else
                error('wrong mode');
            end
            disp('done')
        end%includeDiagNeighborFeatures

        function includeLocalDiagNeighborFeatures(self, designMatrix, mode)
            %Includes feature function information of direct and diagonal 
            %neighboring cells can only be exec. after standardization/rescal.!
            %nc/nf: coarse/fine elements in x/y direction
            disp(strcat('Including nearest + diagonal neighbor feature', ...
                ' function information separately for each cell...'))
            nFeatureFunctionsTotal = size(designMatrix{1}, 2);
            nData = numel(designMatrix);
%             PhiCell = repmat(PhiCell, nTrain, 1);
            
            for n = 1:nData
                %Only assign nonzero values to design matrix for neighboring 
                %elements if neighbor in respective direction exists
                k = 0;
                for i = 1:self.coarseMesh.nEl
                    PhiCell{n}(i, (k*nFeatureFunctionsTotal + 1):...
                        ((k + 1)*nFeatureFunctionsTotal)) =...
                        designMatrix{n}(i, :);
                    self.neighborDictionary((k*nFeatureFunctionsTotal + 1):...
                        ((k + 1)*nFeatureFunctionsTotal), 1) = ...
                        (1:nFeatureFunctionsTotal)'; %feature index
                    self.neighborDictionary((k*nFeatureFunctionsTotal + 1):...
                        ((k + 1)*nFeatureFunctionsTotal), 2) = ...
                        i; %coarse element index
                    self.neighborDictionary((k*nFeatureFunctionsTotal + 1):...
                        ((k + 1)*nFeatureFunctionsTotal), 3) = ...
                        0; %center element
                    k = k + 1;
                    if(mod(i, self.coarseMesh.nElX) ~= 0)
                        %right neighbor of coarse element exists
                        PhiCell{n}(i, (k*nFeatureFunctionsTotal + 1):...
                            ((k + 1)*nFeatureFunctionsTotal)) =...
                            designMatrix{n}(i + 1, :);
                        self.neighborDictionary((k*nFeatureFunctionsTotal+1):...
                            ((k + 1)*nFeatureFunctionsTotal), 1) = ...
                            (1:nFeatureFunctionsTotal)'; %feature index
                        self.neighborDictionary((k*nFeatureFunctionsTotal+1):...
                            ((k + 1)*nFeatureFunctionsTotal), 2) = ...
                            i; %coarse element index
                        self.neighborDictionary((k*nFeatureFunctionsTotal+1):...
                            ((k + 1)*nFeatureFunctionsTotal), 3) = ...
                            1; %right neighbor
                        k = k + 1;
                        
                        if(i <= self.coarseMesh.nElX*...
                                (self.coarseMesh.nElY - 1))
                            %upper right neighbor of coarse element exists
                            PhiCell{n}(i, (k*nFeatureFunctionsTotal + 1):...
                                ((k + 1)*nFeatureFunctionsTotal)) =...
                                designMatrix{n}(i + ...
                                self.coarseMesh.nElX + 1, :);
                            self.neighborDictionary((k*nFeatureFunctionsTotal...
                                + 1):((k + 1)*nFeatureFunctionsTotal), 1) = ...
                                (1:nFeatureFunctionsTotal)'; %feature index
                            self.neighborDictionary((k*nFeatureFunctionsTotal...
                                + 1):((k + 1)*nFeatureFunctionsTotal), 2) = ...
                                i; %coarse element index
                            self.neighborDictionary((k*nFeatureFunctionsTotal...
                                + 1):((k + 1)*nFeatureFunctionsTotal), 3) = ...
                                2; % upper right neighbor
                            k = k + 1;
                        end
                        
                    end
                    
                    
                    if(i <= self.coarseMesh.nElX*...
                            (self.coarseMesh.nElY - 1))
                        %upper neighbor of coarse element exists
                        PhiCell{n}(i, (k*nFeatureFunctionsTotal + 1):...
                            ((k + 1)*nFeatureFunctionsTotal)) =...
                            designMatrix{n}(i + self.coarseMesh.nElX, :);
                        self.neighborDictionary((k*nFeatureFunctionsTotal...
                            + 1):((k + 1)*nFeatureFunctionsTotal), 1) = ...
                            (1:nFeatureFunctionsTotal)'; %feature index
                        self.neighborDictionary((k*nFeatureFunctionsTotal...
                            + 1):((k + 1)*nFeatureFunctionsTotal), 2) = ...
                            i; %coarse element index
                        self.neighborDictionary((k*nFeatureFunctionsTotal...
                            + 1):((k + 1)*nFeatureFunctionsTotal), 3) = ...
                            2; %upper neighbor
                        k = k + 1;
                        
                        if(mod(i - 1, self.coarseMesh.nElX) ~= 0)
                            %upper left neighbor of coarse element exists
                            PhiCell{n}(i, (k*nFeatureFunctionsTotal + 1):...
                                ((k + 1)*nFeatureFunctionsTotal)) =...
                                designMatrix{n}(i +...
                                self.coarseMesh.nElX - 1, :);
                            self.neighborDictionary((k*nFeatureFunctionsTotal...
                                + 1):((k + 1)*nFeatureFunctionsTotal), 1) = ...
                                (1:nFeatureFunctionsTotal)'; %feature index
                            self.neighborDictionary((k*nFeatureFunctionsTotal...
                                + 1):((k + 1)*nFeatureFunctionsTotal), 2) = ...
                                i; %coarse element index
                            self.neighborDictionary((k*nFeatureFunctionsTotal...
                                + 1):((k + 1)*nFeatureFunctionsTotal), 3) = ...
                                4; % upper left neighbor
                            k = k + 1;
                        end
                        
                    end
                    
                    
                    if(mod(i - 1, self.coarseMesh.nElX) ~= 0)
                        %left neighbor of coarse element exists
                        PhiCell{n}(i, (k*nFeatureFunctionsTotal + 1):...
                            ((k + 1)*nFeatureFunctionsTotal)) =...
                            designMatrix{n}(i - 1, :);
                        self.neighborDictionary((k*nFeatureFunctionsTotal...
                            + 1):((k + 1)*nFeatureFunctionsTotal), 1) = ...
                            (1:nFeatureFunctionsTotal)'; %feature index
                        self.neighborDictionary((k*nFeatureFunctionsTotal...
                            + 1):((k + 1)*nFeatureFunctionsTotal), 2) = ...
                            i; %coarse element index
                        self.neighborDictionary((k*nFeatureFunctionsTotal...
                            + 1):((k + 1)*nFeatureFunctionsTotal), 3) = ...
                            3; %left neighbor
                        k = k + 1;
                        
                        if(i > self.coarseMesh.nElX)
                            %lower left neighbor of coarse element exists
                            PhiCell{n}(i, (k*nFeatureFunctionsTotal + 1):...
                                ((k + 1)*nFeatureFunctionsTotal)) =...
                                designMatrix{n}(i - ...
                                self.coarseMesh.nElX - 1, :);
                            self.neighborDictionary((k*nFeatureFunctionsTotal...
                                + 1):((k + 1)*nFeatureFunctionsTotal), 1) = ...
                                (1:nFeatureFunctionsTotal)'; %feature index
                            self.neighborDictionary((k*nFeatureFunctionsTotal...
                                + 1):((k + 1)*nFeatureFunctionsTotal), 2) = ...
                                i; %coarse element index
                            self.neighborDictionary((k*nFeatureFunctionsTotal...
                                + 1):((k + 1)*nFeatureFunctionsTotal), 3) = ...
                                6; % lower left neighbor
                            k = k + 1;
                        end
                        
                    end
                    
                    
                    if(i > self.coarseMesh.nElX)
                        %lower neighbor of coarse element exists
                        PhiCell{n}(i, (k*nFeatureFunctionsTotal+ 1):((k + 1)*...
                            nFeatureFunctionsTotal)) =...
                            designMatrix{n}(i - self.coarseMesh.nElX, :);
                        self.neighborDictionary((k*nFeatureFunctionsTotal...
                            + 1):((k + 1)*nFeatureFunctionsTotal), 1) = ...
                            (1:nFeatureFunctionsTotal)'; %feature index
                        self.neighborDictionary((k*nFeatureFunctionsTotal...
                            + 1):((k + 1)*nFeatureFunctionsTotal), 2) = ...
                            i; %coarse element index
                        self.neighborDictionary((k*nFeatureFunctionsTotal...
                            + 1):((k + 1)*nFeatureFunctionsTotal), 3) = ...
                            4; %lower neighbor
                        k = k + 1;
                        
                        if(mod(i, self.coarseMesh.nElX) ~= 0)
                            %lower right neighbor of coarse element exists
                            PhiCell{n}(i, (k*nFeatureFunctionsTotal + 1):...
                                ((k + 1)*nFeatureFunctionsTotal)) =...
                                designMatrix{n}(i -...
                                self.coarseMesh.nElX + 1, :);
                            self.neighborDictionary((k*nFeatureFunctionsTotal...
                                + 1):((k + 1)*nFeatureFunctionsTotal), 1) = ...
                                (1:nFeatureFunctionsTotal)'; %feature index
                            self.neighborDictionary((k*nFeatureFunctionsTotal...
                                + 1):((k + 1)*nFeatureFunctionsTotal), 2) = ...
                                i; %coarse element index
                            self.neighborDictionary((k*nFeatureFunctionsTotal...
                                + 1):((k + 1)*nFeatureFunctionsTotal), 3) = ...
                                8; % lower right neighbor
                            k = k + 1;
                        end
                        
                    end
                end
            end
            if strcmp(mode, 'train')
                self.designMatrix = PhiCell;
            elseif strcmp(mode, 'test')
                self.testDesignMatrix = PhiCell;
            else
                error('wrong mode');
            end
            disp('done')
        end%includeLocalDiagNeighborFeatures

        function localTheta_c(self, designMatrix, mode)
            %Sets separate coefficients theta_c for each macro-cell in a 
            %single microstructure sample can never be executed before 
            %rescaling/standardization of design Matrix!
            debug = false; %debug mode
            disp(strcat('Using sep. feature coeff. theta_c for each', ...
                'macro-cell in a microstructure...'))
            nFeatureFunctionsTotal = size(designMatrix{1}, 2);
            PhiCell{1} = zeros(self.coarseMesh.nEl,...
                self.coarseMesh.nEl*nFeatureFunctionsTotal);
            nData = numel(designMatrix);
            PhiCell = repmat(PhiCell, nData, 1);
            
            %Reassemble design matrix
            for n = 1:nData
                for i = 1:self.coarseMesh.nEl
                    PhiCell{n}(i, ((i - 1)*nFeatureFunctionsTotal + 1):...
                        (i*nFeatureFunctionsTotal)) = ...
                        designMatrix{n}(i, :);
                end
                PhiCell{n} = sparse(PhiCell{n});
            end
            if debug
                firstDesignMatrixBeforeLocal = designMatrix{1}
                firstDesignMatrixAfterLocal = full(PhiCell{1})
                pause
            end
            if strcmp(mode, 'train')
                self.designMatrix = PhiCell;
            elseif strcmp(mode, 'test')
                self.testDesignMatrix = PhiCell;
            else
                error('wrong mode');
            end
            disp('done')
        end%localTheta_c
        
        function computeSumPhiTPhi(self)
            self.sumPhiTPhi = 0;
            for n = 1:numel(self.designMatrix)
                self.sumPhiTPhi = self.sumPhiTPhi +...
                    self.designMatrix{n}'*self.designMatrix{n};
            end
            if strcmp(self.mode, 'useLocal')
                self.sumPhiTPhi = sparse(self.sumPhiTPhi);
            end
        end
        

        
        %% plot functions
        function p = plotTrainingInput(self, samples, titl)
            %Load microstructures
            samplesTemp = min(samples):max(samples);
            cond = self.trainingDataMatfile.cond(:, samplesTemp);
            samples = samples - min(samples) + 1;
            f = figure;
            self.loadTrainingData;
            xLines = cumsum(self.coarseGridVectorX)*self.nElFX;
            yLines = cumsum(self.coarseGridVectorY)*self.nElFY;
            for i = 1:3
                subplot(1,3,i);
                p(i) = imagesc(reshape(cond(:, samples(i)),...
                    self.fineMesh.nElX, self.fineMesh.nElY));
                grid off;
                axis square;
                xticks({});
                yticks({});
                if nargin > 2
                    %to plot numerical title
                    title(num2str(titl(i)));
                end
                plotGrid = false;
                if plotGrid
                    for x = 1:(numel(self.coarseGridVectorX) - 1)
                        line([xLines(x), xLines(x)], [0, self.nElFX],...
                            'Color', 'w')
                    end
                    for y = 1:(numel(self.coarseGridVectorY) - 1)
                        line([0, self.nElFX], [yLines(y), yLines(y)],...
                            'Color', 'w')
                    end
                end
            end
            
        end

        function [p, im] = plotTrainingOutput(self, samples, titl)
            %Load microstructures
            samplesTemp = min(samples):max(samples);
            uf = self.trainingDataMatfile.uf(:, samplesTemp);
            cond = self.trainingDataMatfile.cond(:, samplesTemp);
            samples = samples - min(samples) + 1;
            min_uf = min(min(uf(:, samples)));
            max_uf = max(max(uf(:, samples)));
            %to use same color scale
            cond = ((min_uf - max_uf)/(min(min(cond)) -...
                max(max(cond))))*cond + max_uf - ...
                ((min_uf - max_uf)/(min(min(cond)) -...
                max(max(cond))))*max(max(cond));
            f1 = figure;
            f2 = figure;
            self.loadTrainingData;
            xLines = cumsum(self.coarseGridVectorX)*self.nElFX;
            yLines = cumsum(self.coarseGridVectorY)*self.nElFY;
            for i = 1:12
                figure(f1)
                subplot(4, 3, i);
                p(i) = surf(reshape(uf(:, samples(i)),...
                    (self.fineMesh.nElX + 1),...
                    (self.fineMesh.nElY + 1)));
                caxis([min_uf, max_uf])
                hold
                im(i) = imagesc(reshape(cond(:, samples(i)),...
                    self.fineMesh.nElX, self.fineMesh.nElY));
                p(i).LineStyle = 'none';
                grid on;
                axis tight;
                box on;
                axis square;
                zlim([min_uf, max_uf])
%                 xlabel('x')
%                 ylabel('y')
                zl = zlabel('$y(\vec{x})$');
                zl.Interpreter = 'latex';
                zl.Rotation = 90;
                zl.FontSize = 26;
                xticklabels({})
                yticklabels({})
                zticklabels({})
%                 xticks({});
%                 yticks({});
                if nargin > 2
                    %to plot numerical title
                    title(num2str(titl(i)));
                end
%                 for x = 1:(numel(obj.coarseGridVectorX) - 1)
%                     line([xLines(x), xLines(x)], [0, obj.nElFX], 'Color', 'w')
%                 end
%                 for y = 1:(numel(obj.coarseGridVectorY) - 1)
%                     line([0, obj.nElFX], [yLines(y), yLines(y)], 'Color', 'w')
%                 end

                figure(f2)
                subplot(4, 3, i);
                [~,p2(i)] = contourf(reshape(uf(:, samples(i)),...
                    (self.fineMesh.nElX + 1),...
                    (self.fineMesh.nElY + 1)), 8);
                caxis([min_uf, max_uf])
                grid off;
                p2(i).LineStyle = 'none';
                xticks({});
                yticks({});
                axis square;
                colorbar
            end
            
            plotSingleSamples = false;
            if plotSingleSamples
                for i = 1:10
                    f = figure;
                    subplot(1,2,1)
                    p(i) = imagesc(reshape(cond(:, samples(i)),...
                        self.fineMesh.nElX, self.fineMesh.nElY));
                    grid off;
                    axis square;
                    xticks({});
                    yticks({});
                    subplot(1,2,2)
                    q(i) = surf(reshape(uf(:, samples(i)),...
                        (self.fineMesh.nElX + 1),...
                        (self.fineMesh.nElY + 1)));
                    q(i).LineStyle = 'none';
                    grid on;
                    box on;
                    axis square;
                    axis tight;
                    xticks([64 128 192]);
                    yticks([64 128 192]);
                    zticks(100:100:800);
                    xticklabels({});
                    yticklabels({});
                    zticklabels({});
                    zlim([0 800])
                    print(f, strcat('~/images/uncecomp17/fineScaleSample',...
                        num2str(i)), '-dpng', '-r300')
                end
            end
            
        end

        function p = plot_p_c_regression(self, features)
            %Plots regressions of single features to the data <X>_q
            totalFeatures = size(self.featureFunctions, 2) + size(self.globalFeatureFunctions, 2);
            totalFeatures

            iter = 1;
            for feature = features
                if(~(totalFeatures ~= 1) || iter == 1)
                    f = figure;
                end

                k = 1;
                if strcmp(self.mode, 'useLocal')
                    mink = Inf*ones(self.coarseMesh.nEl, 1);
                    maxk = -Inf*ones(self.coarseMesh.nEl, 1);
                    for i = 1:self.coarseMesh.nElX
                        for j = 1:self.coarseMesh.nElY
                            subplot(self.coarseMesh.nElX,...
                                self.coarseMesh.nElY, k);
                            for s = 1:self.nTrain
                                yData(k, s) = self.XMean(k, s);
                                for l = 1:totalFeatures
                                    if l~= feature
                                        yData(k, s) = yData(k, s) -...
                                            self.theta_c.theta(totalFeatures*...
                                            (k - 1) + l)*...
                                            self.designMatrix{s}(k,...
                                            totalFeatures*(k - 1) + l);
                                    end
                                end
                                plot(self.designMatrix{s}(k,...
                                    totalFeatures*(k - 1) + feature),...
                                    yData(k, s), 'xb')
                                if(self.designMatrix{s}(k, totalFeatures*...
                                        (k - 1) + feature) < mink(k))
                                    mink(k) = self.designMatrix{s}(k,...
                                        totalFeatures*(k - 1) + feature);
                                elseif(self.designMatrix{s}(k, totalFeatures*...
                                        (k - 1) + feature) > maxk(k))
                                    maxk(k) = self.designMatrix{s}(k,...
                                        totalFeatures*(k - 1) + feature);
                                end
                                hold on;
                            end
                            x = linspace(mink(k), maxk(k), 100);
                            y = 0*x;
                            
                            useOffset = false;
                            if useOffset
                                %it is important that the offset feature
                                %phi(lambda) = 1 is the very first feature
                                y = self.theta_c.theta(...
                                    totalFeatures*(k - 1) + 1) +...
                                    self.theta_c.theta(totalFeatures*...
                                    (k - 1) + feature)*x;
                            else
                                y = self.theta_c.theta(totalFeatures*...
                                    (k - 1) + feature)*x;
                            end
                            plot(x, y);
                            axis tight;
                            axis square;
                            xl = xlabel('Feature function output $\phi_i$');
                            xl.Interpreter = 'latex';
                            yl = ylabel(...
                                '$<X_k> - \sum_{j\neq i} \theta_j \phi_j$');
                            yl.Interpreter = 'latex';
                            k = k + 1;
                        end
                    end
                elseif strcmp(self.mode, 'none')
                    for i = 1:self.coarseMesh.nElX
                        for j = 1:self.coarseMesh.nElY
                            for s = 1:self.nTrain
                                plot(self.designMatrix{s}(k, feature),...
                                    self.XMean(k, s), 'xb')
                                hold on;
                            end
                            k = k + 1;
                        end
                    end
                    
                    x = linspace(min(min(cell2mat(self.designMatrix))),...
                        max(max(cell2mat(self.designMatrix))), 100);
                    y = 0*x;
                    y = self.theta_c.theta(feature)*x;
                    plot(x, y);
                    axis tight;
                end
                iter = iter + 1;
            end
        end
        
        function plotTheta(self, figHandle)
            %Plots the current theta_c
            if isempty(self.thetaArray)
                self.thetaArray = self.theta_c.theta';
            else
                if(size(self.theta_c.theta, 1) > size(self.thetaArray, 2))
                    %New basis function included. Expand array
                    self.thetaArray = [self.thetaArray,...
                        zeros(size(self.thetaArray, 1),...
                        numel(self.theta_c.theta) - size(self.thetaArray, 2))];
                    self.thetaArray = [self.thetaArray; self.theta_c.theta'];
                else
                    self.thetaArray = [self.thetaArray; self.theta_c.theta'];
                end
            end
            if isempty(self.thetaHyperparamArray)
                self.thetaHyperparamArray = self.thetaPriorHyperparam';
            else
                if(size(self.thetaPriorHyperparam, 1) >...
                        size(self.thetaHyperparamArray, 2))
                    %New basis function included. Expand array
                    self.thetaHyperparamArray = [self.thetaHyperparamArray,...
                        zeros(size(self.thetaHyperparamArray, 1),...
                        numel(self.thetaPriorHyperparam) -...
                        size(self.thetaHyperparamArray, 2))];
                    self.thetaHyperparamArray = [self.thetaHyperparamArray;...
                        self.thetaPriorHyperparam'];
                else
                    self.thetaHyperparamArray = [self.thetaHyperparamArray;...
                        self.thetaPriorHyperparam'];
                end
            end
            if isempty(self.sigmaArray)
                self.sigmaArray = diag(self.theta_c.Sigma)';
            else
                self.sigmaArray = [self.sigmaArray; diag(self.theta_c.Sigma)'];
            end
%            figure(figHandle);
            sb1 = subplot(3, 2, 1, 'Parent', figHandle);
            plot(self.thetaArray, 'linewidth', 1, 'Parent', sb1)
            axis(sb1, 'tight');
            sb1.YLim = [(min(self.thetaArray(end, :)) - 1), (max(self.thetaArray(end, :)) + 1)];
            
            sb2 = subplot(3, 2, 2, 'Parent', figHandle);
            bar(self.theta_c.theta, 'linewidth', 1, 'Parent', sb2)
            axis(sb2, 'tight');
            
            sb3 = subplot(3, 2, 3, 'Parent', figHandle);
            semilogy(sqrt(self.sigmaArray), 'linewidth', 1, 'Parent', sb3)
            axis(sb3, 'tight');
            
            sb4 = subplot(3, 2, 4, 'Parent', figHandle);
            if ~self.theta_c.full_Sigma
                imagesc(reshape(diag(sqrt(self.theta_c.Sigma(...
                    1:self.coarseMesh.nEl,...
                    1:self.coarseMesh.nEl))),...
                    self.coarseMesh.nElX,self.coarseMesh.nElY),...
                    'Parent', sb4)
            else
                imagesc(reshape(sqrt(diag(self.theta_c.Sigma)),...
                    self.coarseMesh.nElX,self.coarseMesh.nElY),...
                    'Parent', sb4)
            end
            sb4.Title.String = '$\sigma_k$';
            colorbar('Parent', figHandle);
            sb4.GridLineStyle = 'none';
            axis(sb4, 'square');
            
            sb5 = subplot(3, 2, 5, 'Parent', figHandle);
            semilogy(self.thetaHyperparamArray, 'linewidth', 1, 'Parent', sb5)
            axis(sb5, 'tight');
            
            sb6 = subplot(3, 2, 6, 'Parent', figHandle);
            sb6.YScale = 'log';
            bar(self.thetaPriorHyperparam, 'linewidth', 1, 'Parent', sb6)
            axis(sb6, 'tight');
            drawnow
        end
        
        function plotCurrentState(self, fig, dataOffset)
            %Plots the current modal effective property and the modal
            %reconstruction for 2 -training- samples
            %load finescale conductivity field
            conductivity = self.trainingDataMatfile.cond(:, self.trainingSamples(1:4));
            for i = 1:4
                Lambda_eff_mode = conductivityBackTransform(self.designMatrix{i + dataOffset}*self.theta_c.theta,...
                    self.conductivityTransformation);
                sb1 = subplot(4, 3, 1 + (i - 1)*3, 'Parent', fig);
                imagesc(reshape(Lambda_eff_mode, self.coarseMesh.nElX, self.coarseMesh.nElY), 'Parent', sb1)
                sb1.YDir = 'normal';
                axis(sb1, 'tight');
                axis(sb1, 'square');
                sb1.GridLineStyle = 'none';
                sb1.XTick = [];
                sb1.YTick = [];
                cbp_lambda = colorbar('Parent', fig);
                sb2 = subplot(4, 3, 2 + (i - 1)*3, 'Parent', fig);
                nx = self.nElFX + 1;
                ny = self.nElFY + 1;
                XX = meshgrid(linspace(0, 1, nx));
                [~, YY] = meshgrid(linspace(0, 1, ny));
                conductivityHandle = imagesc(reshape(conductivity(:, i), self.nElFX, self.nElFY), 'Parent', sb2);
                sb2.YDir = 'normal';
                sb2.GridLineStyle = 'none';
                sb2.XTick = [];
                sb2.YTick = [];
                
                cbp_true = colorbar('Parent', fig);
                
                sb3 = subplot(4, 3, 3 + (i - 1)*3, 'Parent', fig);
                
                isotropicDiffusivity = true;
                if isotropicDiffusivity
                    coarseFEMout = heat2d(self.coarseMesh, Lambda_eff_mode);
                else
                    D = zeros(2, 2, self.coarseMesh.nEl);
                    for j = 1:self.coarseMesh.nEl
                        D(:, :, j) =  Lambda_eff_mode(j)*eye(2);
                    end
                    coarseFEMout = heat2d(self.coarseMesh, D);
                end
                
                uc = coarseFEMout.u;
                nxc = self.coarseMesh.nElX + 1;
                nyc = self.coarseMesh.nElY + 1;
                XXc = meshgrid(linspace(0, 1, self.coarseMesh.nElX + 1));
                [~, YYc] = meshgrid(linspace(0, 1, self.coarseMesh.nElY + 1));
                reconstructionHandle = surf(XXc, YYc, reshape(uc, nxc, nyc), 'Parent', sb3);
                reconstructionHandle.LineStyle = 'none';
                reconstructionHandle.FaceColor = 'b';
                hold(sb3, 'on');
                fineDataHandle2 = surf(XX, YY, reshape(self.fineScaleDataOutput(:, i + dataOffset),...
                    nx, ny), 'Parent', sb3);
                fineDataHandle2.LineStyle = 'none';
                hold(sb3, 'off');
                axis(sb3, 'tight');
                %sb3.ZLim = [mean(self.trainingData.P{i + dataOffset}) - ...
                %3*std(self.trainingData.P{i + dataOffset}), ...
                %mean(self.trainingData.P{i + dataOffset}) + ...
                %3*std(self.trainingData.P{i + dataOffset})];
                caxis(sb3, sb3.ZLim);
                axis(sb3, 'square');
                sb3.Box = 'on';
                sb3.BoxStyle = 'full';
                sb3.XTick = [];
                sb3.YTick = [];
                sb3.View = [-125, -5];
                cbp_reconst = colorbar('Parent', fig);
            end
            drawnow
        end

        
        
        function [X, Y, corrX_log_p_cf, corrX, corrY, corrXp_cf, corrpcfX, corrpcfY] = findMeshRefinement(self)
            %Script to sample d_log_p_cf under p_c to find where 
            %to refine mesh next

            uf = self.trainingDataMatfile.uf(:, self.nStart:(self.nStart + self.nTrain - 1));
            
            self.loadTrainedParams;
            theta_cfTemp = self.theta_cf;
            
            %comment this for inclusion of variances S of p_cf
            theta_cfTemp.S = ones(size(theta_cfTemp.S));
            
            theta_cfTemp.Sinv = sparse(1:self.fineMesh.nNodes,...
                1:self.fineMesh.nNodes, 1./theta_cfTemp.S);
            theta_cfTemp.Sinv_vec = 1./theta_cfTemp.S;
            %precomputation to save resources
            theta_cfTemp.WTSinv = theta_cfTemp.W'*theta_cfTemp.Sinv;
            theta_cfTemp.sumLogS = sum(log(theta_cfTemp.S));
            
            %% Compute design matrices
            if isempty(self.designMatrix)
                self.computeDesignMatrix('train');
            end
            
            nSamples = 1000;
            d_log_p_cf_mean = 0;
            log_p_cf_mean = 0;
            p_cfMean = 0;
            p_cfSqMean = 0;
            log_p_cfSqMean = 0;
            d_log_p_cf_sqMean = 0;
            k = 1;
            XsampleMean = 0;
            XsampleSqMean = 0;
            Xlog_p_cf_mean = 0;
            Xp_cfMean = 0;
            for i = self.nStart:(self.nStart + self.nTrain - 1)
                mu_i = self.designMatrix{i}*self.theta_c.theta;
                XsampleMean = ((i - 1)/i)*XsampleMean + (1/i)*mu_i;
                XsampleSqMean = ((i - 1)/i)*XsampleSqMean + (1/i)*mu_i.^2;
                for j = 1:nSamples
                    Xsample = mvnrnd(mu_i, self.theta_c.Sigma)';
                    conductivity = conductivityBackTransform(Xsample, self.conductivityTransformation);
                    [lg_p_cf, d_log_p_cf] = log_p_cf(uf(:, i), self.coarseMesh, conductivity,...
                        theta_cfTemp, self.conductivityTransformation);
                    d_log_p_cf_mean = ((k - 1)/k)*d_log_p_cf_mean + (1/k)*d_log_p_cf;
                    d_log_p_cf_sqMean = ((k - 1)/k)*d_log_p_cf_sqMean + (1/k)*d_log_p_cf.^2;
                    log_p_cf_mean = ((k - 1)/k)*log_p_cf_mean + (1/k)*lg_p_cf;
                    log_p_cfSqMean = ((k - 1)/k)*log_p_cfSqMean + (1/k)*lg_p_cf^2;
                    p_cfMean = ((k - 1)/k)*p_cfMean + (1/k)*exp(lg_p_cf);
                    p_cfSqMean = ((k - 1)/k)*p_cfSqMean + (1/k)*exp(2*lg_p_cf);
                    Xlog_p_cf_mean = ((k - 1)/k)*Xlog_p_cf_mean + (1/k)*Xsample*lg_p_cf;
                    Xp_cfMean = ((k - 1)/k)*Xp_cfMean + (1/k)*Xsample*exp(lg_p_cf);
                    k = k + 1;
                end
            end
            covX_log_p_cf = Xlog_p_cf_mean - XsampleMean*log_p_cf_mean;
            var_log_p_cf = log_p_cfSqMean - log_p_cf_mean^2;
            varX = XsampleSqMean - XsampleMean.^2;
            corrX_log_p_cf = covX_log_p_cf./(sqrt(var_log_p_cf)*sqrt(varX));
            
            covXp_cf = Xp_cfMean - XsampleMean*p_cfMean
            var_p_cf = p_cfSqMean - p_cfMean^2
            corrXp_cf = covXp_cf./(sqrt(var_p_cf)*sqrt(varX))
            
            d_log_p_cf_mean
            d_log_p_cf_var = d_log_p_cf_sqMean - d_log_p_cf_mean.^2
            d_log_p_cf_std = sqrt(d_log_p_cf_var)
            d_log_p_cf_err = d_log_p_cf_std/sqrt(nSamples*self.nTrain)
            d_log_p_cf_sqMean
            load('./data/noPriorSigma')
            noPriorSigma
            log_noPriorSigma = log(noPriorSigma)
            
            disp('Sum of grad squares in x-direction:')
            for i = 1:self.coarseMesh.nElY
                X(i) = sum(d_log_p_cf_sqMean(((i - 1)*self.coarseMesh.nElX + 1):(i*self.coarseMesh.nElX)));
                corrX(i) = sum(abs(corrX_log_p_cf(((i - 1)*self.coarseMesh.nElX + 1):(i*self.coarseMesh.nElX))));
                corrpcfX(i) = sum((corrXp_cf(((i - 1)*self.coarseMesh.nElX + 1):(i*self.coarseMesh.nElX))).^2);
            end
            
            disp('Sum of grad squares in y-direction:')
            for i = 1:self.coarseMesh.nElX
                Y(i) = sum(d_log_p_cf_sqMean(i:self.coarseMesh.nElX:((self.coarseMesh.nElY - 1)*...
                    self.coarseMesh.nElX + i)));
                corrY(i) = sum(abs(corrX_log_p_cf(i:self.coarseMesh.nElX:((self.coarseMesh.nElY - 1)*...
                    self.coarseMesh.nElX + i))));
                corrpcfY(i) = sum((corrXp_cf(i:self.coarseMesh.nElX:((self.coarseMesh.nElY - 1)*...
                    self.coarseMesh.nElX + i))).^2);
            end
        end

        
        
        %% Setter functions
        function setConductivityDistributionParams(self, condDistParams)
            self.conductivityDistributionParams = condDistParams;
            self.generateFineScaleDataPath;
        end

        function setBoundaryConditions(self, boundaryConditions)
            %Coefficients of boundary condition functions 
            %must be given as string
            assert(ischar(boundaryConditions), 'boundaryConditions must be given as string');
            self.boundaryConditions = boundaryConditions;
            self.genBoundaryConditionFunctions;
        end
        
        function pcaComponents = globalPCA(self)
            %Compute PCA on global microstructure - load when file exists
            fprintf('Performing PCA on global microstructure...')
            if exist(strcat(self.fineScaleDataPath, 'globalPCA.mat'), 'file')
                load(strcat(self.fineScaleDataPath, 'globalPCA.mat'));
            else
                if(size(self.trainingDataMatfile.cond, 2) < self.pcaSamples)
                    self.pcaSamples = size(self.trainingDataMatfile.cond, 2);
                    warning('Less samples than specified are available for PCA')
                end
                
                condUnsup = self.trainingDataMatfile.cond(:, 1:self.pcaSamples)';
                pcaComponents = pca(condUnsup, 'NumComponents', self.globalPcaComponents);
                save(strcat(self.fineScaleDataPath, 'globalPCA.mat'), 'pcaComponents');
            end
            
            pltPca = false;
            if pltPca
                figure
                for i = 1:min([size(pcaComponents, 2) 25])
                    subplot(5,5,i)
                    imagesc(reshape(pcaComponents(:, i), 256, 256))
                    axis square
                    grid off
                    colorbar
                    xticks({})
                    yticks({})
                end
                drawnow
            end
            fprintf('   done.\n')
        end
        
        function pcaComponents = localPCA(self)
            %Perform PCA on local macro-cells
            
            fprintf('Performing PCA on every macro-cell...')
            filename = strcat(self.fineScaleDataPath, 'localPCA', num2str(self.coarseGridVectorX),...
                   num2str(self.coarseGridVectorY), '.mat');
            if exist(filename, 'file')
                load(filename);
            else
                if(size(self.trainingDataMatfile.cond, 2) < self.pcaSamples)
                    self.pcaSamples = size(self.trainingDataMatfile.cond, 2);
                    warning('Less samples than specified are available for PCA')
                end
                conductivity_pca = self.trainingDataMatfile.cond(:, 1:self.pcaSamples);
                lambdak = self.get_coarseElementConductivities(conductivity_pca);
                averageMacroCells = true;
                if averageMacroCells
                    iter = 1;
                    for k = 1:self.coarseMesh.nEl
                        lambdakArray = zeros(numel(lambdak{1, k}),...
                            self.pcaSamples*self.coarseMesh.nEl);
                        for n = 1:self.pcaSamples
                            lambdakArray(:, iter) = lambdak{n, k}(:);
                            iter = iter + 1;
                        end
                    end
                    pcaComponents = pca(lambdakArray', 'NumComponents',...
                        self.localPcaComponents);
                else
                    for k = 1:self.coarseMesh.nEl
                        lambdakArray = zeros(numel(lambdak{1, k}),...
                            self.pcaSamples);
                        for n = 1:self.pcaSamples
                            lambdakArray(:, n) = lambdak{n, k}(:);
                        end
                        pcaComponents(:, :, k) = pca(lambdakArray',...
                            'NumComponents', self.localPcaComponents);
                    end
                end
                save(filename, 'pcaComponents');
            end
            pltPca = false;
            if pltPca
                figure
                for i = 1:min([size(pcaComponents, 2) 25])
                    subplot(5,5,i)
                    imagesc(reshape(pcaComponents(:,i), 64, 64))
                    axis square
                    grid off
                    colorbar
                    xticks({})
                    yticks({})
                end
                drawnow
            end
            fprintf('   done.\n')
        end

        function [phi, phiGlobal] = setFeatureFunctions(self)
            %Set up feature function handles; First cell array index is
            %for macro-cell. This allows different features for different
            %macro-cells
            addpath('./featureFunctions')   %Path to feature function library
            log_cutoff = 1e-5;
            phi = {};
            phiGlobal = {};
            %avoid broadcasting overhead in designMatrix
            ct = self.conductivityTransformation;
            
            for k = 1:self.coarseMesh.nEl
                nFeatures = 0;
                %constant bias
                phi{nFeatures + 1} = @(lambda, mask) constant(lambda, mask);
                nFeatures = nFeatures + 1;

%                 phi{k, nFeatures + 1} = @(lambda)...
%                     conductivityTransform(generalizedMean(lambda, -1), ct);
%                 nFeatures = nFeatures + 1;
%                 phi{k, nFeatures + 1} = @(lambda)...
%                     conductivityTransform(generalizedMean(lambda, -.5), ct);
%                 nFeatures = nFeatures + 1;
%                 phi{k, nFeatures + 1} = @(lambda)...
%                     conductivityTransform(generalizedMean(lambda, 0), ct);
%                 nFeatures = nFeatures + 1;
%                 phi{k, nFeatures + 1} = @(lambda)...
%                     conductivityTransform(generalizedMean(lambda, .5), ct);
%                 nFeatures = nFeatures + 1;
                phi{nFeatures + 1} = @(lambda, mask) generalizedMean(lambda, mask, 1);
                nFeatures = nFeatures + 1;

                
%                 phi{k, nFeatures + 1} = @(lambda) gaussLinFilt(lambda, nan, 1);
%                 nFeatures = nFeatures + 1;
%                 phi{k, nFeatures + 1} = @(lambda) gaussLinFilt(lambda, nan, 2);
%                 nFeatures = nFeatures + 1;
%                 phi{k, nFeatures + 1} = @(lambda) gaussLinFilt(lambda, nan, 4);
%                 nFeatures = nFeatures + 1;
%                 phi{k, nFeatures + 1} = @(lambda) gaussLinFilt(lambda, nan, 8);
%                 nFeatures = nFeatures + 1;
%                 phi{k, nFeatures + 1} = @(lambda) gaussLinFilt(lambda, nan, 16);
%                 nFeatures = nFeatures + 1;
% 
%                 phi{k, nFeatures + 1} = @(lambda) std(lambda(:));
%                 nFeatures = nFeatures + 1;
%                 phi{k, nFeatures + 1} = @(lambda) log(std(lambda(:)) +...
%                     log_cutoff);
%                 nFeatures = nFeatures + 1;
%                 
%                 phi{k, nFeatures + 1} = @(lambda) isingEnergy(lambda);
%                 nFeatures = nFeatures + 1;
%                 
%                 phi{k, nFeatures + 1} = @(lambda) generalizedMeanBoundary(...
%                     lambda, -1, 'left');
%                 nFeatures = nFeatures + 1;
%                 phi{k, nFeatures + 1} = @(lambda) generalizedMeanBoundary(...
%                     lambda, -1, 'lower');
%                 nFeatures = nFeatures + 1;
%                 phi{k, nFeatures + 1} = @(lambda) generalizedMeanBoundary(...
%                     lambda, -1, 'right');
%                 nFeatures = nFeatures + 1;
%                 phi{k, nFeatures + 1} = @(lambda) generalizedMeanBoundary(...
%                     lambda, -1, 'upper');
%                 nFeatures = nFeatures + 1;
%                 
%                 phi{k, nFeatures + 1} = @(lambda) generalizedMeanBoundary(...
%                     lambda, 0, 'left');
%                 nFeatures = nFeatures + 1;
%                 phi{k, nFeatures + 1} = @(lambda) generalizedMeanBoundary(...
%                     lambda, 0, 'lower');
%                 nFeatures = nFeatures + 1;
%                 phi{k, nFeatures + 1} = @(lambda) generalizedMeanBoundary(...
%                     lambda, 0, 'right');
%                 nFeatures = nFeatures + 1;
%                 phi{k, nFeatures + 1} = @(lambda) generalizedMeanBoundary(...
%                     lambda, 0, 'upper');
%                 nFeatures = nFeatures + 1;
%                 
%                 phi{k, nFeatures + 1} = @(lambda) generalizedMeanBoundary(...
%                     lambda, 1, 'left');
%                 nFeatures = nFeatures + 1;
%                 phi{k, nFeatures + 1} = @(lambda) generalizedMeanBoundary(...
%                     lambda, 1, 'lower');
%                 nFeatures = nFeatures + 1;
%                 phi{k, nFeatures + 1} = @(lambda) generalizedMeanBoundary(...
%                     lambda, 1, 'right');
%                 nFeatures = nFeatures + 1;
%                 phi{k, nFeatures + 1} = @(lambda) generalizedMeanBoundary(...
%                     lambda, 1, 'upper');
%                 nFeatures = nFeatures + 1;
            end
            
%             nGlobalFeatures = 0;
%             %Unsupervised pretraining: compute PCA components
%             if(self.globalPcaComponents > 0)
%                 globalComponents = self.globalPCA;
%             end
%             %PCA projection
%             for n = 1:self.globalPcaComponents
%                 for k = 1:self.coarseMesh.nEl
%                     phiGlobal{k, nGlobalFeatures + n} = @(lambda)...
%                         globalComponents(:, n)'*lambda(:);
%                 end
%             end
%             
%             %local PCA projection
%             if(self.localPcaComponents > 0)
%                 localComponents = self.localPCA;
%                 for n = 1:self.localPcaComponents
%                     for k = 1:self.coarseMesh.nEl
%                         if(ndims(localComponents) == 3)
%                             %Separate PCA on every macro cell
%                             phi{k, nFeatures + n} = @(lambda)...
%                                 localComponents(:, n, k)'*lambda(:);
%                         else
%                             phi{k, nFeatures + n} = @(lambda)...
%                                 localComponents(:, n)'*lambda(:);
%                         end
%                     end
%                 end
%             end
                        
            self.featureFunctions = phi;
            self.globalFeatureFunctions = phiGlobal;
        end

        function setCoarseGrid(self, coarseGridX, coarseGridY)
            %coarseGridX and coarseGridY are coarse model grid vectors
            self.coarseGridVectorX = coarseGridX;
            self.coarseGridVectorY = coarseGridY;
        end
    end
    
end


















