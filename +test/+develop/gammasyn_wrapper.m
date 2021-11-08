clear; close all;

%% Settings
loaded_system = test.develop.load_system('ABS_I_Anteil_geregelt','combined');

%% pole area parameters
a = 15;
b = 1000;
R = 50;

%% gammasyn options
solver = optimization.solver.Optimizer.FMINCON; %FMINCON;% solver to use
options = optimization.options.OptionFactory.instance.options(solver,...
	'Retries',						1,...											% number of retries
	'Algorithm',					solver.getDefaultAlgorithm(),...				% algorithm of solver, for not builtin solvers name of the solver, e.g. 'snopt' for SNOPT
	'FunctionTolerance',			1E-5,...
	'StepTolerance',				1E-5,...
	'ConstraintTolerance',			1.4e-5,...
	'MaxFunctionEvaluations',		5E3,...
	'MaxIterations',				5E3,...
	'MaxSQPIter',					5E3,...
	'SpecifyObjectiveGradient',		true,...
	'SpecifyConstraintGradient',	true,...
	'CheckGradients',				false,...
	'FunValCheck',					false,...
	'FiniteDifferenceType',			'forward',...
	'Diagnostics',					false,...
	'Display',						'iter-detailed'...
);
objectiveoptions = struct(...
	'usecompiled',				false,...											% indicator, if compiled functions should be used
	'type',						GammaJType.LINEAR,...								% type of pole area weighting in objective function
	'allowvarorder',			false,...											%allow variable state number for different multi models
	'eigenvaluederivative',		GammaEigenvalueDerivativeType.VANDERAA,...
	'errorhandler',				GammaErrorHandler.USER,...
	'errorhandler_function',	@control.design.gamma.test.errorhandler_log,...
	'strategy',					GammaSolutionStrategy(0)...
	);

%% Pole area
weight = {5, repmat([1, 10], loaded_system.mu, 1)};
polearea = {repmat([
	control.design.gamma.area.Circlesquare(R),	control.design.gamma.area.Hyperbolasquare(a, b),...
	control.design.gamma.area.Imag(1, a) %Imag(1, a)
], loaded_system.mu, 1), repmat([
	control.design.gamma.area.Circlesquare(R/2, -R/2, 0), control.design.gamma.area.Imag(1, a)
], loaded_system.mu, 1)};
polearea = polearea{1};
weight = weight{1};

%% Transform Systems
KDF_0 = loaded_system.KDF_0{1};

for jj = 1:loaded_system.mu
	system(jj, 1) = struct(...
		'A',	loaded_system.Ai{jj},...
		'B',	loaded_system.Bi{jj},...
		'C',	eye(6)... %[loaded_system.C1; loaded_system.C2i{jj}]...
	);
end

%% K_fixed for ABS

K_fixed_A = zeros(2, 6, 2);
K_fixed_B = zeros(2, 1);
K_fixed_A(1, 6, 1) = 1;
K_fixed_A(2, 5, 2) = 1;
K_fixed = {K_fixed_A, K_fixed_B};
K_fixed = [];

%% gammasyn_couplingcontrol
[Kopt, Jopt, information] = control.design.gamma.gammasyn(system, polearea, weight, K_fixed, KDF_0, options, objectiveoptions);
if isempty(Kopt)
	return;
end

Kopt
information

%% Analysis
for jj = 1:loaded_system.mu
	system(jj, 1).C = [loaded_system.C1; loaded_system.C2i{jj}];
end

T = 2;
if ~isempty(Kopt)
	if 1
		gain_ratio = test.develop.analyze_results(system, Kopt, T, [], 0);
		minimal_deveation = 100/min(abs(gain_ratio))
	end
end