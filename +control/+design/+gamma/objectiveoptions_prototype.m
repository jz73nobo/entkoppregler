function [objectiveoption_prototype] = objectiveoptions_prototype()
	%OBJECTIVEOPTIONS_PROTOTYPE return prototype structure for objective options
	%	Output:
	%		objectiveoption_prototype:	prototype for objective options
	% TODO: add deprecation warning
	coupling_prototype = control.design.gamma.objectiveoptions_prototype_coupling();
	objectiveoption_prototype = struct(...
		'usecompiled',			false,...
		'numthreads',			uint32(configuration.matlab.numthreads()),...
		'type',					GammaJType.ZERO,...
		'weight',				0,...
		'allowvarorder',		false,...
		'eigenvaluederivative',	GammaEigenvalueDerivativeType.getDefaultValue(),...
		'eigenvaluefilter',		GammaEigenvalueFilterType.getDefaultValue(),...
		'eigenvalueignoreinf',	false,...
		'couplingcontrol',		coupling_prototype,...
		'objective',			struct(...
			'preventNaN',		false,...
			'kreisselmeier',	struct(...
				'rho',				20,...
				'max',				(0)...
			),...
			'lyapunov',			struct(...
				'Q',				[]...
			),...
			'normgain',			struct(...
				'R',	[],...
				'K',	[],...
				'F',	[]...
			)...
		)...
	);
end