function [gain_ratio, gains_all, poles_all, F] = analyze_results(systems, K_opt, polearea, number_couplingconditions)
	gain_ratio = NaN;
	gains_all = NaN;
	if nargin <= 3
		polearea = [];
	end
	if isa(K_opt, 'double')
		R = K_opt;
		F = eye(size(R, 1));
	else
		R = K_opt{1};
		F = K_opt{3};
	end
	couplingcontrol = number_couplingconditions > 0;

	number_models = length(systems);
	number_states = size(systems(1).A, 2);
	number_controls = size(systems(1).B, 2);
	number_references = size(systems(1).C_ref, 1);

	if couplingcontrol
		F1 = F(:, 1:number_controls - number_couplingconditions);
	else
		F1 = F;
	end
	poles_all = cell(number_models, 1);
	all_E_matrices = struct2cell(systems).';
	all_E_matrices = all_E_matrices(:, 1);
	all_E_matrices = blkdiag(all_E_matrices{:});
	if rank(all_E_matrices) ~= number_models*number_states
		return;
	end

	%% calculate optimal prefilter gain to match stationary accuracy in quadratically optimal sense
	N_R = cell(number_models, 1);
	N_L = cell(number_models, 1);
	parfor ii = 1:number_models
		E = systems(ii).E;
		A = E\systems(ii).A;
		B = E\systems(ii).B;
		C = systems(ii).C_ref;
		if isfield(systems(ii), 'D_ref')
			D = systems(ii).D_ref;
		else
			D = zeros(number_references, number_controls);
		end

		if couplingcontrol
			C = C(1:end-number_couplingconditions, :);
			D = D(1:end-number_couplingconditions, :);
		end
		N_R{ii} = (C-D*R)/(-A + B*R)*B*F1 + D*F1;
		N_L{ii} = eye(number_controls - number_couplingconditions);
	end
	N_R = cat(1, N_R{:});
	N_L = cat(1, N_L{:});
	Q_F  = N_R\N_L;
	if couplingcontrol
		F = [F1*Q_F, F(:, end - number_couplingconditions + 1:end)];
	else
		F = F1*Q_F;
	end
	fprintf('\nThe prefilter for quadratically optimal stationary accuracy is:\n')
	disp(F);

	%% calculate poles
	ss_all = cell(number_models, 1);
	gain_ratio = zeros(number_states, 1);
	gains_all = cell(number_models, 1);
	parfor ii = 1:number_models
		E = systems(ii).E;
		A = E\systems(ii).A;
		B = E\systems(ii).B;
		C = systems(ii).C_ref;
		if isfield(systems(ii), 'D_ref')
			D = systems(ii).D_ref;
		else
			D = zeros(number_references, number_controls);
		end
		ss_all{ii} = ss(A - B*R, B*F, C-D*R, D*F);
		poles = pole(ss_all{ii});
		if any(poles > 0)
			warning('System %d unstable!',ii);
		end
		gains_all{ii} = dcgain(ss_all{ii});
		gain_ratio(ii) = gains_all{ii}(1,1)/gains_all{ii}(2,1);
		poles_all{ii} = poles;
	end
	poles_all = cat(1, poles_all{:});
	T_max = max(1./abs(poles_all));

	%% simulation specification
	T = max([5, round(T_max*3.5)*10]); % round 35*T_max to nearest 10, but limit to at least 5
	steps = 1000*T;
	t = linspace(0,T,steps+1);
	w = zeros(number_controls,steps+1);
	w(1,t<T/3) = 1;
	w(1,t>=T/3 & t<2*T/3) = 2;
	w(1,t>=2*T/3 & t<3*T/4) = 1.5;
	w(1,t>=3*T/4) = 0;
	number_references = number_controls - number_couplingconditions;
	for ii = 2:number_references
		shift = ceil(ii/number_references*steps);
		w(ii, :) = ii*[w(1, shift + 1:end), w(1, 1:shift)];
	end
	x0 = linspace(-0.7, 0.7, number_states);

	%% plot
	figure();
	red = linspace(1, 0, number_controls);
	green = linspace(0, 0, number_controls);
	blue = linspace(0, 1, number_controls);
	colors = [red.', green.', blue.'];

	ax_out = subplot(6,2,1:2:5);
	plot(t, w, 'k');
	xlim([0, T]);
	set(gca,'XTickLabel',[]);
	set(gca, 'ColorOrder', colors);
	title('Outputs')
	hold on;
	grid on;

	ax_in = subplot(6,2,7:2:11);
	xlim([0, T]);
	set(gca, 'ColorOrder', colors);
	xlabel('t/s');
	title('Inputs')
	hold on;
	grid on;

	subplot(6,2,2:2:12);
	title('Eigenvalues');
	hold on;
	grid on;

	plots_y = cell(number_models, number_controls);
	plots_u = cell(number_models, number_controls);

	parfor ii = 1:number_models
		sys = ss_all{ii}
		[y,ts,x] = lsim(sys,w,t,x0);
		u = (-R*x.'+F*w).';
		for jj = 1:number_controls
			subplot(6,2,1:2:5);
			plots_y{ii, jj} = plot(ts, y(:, jj));
			subplot(6,2,7:2:11);
			plots_u{ii, jj} = plot(ts, u(:, jj));
		end
	end
	subplot(6,2,2:2:12);
	scatter(real(poles_all).',imag(poles_all).','x','MarkerEdgeColor','k');

	% plot pole area
	if ~isempty(polearea)
		if iscell(polearea) && ~isempty(polearea{1})
			area = polearea{1}.unique(polearea{1});
		else
			area = polearea.unique(polearea);
		end
		if ~isempty(area)
			hold('on');
			area.plotborder();
		end
	end

	% legends
	legend_strings_y = cell(number_controls, 1);
	legend_strings_u = cell(number_controls, 1);
	for ii = 1:number_controls
		legend_strings_y{ii, 1} = sprintf('y_%d', ii);
		legend_strings_u{ii, 1} = sprintf('u_%d', ii);
	end
	subplot(6,2,1:2:5);
	lgd_y = legend([plots_y{1, :}], 'Location', 'southeast', 'Orientation', 'horizontal');
	lgd_y.String = legend_strings_y;
	subplot(6,2,7:2:11);
	lgd_u = legend([plots_u{1, :}]);
	lgd_u.String = legend_strings_u;
	linkaxes([ax_out ax_in],'x');
end