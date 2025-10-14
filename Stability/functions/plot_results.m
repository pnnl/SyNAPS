%% plotting

function [sample_set, assess_stability] = plot_results(q,q0,Aq,constr_func,plot_opt)

% read plot options
plot_vars = plot_opt.varidx;
switch plot_opt.valtyp
    case 'abs'
        plot_vals = plot_opt.ranges;
    case 'rel'
        plot_vals = kron(q0(plot_vars),ones(size(plot_vars))) + plot_opt.ranges;
    otherwise
        error('unknown varidx specified in plot_opt!')
end
npoints = plot_opt.npoint;

omit_vars = setdiff(1:length(q),plot_vars);

plot_vars_label = convert_into_nice_labels(q(plot_vars));

switch plot_opt.figtyp
    case 'verify'

        % analytical estimate
        figure;
        plot(-1e3,-1e3,'go','MarkerSize',8,'LineWidth',2); hold on
        plot(-1e3,-1e3,'bx','MarkerSize',8,'LineWidth',2)
        plot(-1e3,-1e3,'k+','MarkerSize',8,'LineWidth',2)
        plot(-1e3,-1e3,'ro','MarkerSize',8,'LineWidth',2)
        legend('(stable, stable)','(--,stable)','(--,unstable)',...
            '(stable, unstable)','AutoUpdate','off','Orientation','horizontal')        

        % set plot axes' limits
        xlim(plot_vals(1:2))
        ylim(plot_vals(3:4))
        xlabel(plot_vars_label{1});%char(q(plot_vars(1))))
        ylabel(plot_vars_label{2});%char(q(plot_vars(2))))

        switch plot_opt.algo
            case 'oneshot'
                % plot contour boundary
                if length(constr_func)<=32
                    constr_sub = vpa(subs(constr_func,q(omit_vars),q0(omit_vars)));
                    fcontour(eval(plot_func_script(q,plot_vars,constr_sub)),plot_vals,'LevelList',0); hold on
                else
                    fprintf('\n*** [Plotting] Skipping contour plot since >32 constraints! ***\n')
                end
                % numerical validation

                sample_set = NaN(npoints,length(q));
                sample_set(:,omit_vars) = kron(ones(npoints,1),q0(omit_vars));

                % rng default % for reproducibility
                sample_norm = lhsdesign(npoints,length(plot_vars));

                for ip = 1:length(plot_vars)
                    min_val = plot_vals(2*ip-1);
                    max_val = plot_vals(2*ip);
                    sample_set(:,plot_vars(ip)) = min_val + (max_val-min_val)*sample_norm(:,ip);
                end

                assess_stability = NaN(npoints,2); tic
                for in = 1:npoints
                    est_stability = min(double(subs(constr_func,q,sample_set(in,:)))); % analytical
                    act_stability = -max(real(eig(double(subs(Aq,q,sample_set(in,:)))))); % numerical

                    assess_stability(in,:) = [est_stability act_stability];

                    % now plot the sample points with color
                    if assess_stability(in,:)>-plot_opt.tol == [0 0]
                        plot_color = 'k+';
                    elseif assess_stability(in,:)>-plot_opt.tol == [0 1]
                        plot_color = 'bx';
                    elseif assess_stability(in,:)>-plot_opt.tol == [1 0]
                        plot_color = 'ro';
                    else
                        plot_color = 'go';
                    end

                    plot(sample_set(in,plot_vars(1)),sample_set(in,plot_vars(2)),plot_color,'MarkerSize',8,'LineWidth',2)

                    if mod(in*100/npoints,10)==0
                        fprintf('\nCompleted (%4d/%4d) %d%%\t(t=%.2f)\n',in,npoints,floor(in*100/npoints),toc)
                    end
                end

                % plot the supplied point (q0)
                plot(q0(plot_vars(1)),q0(plot_vars(2)),'g+','MarkerSize',8,'LineWidth',2)

            otherwise
                avar = setdiff(symvar(constr_func),q);

                % plot contour boundary

                constr_sub = vpa(subs(constr_func(1:end-length(avar)-1),q(omit_vars),q0(omit_vars)));

                num_a_points = 100;
                a_sample_pts = lhsdesign(num_a_points,length(avar));
                a_sample_pts = [a_sample_pts(sum(a_sample_pts,2)<=1,:); zeros(1,length(avar))];
                for ia = 1:length(avar)
                    aval_i = zeros(size(avar));
                    aval_i(ia) = 1;
                    a_sample_pts = [a_sample_pts; aval_i];
                end

                for iSample = 1:size(a_sample_pts,1)
                    constr_i = vpa(subs(constr_sub,avar,a_sample_pts(iSample,:)));
                    fcontour(eval(plot_func_script(q,plot_vars,constr_i)),plot_vals,'LevelList',0); hold on
                end

                % numerical validation

                sample_set = NaN(npoints,length(q));
                sample_set(:,omit_vars) = kron(ones(npoints,1),q0(omit_vars));

                % rng default % for reproducibility
                sample_norm = lhsdesign(npoints,length(plot_vars));

                for ip = 1:length(plot_vars)
                    min_val = plot_vals(2*ip-1);
                    max_val = plot_vals(2*ip);
                    sample_set(:,plot_vars(ip)) = min_val + (max_val-min_val)*sample_norm(:,ip);
                end

                assess_stability = NaN(npoints,2); tic
                for in = 1:npoints
                    % solve an optimization
                    constr_sub = vpa(simplify(subs(constr_func,q,sample_set(in,:)))); % contains 'a'
                    constr_sub = vpa(constr_sub./max(1,abs(double(subs(constr_sub,avar,zeros(size(avar)))))));

                    cvx_begin quiet; variable cmin;
                    for ia = 1:length(avar)
                        eval(sprintf('variable %s nonnegative',char(avar(ia))));
                    end
                    subject to;
                    for ic = 1:length(constr_sub)-3*length(q)-1
                        eval(sprintf('%s >= cmin;',char(constr_sub(ic))));
                    end
                    for ic = length(constr_sub)-3*length(q):length(constr_sub)
                        eval(sprintf('%s >= 0;',char(constr_sub(ic))));
                    end
                    minimize -cmin; cvx_end

                    est_stability = cmin; % analytical

                    act_stability = -max(real(eig(double(subs(Aq,q,sample_set(in,:)))))); % numerical

                    assess_stability(in,:) = [est_stability act_stability];

                    % now plot the sample points with color
                    if assess_stability(in,:)>-plot_opt.tol == [0 0]
                        plot_color = 'k+';
                    elseif assess_stability(in,:)>-plot_opt.tol == [0 1]
                        plot_color = 'bx';
                    elseif assess_stability(in,:)>-plot_opt.tol == [1 0]
                        plot_color = 'ro';
                    else
                        plot_color = 'go';
                    end
                    
                    plot(sample_set(in,plot_vars(1)),sample_set(in,plot_vars(2)),plot_color,'MarkerSize',8,'LineWidth',2)

                    if mod(in*100/npoints,10)==0
                        fprintf('\nCompleted (%4d/%4d) %d%%\t(t=%.2f)\n',in,npoints,floor(in*100/npoints),toc)
                    end
                end
        end

    case 'demo'
        % set plot axes' limits
        xlim(plot_vals(1:2))
        ylim(plot_vals(3:4))
        xlabel(plot_vars_label{1});%char(q(plot_vars(1))))
        ylabel(plot_vars_label{1});%char(q(plot_vars(2))))

        switch plot_opt.algo
            case 'oneshot'
                % plot contour boundary

                constr_sub = vpa(subs(constr_func,q(omit_vars),q0(omit_vars)));

                fcontour(eval(plot_func_script(q,plot_vars,constr_sub)),plot_vals,...
                    'LevelList',0,'LineColor','b','LineStyle','--','LineWidth',2); hold on

                % plot the supplied point (q0)
                plot(q0(plot_vars(1)),q0(plot_vars(2)),'b+','MarkerSize',8,'LineWidth',2)

            otherwise
                error('plot_opt.algo NOT suitable for DEMO plot!')
        end

        sample_set          = [];
        assess_stability    = [];

    otherwise
        error('Unknown plot_opt.figtyp encountered!')
end

end

%% plotting scipt for contours

function plot_func = plot_func_script(q,plot_vars,constr_sub)

plot_func1 = sprintf('@(%s)',q(plot_vars(1)));
for ip = 2:length(plot_vars)
    plot_func1 = strcat(plot_func1(1:end-1), sprintf(',%s)',q(plot_vars(ip))));
end

plot_func2 = sprintf('%s',constr_sub(1));
for ic = 2:length(constr_sub)
    plot_func2 = sprintf('min(%s,%s)',constr_sub(ic),plot_func2);
end

% plot_func2 = sprintf('min([ (%s) ])',constr_sub(1));
% for ic = 2:length(constr_sub)
%     plot_func2 = strcat(plot_func2(1:end-2), sprintf('(%s) ])',constr_sub(ic)));
% end

plot_func = strcat(plot_func1, plot_func2);

end



%% function to make the labels look nice
function plot_vars_label = convert_into_nice_labels(qvars)

plot_vars_label = cell(size(qvars));

for ivar = 1:length(qvars)
    split_vars = regexp(char(qvars(ivar)),'\_','split');
    plot_vars_label{ivar} = split_vars{1};
    for jvar = 2:length(split_vars)
        plot_vars_label{ivar} = strcat(plot_vars_label{ivar}, ...
            '_{',split_vars{jvar},'}');
    end
end

end