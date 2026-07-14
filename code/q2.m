% =========================================================================
% 华中杯 A题 问题二：环保政策约束下的双目标调度优化
% 模型: TD-HVRPTW + 绿色区限行 + 时段分流 + epsilon-约束法
% 算法: 自适应大邻域搜索 (ALNS) + Pareto前沿膝点分析
% =========================================================================
clear; clc; close all; warning off; rng(1);

%% 1. 数据预处理
disp('===== 问题二：绿色政策调度 + 双目标优化引擎启动 =====');

dist_matrix = readmatrix('距离矩阵.xlsx');
dist_matrix(isnan(dist_matrix)) = 0;
dist_matrix = max(dist_matrix, dist_matrix');

coords_table = readtable('客户坐标信息.xlsx');
coords = [coords_table{:,2}, coords_table{:,3}, coords_table{:,4}];
n_cust = size(coords,1) - 1;

orders = readtable('订单信息.xlsx');
weight_col = fillmissing(orders{:,2}, 'constant', 0);
volume_col = fillmissing(orders{:,3}, 'constant', 0);
target_col = orders{:,4};

demand_w = zeros(n_cust,1);
demand_v = zeros(n_cust,1);
for i = 1:n_cust
    idx = (target_col == i);
    demand_w(i) = sum(weight_col(idx));
    demand_v(i) = sum(volume_col(idx));
end

tw_table = readtable('时间窗.xlsx');
tw = zeros(n_cust, 2);
for i = 1:n_cust
    s = char(tw_table{i,2}); p = sscanf(s,'%d:%d'); tw(i,1) = p(1)*60+p(2);
    s = char(tw_table{i,3}); p = sscanf(s,'%d:%d'); tw(i,2) = p(1)*60+p(2);
end

fprintf('数据加载完成: %d 客户, %d 订单\n', n_cust, height(orders));

%% 1.5 需求拆分 — 超容量客户拆为多个虚拟客户
MAX_CAP_W = 3000; MAX_CAP_V = 15.0;
orig_n_cust = n_cust;
orig_dw = demand_w; orig_dv = demand_v; orig_tw = tw;
orig_coords = coords; orig_dist = dist_matrix;

new_dw = []; new_dv = []; new_tw = []; new_coords_ext = [coords(1,:)];
parent_map = [];
for i = 1:orig_n_cust
    w = orig_dw(i); v = orig_dv(i);
    n_split = max(ceil(w / MAX_CAP_W), ceil(v / MAX_CAP_V));
    sw_each = w / n_split;
    sv_each = v / n_split;
    for s = 1:n_split
        new_dw(end+1,1) = sw_each; new_dv(end+1,1) = sv_each;
        new_tw(end+1,:) = orig_tw(i,:);
        new_coords_ext(end+1,:) = orig_coords(i+1,:);
        parent_map(end+1) = i;
    end
end
BIG_ONLY_W = 1500;
BIG_ONLY_V = 10.8;
LARGE_STOCK = 70;
is_big_only = (new_dw > BIG_ONLY_W) | (new_dv > BIG_ONLY_V);
while sum(is_big_only) > LARGE_STOCK
    idxs = find(is_big_only);
    [~, p] = min(new_dw(idxs) / MAX_CAP_W + new_dv(idxs) / MAX_CAP_V);
    idx = idxs(p);
    new_dw(idx) = new_dw(idx) / 2;
    new_dv(idx) = new_dv(idx) / 2;
    new_dw(end+1,1) = new_dw(idx); new_dv(end+1,1) = new_dv(idx);
    new_tw(end+1,:) = new_tw(idx,:);
    new_coords_ext(end+1,:) = new_coords_ext(idx+1,:);
    parent_map(end+1) = parent_map(idx);
    is_big_only = (new_dw > BIG_ONLY_W) | (new_dv > BIG_ONLY_V);
end
n_cust = length(new_dw);
demand_w = new_dw; demand_v = new_dv; tw = new_tw; coords = new_coords_ext;
new_dist = zeros(n_cust+1);
for i = 0:n_cust
    pi = 0; if i>0, pi = parent_map(i); end
    for j = 0:n_cust
        pj = 0; if j>0, pj = parent_map(j); end
        new_dist(i+1, j+1) = orig_dist(pi+1, pj+1);
    end
end
dist_matrix = new_dist;
fprintf('需求拆分: %d -> %d 虚拟客户\n', orig_n_cust, n_cust);

%% 2. 绿色配送区识别
cx = coords(2:end, 2);
cy = coords(2:end, 3);
dist_to_center = sqrt(cx.^2 + cy.^2);
R_GREEN = 10;
is_green = dist_to_center <= R_GREEN;
green_custs = find(is_green);
normal_custs = find(~is_green);
fprintf('绿色区客户: %d 个, 常规区客户: %d 个\n', length(green_custs), length(normal_custs));

%% 3. 车型定义与政策参数
% [载重kg, 容积m3, 启动成本, 满载修正, 碳排系数kg/L或kg/kWh, 能源单价, 是否燃油, 库存]
VT = [
    3000, 13.5, 400, 0.40, 2.547, 7.61, 1, 60;
    1500, 10.8, 400, 0.40, 2.547, 7.61, 1, 50;
    1250,  6.5, 400, 0.40, 2.547, 7.61, 1, 50;
    3000, 15.0, 400, 0.35, 0.501, 1.64, 0, 10;
    1250,  8.5, 400, 0.35, 0.501, 1.64, 0, 15;
];
n_vtype = size(VT,1);

CARBON_PRICE = 0.65;
EARLY_PEN    = 20/60;
LATE_PEN     = 50/60;
SERVICE_TIME = 20;
T_START      = 480;          % 8:00
T_FORBID_END = 960;          % 16:00 限行结束

%% 4. ALNS 参数
ALNS_ITER    = 2000;
N_DESTROY    = 3;
N_REPAIR     = 2;
SIGMA        = [33, 9, 3];
REACT_FACTOR = 0.1;
DESTROY_RATE = [0.1, 0.4];
SA_T0        = 1000;
SA_COOL      = 0.9995;

%% 5. epsilon-约束法双目标扫描
N_EPSILON = 8;

disp('>>> 锚点1: 求解成本最优解...');
[sol_cost_opt, cost_opt, det_cost, carb_cost_opt, hist_cost] = ...
    run_alns_q2(n_cust, demand_w, demand_v, tw, dist_matrix, coords, VT, ...
    CARBON_PRICE, EARLY_PEN, LATE_PEN, SERVICE_TIME, T_START, ...
    is_green, T_FORBID_END, ALNS_ITER, SIGMA, REACT_FACTOR, DESTROY_RATE, SA_T0, SA_COOL, ...
    'cost', inf);
fprintf('  成本最优: 总成本=%.2f 元, 碳排=%.2f kg\n', cost_opt, carb_cost_opt);

disp('>>> 锚点2: 求解碳排最优解...');
[sol_carb_opt, cost_carb_opt, det_carb, carb_min, hist_carb] = ...
    run_alns_q2(n_cust, demand_w, demand_v, tw, dist_matrix, coords, VT, ...
    CARBON_PRICE, EARLY_PEN, LATE_PEN, SERVICE_TIME, T_START, ...
    is_green, T_FORBID_END, ALNS_ITER, SIGMA, REACT_FACTOR, DESTROY_RATE, SA_T0, SA_COOL, ...
    'carbon', inf);
fprintf('  碳排最优: 总成本=%.2f 元, 碳排=%.2f kg\n', cost_carb_opt, carb_min);

pareto_costs   = [cost_opt; cost_carb_opt];
pareto_carbons = [carb_cost_opt; carb_min];
pareto_sols    = {sol_cost_opt; sol_carb_opt};
pareto_details = [det_cost; det_carb];

eps_vals = linspace(carb_min, carb_cost_opt, N_EPSILON+2);
eps_vals = eps_vals(2:end-1);

prev_sol = sol_cost_opt;
for ep_i = 1:N_EPSILON
    fprintf('>>> Pareto扫描 %d/%d: 碳排上界=%.1f kg\n', ep_i, N_EPSILON, eps_vals(ep_i));
    [s_ep, c_ep, d_ep, cb_ep, ~] = ...
        run_alns_q2(n_cust, demand_w, demand_v, tw, dist_matrix, coords, VT, ...
        CARBON_PRICE, EARLY_PEN, LATE_PEN, SERVICE_TIME, T_START, ...
        is_green, T_FORBID_END, ALNS_ITER, SIGMA, REACT_FACTOR, DESTROY_RATE, SA_T0, SA_COOL, ...
        'cost', eps_vals(ep_i), prev_sol);
    pareto_costs(end+1)   = c_ep;
    pareto_carbons(end+1) = cb_ep;
    pareto_sols{end+1}    = s_ep;
    pareto_details(end+1,:) = d_ep;
    prev_sol = s_ep;
    fprintf('  => 成本=%.2f, 碳排=%.2f\n', c_ep, cb_ep);
end

%% 6. Pareto前沿提取与膝点检测
[pareto_carbons_s, si] = sort(pareto_carbons);
pareto_costs_s   = pareto_costs(si);
pareto_details_s = pareto_details(si,:);
pareto_sols_s    = pareto_sols(si);

pf_mask = true(size(pareto_costs_s));
for i = 1:length(pareto_costs_s)
    for j = 1:length(pareto_costs_s)
        if i~=j && pareto_costs_s(j)<=pareto_costs_s(i) && pareto_carbons_s(j)<=pareto_carbons_s(i)
            if pareto_costs_s(j)<pareto_costs_s(i) || pareto_carbons_s(j)<pareto_carbons_s(i)
                pf_mask(i) = false; break;
            end
        end
    end
end
pf_costs   = pareto_costs_s(pf_mask);
pf_carbons = pareto_carbons_s(pf_mask);
pf_details = pareto_details_s(pf_mask,:);
pf_sols    = pareto_sols_s(pf_mask);

knee_idx = 1;
if length(pf_costs) >= 3
    cn = (pf_costs - min(pf_costs)) / max(max(pf_costs)-min(pf_costs), 1);
    en = (pf_carbons - min(pf_carbons)) / max(max(pf_carbons)-min(pf_carbons), 1);
    p1 = [en(1), cn(1)]; p2 = [en(end), cn(end)];
    lv = p2 - p1; ll = norm(lv);
    md = 0;
    for i = 2:length(pf_costs)-1
        pt = [en(i), cn(i)];
        d = abs(det([lv; pt-p1])) / ll;
        if d > md, md = d; knee_idx = i; end
    end
end

best_sol    = pf_sols{knee_idx};
best_cost   = pf_costs(knee_idx);
best_carbon = pf_carbons(knee_idx);
best_detail = pf_details(knee_idx,:);

%% 7. 结果统计
fprintf('\n============ 问题二最终结果 (膝点方案) ============\n');
fprintf('最优总成本: %.2f 元\n', best_cost);
fprintf('碳排放量: %.2f kg\n', best_carbon);
fprintf('成本构成: 固定=%.1f 能耗=%.1f 碳排=%.1f 早到=%.1f 晚到=%.1f\n', ...
    best_detail(1), best_detail(2), best_detail(3), best_detail(4), best_detail(5));
fprintf('派车总数: %d 辆\n', length(best_sol));

vtype_count = zeros(1, n_vtype);
for r = 1:length(best_sol)
    vtype_count(best_sol(r).vtype) = vtype_count(best_sol(r).vtype)+1;
end
vt_names = {'燃油大(3t)','燃油中(1.5t)','燃油小(1.25t)','电动大(3t)','电动小(1.25t)'};
for v = 1:n_vtype
    if vtype_count(v)>0, fprintf('  %s: %d 辆\n', vt_names{v}, vtype_count(v)); end
end

green_by_ev = 0; green_by_fuel = 0;
for r = 1:length(best_sol)
    vt = best_sol(r).vtype;
    for c = best_sol(r).path
        if is_green(c)
            if VT(vt,7)==0, green_by_ev=green_by_ev+1; else, green_by_fuel=green_by_fuel+1; end
        end
    end
end
fprintf('时段分流统计: 电车送绿色区 %d 单, 燃油16:00后补位 %d 单\n', green_by_ev, green_by_fuel);
fprintf('================================================\n');

export_q2_excel('问题二_环保政策调度方案.xlsx', best_sol, best_cost, best_detail, best_carbon, ...
    pareto_costs, pareto_carbons, VT, demand_w, demand_v, tw, dist_matrix, parent_map, ...
    is_green, T_FORBID_END, T_START, SERVICE_TIME, CARBON_PRICE, vtype_count, ...
    green_by_ev, green_by_fuel, vt_names);

%% 8. 可视化

% 图1: Pareto前沿与膝点
figure('Color','w','Position',[50 100 720 500]); hold on; grid on;
scatter(pareto_carbons, pareto_costs, 50, [.7 .7 .7], 'filled', 'DisplayName','候选解');
plot(pf_carbons, pf_costs, '-o', 'Color',[0 0.45 0.74], 'LineWidth',2, ...
    'MarkerFaceColor',[0 0.45 0.74], 'MarkerSize',7, 'DisplayName','Pareto前沿');
scatter(pf_carbons(knee_idx), pf_costs(knee_idx), 200, 'rp', 'filled', ...
    'MarkerEdgeColor','k', 'DisplayName', sprintf('膝点(%.0f元,%.0fkg)', best_cost, best_carbon));
text(best_carbon+2, best_cost, sprintf(' 膝点\n 成本=%.0f元\n 碳排=%.0fkg', best_cost, best_carbon), ...
    'FontSize',10, 'FontWeight','bold', 'Color','r');
xlabel('碳排放量 (kg)'); ylabel('总配送成本 (元)');
title('\epsilon-约束法 Pareto前沿与膝点','FontWeight','bold','FontSize',13);
legend('Location','northeast'); set(gca,'FontName','Microsoft YaHei');

% 图2: 绿色区路径拓扑 (时段分流)
figure('Color','w','Position',[100 80 900 700]); hold on; grid on;
theta = linspace(0,2*pi,300);
fill(R_GREEN*cos(theta), R_GREEN*sin(theta), [0.92 0.98 0.92], ...
    'EdgeColor',[0 0.6 0], 'LineStyle','--', 'LineWidth',1.5, 'DisplayName','10km绿色配送区');
scatter(cx(is_green), cy(is_green), 60, 'o', 'MarkerEdgeColor',[0 0.6 0], ...
    'LineWidth',1.3, 'DisplayName','绿色区客户');
scatter(cx(~is_green), cy(~is_green), 30, 'k', '.', 'DisplayName','常规区客户');
plot(0, 0, 'rp', 'MarkerSize',18, 'MarkerFaceColor','y', 'DisplayName','配送中心');
clrs = lines(length(best_sol));
for r = 1:length(best_sol)
    vt = best_sol(r).vtype; path = best_sol(r).path;
    seq = [0, path, 0] + 1;
    X = coords(seq,2); Y = coords(seq,3);
    if VT(vt,7)==1
        plot(X, Y, '--', 'Color', clrs(r,:)*0.7+0.3, 'LineWidth',1);
    else
        plot(X, Y, '-', 'Color', clrs(r,:), 'LineWidth',1.5);
    end
end
h1 = plot(nan,nan,'-','Color',[0.1 0.7 0.2],'LineWidth',2);
h2 = plot(nan,nan,'--','Color',[0.8 0.3 0.1],'LineWidth',1.5);
legend([h1,h2], {'新能源车路径','燃油车路径'}, 'Location','northeastoutside');
xlabel('X (km)'); ylabel('Y (km)'); axis equal;
title('环保政策下时段分流调度拓扑','FontWeight','bold','FontSize',13);
set(gca,'FontName','Microsoft YaHei');

% 图3: 成本构成饼图
figure('Color','w','Position',[820 100 550 450]);
labels = {sprintf('固定启动\n%.2f元',best_detail(1)), sprintf('能耗成本\n%.2f元',best_detail(2)), ...
    sprintf('碳排成本\n%.2f元',best_detail(3)), sprintf('早到等待\n%.2f元',best_detail(4)), ...
    sprintf('晚到惩罚\n%.2f元',best_detail(5))};
pie(max(best_detail,0.01), [0 0 0 0.1 0.15], labels);
title(sprintf('膝点方案成本结构 (总计%.2f元)',best_cost),'FontWeight','bold','FontSize',13);
set(gca,'FontName','Microsoft YaHei');

% 图4: 车型使用柱状图
figure('Color','w','Position',[500 200 650 420]);
bar(vtype_count, 'FaceColor','flat', ...
    'CData',[.85 .33 .1;.93 .69 .13;1 .85 .2;.2 .7 .3;.1 .5 .6]);
set(gca,'XTickLabel',vt_names,'FontName','Microsoft YaHei');
ylabel('使用数量 (辆)');
title('膝点方案各车型使用情况','FontWeight','bold','FontSize',13);
grid on;

% 图5: ALNS收敛曲线 (成本最优求解过程)
figure('Color','w','Position',[50 50 700 420]); hold on; grid on;
plot(1:length(hist_cost), hist_cost, 'Color',[0 0.45 0.74], 'LineWidth',1.8);
xlabel('迭代次数'); ylabel('最优成本 (元)');
title('ALNS收敛曲线 (成本最优锚点)','FontWeight','bold','FontSize',13);
set(gca,'FontName','Microsoft YaHei');

% 图6: 灵敏度分析 — 碳价对车队结构的影响
cp_range = [0.3, 0.5, 0.65, 0.8, 1.0, 1.5];
fuel_ratio = zeros(size(cp_range));
ev_ratio   = zeros(size(cp_range));
for ci = 1:length(cp_range)
    [~,~,~,vt_cnt_sa] = eval_solution_q2( ...
        best_sol, dist_matrix, VT, demand_w, tw, cp_range(ci), EARLY_PEN, LATE_PEN, ...
        SERVICE_TIME, T_START, is_green, T_FORBID_END);
    fuel_ratio(ci) = sum(vt_cnt_sa(1:3));
    ev_ratio(ci)   = sum(vt_cnt_sa(4:5));
end
figure('Color','w','Position',[300 150 650 420]);
bar(categorical(arrayfun(@(x) sprintf('%.2f',x), cp_range, 'Uni',0)), ...
    [fuel_ratio; ev_ratio]', 'stacked');
legend({'燃油车','新能源车'}, 'Location','northwest');
xlabel('碳排放单价 (元/kg)'); ylabel('车辆数 (辆)');
title('碳价灵敏度分析 — 车队结构变化','FontWeight','bold','FontSize',13);
set(gca,'FontName','Microsoft YaHei'); grid on;

disp('>>> 全部图表已生成。');

%% ==================== 辅助函数 ====================

function export_q2_excel(file, sol, best_cost, best_detail, best_carbon, pareto_costs, pareto_carbons, VT, dw, dv, tw, D, parent_map, is_green, t_forbid_end, t0, svc, carbon_price, vtype_count, green_by_ev, green_by_fuel, vt_names)
    if exist(file,'file'), delete(file); end
    rows = {};
    for r = 1:length(sol)
        vt = sol(r).vtype; path = sol(r).path; cfg = VT(vt,:);
        path_orig = arrayfun(@(x) parent_map(x), path);
        path_str = '0';
        for p = path_orig, path_str = [path_str, sprintf('->%d', p)]; end
        path_str = [path_str, '->0'];
        [rc, det, carb, dist_sum, ret_t] = route_metrics_q2(path, vt, D, VT, dw, tw, svc, t0, carbon_price, is_green, t_forbid_end);
        green_cnt = sum(is_green(path));
        fuel_green_forbid = 0;
        if cfg(7)==1
            fuel_green_forbid = has_fuel_green_forbid_q2(path, vt, D, VT, tw, svc, t0, is_green, t_forbid_end);
        end
        start_t = t0;
        if cfg(7)==1 && any(is_green(path)), start_t = max(start_t, t_forbid_end); end
        rows(end+1,:) = {r, vt_names{vt}, cfg(7)==1, round(sum(dw(path))), round(cfg(1)), ...
            round(sum(dv(path))*10)/10, round(cfg(2)*10)/10, length(path), green_cnt, ...
            fuel_green_forbid, path_str, fmt_time_q2(start_t), fmt_time_q2(ret_t), round(dist_sum*10)/10, ...
            round(det(1)*100)/100, round(det(2)*100)/100, round(det(3)*100)/100, ...
            round(carb*100)/100, round(det(4)*100)/100, round(det(5)*100)/100, round(rc*100)/100};
    end
    h = {'车辆编号','车型','是否燃油','载货量kg','额定载重kg','载货体积m3','额定容积m3','服务客户数','绿色区客户数','是否违反限行','行驶路径','出发时间','返回时间','总距离km','固定成本元','能耗成本元','碳排成本元','碳排放量kg','早到惩罚元','晚到惩罚元','该车总成本元'};
    writetable(cell2table(rows,'VariableNames',h), file, 'Sheet','车辆使用方案与成本');

    detail = {}; row = 0;
    for r = 1:length(sol)
        vt = sol(r).vtype; path = sol(r).path; cfg = VT(vt,:);
        curr = 0; t = t0; load_now = sum(dw(path));
        if cfg(7)==1 && any(is_green(path)), t = max(t, t_forbid_end); end
        for k = 1:length(path)
            c = path(k); d = D(curr+1,c+1); dt = travel_time_tv(t,d); arr = t + dt;
            va = max(5, min(60, d/max(dt/60,0.001)));
            wait = 0; late = 0; status = '准时'; st = arr;
            if tw(c,1)>0 && arr<tw(c,1), wait=tw(c,1)-arr; st=tw(c,1); status=sprintf('早到等%.0fmin',wait); end
            if tw(c,2)>0 && arr>tw(c,2), late=arr-tw(c,2); status=sprintf('迟到%.0fmin',late); end
            forbid = is_green(c) && cfg(7)==1 && arr < t_forbid_end;
            if cfg(7)==1, bc=(0.0025*va^2-0.2554*va+31.75)/100; else, bc=(0.0014*va^2-0.12*va+36.19)/100; end
            ac = bc*d*(1+cfg(4)*min(load_now/cfg(1),1.5));
            row = row+1;
            detail(row,:) = {r, vt_names{vt}, k, parent_map(c), is_green(c), forbid, round(d*10)/10, round(va*10)/10, fmt_time_q2(arr), sprintf('%s-%s',fmt_time_q2(tw(c,1)),fmt_time_q2(tw(c,2))), status, round(wait), round(late), fmt_time_q2(st), fmt_time_q2(st+svc), round(dw(c)), round(load_now-dw(c)), round(ac*cfg(6)*100)/100, round(ac*cfg(5)*100)/100};
            t = st + svc; curr = c; load_now = load_now - dw(c);
        end
        db = D(curr+1,1); dtb = travel_time_tv(t,db); ret = t + dtb;
        row = row+1;
        detail(row,:) = {r, vt_names{vt}, length(path)+1, 0, false, false, round(db*10)/10, round((db/max(dtb/60,0.001))*10)/10, fmt_time_q2(ret), '-', '返回仓库', 0, 0, '-', fmt_time_q2(ret), 0, 0, 0, 0};
    end
    h2 = {'车辆编号','车型','访问序号','客户编号','是否绿色区','是否违反限行','本段距离km','本段均速kmh','到达时间','客户时间窗','时间窗状态','等待min','迟到min','开始服务','离开时间','配送量kg','车上剩余kg','本段能耗成本元','本段碳排放kg'};
    writetable(cell2table(detail,'VariableNames',h2), file, 'Sheet','逐站点行程明细');

    q1_cost = NaN; q1_carbon = NaN; q1_veh = NaN;
    if exist('问题一_完整调度方案.xlsx','file')
        try
            q1 = readtable('问题一_完整调度方案.xlsx','Sheet','总体成本汇总','VariableNamingRule','preserve');
            q1_cost = q1{strcmp(q1{:,1},'总配送成本(元)'),2};
            q1_carbon = q1{strcmp(q1{:,1},'碳排放总量(kg)'),2};
            q1_veh = q1{strcmp(q1{:,1},'派车总数(辆)'),2};
        catch
        end
    end
    summary = {'政策后总配送成本(元)',round(best_cost*100)/100; '政策后碳排放总量(kg)',round(best_carbon*100)/100; '政策后派车总数(辆)',length(sol); '问题一总配送成本(元)',q1_cost; '问题一碳排放总量(kg)',q1_carbon; '问题一派车总数(辆)',q1_veh; '成本变化(元)',round((best_cost-q1_cost)*100)/100; '碳排变化(kg)',round((best_carbon-q1_carbon)*100)/100; '车辆数变化(辆)',length(sol)-q1_veh; '电车送绿色区客户数',green_by_ev; '燃油车16点后补位绿色区客户数',green_by_fuel; '燃油大(辆)',vtype_count(1); '燃油中(辆)',vtype_count(2); '燃油小(辆)',vtype_count(3); '电动大(辆)',vtype_count(4); '电动小(辆)',vtype_count(5); '固定启动成本(元)',round(best_detail(1)*100)/100; '能耗成本(元)',round(best_detail(2)*100)/100; '碳排成本(元)',round(best_detail(3)*100)/100; '早到等待惩罚(元)',round(best_detail(4)*100)/100; '晚到惩罚(元)',round(best_detail(5)*100)/100};
    writetable(cell2table(summary,'VariableNames',{'指标','数值'}), file, 'Sheet','政策影响与成本汇总');

    ptab = [(1:length(pareto_costs))', pareto_costs(:), pareto_carbons(:)];
    writetable(array2table(ptab,'VariableNames',{'候选序号','总成本元','碳排放kg'}), file, 'Sheet','Pareto候选解');
    fprintf('>>> 已导出: %s\n', file);
end

function [total, det, carbon, dist_sum, ret_t] = route_metrics_q2(path, vt, D, VT, dw, tw, svc, t0, cp, is_green, t_forbid_end)
    cfg = VT(vt,:); det = [cfg(3),0,0,0,0]; carbon = 0; dist_sum = 0; curr = 0; t = t0; load_now = sum(dw(path));
    if cfg(7)==1 && any(is_green(path)), t = max(t, t_forbid_end); end
    for k = 1:length(path)
        c = path(k); d = D(curr+1,c+1); dist_sum = dist_sum + d; dt = travel_time_tv(t,d); va = max(5,min(60,d/max(dt/60,0.001)));
        if cfg(7)==1, bc=(0.0025*va^2-0.2554*va+31.75)/100; else, bc=(0.0014*va^2-0.12*va+36.19)/100; end
        ac = bc*d*(1+cfg(4)*min(load_now/cfg(1),1.5)); det(2)=det(2)+ac*cfg(6); segc=ac*cfg(5); det(3)=det(3)+segc*cp; carbon=carbon+segc;
        t = t + dt; if tw(c,1)>0 && t<tw(c,1), det(4)=det(4)+(tw(c,1)-t)*20/60; t=tw(c,1); end; if tw(c,2)>0 && t>tw(c,2), det(5)=det(5)+(t-tw(c,2))*50/60; end
        t = t + svc; curr = c; load_now = load_now - dw(c);
    end
    db = D(curr+1,1); dist_sum = dist_sum + db; dtb = travel_time_tv(t,db); va=max(5,min(60,db/max(dtb/60,0.001)));
    if cfg(7)==1, bc=(0.0025*va^2-0.2554*va+31.75)/100; else, bc=(0.0014*va^2-0.12*va+36.19)/100; end
    det(2)=det(2)+bc*db*cfg(6); segc=bc*db*cfg(5); det(3)=det(3)+segc*cp; carbon=carbon+segc; ret_t=t+dtb; total=sum(det);
end

function yes = has_fuel_green_forbid_q2(path, vt, D, VT, tw, svc, t0, is_green, t_forbid_end)
    yes = false; curr=0; t=t0;
    if VT(vt,7)==1 && any(is_green(path)), t = max(t, t_forbid_end); end
    for k=1:length(path)
        c=path(k); t=t+travel_time_tv(t,D(curr+1,c+1));
        if is_green(c) && VT(vt,7)==1 && t<t_forbid_end, yes=true; return; end
        if tw(c,1)>0 && t<tw(c,1), t=tw(c,1); end
        t=t+svc; curr=c;
    end
end

function s = fmt_time_q2(t)
    s = sprintf('%d:%02d', floor(t/60), round(mod(t,60)));
end

function [best_sol, best_cost, best_detail, best_carbon, best_hist, vt_cnt] = ...
    run_alns_q2(n_cust, dw, dv, tw, D, coords, VT, cp, ep, lp, svc, t0, ...
    is_green, t_forbid_end, max_iter, SIGMA, react, dest_rate, sa_t0, sa_cool, ...
    obj_mode, carbon_limit, warm_start)

    if nargin >= 23 && ~isempty(warm_start)
        sol = warm_start;
    else
        sol = greedy_init_q2(n_cust, dw, dv, tw, D, coords, VT, t0, svc, is_green, t_forbid_end);
    end
    [sol_cost, sol_det, sol_carb] = eval_solution_q2(sol, D, VT, dw, tw, cp, ep, lp, svc, t0, is_green, t_forbid_end);

    if strcmp(obj_mode, 'carbon')
        sol_obj = sol_carb;
    else
        sol_obj = sol_cost;
    end

    best_sol = sol; best_cost = sol_cost; best_detail = sol_det;
    best_carbon = sol_carb; best_obj = sol_obj;
    curr_sol = sol; curr_obj = sol_obj;

    w_d = ones(1,3); w_r = ones(1,2);
    sc_d = zeros(1,3); sc_r = zeros(1,2);
    ct_d = zeros(1,3); ct_r = zeros(1,2);
    sa_t = sa_t0;

    best_hist = zeros(max_iter, 1);

    for iter = 1:max_iter
        d_op = roulette_sel(w_d);
        r_op = roulette_sel(w_r);
        n_rm = max(1, round(n_cust * (dest_rate(1) + rand()*(dest_rate(2)-dest_rate(1)))));

        [removed, partial] = destroy_q2(curr_sol, d_op, n_rm, D, dw, tw);
        try
            new_sol = repair_q2(partial, removed, r_op, D, VT, dw, dv, tw, cp, ep, lp, svc, t0, is_green, t_forbid_end);
        catch
            reward = 0;
            sc_d(d_op) = sc_d(d_op)+reward; ct_d(d_op) = ct_d(d_op)+1;
            sc_r(r_op) = sc_r(r_op)+reward; ct_r(r_op) = ct_r(r_op)+1;
            sa_t = sa_t * sa_cool;
            best_hist(iter) = best_cost;
            continue;
        end
        [new_cost, new_det, new_carb] = eval_solution_q2(new_sol, D, VT, dw, tw, cp, ep, lp, svc, t0, is_green, t_forbid_end);

        % 库存与容量约束检查
        feasible_ok = check_solution_feasible_q2(new_sol, VT, dw, dv);

        if strcmp(obj_mode, 'carbon'), new_obj = new_carb; else, new_obj = new_cost; end

        if new_carb > carbon_limit
            excess = new_carb - carbon_limit;
            new_obj = new_obj + 500 * excess + 10 * excess^2;
        end

        reward = 0;
        if feasible_ok && new_obj < best_obj && new_carb <= carbon_limit * 1.02
            best_sol = new_sol; best_cost = new_cost; best_detail = new_det;
            best_carbon = new_carb; best_obj = new_obj;
            curr_sol = new_sol; curr_obj = new_obj;
            reward = SIGMA(1);
        elseif feasible_ok && new_obj < curr_obj
            curr_sol = new_sol; curr_obj = new_obj;
            reward = SIGMA(2);
        elseif feasible_ok && rand() < exp(-(new_obj - curr_obj)/sa_t)
            curr_sol = new_sol; curr_obj = new_obj;
            reward = SIGMA(3);
        end

        sc_d(d_op) = sc_d(d_op)+reward; ct_d(d_op) = ct_d(d_op)+1;
        sc_r(r_op) = sc_r(r_op)+reward; ct_r(r_op) = ct_r(r_op)+1;
        sa_t = sa_t * sa_cool;
        best_hist(iter) = best_cost;

        if mod(iter,100)==0
            for d=1:3, if ct_d(d)>0, w_d(d)=w_d(d)*(1-react)+react*sc_d(d)/ct_d(d); end; end
            for r=1:2, if ct_r(r)>0, w_r(r)=w_r(r)*(1-react)+react*sc_r(r)/ct_r(r); end; end
            w_d = max(w_d,0.05); w_r = max(w_r,0.05);
            sc_d(:)=0; ct_d(:)=0; sc_r(:)=0; ct_r(:)=0;
        end

        if mod(iter,500)==0
            fprintf('    iter %d/%d | best_cost=%.1f | best_carb=%.1f\n', iter, max_iter, best_cost, best_carbon);
        end
    end

    vt_cnt = zeros(1, size(VT,1));
    for r=1:length(best_sol), vt_cnt(best_sol(r).vtype)=vt_cnt(best_sol(r).vtype)+1; end
end

function sol = greedy_init_q2(n_cust, dw, dv, tw, D, coords, VT, t0, svc, is_green, t_forbid_end)
    sol = struct('path',{},'vtype',{});
    stock = VT(:,8)';
    unvisited = 1:n_cust;

    while ~isempty(unvisited)
        avail_types = find(stock > 0);
        if isempty(avail_types)
            error('车辆库存不足，无法为所有客户分配车辆');
        end

        feasible_count = zeros(size(unvisited));
        demand_score = zeros(size(unvisited));
        green_score = zeros(size(unvisited));
        for ii = 1:length(unvisited)
            c = unvisited(ii);
            feasible_count(ii) = sum(dw(c) <= VT(avail_types,1) & dv(c) <= VT(avail_types,2));
            demand_score(ii) = dw(c)/max(VT(:,1)) + dv(c)/max(VT(:,2));
            green_score(ii) = is_green(c);
        end
        if any(feasible_count == 0)
            bad = unvisited(find(feasible_count==0,1));
            error('车辆库存不足或剩余车辆容量不足，客户%d无法分配: %.1fkg, %.2fm3', bad, dw(bad), dv(bad));
        end
        [~, ord] = sortrows([feasible_count(:), -green_score(:), -demand_score(:)], [1, 2, 3]);
        seed = unvisited(ord(1));

        candidate_vt = avail_types(dw(seed) <= VT(avail_types,1) & dv(seed) <= VT(avail_types,2));
        if is_green(seed)
            ev = candidate_vt(VT(candidate_vt,7)==0);
            if ~isempty(ev), candidate_vt = ev; end
        end
        [~, cap_ord] = sortrows([VT(candidate_vt,1), VT(candidate_vt,2), VT(candidate_vt,7)], [1, 2, 3]);
        best_vt = candidate_vt(cap_ord(1));

        path = seed;
        load_w = dw(seed); load_v = dv(seed); curr = seed; t = t0;
        dt0 = travel_time_tv(t, D(1, seed+1));
        t = t + dt0;
        if tw(seed,1)>0 && t<tw(seed,1), t=tw(seed,1); end
        t = t + svc;
        remain = setdiff(unvisited, seed);

        while true
            best_next = 0; best_inc = inf;
            for j = 1:length(remain)
                c = remain(j);
                if load_w+dw(c)>VT(best_vt,1) || load_v+dv(c)>VT(best_vt,2), continue; end
                if is_green(c) && VT(best_vt,7)==1
                    arr_est = t + travel_time_tv(t, D(curr+1,c+1));
                    if arr_est < t_forbid_end, continue; end
                end
                dt = travel_time_tv(t, D(curr+1,c+1));
                arr = t + dt;
                inc = D(curr+1,c+1);
                if tw(c,1)>0 && arr<tw(c,1), inc=inc+(tw(c,1)-arr)*0.5; end
                if tw(c,2)>0 && arr>tw(c,2), inc=inc+(arr-tw(c,2))*2; end
                if inc < best_inc, best_inc=inc; best_next=j; end
            end
            if best_next==0, break; end
            c = remain(best_next);
            dt = travel_time_tv(t, D(curr+1,c+1));
            t = t + dt;
            if tw(c,1)>0 && t<tw(c,1), t=tw(c,1); end
            t = t + svc;
            path = [path, c]; load_w=load_w+dw(c); load_v=load_v+dv(c);
            curr = c; remain(best_next) = [];
        end

        sol(end+1).path = path;
        sol(end).vtype = best_vt;
        stock(best_vt) = stock(best_vt) - 1;
        unvisited = setdiff(unvisited, path);
    end
end

function [total, detail, carbon_kg, vt_cnt] = eval_solution_q2(sol, D, VT, dw, tw, cp, ep, lp, svc, t0, is_green, t_forbid_end)
    c_fix=0; c_eng=0; c_carb=0; c_early=0; c_late=0; carbon_kg=0;
    n_vt = size(VT,1);
    vt_cnt = zeros(1, n_vt);
    penalty_green = 0;

    for r = 1:length(sol)
        vt = sol(r).vtype; path = sol(r).path;
        cfg = VT(vt,:);
        c_fix = c_fix + cfg(3);
        vt_cnt(vt) = vt_cnt(vt) + 1;
        curr = 0; t = t0;
        if cfg(7)==1 && any(is_green(path))
            t = max(t, t_forbid_end);
        end
        total_load = sum(dw(path)); cur_load = total_load;

        for k = 1:length(path)
            c = path(k);
            dist = D(curr+1, c+1);
            dt = travel_time_tv(t, dist);
            v_avg = dist / max(dt/60, 0.001);
            v_avg = max(5, min(60, v_avg));

            if cfg(7)==1
                base_cons = (0.0025*v_avg^2 - 0.2554*v_avg + 31.75)/100;
            else
                base_cons = (0.0014*v_avg^2 - 0.12*v_avg + 36.19)/100;
            end
            lr = min(cur_load/cfg(1), 1.5);
            ac = base_cons * dist * (1 + cfg(4)*lr);

            c_eng = c_eng + ac * cfg(6);
            seg_carbon = ac * cfg(5);
            c_carb = c_carb + seg_carbon * cp;
            carbon_kg = carbon_kg + seg_carbon;

            t = t + dt;

            if is_green(c) && cfg(7)==1 && t < t_forbid_end
                penalty_green = penalty_green + 5000;
            end

            if tw(c,1)>0 && t<tw(c,1)
                c_early = c_early + (tw(c,1)-t)*ep;
                t = tw(c,1);
            end
            if tw(c,2)>0 && t>tw(c,2)
                c_late = c_late + (t-tw(c,2))*lp;
            end
            t = t + svc; curr = c; cur_load = cur_load - dw(c);
        end

        dist_back = D(curr+1, 1);
        dt_b = travel_time_tv(t, dist_back);
        vb = dist_back / max(dt_b/60, 0.001); vb = max(5, min(60, vb));
        if cfg(7)==1, bcb=(0.0025*vb^2-0.2554*vb+31.75)/100;
        else, bcb=(0.0014*vb^2-0.12*vb+36.19)/100; end
        c_eng = c_eng + bcb*dist_back*cfg(6);
        seg_cb = bcb*dist_back*cfg(5);
        c_carb = c_carb + seg_cb*cp;
        carbon_kg = carbon_kg + seg_cb;
    end

    detail = [c_fix, c_eng, c_carb, c_early, c_late];
    total = sum(detail) + penalty_green;
end

function ok = check_solution_feasible_q2(sol, VT, dw, dv)
    ok = true;
    used = zeros(1, size(VT,1));
    for r = 1:length(sol)
        vt = sol(r).vtype;
        used(vt) = used(vt) + 1;
        if used(vt) > VT(vt,8)
            ok = false; return;
        end
        if sum(dw(sol(r).path)) > VT(vt,1) + 1e-6 || sum(dv(sol(r).path)) > VT(vt,2) + 1e-6
            ok = false; return;
        end
    end
end

function [cost, detail] = eval_single_route_q2(path, vt, D, VT, dw, tw, cp, ep, lp, svc, t0, is_green, t_forbid_end)
    cfg = VT(vt,:);
    c_eng=0; c_carb=0; c_early=0; c_late=0;
    curr=0; t=t0;
    if cfg(7)==1 && any(is_green(path))
        t = max(t, t_forbid_end);
    end
    total_load = sum(dw(path)); cur_load = total_load;
    pen_g = 0;

    for k=1:length(path)
        c = path(k); dist = D(curr+1, c+1);
        dt = travel_time_tv(t, dist);
        va = dist/max(dt/60,0.001); va = max(5,min(60,va));
        if cfg(7)==1, bc=(0.0025*va^2-0.2554*va+31.75)/100;
        else, bc=(0.0014*va^2-0.12*va+36.19)/100; end
        lr = min(cur_load/cfg(1), 1.5);
        ac = bc*dist*(1+cfg(4)*lr);
        c_eng = c_eng + ac*cfg(6);
        c_carb = c_carb + ac*cfg(5)*cp;
        t = t+dt;
        if is_green(c) && cfg(7)==1 && t < t_forbid_end
            pen_g = pen_g + 5000;
        end
        if tw(c,1)>0 && t<tw(c,1), c_early=c_early+(tw(c,1)-t)*ep; t=tw(c,1); end
        if tw(c,2)>0 && t>tw(c,2), c_late=c_late+(t-tw(c,2))*lp; end
        t=t+svc; curr=c; cur_load=cur_load-dw(c);
    end
    db=D(curr+1,1); dtb=travel_time_tv(t,db);
    vb=db/max(dtb/60,0.001); vb=max(5,min(60,vb));
    if cfg(7)==1, bcb=(0.0025*vb^2-0.2554*vb+31.75)/100;
    else, bcb=(0.0014*vb^2-0.12*vb+36.19)/100; end
    c_eng=c_eng+bcb*db*cfg(6); c_carb=c_carb+bcb*db*cfg(5)*cp;
    detail = [cfg(3), c_eng, c_carb, c_early, c_late];
    cost = sum(detail) + pen_g;
end

function [removed, partial] = destroy_q2(sol, op, n_rm, D, dw, tw)
    all_c = [];
    for r=1:length(sol), all_c = [all_c, sol(r).path]; end
    switch op
        case 1
            idx = randperm(length(all_c), min(n_rm, length(all_c)));
            removed = all_c(idx);
        case 2
            scores = zeros(size(all_c));
            for i=1:length(all_c)
                c = all_c(i);
                scores(i) = dw(c) / max(tw(c,2)-tw(c,1), 1);
            end
            [~, idx] = sort(scores, 'descend');
            removed = all_c(idx(1:min(n_rm, length(idx))));
        case 3
            seed = all_c(randi(length(all_c)));
            dists = zeros(size(all_c));
            for i=1:length(all_c)
                c = all_c(i);
                dists(i) = D(seed+1, c+1) + abs(dw(seed)-dw(c))*0.1;
            end
            [~, idx] = sort(dists);
            removed = all_c(idx(1:min(n_rm, length(idx))));
    end
    partial = sol;
    for r = 1:length(partial)
        partial(r).path = setdiff(partial(r).path, removed, 'stable');
    end
    partial = partial(arrayfun(@(s) ~isempty(s.path), partial));
end

function sol = repair_q2(partial, removed, op, D, VT, dw, dv, tw, cp, ep, lp, svc, t0, is_green, t_forbid_end)
    sol = partial;
    stock_used = zeros(1, size(VT,1));
    for r=1:length(sol), stock_used(sol(r).vtype) = stock_used(sol(r).vtype)+1; end

    to_ins = removed(randperm(length(removed)));
    for ii = 1:length(to_ins)
        c = to_ins(ii);
        best_inc = inf; best_r = 0; best_pos = 0;

        for r = 1:length(sol)
            vt = sol(r).vtype; path = sol(r).path;
            cw = sum(dw(path)); cv = sum(dv(path));
            if cw+dw(c)>VT(vt,1) || cv+dv(c)>VT(vt,2), continue; end

            if is_green(c) && VT(vt,7)==1
                arr_est = estimate_arrival(path, 0, D, t0, svc, tw);
                if arr_est < t_forbid_end, continue; end
            end

            for pos = 1:length(path)+1
                np = [path(1:pos-1), c, path(pos:end)];
                [nc,~] = eval_single_route_q2(np, vt, D, VT, dw, tw, cp, ep, lp, svc, t0, is_green, t_forbid_end);
                [oc,~] = eval_single_route_q2(path, vt, D, VT, dw, tw, cp, ep, lp, svc, t0, is_green, t_forbid_end);
                inc = nc - oc;
                if op==2, inc = inc - 0.1*rand()*abs(inc); end
                if inc < best_inc, best_inc=inc; best_r=r; best_pos=pos; end
            end
        end

        if best_r > 0 && best_inc < 1e6
            sol(best_r).path = [sol(best_r).path(1:best_pos-1), c, sol(best_r).path(best_pos:end)];
        else
            vt = pick_vehicle_q2(c, dw(c), dv(c), VT, stock_used, is_green, t_forbid_end);
            sol(end+1).path = c;
            sol(end).vtype = vt;
            stock_used(vt) = stock_used(vt) + 1;
        end
    end
end

function vt = pick_vehicle_q2(c, w, v, VT, used, is_green, t_forbid_end)
    if is_green(c)
        order = [5, 4];
    else
        order = [5, 4, 3, 2, 1];
    end
    for k = order
        if w<=VT(k,1) && v<=VT(k,2) && used(k)<VT(k,8)
            vt = k; return;
        end
    end
    for k = [3, 2, 1]
        if w<=VT(k,1) && v<=VT(k,2) && used(k)<VT(k,8)
            vt = k; return;
        end
    end
    error('车辆库存不足或剩余车辆容量不足，无法分配 %.1fkg, %.2fm3 的需求', w, v);
end

function arr = estimate_arrival(path, pos_idx, D, t0, svc, tw)
    t = t0; curr = 0;
    n = min(pos_idx, length(path));
    for k = 1:n
        c = path(k);
        dt = travel_time_tv(t, D(curr+1, c+1));
        t = t + dt;
        if tw(c,1)>0 && t<tw(c,1), t=tw(c,1); end
        t = t + svc;
        curr = c;
    end
    arr = t;
end

function idx = roulette_sel(w)
    cw = cumsum(w)/sum(w);
    idx = find(cw >= rand(), 1);
    if isempty(idx), idx = length(w); end
end

function tt = travel_time_tv(tStart, dist)
    ct = tStart; rd = dist;
    while rd > 0.01
        if (ct>=540 && ct<600) || (ct>=780 && ct<900), v=55.3/60;
        elseif (ct>=480 && ct<540) || (ct>=690 && ct<780), v=9.8/60;
        else, v=35.4/60; end
        nxt = [480,540,600,690,780,900,1020,1440];
        ni = find(nxt > ct, 1);
        if isempty(ni), nt=1440; else, nt=nxt(ni); end
        can_go = v*(nt-ct);
        if rd <= can_go, ct=ct+rd/v; rd=0;
        else, rd=rd-can_go; ct=nt; end
        if ct>=1440, break; end
    end
    tt = max(0.1, ct-tStart);
end
