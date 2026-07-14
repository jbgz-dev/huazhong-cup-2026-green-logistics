% =========================================================================
% 华中杯 A题 问题三：动态事件下的实时车辆调度策略
% 模型: 事件触发局部动态重调度 (Event-Triggered Local Rescheduling)
% 算法: 锁定未受影响车辆 + 受影响路径ALNS局部重优化
% =========================================================================
clear; clc; close all; warning off; rng(1);

%% 1. 数据预处理 (与问题二一致)
disp('===== 问题三：动态事件触发局部重调度引擎启动 =====');

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

%% 1.5 需求拆分
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

cx = coords(2:end, 2);
cy = coords(2:end, 3);
dist_to_center = sqrt(cx.^2 + cy.^2);
R_GREEN = 10;
is_green = dist_to_center <= R_GREEN;

fprintf('数据加载完成: %d 客户\n', n_cust);

%% 2. 车型与参数定义
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
T_START      = 480;
T_FORBID_END = 960;

%% 3. 生成问题二基准静态方案 (作为动态调度的初始状态)
disp('>>> 生成静态基准方案 (问题二最优解)...');
base_sol = greedy_init_q2(n_cust, demand_w, demand_v, tw, dist_matrix, coords, VT, ...
    T_START, SERVICE_TIME, is_green, T_FORBID_END);

ALNS_ITER_BASE = 1500;
SA_T0 = 1000; SA_COOL = 0.9995;
SIGMA = [33, 9, 3]; REACT_FACTOR = 0.1;
DESTROY_RATE = [0.1, 0.4];

[base_sol, base_cost, base_det, base_carb] = run_alns_local( ...
    base_sol, n_cust, demand_w, demand_v, tw, dist_matrix, VT, ...
    CARBON_PRICE, EARLY_PEN, LATE_PEN, SERVICE_TIME, T_START, ...
    is_green, T_FORBID_END, ALNS_ITER_BASE, SIGMA, REACT_FACTOR, ...
    DESTROY_RATE, SA_T0, SA_COOL, 1:n_cust);

fprintf('静态基准: 总成本=%.2f 元, 碳排=%.2f kg, %d 辆车\n', ...
    base_cost, base_carb, length(base_sol));

%% 4. 定义动态突发事件场景
% 事件触发时刻: 12:00 (720min) — 车辆已执行部分配送
T_EVENT = 720;

% 原始客户编号 -> 虚拟客户编号映射
orig2virt = cell(orig_n_cust, 1);
for vi = 1:n_cust
    oi = parent_map(vi);
    orig2virt{oi} = [orig2virt{oi}, vi];
end

% 场景1: 订单取消 — 客户15、客户32取消订单
event1.type = 'cancel';
event1.name = '订单取消';
event1.custs = [orig2virt{15}, orig2virt{32}];
event1.custs_orig = [15, 32];
event1.time = T_EVENT;

% 场景2: 新增紧急订单 — 绿色区内2个新客户(复用绿色区有订单客户的坐标模拟)
event2.type = 'new_order';
event2.name = '新增紧急订单';
green_orig = find(sqrt(orig_coords(2:end,2).^2 + orig_coords(2:end,3).^2) <= 10 & ~cellfun(@isempty, orig2virt));
event2_orig = green_orig(1:min(2, length(green_orig)))';
event2.custs = arrayfun(@(x) orig2virt{x}(1), event2_orig);
event2.custs_orig = event2_orig;
event2.new_tw = [750, 900; 780, 960];
event2.new_dw = [200; 350];
event2.new_dv = [1.5; 2.0];
event2.time = T_EVENT;

% 场景3: 配送地址变更 — 客户40地址变更(距离修正)
event3.type = 'address_change';
event3.name = '配送地址变更';
event3.custs = orig2virt{40};
event3.custs_orig = [40];
event3.time = T_EVENT;

% 场景4: 时间窗调整 — 客户10、客户55时间窗收紧
event4.type = 'tw_change';
event4.name = '时间窗收紧';
event4.custs = [orig2virt{10}, orig2virt{55}];
event4.custs_orig = [10, 55];
n_v10 = length(orig2virt{10}); n_v55 = length(orig2virt{55});
event4.new_tw = [repmat([720, 780], n_v10, 1); repmat([750, 840], n_v55, 1)];
event4.time = T_EVENT;

events = {event1, event2, event3, event4};
n_events = length(events);

%% 5. 模拟车辆执行状态 (截至事件触发时刻)
disp('>>> 模拟车辆执行状态至事件触发时刻...');

[route_status, delivered, in_transit] = simulate_execution( ...
    base_sol, dist_matrix, tw, demand_w, T_START, SERVICE_TIME, T_EVENT);

n_delivered = sum(cellfun(@length, {route_status.delivered}));
n_pending = sum(cellfun(@length, {route_status.pending}));
fprintf('事件触发时刻 %.0f min: 已送达 %d 单, 待配送 %d 单\n', ...
    T_EVENT, n_delivered, n_pending);

%% 6. 逐场景执行局部动态重调度
results = struct('name',{},'cost',{},'carbon',{},'detail',{},'n_vehicles',{},...
    'n_affected',{},'sol',{},'delta_cost',{},'delta_carbon',{},...
    'tw',{},'dw',{},'dv',{},'dist',{});

ALNS_ITER_LOCAL = 800;

for ei = 1:n_events
    evt = events{ei};
    fprintf('\n>>> 场景%d: %s (触发时刻=%.0fmin)\n', ei, evt.name, evt.time);

    [affected_routes, modified_tw, modified_dw, modified_dv, modified_dist, active_custs] = ...
        apply_event(evt, base_sol, route_status, tw, demand_w, demand_v, ...
        dist_matrix, n_cust, is_green);

    locked_sol = [];
    replan_custs = [];
    for r = 1:length(base_sol)
        if ~ismember(r, affected_routes)
            locked_sol = [locked_sol, base_sol(r)];
        else
            replan_custs = [replan_custs, route_status(r).pending];
        end
    end

    if strcmp(evt.type, 'cancel')
        replan_custs = setdiff(replan_custs, evt.custs);
    end

    if strcmp(evt.type, 'new_order')
        for nc = 1:length(evt.custs)
            c = evt.custs(nc);
            if ~ismember(c, replan_custs)
                replan_custs = [replan_custs, c];
            end
            modified_tw(c,:) = evt.new_tw(nc,:);
            modified_dw(c) = evt.new_dw(nc);
            modified_dv(c) = evt.new_dv(nc);
        end
    end

    if strcmp(evt.type, 'tw_change')
        for tc = 1:length(evt.custs)
            c = evt.custs(tc);
            modified_tw(c,:) = evt.new_tw(tc,:);
        end
    end

    replan_custs = unique(replan_custs);
    fprintf('  受影响路径: %d 条, 待重规划客户: %d 个\n', ...
        length(affected_routes), length(replan_custs));

    if isempty(replan_custs)
        replan_sol = [];
    else
        replan_init = greedy_init_partial(replan_custs, modified_dw, modified_dv, ...
            modified_tw, modified_dist, VT, evt.time, SERVICE_TIME, is_green, T_FORBID_END);

        [replan_sol, ~, ~, ~] = run_alns_local( ...
            replan_init, n_cust, modified_dw, modified_dv, modified_tw, modified_dist, VT, ...
            CARBON_PRICE, EARLY_PEN, LATE_PEN, SERVICE_TIME, evt.time, ...
            is_green, T_FORBID_END, ALNS_ITER_LOCAL, SIGMA, REACT_FACTOR, ...
            DESTROY_RATE, SA_T0, SA_COOL, replan_custs);
    end

    final_sol = [locked_sol, replan_sol];
    [final_cost, final_det, final_carb] = eval_solution_q3( ...
        final_sol, modified_dist, VT, modified_dw, modified_tw, ...
        CARBON_PRICE, EARLY_PEN, LATE_PEN, SERVICE_TIME, T_START, is_green, T_FORBID_END);

    results(ei).name = evt.name;
    results(ei).cost = final_cost;
    results(ei).carbon = final_carb;
    results(ei).detail = final_det;
    results(ei).n_vehicles = length(final_sol);
    results(ei).n_affected = length(affected_routes);
    results(ei).sol = final_sol;
    results(ei).delta_cost = final_cost - base_cost;
    results(ei).delta_carbon = final_carb - base_carb;
    results(ei).tw = modified_tw;
    results(ei).dw = modified_dw;
    results(ei).dv = modified_dv;
    results(ei).dist = modified_dist;

    fprintf('  重调度结果: 成本=%.2f (Δ=%+.2f), 碳排=%.2f (Δ=%+.2f), %d辆车\n', ...
        final_cost, results(ei).delta_cost, final_carb, results(ei).delta_carbon, ...
        length(final_sol));
end

%% 7. 全局重规划对比 (方案一基准: 所有待配送订单全部重新规划)
disp('');
disp('>>> 全局重规划对比 (Global Rescheduling)...');
global_results = struct('name',{},'cost',{},'carbon',{});

for ei = 1:n_events
    evt = events{ei};
    all_pending = [];
    for r = 1:length(base_sol)
        all_pending = [all_pending, route_status(r).pending];
    end

    mod_tw = tw; mod_dw = demand_w; mod_dv = demand_v; mod_dist = dist_matrix;

    if strcmp(evt.type, 'cancel')
        all_pending = setdiff(all_pending, evt.custs);
    elseif strcmp(evt.type, 'new_order')
        for nc = 1:length(evt.custs)
            c = evt.custs(nc);
            if ~ismember(c, all_pending), all_pending = [all_pending, c]; end
            mod_tw(c,:) = evt.new_tw(nc,:);
            mod_dw(c) = evt.new_dw(nc);
            mod_dv(c) = evt.new_dv(nc);
        end
    elseif strcmp(evt.type, 'tw_change')
        for tc = 1:length(evt.custs)
            mod_tw(evt.custs(tc),:) = evt.new_tw(tc,:);
        end
    elseif strcmp(evt.type, 'address_change')
        for ci = 1:length(evt.custs)
            c = evt.custs(ci);
            mod_dist(c+1,:) = mod_dist(c+1,:) * 1.3;
            mod_dist(:,c+1) = mod_dist(:,c+1) * 1.3;
        end
    end

    all_pending = unique(all_pending);
    g_init = greedy_init_partial(all_pending, mod_dw, mod_dv, mod_tw, mod_dist, ...
        VT, evt.time, SERVICE_TIME, is_green, T_FORBID_END);
    [g_sol, g_cost, ~, g_carb] = run_alns_local(g_init, n_cust, mod_dw, mod_dv, ...
        mod_tw, mod_dist, VT, CARBON_PRICE, EARLY_PEN, LATE_PEN, SERVICE_TIME, ...
        evt.time, is_green, T_FORBID_END, ALNS_ITER_LOCAL, SIGMA, REACT_FACTOR, ...
        DESTROY_RATE, SA_T0, SA_COOL, all_pending);

    locked_cost = 0; locked_carb = 0;
    for r = 1:length(base_sol)
        delivered_path = route_status(r).delivered;
        if ~isempty(delivered_path)
            [rc, rd, rcb] = eval_solution_q3( ...
                struct('path',delivered_path,'vtype',base_sol(r).vtype), ...
                dist_matrix, VT, demand_w, tw, CARBON_PRICE, EARLY_PEN, LATE_PEN, ...
                SERVICE_TIME, T_START, is_green, T_FORBID_END);
            locked_cost = locked_cost + rc;
            locked_carb = locked_carb + rcb;
        end
    end

    global_results(ei).name = evt.name;
    global_results(ei).cost = g_cost + locked_cost;
    global_results(ei).carbon = g_carb + locked_carb;
    global_results(ei).n_vehicles = length(g_sol) + length(base_sol);
end

%% 8. 结果汇总输出
fprintf('\n\n============ 问题三最终结果汇总 ============\n');
fprintf('静态基准方案: 成本=%.2f 元, 碳排=%.2f kg\n\n', base_cost, base_carb);
fprintf('%-14s | %-10s | %-10s | %-8s | %-10s | %-10s | %-10s\n', ...
    '场景', '局部成本', '局部碳排', '车辆数', '全局成本', '局部Δ成本', '全局Δ成本');
fprintf('%s\n', repmat('-', 1, 90));
for ei = 1:n_events
    fprintf('%-12s | %10.1f | %10.1f | %8d | %10.1f | %+10.1f | %+10.1f\n', ...
        results(ei).name, results(ei).cost, results(ei).carbon, results(ei).n_vehicles, ...
        global_results(ei).cost, results(ei).delta_cost, ...
        global_results(ei).cost - base_cost);
end
fprintf('================================================\n');
fprintf('局部重调度平均额外成本: %.2f 元\n', mean([results.delta_cost]));
fprintf('全局重规划平均额外成本: %.2f 元\n', mean([global_results.cost]) - base_cost);
fprintf('局部策略稳定性优势: 平均仅影响 %.1f 条路径\n', mean([results.n_affected]));

export_q3_excel('问题三_动态事件调度方案.xlsx', base_sol, base_cost, base_det, base_carb, ...
    results, global_results, VT, demand_w, demand_v, tw, dist_matrix, parent_map, ...
    is_green, T_FORBID_END, T_START, SERVICE_TIME, CARBON_PRICE);

%% 9. 可视化

% 图1: 局部重调度 vs 全局重规划 成本对比
figure('Color','w','Position',[50 100 750 480]); hold on; grid on;
scene_names = {results.name};
local_costs = [results.cost];
global_costs = [global_results.cost];
x = 1:n_events;
bar_data = [local_costs; global_costs]';
b = bar(x, bar_data, 'grouped');
b(1).FaceColor = [0.2 0.6 0.85];
b(2).FaceColor = [0.85 0.33 0.1];
yline(base_cost, '--k', sprintf('静态基准 %.0f元', base_cost), 'LineWidth',1.5, 'FontSize',10);
set(gca, 'XTickLabel', scene_names, 'FontName','Microsoft YaHei');
ylabel('总配送成本 (元)'); xlabel('突发事件场景');
title('图1: 局部重调度 vs 全局重规划成本对比','FontWeight','bold','FontSize',13);
legend({'局部重调度','全局重规划','静态基准'}, 'Location','northwest');

% 图2: 各场景额外成本增量对比
figure('Color','w','Position',[100 80 700 450]); hold on; grid on;
delta_local = [results.delta_cost];
delta_global = [global_results.cost] - base_cost;
bar_delta = [delta_local; delta_global]';
b2 = bar(x, bar_delta, 'grouped');
b2(1).FaceColor = [0.2 0.7 0.3];
b2(2).FaceColor = [0.9 0.5 0.1];
set(gca, 'XTickLabel', scene_names, 'FontName','Microsoft YaHei');
ylabel('额外成本增量 (元)'); xlabel('突发事件场景');
title('图2: 动态调度额外成本评估','FontWeight','bold','FontSize',13);
legend({'局部重调度Δ','全局重规划Δ'}, 'Location','northwest');

% 图3: 场景1(订单取消)路径演变 — 调度前 vs 调度后
figure('Color','w','Position',[150 60 1200 500]);
subplot(1,2,1); hold on; grid on;
theta = linspace(0,2*pi,300);
fill(R_GREEN*cos(theta), R_GREEN*sin(theta), [0.92 0.98 0.92], ...
    'EdgeColor',[0 0.6 0], 'LineStyle','--', 'LineWidth',1.2);
scatter(cx, cy, 25, 'b', 'filled');
plot(0, 0, 'rp', 'MarkerSize',15, 'MarkerFaceColor','y');
clrs = lines(length(base_sol));
for r = 1:length(base_sol)
    seq = [0, base_sol(r).path, 0] + 1;
    plot(coords(seq,2), coords(seq,3), '-', 'Color', clrs(r,:), 'LineWidth',1);
end
if ~isempty(events{1}.custs)
    scatter(cx(events{1}.custs), cy(events{1}.custs), 120, 'rx', 'LineWidth',2.5);
end
title('调度前 (静态方案)','FontWeight','bold','FontSize',12);
xlabel('X (km)'); ylabel('Y (km)'); axis equal;
set(gca,'FontName','Microsoft YaHei');

subplot(1,2,2); hold on; grid on;
fill(R_GREEN*cos(theta), R_GREEN*sin(theta), [0.92 0.98 0.92], ...
    'EdgeColor',[0 0.6 0], 'LineStyle','--', 'LineWidth',1.2);
scatter(cx, cy, 25, 'b', 'filled');
plot(0, 0, 'rp', 'MarkerSize',15, 'MarkerFaceColor','y');
sol1 = results(1).sol;
clrs2 = lines(length(sol1));
for r = 1:length(sol1)
    seq = [0, sol1(r).path, 0] + 1;
    plot(coords(seq,2), coords(seq,3), '-', 'Color', clrs2(r,:), 'LineWidth',1);
end
if ~isempty(events{1}.custs)
    scatter(cx(events{1}.custs), cy(events{1}.custs), 120, 'rx', 'LineWidth',2.5);
end
title('调度后 (局部重调度)','FontWeight','bold','FontSize',12);
xlabel('X (km)'); ylabel('Y (km)'); axis equal;
set(gca,'FontName','Microsoft YaHei');
sgtitle('图3: 订单取消场景路径演变','FontWeight','bold','FontSize',14);

% 图4: 碳排放对比
figure('Color','w','Position',[200 100 700 450]); hold on; grid on;
local_carbs = [results.carbon];
global_carbs = [global_results.carbon];
bar_carb = [local_carbs; global_carbs]';
b3 = bar(x, bar_carb, 'grouped');
b3(1).FaceColor = [0.2 0.7 0.3];
b3(2).FaceColor = [0.6 0.2 0.2];
yline(base_carb, '--k', sprintf('基准 %.0fkg', base_carb), 'LineWidth',1.5, 'FontSize',10);
set(gca, 'XTickLabel', scene_names, 'FontName','Microsoft YaHei');
ylabel('碳排放量 (kg)'); xlabel('突发事件场景');
title('图4: 各场景碳排放对比','FontWeight','bold','FontSize',13);
legend({'局部重调度','全局重规划','静态基准'}, 'Location','northwest');

% 图5: 受影响路径数量与扰动度
figure('Color','w','Position',[250 120 650 420]); hold on; grid on;
n_aff = [results.n_affected];
n_total = length(base_sol) * ones(1, n_events);
bar([n_aff; n_total-n_aff]', 'stacked');
set(gca, 'XTickLabel', scene_names, 'FontName','Microsoft YaHei');
ylabel('路径数量'); xlabel('突发事件场景');
title('图5: 局部重调度扰动范围分析','FontWeight','bold','FontSize',13);
legend({'受影响路径','锁定路径'}, 'Location','northwest');
colormap([0.85 0.33 0.1; 0.7 0.85 0.7]);

% 图6: 成本构成堆叠对比 (基准 vs 4场景)
figure('Color','w','Position',[300 80 800 480]);
all_details = [base_det; cat(1, results.detail)];
bar_labels = ['静态基准', scene_names];
b4 = bar(all_details, 'stacked');
colors = [.85 .33 .1; .93 .69 .13; .5 .5 .5; .2 .7 .3; .8 .2 .2];
for k = 1:5, b4(k).FaceColor = colors(k,:); end
set(gca, 'XTickLabel', bar_labels, 'FontName','Microsoft YaHei');
ylabel('成本 (元)');
title('图6: 各场景成本构成对比','FontWeight','bold','FontSize',13);
legend({'固定启动','能耗','碳排','早到','晚到'}, 'Location','northeast');
grid on;

disp('>>> 全部图表已生成。');

%% ==================== 辅助函数 ====================

function export_q3_excel(file, base_sol, base_cost, base_det, base_carb, results, global_results, VT, dw0, dv0, tw0, D0, parent_map, is_green, t_forbid_end, t0, svc, cp)
    if exist(file,'file'), delete(file); end
    vt_names = {'燃油大(3t)','燃油中(1.5t)','燃油小(1.25t)','电动大(3t)','电动小(1.25t)'};
    summary = {'静态基准成本(元)',round(base_cost*100)/100; '静态基准碳排(kg)',round(base_carb*100)/100; '静态基准车辆数',length(base_sol); '固定启动成本(元)',round(base_det(1)*100)/100; '能耗成本(元)',round(base_det(2)*100)/100; '碳排成本(元)',round(base_det(3)*100)/100; '早到等待惩罚(元)',round(base_det(4)*100)/100; '晚到惩罚(元)',round(base_det(5)*100)/100};
    writetable(cell2table(summary,'VariableNames',{'指标','数值'}), file, 'Sheet','静态基准汇总');

    cmp = {};
    for i = 1:length(results)
        cmp(end+1,:) = {results(i).name, round(results(i).cost*100)/100, round(results(i).carbon*100)/100, results(i).n_vehicles, results(i).n_affected, round(results(i).delta_cost*100)/100, round(results(i).delta_carbon*100)/100, round(global_results(i).cost*100)/100, round(global_results(i).carbon*100)/100, round((global_results(i).cost-base_cost)*100)/100};
    end
    writetable(cell2table(cmp,'VariableNames',{'场景','局部重调度成本元','局部重调度碳排kg','局部车辆数','受影响路径数','局部成本变化元','局部碳排变化kg','全局重规划成本元','全局重规划碳排kg','全局成本变化元'}), file, 'Sheet','动态事件影响对比');

    veh = {}; detail = {};
    write_solution_rows('静态基准', base_sol, dw0, dv0, tw0, D0);
    for i = 1:length(results)
        write_solution_rows(results(i).name, results(i).sol, results(i).dw, results(i).dv, results(i).tw, results(i).dist);
    end
    writetable(cell2table(veh,'VariableNames',{'场景','车辆编号','车型','是否燃油','载货量kg','额定载重kg','载货体积m3','额定容积m3','服务客户数','绿色区客户数','行驶路径','出发时间','返回时间','总距离km','固定成本元','能耗成本元','碳排成本元','碳排放量kg','早到惩罚元','晚到惩罚元','该车总成本元'}), file, 'Sheet','各场景车辆使用方案');
    writetable(cell2table(detail,'VariableNames',{'场景','车辆编号','车型','访问序号','客户编号','是否绿色区','是否违反限行','本段距离km','到达时间','客户时间窗','时间窗状态','等待min','迟到min','开始服务','离开时间','配送量kg','车上剩余kg'}), file, 'Sheet','各场景逐站明细');
    fprintf('>>> 已导出: %s\n', file);

    function write_solution_rows(scene, sol, dw, dv, tw, D)
        for r = 1:length(sol)
            vt = sol(r).vtype; path = sol(r).path; cfg = VT(vt,:);
            path_orig = arrayfun(@(x) parent_map(x), path);
            path_str = '0'; for p = path_orig, path_str = [path_str, sprintf('->%d', p)]; end; path_str = [path_str,'->0'];
            [rc, det, carb, dist_sum, ret_t] = route_metrics_q3_export(path, vt, D, VT, dw, tw, svc, t0, cp);
            veh(end+1,:) = {scene, r, vt_names{vt}, cfg(7)==1, round(sum(dw(path))), round(cfg(1)), round(sum(dv(path))*10)/10, round(cfg(2)*10)/10, length(path), sum(is_green(path)), path_str, fmt_time_q3(t0), fmt_time_q3(ret_t), round(dist_sum*10)/10, round(det(1)*100)/100, round(det(2)*100)/100, round(det(3)*100)/100, round(carb*100)/100, round(det(4)*100)/100, round(det(5)*100)/100, round(rc*100)/100};
            curr = 0; t = t0; load_now = sum(dw(path));
            for k = 1:length(path)
                c = path(k); d = D(curr+1,c+1); dt = travel_time_tv(t,d); arr = t+dt; wait=0; late=0; st=arr; status='准时';
                if tw(c,1)>0 && arr<tw(c,1), wait=tw(c,1)-arr; st=tw(c,1); status=sprintf('早到等%.0fmin',wait); end
                if tw(c,2)>0 && arr>tw(c,2), late=arr-tw(c,2); status=sprintf('迟到%.0fmin',late); end
                forbid = is_green(c) && cfg(7)==1 && arr<t_forbid_end;
                detail(end+1,:) = {scene, r, vt_names{vt}, k, parent_map(c), is_green(c), forbid, round(d*10)/10, fmt_time_q3(arr), sprintf('%s-%s',fmt_time_q3(tw(c,1)),fmt_time_q3(tw(c,2))), status, round(wait), round(late), fmt_time_q3(st), fmt_time_q3(st+svc), round(dw(c)), round(load_now-dw(c))};
                t=st+svc; curr=c; load_now=load_now-dw(c);
            end
            db=D(curr+1,1); ret=t+travel_time_tv(t,db);
            detail(end+1,:) = {scene, r, vt_names{vt}, length(path)+1, 0, false, false, round(db*10)/10, fmt_time_q3(ret), '-', '返回仓库', 0, 0, '-', fmt_time_q3(ret), 0, 0};
        end
    end
end

function [total, det, carbon, dist_sum, ret_t] = route_metrics_q3_export(path, vt, D, VT, dw, tw, svc, t0, cp)
    cfg = VT(vt,:); det = [cfg(3),0,0,0,0]; carbon=0; dist_sum=0; curr=0; t=t0; load_now=sum(dw(path));
    for k=1:length(path)
        c=path(k); d=D(curr+1,c+1); dist_sum=dist_sum+d; dt=travel_time_tv(t,d); va=max(5,min(60,d/max(dt/60,0.001)));
        if cfg(7)==1, bc=(0.0025*va^2-0.2554*va+31.75)/100; else, bc=(0.0014*va^2-0.12*va+36.19)/100; end
        ac=bc*d*(1+cfg(4)*min(load_now/cfg(1),1.5)); det(2)=det(2)+ac*cfg(6); seg=ac*cfg(5); det(3)=det(3)+seg*cp; carbon=carbon+seg;
        t=t+dt; if tw(c,1)>0 && t<tw(c,1), det(4)=det(4)+(tw(c,1)-t)*20/60; t=tw(c,1); end; if tw(c,2)>0 && t>tw(c,2), det(5)=det(5)+(t-tw(c,2))*50/60; end
        t=t+svc; curr=c; load_now=load_now-dw(c);
    end
    db=D(curr+1,1); dist_sum=dist_sum+db; dtb=travel_time_tv(t,db); va=max(5,min(60,db/max(dtb/60,0.001)));
    if cfg(7)==1, bc=(0.0025*va^2-0.2554*va+31.75)/100; else, bc=(0.0014*va^2-0.12*va+36.19)/100; end
    det(2)=det(2)+bc*db*cfg(6); seg=bc*db*cfg(5); det(3)=det(3)+seg*cp; carbon=carbon+seg; ret_t=t+dtb; total=sum(det);
end

function s = fmt_time_q3(t)
    s = sprintf('%d:%02d', floor(t/60), round(mod(t,60)));
end

function [status, delivered_all, in_transit_all] = simulate_execution( ...
        sol, D, tw, dw, t0, svc, t_event)
    n_routes = length(sol);
    status = struct('delivered',{},'pending',{},'current_node',{},'current_time',{});
    delivered_all = [];
    in_transit_all = [];

    for r = 1:n_routes
        path = sol(r).path;
        curr = 0; t = t0;
        del = []; pend = path;

        for k = 1:length(path)
            c = path(k);
            dt = travel_time_tv(t, D(curr+1, c+1));
            arr = t + dt;
            if arr > t_event
                pend = path(k:end);
                break;
            end
            t = arr;
            if tw(c,1)>0 && t < tw(c,1), t = tw(c,1); end
            t = t + svc;
            del = [del, c];
            curr = c;
            if k == length(path), pend = []; end
        end

        status(r).delivered = del;
        status(r).pending = setdiff(path, del, 'stable');
        status(r).current_node = curr;
        status(r).current_time = min(t, t_event);
        delivered_all = [delivered_all, del];
    end
    in_transit_all = [];
    for r = 1:n_routes
        in_transit_all = [in_transit_all, status(r).pending];
    end
end

function [affected, mod_tw, mod_dw, mod_dv, mod_dist, active] = ...
        apply_event(evt, sol, status, tw, dw, dv, D, n_cust, is_green)
    mod_tw = tw; mod_dw = dw; mod_dv = dv; mod_dist = D;
    affected = [];
    active = 1:n_cust;

    switch evt.type
        case 'cancel'
            for r = 1:length(sol)
                if any(ismember(evt.custs, status(r).pending))
                    affected = [affected, r];
                end
            end
        case 'new_order'
            for r = 1:length(sol)
                path = sol(r).path;
                for c = evt.custs
                    if any(ismember(c, path))
                        affected = [affected, r];
                    end
                end
            end
            if isempty(affected) && ~isempty(sol)
                [~, mi] = min(arrayfun(@(s) length(s.pending), status));
                affected = mi;
            end
        case 'address_change'
            for ci = 1:length(evt.custs)
                c = evt.custs(ci);
                mod_dist(c+1, :) = mod_dist(c+1, :) * 1.3;
                mod_dist(:, c+1) = mod_dist(:, c+1) * 1.3;
            end
            for r = 1:length(sol)
                if any(ismember(evt.custs, status(r).pending))
                    affected = [affected, r];
                end
            end
        case 'tw_change'
            for r = 1:length(sol)
                if any(ismember(evt.custs, status(r).pending))
                    affected = [affected, r];
                end
            end
    end
    affected = unique(affected);
end

function sol = greedy_init_partial(custs, dw, dv, tw, D, VT, t0, svc, is_green, t_forbid_end)
    sol = struct('path',{},'vtype',{});
    stock = VT(:,8)';
    unvisited = custs;

    while ~isempty(unvisited)
        best_vt = 0; best_path = []; best_uc = inf;
        for vt = [5, 4, 3, 2, 1]
            if stock(vt) <= 0, continue; end
            path = []; lw = 0; lv = 0; curr = 0; t = t0;
            remain = unvisited;
            while true
                bn = 0; bi = inf;
                for j = 1:length(remain)
                    c = remain(j);
                    if lw+dw(c)>VT(vt,1) || lv+dv(c)>VT(vt,2), continue; end
                    if is_green(c) && VT(vt,7)==1
                        arr_est = t + travel_time_tv(t, D(curr+1,c+1));
                        if arr_est < t_forbid_end, continue; end
                    end
                    dt = travel_time_tv(t, D(curr+1,c+1));
                    arr = t + dt;
                    inc = D(curr+1,c+1);
                    if tw(c,1)>0 && arr<tw(c,1), inc=inc+(tw(c,1)-arr)*0.5; end
                    if tw(c,2)>0 && arr>tw(c,2), inc=inc+(arr-tw(c,2))*2; end
                    if inc < bi, bi=inc; bn=j; end
                end
                if bn==0, break; end
                c = remain(bn);
                dt = travel_time_tv(t, D(curr+1,c+1));
                t = t + dt;
                if tw(c,1)>0 && t<tw(c,1), t=tw(c,1); end
                t = t + svc;
                path = [path, c]; lw=lw+dw(c); lv=lv+dv(c);
                curr = c; remain(bn) = [];
            end
            if ~isempty(path)
                [rc,~] = eval_single_route_q3(path, vt, D, VT, dw, tw, ...
                    0.65, 20/60, 50/60, svc, t0, is_green, t_forbid_end);
                uc = rc / length(path);
                if uc < best_uc, best_uc=uc; best_vt=vt; best_path=path; end
            end
        end
        if isempty(best_path)
            c0 = unvisited(1);
            feasible_vt = find(stock > 0 & dw(c0) <= VT(:,1)' & dv(c0) <= VT(:,2)', 1);
            if isempty(feasible_vt)
                error('车辆库存不足或剩余车辆容量不足，客户%d无法分配: %.1fkg, %.2fm3', c0, dw(c0), dv(c0));
            end
            best_vt = feasible_vt; best_path = c0;
        end
        sol(end+1).path = best_path;
        sol(end).vtype = best_vt;
        stock(best_vt) = stock(best_vt) - 1;
        unvisited = setdiff(unvisited, best_path);
    end
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
                    ae = t + travel_time_tv(t, D(curr+1,c+1));
                    if ae < t_forbid_end, continue; end
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

function [total, detail, carbon_kg] = eval_solution_q3(sol, D, VT, dw, tw, cp, ep, lp, svc, t0, is_green, t_forbid_end)
    c_fix=0; c_eng=0; c_carb=0; c_early=0; c_late=0; carbon_kg=0;
    pen_g = 0;
    for r = 1:length(sol)
        vt = sol(r).vtype; path = sol(r).path;
        cfg = VT(vt,:);
        c_fix = c_fix + cfg(3);
        curr = 0; t = t0;
        total_load = sum(dw(path)); cur_load = total_load;
        for k = 1:length(path)
            c = path(k);
            dist = D(curr+1, c+1);
            dt = travel_time_tv(t, dist);
            va = dist / max(dt/60, 0.001); va = max(5, min(60, va));
            if cfg(7)==1, bc=(0.0025*va^2-0.2554*va+31.75)/100;
            else, bc=(0.0014*va^2-0.12*va+36.19)/100; end
            lr = min(cur_load/cfg(1), 1.5);
            ac = bc * dist * (1 + cfg(4)*lr);
            c_eng = c_eng + ac * cfg(6);
            seg_c = ac * cfg(5);
            c_carb = c_carb + seg_c * cp;
            carbon_kg = carbon_kg + seg_c;
            t = t + dt;
            if is_green(c) && cfg(7)==1 && t < t_forbid_end
                pen_g = pen_g + 5000;
            end
            if tw(c,1)>0 && t<tw(c,1), c_early=c_early+(tw(c,1)-t)*ep; t=tw(c,1); end
            if tw(c,2)>0 && t>tw(c,2), c_late=c_late+(t-tw(c,2))*lp; end
            t = t + svc; curr = c; cur_load = cur_load - dw(c);
        end
        db = D(curr+1, 1); dtb = travel_time_tv(t, db);
        vb = db/max(dtb/60,0.001); vb = max(5,min(60,vb));
        if cfg(7)==1, bcb=(0.0025*vb^2-0.2554*vb+31.75)/100;
        else, bcb=(0.0014*vb^2-0.12*vb+36.19)/100; end
        c_eng = c_eng + bcb*db*cfg(6);
        seg_cb = bcb*db*cfg(5);
        c_carb = c_carb + seg_cb*cp;
        carbon_kg = carbon_kg + seg_cb;
    end
    detail = [c_fix, c_eng, c_carb, c_early, c_late];
    total = sum(detail) + pen_g;
end

function ok = check_solution_feasible_q3(sol, VT, dw, dv)
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

function [cost, detail] = eval_single_route_q3(path, vt, D, VT, dw, tw, cp, ep, lp, svc, t0, is_green, t_forbid_end)
    cfg = VT(vt,:);
    c_eng=0; c_carb=0; c_early=0; c_late=0;
    curr=0; t=t0;
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
        if is_green(c) && cfg(7)==1 && t < t_forbid_end, pen_g=pen_g+5000; end
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

function [best_sol, best_cost, best_det, best_carb] = run_alns_local( ...
        sol, n_cust, dw, dv, tw, D, VT, cp, ep, lp, svc, t0, ...
        is_green, t_forbid_end, max_iter, SIGMA, react, d_rate, sa_t0, sa_cool, active_custs)

    [cur_cost, cur_det, cur_carb] = eval_solution_q3(sol, D, VT, dw, tw, cp, ep, lp, svc, t0, is_green, t_forbid_end);
    best_sol = sol; best_cost = cur_cost; best_det = cur_det; best_carb = cur_carb;
    cur_sol = sol;

    N_D = 3; N_R = 2;
    w_d = ones(1,N_D)/N_D; w_r = ones(1,N_R)/N_R;
    score_d = zeros(1,N_D); score_r = zeros(1,N_R);
    cnt_d = zeros(1,N_D); cnt_r = zeros(1,N_R);
    sa_temp = sa_t0;

    all_c = [];
    for r=1:length(sol), all_c = [all_c, sol(r).path]; end
    n_active = length(all_c);

    for it = 1:max_iter
        d_op = roulette_sel(w_d);
        r_op = roulette_sel(w_r);
        n_rm = randi([max(1,round(d_rate(1)*n_active)), max(2,round(d_rate(2)*n_active))]);

        [removed, partial] = destroy_q3(cur_sol, d_op, n_rm, D, dw, tw);
        try
            new_sol = repair_q3(partial, removed, r_op, D, VT, dw, dv, tw, cp, ep, lp, svc, t0, is_green, t_forbid_end);
        catch
            reward = 0;
            score_d(d_op) = score_d(d_op) + reward;
            score_r(r_op) = score_r(r_op) + reward;
            cnt_d(d_op) = cnt_d(d_op) + 1;
            cnt_r(r_op) = cnt_r(r_op) + 1;
            sa_temp = sa_temp * sa_cool;
            continue;
        end

        [new_cost, new_det, new_carb] = eval_solution_q3(new_sol, D, VT, dw, tw, cp, ep, lp, svc, t0, is_green, t_forbid_end);

        % 库存与容量约束检查
        feasible_ok = check_solution_feasible_q3(new_sol, VT, dw, dv);

        reward = 0;
        if feasible_ok && new_cost < best_cost
            best_sol = new_sol; best_cost = new_cost; best_det = new_det; best_carb = new_carb;
            cur_sol = new_sol; cur_cost = new_cost;
            reward = SIGMA(1);
        elseif feasible_ok && new_cost < cur_cost
            cur_sol = new_sol; cur_cost = new_cost;
            reward = SIGMA(2);
        elseif feasible_ok && rand() < exp(-(new_cost-cur_cost)/max(sa_temp,1))
            cur_sol = new_sol; cur_cost = new_cost;
            reward = SIGMA(3);
        end

        score_d(d_op) = score_d(d_op) + reward;
        score_r(r_op) = score_r(r_op) + reward;
        cnt_d(d_op) = cnt_d(d_op) + 1;
        cnt_r(r_op) = cnt_r(r_op) + 1;
        sa_temp = sa_temp * sa_cool;

        if mod(it, 100) == 0
            for k=1:N_D
                if cnt_d(k)>0
                    w_d(k) = w_d(k)*(1-react) + react*score_d(k)/cnt_d(k);
                end
            end
            for k=1:N_R
                if cnt_r(k)>0
                    w_r(k) = w_r(k)*(1-react) + react*score_r(k)/cnt_r(k);
                end
            end
            w_d = max(w_d, 0.05); w_d = w_d/sum(w_d);
            w_r = max(w_r, 0.05); w_r = w_r/sum(w_r);
            score_d(:)=0; score_r(:)=0; cnt_d(:)=0; cnt_r(:)=0;
        end
    end
end

function [removed, partial] = destroy_q3(sol, op, n_rm, D, dw, tw)
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

function sol = repair_q3(partial, removed, op, D, VT, dw, dv, tw, cp, ep, lp, svc, t0, is_green, t_forbid_end)
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
                ae = estimate_arr(path, 0, D, t0, svc, tw);
                if ae < t_forbid_end, continue; end
            end
            for pos = 1:length(path)+1
                np = [path(1:pos-1), c, path(pos:end)];
                [nc,~] = eval_single_route_q3(np, vt, D, VT, dw, tw, cp, ep, lp, svc, t0, is_green, t_forbid_end);
                [oc,~] = eval_single_route_q3(path, vt, D, VT, dw, tw, cp, ep, lp, svc, t0, is_green, t_forbid_end);
                inc = nc - oc;
                if op==2, inc = inc - 0.1*rand()*abs(inc); end
                if inc < best_inc, best_inc=inc; best_r=r; best_pos=pos; end
            end
        end
        if best_r > 0 && best_inc < 1e6
            sol(best_r).path = [sol(best_r).path(1:best_pos-1), c, sol(best_r).path(best_pos:end)];
        else
            vt = pick_veh(c, dw(c), dv(c), VT, stock_used, is_green);
            sol(end+1).path = c;
            sol(end).vtype = vt;
            stock_used(vt) = stock_used(vt) + 1;
        end
    end
end

function vt = pick_veh(c, w, v, VT, used, is_green)
    if is_green(c), order = [5, 4]; else, order = [5, 4, 3, 2, 1]; end
    for k = order
        if w<=VT(k,1) && v<=VT(k,2) && used(k)<VT(k,8), vt = k; return; end
    end
    for k = [3, 2, 1]
        if w<=VT(k,1) && v<=VT(k,2) && used(k)<VT(k,8), vt = k; return; end
    end
    error('车辆库存不足或剩余车辆容量不足，无法分配 %.1fkg, %.2fm3 的需求', w, v);
end

function arr = estimate_arr(path, pos_idx, D, t0, svc, tw)
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
