% =========================================================================
% 华中杯 A题 问题一：静态调度 (TD-HVRPTW + ALNS)
% 时变软时间窗异构车辆路径模型 + 自适应大邻域搜索算法
% =========================================================================
clear; clc; close all; warning off; rng(1);

%% 1. 数据预处理
disp('===== 问题一：TD-HVRPTW + ALNS 求解引擎启动 =====');

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
MAX_CAP_W = 3000;  % 最大车辆载重
MAX_CAP_V = 15.0;  % 最大车辆容积
orig_n_cust = n_cust;
orig_demand_w = demand_w;
orig_demand_v = demand_v;
orig_tw = tw;
orig_coords = coords;
orig_dist = dist_matrix;

new_dw = []; new_dv = []; new_tw = []; new_coords_ext = [coords(1,:)];
parent_map = [];

for i = 1:orig_n_cust
    w = orig_demand_w(i); v = orig_demand_v(i);
    n_split = max(ceil(w / MAX_CAP_W), ceil(v / MAX_CAP_V));
    sw_each = w / n_split;
    sv_each = v / n_split;
    for s = 1:n_split
        new_dw(end+1,1) = sw_each;
        new_dv(end+1,1) = sv_each;
        new_tw(end+1,:) = orig_tw(i,:);
        new_coords_ext(end+1,:) = orig_coords(i+1,:);
        parent_map(end+1) = i;
    end
end

% 二次拆分：控制只能由3吨大车服务的虚拟客户数量，避免大车库存不足
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
    new_dw(end+1,1) = new_dw(idx);
    new_dv(end+1,1) = new_dv(idx);
    new_tw(end+1,:) = new_tw(idx,:);
    new_coords_ext(end+1,:) = new_coords_ext(idx+1,:);
    parent_map(end+1) = parent_map(idx);
    is_big_only = (new_dw > BIG_ONLY_W) | (new_dv > BIG_ONLY_V);
end

n_cust = length(new_dw);
demand_w = new_dw;
demand_v = new_dv;
tw = new_tw;
coords = new_coords_ext;

new_dist = zeros(n_cust+1);
for i = 0:n_cust
    pi = 0; if i>0, pi = parent_map(i); end
    for j = 0:n_cust
        pj = 0; if j>0, pj = parent_map(j); end
        new_dist(i+1, j+1) = orig_dist(pi+1, pj+1);
    end
end
dist_matrix = new_dist;

fprintf('需求拆分: %d 原始客户 -> %d 虚拟客户 (拆分了 %d 个超容量客户)\n', ...
    orig_n_cust, n_cust, n_cust - orig_n_cust);
fprintf('虚拟客户需求分布: 重量 min=%.0f max=%.0f mean=%.0f | 体积 min=%.1f max=%.1f mean=%.1f\n', ...
    min(new_dw), max(new_dw), mean(new_dw), min(new_dv), max(new_dv), mean(new_dv));
fprintf('总需求: 重量=%.0fkg 体积=%.1fm3 | 大车(3t)理论最少=%.0f辆\n', ...
    sum(new_dw), sum(new_dv), ceil(sum(new_dw)/3000));

%% 2. 车型定义 (题目原始参数)
% [载重kg, 容积m3, 启动成本, 满载修正, 碳排转换系数, 能源单价, 是否燃油, 库存]
VT = [
    3000, 13.5, 400, 0.40, 2.547, 7.61, 1, 60;   % 燃油大
    1500, 10.8, 400, 0.40, 2.547, 7.61, 1, 50;   % 燃油中
    1250,  6.5, 400, 0.40, 2.547, 7.61, 1, 50;   % 燃油小
    3000, 15.0, 400, 0.35, 0.501, 1.64, 0, 10;   % 电动大
    1250,  8.5, 400, 0.35, 0.501, 1.64, 0, 15;   % 电动小
];
n_vtype = size(VT,1);

% 成本参数
CARBON_PRICE = 0.65;       % 碳排放成本 (元/kg)
EARLY_PEN    = 20/60;      % 早到惩罚 (元/分钟) = 20元/小时
LATE_PEN     = 50/60;      % 晚到惩罚 (元/分钟) = 50元/小时
SERVICE_TIME = 20;         % 服务时间 (分钟)
T_START      = 480;        % 出发时刻 8:00 = 480min

%% 3. ALNS 参数
ALNS_ITER    = 3000;       % 总迭代次数
N_DESTROY    = 3;          % 破坏算子数
N_REPAIR     = 2;          % 修复算子数
SIGMA        = [33, 9, 3]; % 奖励分数: 全局最优/当前改进/接受
REACT_FACTOR = 0.1;        % 权重反应因子
DESTROY_RATE = [0.1, 0.4]; % 移除比例范围
SA_T0        = 1000;       % SA初始温度
SA_COOL      = 0.9995;     % SA冷却率

%% 4. 构造初始解 (贪婪插入)
disp('>>> 构造初始可行解...');
sol = greedy_init(n_cust, demand_w, demand_v, tw, dist_matrix, VT, T_START, SERVICE_TIME);
[sol_cost, sol_detail, ~] = eval_solution(sol, dist_matrix, VT, demand_w, tw, ...
    CARBON_PRICE, EARLY_PEN, LATE_PEN, SERVICE_TIME, T_START);
cust_per_veh = arrayfun(@(s) length(s.path), sol);
fprintf('初始解: %.2f 元, %d 辆车, 平均每车%.1f客户, 最多%d客户\n', ...
    sol_cost, length(sol), mean(cust_per_veh), max(cust_per_veh));
vt_init = zeros(1,n_vtype);
for rr=1:length(sol), vt_init(sol(rr).vtype)=vt_init(sol(rr).vtype)+1; end
fprintf('初始车型分布: 燃油大=%d/%d 燃油中=%d/%d 燃油小=%d/%d 电动大=%d/%d 电动小=%d/%d\n', ...
    vt_init(1),VT(1,8), vt_init(2),VT(2,8), vt_init(3),VT(3,8), vt_init(4),VT(4,8), vt_init(5),VT(5,8));

%% 5. ALNS 主循环
disp('>>> ALNS 求解中...');
[~, ~, init_carbon] = eval_solution(sol, dist_matrix, VT, demand_w, tw, ...
    CARBON_PRICE, EARLY_PEN, LATE_PEN, SERVICE_TIME, T_START);
best_sol = sol; best_cost = sol_cost; best_detail = sol_detail; best_carbon = init_carbon;
curr_sol = sol; curr_cost = sol_cost;

w_destroy = ones(1, N_DESTROY); w_repair = ones(1, N_REPAIR);
score_d = zeros(1, N_DESTROY); score_r = zeros(1, N_REPAIR);
count_d = zeros(1, N_DESTROY); count_r = zeros(1, N_REPAIR);
sa_temp = SA_T0;

cost_history = zeros(ALNS_ITER, 1);
best_history = zeros(ALNS_ITER, 1);

for iter = 1:ALNS_ITER
    % 自适应选择算子
    d_op = roulette_select(w_destroy);
    r_op = roulette_select(w_repair);

    % 确定移除数量
    n_remove = max(1, round(n_cust * (DESTROY_RATE(1) + rand()*(DESTROY_RATE(2)-DESTROY_RATE(1)))));

    % 破坏
    [removed, partial] = destroy(curr_sol, d_op, n_remove, dist_matrix, demand_w, tw);

    % 修复
    try
        new_sol = repair(partial, removed, r_op, dist_matrix, VT, demand_w, demand_v, tw, ...
            CARBON_PRICE, EARLY_PEN, LATE_PEN, SERVICE_TIME, T_START);
    catch
        reward = 0;
        score_d(d_op) = score_d(d_op) + reward;
        count_d(d_op) = count_d(d_op) + 1;
        score_r(r_op) = score_r(r_op) + reward;
        count_r(r_op) = count_r(r_op) + 1;
        sa_temp = sa_temp * SA_COOL;
        cost_history(iter) = curr_cost;
        best_history(iter) = best_cost;
        continue;
    end

    [new_cost, new_detail, new_carbon] = eval_solution(new_sol, dist_matrix, VT, demand_w, tw, ...
        CARBON_PRICE, EARLY_PEN, LATE_PEN, SERVICE_TIME, T_START);

    % 库存与容量约束检查
    feasible_ok = check_solution_feasible(new_sol, VT, demand_w, demand_v);

    % 评价与接受
    reward = 0;
    if feasible_ok && new_cost < best_cost
        best_sol = new_sol; best_cost = new_cost; best_detail = new_detail; best_carbon = new_carbon;
        curr_sol = new_sol; curr_cost = new_cost;
        reward = SIGMA(1);
    elseif feasible_ok && new_cost < curr_cost
        curr_sol = new_sol; curr_cost = new_cost;
        reward = SIGMA(2);
    elseif feasible_ok && rand() < exp(-(new_cost - curr_cost)/sa_temp)
        curr_sol = new_sol; curr_cost = new_cost;
        reward = SIGMA(3);
    end

    score_d(d_op) = score_d(d_op) + reward;
    count_d(d_op) = count_d(d_op) + 1;
    score_r(r_op) = score_r(r_op) + reward;
    count_r(r_op) = count_r(r_op) + 1;

    sa_temp = sa_temp * SA_COOL;
    cost_history(iter) = curr_cost;
    best_history(iter) = best_cost;

    % 每100次更新权重
    if mod(iter, 100) == 0
        for d = 1:N_DESTROY
            if count_d(d) > 0
                w_destroy(d) = w_destroy(d)*(1-REACT_FACTOR) + REACT_FACTOR*score_d(d)/count_d(d);
            end
        end
        for r = 1:N_REPAIR
            if count_r(r) > 0
                w_repair(r) = w_repair(r)*(1-REACT_FACTOR) + REACT_FACTOR*score_r(r)/count_r(r);
            end
        end
        w_destroy = max(w_destroy, 0.05); w_repair = max(w_repair, 0.05);
        score_d(:)=0; count_d(:)=0; score_r(:)=0; count_r(:)=0;
    end

    if mod(iter, 500) == 0
        fprintf('  迭代 %d/%d | 当前: %.2f | 最优: %.2f | 车辆: %d\n', ...
            iter, ALNS_ITER, curr_cost, best_cost, length(best_sol));
    end
end

disp('===== ALNS 求解完毕 =====');

%% 6. 结果统计与输出
fprintf('\n============ 问题一最终结果 ============\n');
fprintf('最优总成本: %.2f 元\n', best_cost);
fprintf('碳排放总量: %.2f kg\n', best_carbon);
fprintf('派车总数: %d 辆\n', length(best_sol));
fprintf('成本构成: 固定=%.1f  能耗=%.1f  碳排=%.1f  早到=%.1f  晚到=%.1f\n', ...
    best_detail(1), best_detail(2), best_detail(3), best_detail(4), best_detail(5));

vtype_count = zeros(1, n_vtype);
max_over_w = 0; max_over_v = 0;
for r = 1:length(best_sol)
    vt = best_sol(r).vtype;
    vtype_count(vt) = vtype_count(vt)+1;
    max_over_w = max(max_over_w, sum(demand_w(best_sol(r).path)) - VT(vt,1));
    max_over_v = max(max_over_v, sum(demand_v(best_sol(r).path)) - VT(vt,2));
end
if max_over_w > 1e-6 || max_over_v > 1e-6
    error('最终方案存在单车超载: 最大超重 %.2fkg, 最大超体积 %.2fm3', max_over_w, max_over_v);
end
vt_names = {'燃油大(3t)','燃油中(1.5t)','燃油小(1.25t)','电动大(3t)','电动小(1.25t)'};
for v = 1:n_vtype
    if vtype_count(v)>0
        flag = '';
        if vtype_count(v) > VT(v,8), flag = ' *** 超出库存! ***'; end
        fprintf('  %s: %d/%d 辆%s\n', vt_names{v}, vtype_count(v), VT(v,8), flag);
    end
end
fprintf('========================================\n');

%% 7. 可视化
% 图1: ALNS收敛曲线
figure('Color','w','Position',[100 100 700 450]); hold on; grid on;
plot(1:ALNS_ITER, cost_history, 'Color',[0.85 0.33 0.10], 'LineWidth',1);
plot(1:ALNS_ITER, best_history, 'Color',[0 0.45 0.74], 'LineWidth',2);
xlabel('迭代次数'); ylabel('总成本 (元)');
title('ALNS算法收敛曲线','FontWeight','bold','FontSize',13);
legend('当前解','历史最优','Location','northeast');
set(gca,'FontName','Microsoft YaHei');

% 图2: 成本构成饼图
figure('Color','w','Position',[850 100 550 450]);
labels = {sprintf('固定启动\n%.2f元',best_detail(1)), sprintf('能耗成本\n%.2f元',best_detail(2)), ...
    sprintf('碳排成本\n%.2f元',best_detail(3)), sprintf('早到等待\n%.2f元',best_detail(4)), ...
    sprintf('晚到惩罚\n%.2f元',best_detail(5))};
pie(max(best_detail,0.01), [0 0 0 0.1 0.15], labels);
title(sprintf('最优方案成本结构 (总计%.2f元)',best_cost),'FontWeight','bold','FontSize',13);
set(gca,'FontName','Microsoft YaHei');

% 图3: 路径拓扑图
figure('Color','w','Position',[200 150 900 700]); hold on; grid on;
scatter(coords(2:end,2), coords(2:end,3), 50, 'b', 'filled', 'MarkerEdgeColor','k');
scatter(coords(1,2), coords(1,3), 300, 'p', 'MarkerFaceColor','r', 'MarkerEdgeColor','k');
clrs = lines(length(best_sol));
for r = 1:length(best_sol)
    seq = [0, best_sol(r).path, 0] + 1;
    X = coords(seq,2); Y = coords(seq,3);
    if VT(best_sol(r).vtype,7)==1, ls='-'; else, ls='--'; end
    plot(X, Y, ls, 'Color', clrs(r,:), 'LineWidth', 1.3);
end
for i = 2:2:size(coords,1)
    text(coords(i,2)+0.3, coords(i,3)+0.3, num2str(parent_map(i-1)), 'FontSize',7, 'Color',[.4 .4 .4]);
end
text(coords(1,2)+0.8, coords(1,3)+0.8, '配送中心', 'FontSize',11, 'FontWeight','bold', 'Color','r');
xlabel('X (km)'); ylabel('Y (km)');
title('TD-HVRPTW最优路径拓扑 (实线燃油/虚线电动)','FontWeight','bold','FontSize',13);
set(gca,'FontName','Microsoft YaHei');

% 图4: 车型使用柱状图
figure('Color','w','Position',[500 200 600 400]);
bar(vtype_count, 'FaceColor','flat','CData',[.85 .33 .1;.93 .69 .13;1 .85 .2;.2 .7 .3;.1 .5 .6]);
set(gca,'XTickLabel',vt_names,'FontName','Microsoft YaHei');
ylabel('使用数量 (辆)'); title('各车型使用情况','FontWeight','bold','FontSize',13);
grid on;

%% 8. 完整调度方案导出Excel (车辆使用+行驶路径+到达时间+成本构成)
fprintf('\n============ 生成完整调度方案Excel ============\n');
vt_names_short = {'燃油大','燃油中','燃油小','电动大','电动小'};
vt_names_full = {'燃油大(3t)','燃油中(1.5t)','燃油小(1.25t)','电动大(3t)','电动小(1.25t)'};
excel_file = '问题一_完整调度方案.xlsx';
if exist(excel_file,'file'), delete(excel_file); end

% ---- Sheet1: 车辆使用方案与成本汇总 (每辆车一行) ----
veh_summary = {};
for r = 1:length(best_sol)
    vt = best_sol(r).vtype; path = best_sol(r).path; cfg = VT(vt,:);
    path_orig = arrayfun(@(x) parent_map(x), path);
    path_str = '0';
    for pp = path_orig, path_str = [path_str, sprintf('->%d', pp)]; end
    path_str = [path_str, '->0'];
    curr = 0; t = T_START;
    tl = sum(demand_w(path)); tv = sum(demand_v(path)); cl = tl;
    rf = cfg(3); re = 0; rc = 0; rck = 0; rea = 0; rla = 0; td = 0;
    for k = 1:length(path)
        c = path(k); d = dist_matrix(curr+1,c+1); td = td+d;
        dt = travel_time_tv(t,d); va = d/max(dt/60,0.001); va = max(5,min(60,va));
        if cfg(7)==1, bc=(0.0025*va^2-0.2554*va+31.75)/100;
        else, bc=(0.0014*va^2-0.12*va+36.19)/100; end
        lr = min(cl/cfg(1),1.5); ac = bc*d*(1+cfg(4)*lr);
        re = re+ac*cfg(6); sc = ac*cfg(5); rc = rc+sc*CARBON_PRICE; rck = rck+sc;
        t = t+dt;
        if tw(c,1)>0&&t<tw(c,1), rea=rea+(tw(c,1)-t)*EARLY_PEN; t=tw(c,1); end
        if tw(c,2)>0&&t>tw(c,2), rla=rla+(t-tw(c,2))*LATE_PEN; end
        t = t+SERVICE_TIME; curr = c; cl = cl-demand_w(c);
    end
    db = dist_matrix(curr+1,1); td = td+db;
    dtb = travel_time_tv(t,db); vb = db/max(dtb/60,0.001); vb = max(5,min(60,vb));
    if cfg(7)==1, bcb=(0.0025*vb^2-0.2554*vb+31.75)/100;
    else, bcb=(0.0014*vb^2-0.12*vb+36.19)/100; end
    re = re+bcb*db*cfg(6); scb = bcb*db*cfg(5); rc = rc+scb*CARBON_PRICE; rck = rck+scb;
    ret_t = t+dtb;
    ds = sprintf('%d:%02d',floor(T_START/60),round(mod(T_START,60)));
    rs = sprintf('%d:%02d',floor(ret_t/60),round(mod(ret_t,60)));
    rt = rf+re+rc+rea+rla;
    veh_summary(r,:) = {r, vt_names_full{vt}, cfg(7)==1, round(tl), round(cfg(1)), ...
        round(tv*10)/10, round(cfg(2)*10)/10, length(path), path_str, ds, rs, ...
        round(td*10)/10, round(rf*100)/100, round(re*100)/100, round(rc*100)/100, ...
        round(rck*100)/100, round(rea*100)/100, round(rla*100)/100, round(rt*100)/100};
end
h1 = {'车辆编号','车型','是否燃油','载货量kg','额定载重kg','载货体积m3','额定容积m3',...
    '服务客户数','行驶路径','出发时间','返回时间','总距离km',...
    '固定成本元','能耗成本元','碳排成本元','碳排放量kg','早到惩罚元','晚到惩罚元','该车总成本元'};
writetable(cell2table(veh_summary,'VariableNames',h1), excel_file, 'Sheet','车辆使用方案与成本');

% ---- Sheet2: 逐站点行程明细 ----
detail_data = {}; row = 0;
for r = 1:length(best_sol)
    vt = best_sol(r).vtype; path = best_sol(r).path; cfg = VT(vt,:);
    curr = 0; t = T_START; cl = sum(demand_w(path));
    for k = 1:length(path)
        c = path(k); co = parent_map(c);
        d = dist_matrix(curr+1,c+1); dt = travel_time_tv(t,d); arr = t+dt;
        va = d/max(dt/60,0.001); va = max(5,min(60,va));
        wm = 0; lm = 0; ss = '准时';
        if tw(c,1)>0&&arr<tw(c,1), wm=tw(c,1)-arr; ss=sprintf('早到等%.0fmin',wm); st=tw(c,1);
        elseif tw(c,2)>0&&arr>tw(c,2), lm=arr-tw(c,2); ss=sprintf('迟到%.0fmin',lm); st=arr;
        else, st=arr; end
        lv = st+SERVICE_TIME; dw_c = demand_w(c); cl = cl-dw_c;
        if cfg(7)==1, bc=(0.0025*va^2-0.2554*va+31.75)/100;
        else, bc=(0.0014*va^2-0.12*va+36.19)/100; end
        lr = min((cl+dw_c)/cfg(1),1.5); ac = bc*d*(1+cfg(4)*lr);
        sf = ac*cfg(6); sk = ac*cfg(5);
        as = sprintf('%d:%02d',floor(arr/60),round(mod(arr,60)));
        sts = sprintf('%d:%02d',floor(st/60),round(mod(st,60)));
        ls = sprintf('%d:%02d',floor(lv/60),round(mod(lv,60)));
        tws = sprintf('%d:%02d-%d:%02d',floor(tw(c,1)/60),round(mod(tw(c,1),60)),...
            floor(tw(c,2)/60),round(mod(tw(c,2),60)));
        row = row+1;
        detail_data(row,:) = {r, vt_names_full{vt}, k, co, round(d*10)/10, round(va*10)/10, ...
            as, tws, ss, round(wm), round(lm), sts, ls, round(dw_c), round(cl), ...
            round(sf*100)/100, round(sk*100)/100};
        t = lv; curr = c;
    end
    db = dist_matrix(curr+1,1); dtb = travel_time_tv(t,db); ret = t+dtb;
    vb = db/max(dtb/60,0.001); vb = max(5,min(60,vb));
    if cfg(7)==1, bcb=(0.0025*vb^2-0.2554*vb+31.75)/100;
    else, bcb=(0.0014*vb^2-0.12*vb+36.19)/100; end
    rf2 = bcb*db*cfg(6); rk2 = bcb*db*cfg(5);
    rts = sprintf('%d:%02d',floor(ret/60),round(mod(ret,60)));
    row = row+1;
    detail_data(row,:) = {r, vt_names_full{vt}, length(path)+1, 0, round(db*10)/10, ...
        round(vb*10)/10, rts, '-', '返回仓库', 0, 0, '-', rts, 0, 0, ...
        round(rf2*100)/100, round(rk2*100)/100};
end
h2 = {'车辆编号','车型','访问序号','客户编号','本段距离km','本段均速kmh',...
    '到达时间','客户时间窗','时间窗状态','等待min','迟到min',...
    '开始服务','离开时间','配送量kg','车上剩余kg','本段能耗成本元','本段碳排放kg'};
writetable(cell2table(detail_data,'VariableNames',h2), excel_file, 'Sheet','逐站点行程明细');

% ---- Sheet3: 总体成本汇总 ----
cs = {'总配送成本(元)',round(best_cost*100)/100;
    '固定启动成本(元)',round(best_detail(1)*100)/100;
    '能耗成本(元)',round(best_detail(2)*100)/100;
    '碳排放成本(元)',round(best_detail(3)*100)/100;
    '早到等待惩罚(元)',round(best_detail(4)*100)/100;
    '晚到惩罚(元)',round(best_detail(5)*100)/100;
    '碳排放总量(kg)',round(best_carbon*100)/100;
    '派车总数(辆)',length(best_sol);
    '燃油大(辆)',vtype_count(1);'燃油中(辆)',vtype_count(2);
    '燃油小(辆)',vtype_count(3);'电动大(辆)',vtype_count(4);'电动小(辆)',vtype_count(5)};
writetable(cell2table(cs,'VariableNames',{'指标','数值'}), excel_file, 'Sheet','总体成本汇总');

fprintf('>>> 已导出: %s (3个Sheet, %d行明细)\n', excel_file, row);

%% 9. 图5: 代表性路径甘特图 (前5辆车的时间线)
n_show = min(5, length(best_sol));
figure('Color','w','Position',[100 50 1100 500]); hold on;

colors_tw = [0.85 0.93 0.85; 0.93 0.85 0.85; 0.85 0.85 0.95];
bar_h = 0.35;

for ri = 1:n_show
    r = ri;
    vt = best_sol(r).vtype;
    path = best_sol(r).path;
    y_base = n_show - ri + 1;

    curr = 0; t = T_START;
    for k = 1:length(path)
        c = path(k);
        dt = travel_time_tv(t, dist_matrix(curr+1, c+1));
        arr_t = t + dt;

        fill([tw(c,1) tw(c,2) tw(c,2) tw(c,1)], ...
            [y_base-bar_h y_base-bar_h y_base+bar_h y_base+bar_h], ...
            [0.9 0.95 0.9], 'EdgeColor','none', 'FaceAlpha',0.4);

        if arr_t < tw(c,1)
            fill([arr_t tw(c,1) tw(c,1) arr_t], ...
                [y_base-0.15 y_base-0.15 y_base+0.15 y_base+0.15], ...
                [1 0.85 0.4], 'EdgeColor','none');
            start_t = tw(c,1);
        elseif arr_t > tw(c,2)
            plot([tw(c,2) arr_t], [y_base y_base], 'r-', 'LineWidth',2.5);
            start_t = arr_t;
        else
            start_t = arr_t;
        end

        fill([start_t start_t+SERVICE_TIME start_t+SERVICE_TIME start_t], ...
            [y_base-0.2 y_base-0.2 y_base+0.2 y_base+0.2], ...
            [0.2 0.6 0.85], 'EdgeColor',[0.1 0.3 0.5]);

        if k <= 6
            text(start_t+SERVICE_TIME/2, y_base, num2str(parent_map(c)), ...
                'HorizontalAlignment','center', 'FontSize',7, 'Color','w', 'FontWeight','bold');
        end

        t = start_t + SERVICE_TIME;
        curr = c;
    end
end

set(gca, 'YTick', 1:n_show, 'YTickLabel', flip(arrayfun(@(i) ...
    sprintf('车辆%d(%s)', i, vt_names_short{best_sol(i).vtype}), 1:n_show, 'Uni',0)));
xlabel('时间 (分钟, 480=8:00)');
title('图5: 代表性车辆执行时间线 (蓝=服务, 黄=等待, 红=迟到, 绿底=时间窗)', ...
    'FontWeight','bold','FontSize',12);
set(gca,'FontName','Microsoft YaHei');
xlim([460 1100]); grid on;

xt = 480:60:1080;
xt_labels = arrayfun(@(m) sprintf('%d:00', floor(m/60)), xt, 'Uni',0);
set(gca, 'XTick', xt, 'XTickLabel', xt_labels);

%% ==================== 辅助函数 ====================

function sol = greedy_init(n_cust, dw, dv, tw, D, VT, t0, svc)
    sol = struct('path',{},'vtype',{});
    unvisited = 1:n_cust;
    stock = VT(:,8)';

    while ~isempty(unvisited)
        avail_types = find(stock > 0);
        if isempty(avail_types)
            error('车辆库存不足，无法为所有客户分配车辆');
        end

        feasible_count = zeros(size(unvisited));
        demand_score = zeros(size(unvisited));
        for ii = 1:length(unvisited)
            c = unvisited(ii);
            feasible_count(ii) = sum(dw(c) <= VT(avail_types,1) & dv(c) <= VT(avail_types,2));
            demand_score(ii) = dw(c)/max(VT(:,1)) + dv(c)/max(VT(:,2));
        end
        if any(feasible_count == 0)
            bad = unvisited(find(feasible_count==0,1));
            error('车辆库存不足或剩余车辆容量不足，客户%d无法分配: %.1fkg, %.2fm3', bad, dw(bad), dv(bad));
        end
        [~, ord] = sortrows([feasible_count(:), -demand_score(:)], [1, 2]);
        seed = unvisited(ord(1));

        candidate_vt = avail_types(dw(seed) <= VT(avail_types,1) & dv(seed) <= VT(avail_types,2));
        [~, cap_ord] = sortrows([VT(candidate_vt,1), VT(candidate_vt,2)], [1, 2]);
        best_vt = candidate_vt(cap_ord(1));

        path = seed;
        load_w = dw(seed); load_v = dv(seed); curr = seed; t = t0;
        dt0 = travel_time_tv(t, D(1, seed+1));
        t = t + dt0;
        if tw(seed,1)>0 && t < tw(seed,1), t = tw(seed,1); end
        t = t + svc;
        remain = setdiff(unvisited, seed);

        while true
            best_next = 0; best_inc = inf;
            for j = 1:length(remain)
                c = remain(j);
                if load_w + dw(c) > VT(best_vt,1) || load_v + dv(c) > VT(best_vt,2), continue; end
                dt = travel_time_tv(t, D(curr+1, c+1));
                arr = t + dt;
                inc = D(curr+1, c+1);
                if tw(c,1)>0 && arr < tw(c,1), inc = inc + (tw(c,1)-arr)*0.3; end
                if tw(c,2)>0 && arr > tw(c,2), inc = inc + (arr-tw(c,2))*3; end
                if inc < best_inc, best_inc = inc; best_next = j; end
            end
            if best_next == 0, break; end
            c = remain(best_next);
            dt = travel_time_tv(t, D(curr+1, c+1));
            t = t + dt;
            if tw(c,1)>0 && t < tw(c,1), t = tw(c,1); end
            t = t + svc;
            path = [path, c]; load_w = load_w + dw(c); load_v = load_v + dv(c);
            curr = c; remain(best_next) = [];
        end

        sol(end+1).path = path;
        sol(end).vtype = best_vt;
        stock(best_vt) = stock(best_vt) - 1;
        unvisited = setdiff(unvisited, path);
    end
end

function [total, detail, carbon_kg] = eval_solution(sol, D, VT, dw, tw, cp, ep, lp, svc, t0)
    c_fix=0; c_eng=0; c_carb=0; c_early=0; c_late=0; carbon_kg=0;
    for r = 1:length(sol)
        vt = sol(r).vtype; path = sol(r).path;
        cfg = VT(vt,:);
        c_fix = c_fix + cfg(3);
        curr = 0; t = t0;
        total_load = sum(dw(path));
        cur_load = total_load;
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
            load_ratio = min(cur_load / cfg(1), 1.5);
            actual_cons = base_cons * dist * (1 + cfg(4)*load_ratio);
            c_eng = c_eng + actual_cons * cfg(6);
            seg_carbon = actual_cons * cfg(5);
            c_carb = c_carb + seg_carbon * cp;
            carbon_kg = carbon_kg + seg_carbon;
            t = t + dt;
            if tw(c,1)>0 && t < tw(c,1)
                c_early = c_early + (tw(c,1)-t)*ep;
                t = tw(c,1);
            end
            if tw(c,2)>0 && t > tw(c,2)
                c_late = c_late + (t-tw(c,2))*lp;
            end
            t = t + svc;
            curr = c;
            cur_load = cur_load - dw(c);
        end
        dist_back = D(curr+1, 1);
        dt_back = travel_time_tv(t, dist_back);
        v_back = dist_back / max(dt_back/60, 0.001);
        v_back = max(5, min(60, v_back));
        if cfg(7)==1
            bc = (0.0025*v_back^2 - 0.2554*v_back + 31.75)/100;
        else
            bc = (0.0014*v_back^2 - 0.12*v_back + 36.19)/100;
        end
        seg_cb = bc*dist_back*cfg(5);
        c_eng = c_eng + bc*dist_back*cfg(6);
        c_carb = c_carb + seg_cb*cp;
        carbon_kg = carbon_kg + seg_cb;
    end
    detail = [c_fix, c_eng, c_carb, c_early, c_late];
    total = sum(detail);
end

function [removed, partial] = destroy(sol, op, n_rm, D, dw, tw)
    all_custs = [];
    for r=1:length(sol), all_custs = [all_custs, sol(r).path]; end
    switch op
        case 1 % 随机移除
            idx = randperm(length(all_custs), min(n_rm, length(all_custs)));
            removed = all_custs(idx);
        case 2 % 最差移除 (按时间窗紧迫度)
            scores = zeros(size(all_custs));
            for i=1:length(all_custs)
                c = all_custs(i);
                scores(i) = dw(c) / max(tw(c,2)-tw(c,1), 1);
            end
            [~, idx] = sort(scores, 'descend');
            removed = all_custs(idx(1:min(n_rm, length(idx))));
        case 3 % Shaw移除 (相似性)
            seed = all_custs(randi(length(all_custs)));
            dists = zeros(size(all_custs));
            for i=1:length(all_custs)
                c = all_custs(i);
                dists(i) = D(seed+1, c+1) + abs(dw(seed)-dw(c))*0.1;
            end
            [~, idx] = sort(dists);
            removed = all_custs(idx(1:min(n_rm, length(idx))));
    end
    partial = sol;
    for r = 1:length(partial)
        partial(r).path = setdiff(partial(r).path, removed, 'stable');
    end
    partial = partial(arrayfun(@(s) ~isempty(s.path), partial));
end

function sol = repair(partial, removed, op, D, VT, dw, dv, tw, cp, ep, lp, svc, t0)
    sol = partial;
    stock_used = zeros(1, size(VT,1));
    for r=1:length(sol), stock_used(sol(r).vtype) = stock_used(sol(r).vtype)+1; end

    to_insert = removed(randperm(length(removed)));
    for ii = 1:length(to_insert)
        c = to_insert(ii);
        best_cost_inc = inf; best_r = 0; best_pos = 0;

        for r = 1:length(sol)
            vt = sol(r).vtype;
            path = sol(r).path;
            cur_w = sum(dw(path)); cur_v = sum(dv(path));
            if cur_w + dw(c) > VT(vt,1) || cur_v + dv(c) > VT(vt,2), continue; end

            for pos = 1:length(path)+1
                new_path = [path(1:pos-1), c, path(pos:end)];
                [nc,~] = eval_single_route(new_path, vt, D, VT, dw, tw, cp, ep, lp, svc, t0);
                [oc,~] = eval_single_route(path, vt, D, VT, dw, tw, cp, ep, lp, svc, t0);
                inc = nc - oc;
                if op == 2 % regret考虑次优
                    inc = inc - 0.1*rand()*abs(inc);
                end
                if inc < best_cost_inc
                    best_cost_inc = inc; best_r = r; best_pos = pos;
                end
            end
        end

        if best_r > 0 && best_cost_inc < 1e8
            sol(best_r).path = [sol(best_r).path(1:best_pos-1), c, sol(best_r).path(best_pos:end)];
        else
            vt = pick_vehicle(dw(c), dv(c), VT, stock_used);
            sol(end+1).path = c;
            sol(end).vtype = vt;
            stock_used(vt) = stock_used(vt) + 1;
        end
    end

    % 路径合并: 尝试将短路径并入其他路径以节省固定成本
    merged = true;
    while merged
        merged = false;
        for i = length(sol):-1:1
            if length(sol(i).path) > 3, continue; end
            for j = 1:length(sol)
                if j == i, continue; end
                vt_j = sol(j).vtype;
                w_j = sum(dw(sol(j).path)); v_j = sum(dv(sol(j).path));
                w_i = sum(dw(sol(i).path)); v_i = sum(dv(sol(i).path));
                if w_j + w_i > VT(vt_j,1) || v_j + v_i > VT(vt_j,2), continue; end
                merged_path = [sol(j).path, sol(i).path];
                [mc,~] = eval_single_route(merged_path, vt_j, D, VT, dw, tw, cp, ep, lp, svc, t0);
                [oc_j,~] = eval_single_route(sol(j).path, vt_j, D, VT, dw, tw, cp, ep, lp, svc, t0);
                [oc_i,~] = eval_single_route(sol(i).path, sol(i).vtype, D, VT, dw, tw, cp, ep, lp, svc, t0);
                if mc < oc_j + oc_i
                    sol(j).path = merged_path;
                    sol(i) = [];
                    merged = true;
                    break;
                end
            end
            if merged, break; end
        end
    end
end

function vt = pick_vehicle(w, v, VT, used)
    for k = [5, 4, 3, 2, 1]
        if w <= VT(k,1) && v <= VT(k,2) && used(k) < VT(k,8)
            vt = k; return;
        end
    end
    error('车辆库存不足或剩余车辆容量不足，无法分配 %.1fkg, %.2fm3 的需求', w, v);
end

function ok = check_solution_feasible(sol, VT, dw, dv)
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

function [cost, detail] = eval_single_route(path, vt, D, VT, dw, tw, cp, ep, lp, svc, t0)
    cfg = VT(vt,:);
    c_eng=0; c_carb=0; c_early=0; c_late=0;
    curr=0; t=t0;
    total_load = sum(dw(path)); cur_load = total_load;
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
    cost = sum(detail);
end

function idx = roulette_select(w)
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
