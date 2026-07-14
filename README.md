# Huazhong Cup 2026 Green Logistics Optimization

第十八届华中杯大学生数学建模挑战赛 A 题作品：**静态规划-限行合规-动态扰动：低碳配送优化研究**。

本项目围绕城市绿色物流配送调度问题，构建了考虑时变路况、异构车队、软时间窗、绿色配送区限行和动态订单扰动的车辆路径优化模型，并使用改进自适应大邻域搜索算法进行求解。

## Project Highlights

- **Problem 1:** TD-HVRPTW 静态基准调度模型，综合考虑固定成本、动态能耗、碳排放成本与时间窗惩罚。
- **Problem 2:** 绿色限行政策下的双目标优化模型，使用 epsilon-约束法构建成本与碳排放的 Pareto 前沿。
- **Problem 3:** 事件触发式局部弹性重调度策略，用于订单取消、路径扰动等动态场景。
- **Algorithm:** Adaptive Large Neighborhood Search (ALNS) with domain-specific destroy and repair operators.

## Repository Structure

```text
.
├── paper/      # Final paper PDF
├── latex/      # LaTeX source and figures used in the paper
├── code/       # MATLAB scripts for Problem 1, 2 and 3
├── data/       # Input data from the contest attachment
└── results/    # Complete vehicle scheduling result tables
```

## Files

- `paper/A202611901279.pdf`: final submitted paper.
- `latex/main.tex`: main LaTeX source file.
- `code/q1.m`: static TD-HVRPTW scheduling model and ALNS solver.
- `code/q2.m`: green policy constrained bi-objective scheduling model.
- `code/q3.m`: event-triggered local rescheduling model.
- `data/`: customer coordinates, time windows, orders and distance matrix.
- `results/`: complete scheduling plans for Problem 1 and Problem 2.

## How to Run

The code is written in MATLAB. Recommended environment:

- MATLAB R2021a or later
- Statistics and Machine Learning Toolbox is recommended for table/data utilities

Run each problem independently from the `code/` directory. The scripts expect the Excel data files in the current MATLAB path. If needed, copy the files under `data/` into the MATLAB working directory or add `data/` to the MATLAB path.

```matlab
run('q1.m')
run('q2.m')
run('q3.m')
```

## Notes

This repository is for academic portfolio display and learning exchange. Personal submission materials, commitment forms, temporary files and duplicate archives are intentionally excluded.
