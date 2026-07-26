# MetabolicDesigner Pro (MDP)

基于约束的代谢网络 in silico 菌株设计与仿真平台。

## 功能

- **FSEOF** — 通量扫描，识别与目标产物耦合的反应（OE/KD 候选）
- **OptKnock** — 基因敲除策略，维持最低生长并提升目标产物
- **optForce** — 强制通量改变策略（OE/KD/KO）

## 快速开始

```matlab
addpath(genpath('src'));

% 单算法调试
model = loadModel('input_models/iHM1533/iHM1533_heparosan.mat');
result = fseof(model, 'EX_heparosan_e', ...
    'BIOMASS_EcN_iHM1533_core_59p80M', ...
    'Iterations', 21, 'Coefficient', 0.99);

% 批量设计
spec = { ...
    struct('algo','fseof','params',{{'Iterations',21,'Coefficient',0.99}}), ...
    struct('algo','optknock','params',{{'MaxCandidates',150,'NumDel',3}}), ...
    struct('algo','optforce','params',{{'K',3,'NSets',2}}) ...
};
combined = strainDesign( ...
    'input_models/iHM1533/iHM1533_heparosan.mat', ...
    'EX_heparosan_e', ...
    'BIOMASS_EcN_iHM1533_core_59p80M', ...
    spec);
```

## 目录结构

```
src/
├── loadModel.m            加载模型 + ec-model 嗅探
├── strainDesign.m         批量入口 + workspace 管理
└── algorithms/
    ├── fseof.m            FSEOF 算法
    ├── optknock.m         OptKnock 算法
    └── optforce.m        optForce 算法
```

## 依赖

- MATLAB R2019b+
- COBRA Toolbox（FBA、OptKnock、optForce）
- RAVEN Toolbox（可选，FSEOF 候选检测）

## Workspace

运行结果保存在 `workspaces/` 目录：

```
workspaces/{模型名}_{目标反应}/
├── combined_result.mat    所有算法汇总
├── FSEOF/result.mat       FSEOF 结果
├── OptKnock/result.mat    OptKnock 结果
└── optForce/result.mat    optForce 结果
```

重跑某算法只覆盖该算法结果，其他保留。

## 算法参数

| 算法 | 参数 | 默认值 | 说明 |
|------|------|--------|------|
| FSEOF | Iterations | 10 | 扫描步数 |
| FSEOF | Coefficient | 0.9 | 目标最大通量比例（< 1） |
| OptKnock | MaxCandidates | 200 | 最大候选解数 |
| OptKnock | NumDel | 5 | 最大敲除数 |
| OptKnock | MinGrowthFraction | 0.1 | 最低生长比例 |
| OptKnock | VMax | 1000 | 最大通量边界 |
| optForce | K | 2 | 每组反应数 |
| optForce | NSets | 1 | 搜索组数 |
| optForce | MaxCandidates | 500 | 最大候选反应数 |
