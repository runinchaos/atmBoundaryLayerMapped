# atmBoundaryLayerMapped Suite

大气边界层映射边界条件套件 - 支持时变映射数据的完整ABL边界条件实现

## 概述

本套件基于OpenFOAM 2112的`atmBoundaryLayer`家族，增加了**时变映射数据支持**（time-varying mapped data），用于从外部数据源（如实测数据、WRF模式输出等）驱动大气边界层模拟。

## 架构

```
atmBoundaryLayerMapped (基类)
├── atmBoundaryLayerMappedVelocity (速度场，支持映射)
├── atmBoundaryLayerMappedK (湍动能，基于U计算)
├── atmBoundaryLayerMappedEpsilon (耗散率，基于U计算)
└── atmBoundaryLayerMappedOmega (比耗散率，基于U计算)
```

## 核心特性

### 1. Velocity - 时变映射支持
- **类型**: `atmBoundaryLayerMappedVelocity`
- **功能**: 从`boundaryData`读取时变速度数据
- **ABL缩放**: 对mapped数据应用对数律缩放
- **配置**:
```cpp
inlet
{
    type            atmBoundaryLayerMappedVelocity;
    zDir            (0 0 1);
    Zref            10.0;
    z0              uniform 0.1;
    d               uniform 0.0;
    
    // 可选：ABL参数
    kappa           0.41;
    Cmu             0.09;
    
    // 映射数据源在 constant/boundaryData/<patchName>/
}
```

### 2. K/Epsilon/Omega - 与Velocity保持一致
- **自动检测**: 检查同patch的U是否使用`atmBoundaryLayerMappedVelocity`
- **两种模式**:
  - **Mapped模式**: 如果U是mapped的，根据**实际U值**计算u*，再用u*计算k/ε/ω
  - **理论模式**: 如果U不是mapped的，使用ABL理论公式（基于Uref）

#### 计算逻辑
```cpp
if (U使用atmBoundaryLayerMappedVelocity) {
    // 从mapped U计算实际的摩擦速度
    u* = kappa * U_mapped / ln((z - d + z0)/z0)
    
    // 用实际u*计算k/ε/ω
    k = (u*)^2 / sqrt(Cμ) * sqrt(C1*ln(...) + C2)
    ε = (u*)^3 / (kappa*z) * sqrt(C1*ln(...) + C2)
    ω = u* / (kappa*sqrt(Cμ)*z)
} else {
    // 使用理论公式（基于Uref）
    k = k(patch().Cf())
    ε = epsilon(patch().Cf())
    ω = omega(patch().Cf())
}
```

## 文件结构

```
2_DEV/atmBoundaryLayerMapped/
├── atmBoundaryLayerMapped/
│   ├── atmBoundaryLayerMapped.H          # 基类声明
│   └── atmBoundaryLayerMapped.C          # 基类实现
├── atmBoundaryLayerMappedVelocity/
│   ├── atmBoundaryLayerMappedVelocityFvPatchVectorField.H
│   └── atmBoundaryLayerMappedVelocityFvPatchVectorField.C
├── atmBoundaryLayerMappedK/
│   ├── atmBoundaryLayerMappedKFvPatchScalarField.H
│   └── atmBoundaryLayerMappedKFvPatchScalarField.C
├── atmBoundaryLayerMappedEpsilon/
│   ├── atmBoundaryLayerMappedEpsilonFvPatchScalarField.H
│   └── atmBoundaryLayerMappedEpsilonFvPatchScalarField.C
├── atmBoundaryLayerMappedOmega/
│   ├── atmBoundaryLayerMappedOmegaFvPatchScalarField.H
│   └── atmBoundaryLayerMappedOmegaFvPatchScalarField.C
└── Make/
    ├── files
    └── options
```

## 编译

```bash
cd 2_DEV/atmBoundaryLayerMapped
wmake
```

库文件将生成在: `$FOAM_USER_LIBBIN/libatmBoundaryLayerMapped.so`

## 使用示例

### 1. 准备映射数据

在`constant/boundaryData/East/`目录下：
```
East/
├── points              # 采样点坐标
└── 0/
    └── U               # 时速度数据
```

U文件格式示例：
```cpp
// Data on points
3
(
(10 0 0)
(12 0 0)
(15 0 0)
)
```

### 2. 配置U边界条件

`0/U`:
```cpp
inlet
{
    type            atmBoundaryLayerMappedVelocity;
    zDir            (0 0 1);
    Zref            10.0;      // 参考高度
    z0              uniform 0.1;  // 粗糙度长度
    d               uniform 0.0;  // 位移高度
    
    // 自动读取 constant/boundaryData/inlet/0/U
}
```

### 3. 配置K/Epsilon/Omega

**无需特殊配置**，自动检测并使用mapped U：

`0/k`:
```cpp
inlet
{
    type            atmBoundaryLayerMappedK;
    zDir            (0 0 1);
    Zref            10.0;
    z0              uniform 0.1;
    d               uniform 0.0;
    
    // 自动检测U是否为mapped，并据此计算
}
```

`0/epsilon`和`0/omega`类似。

### 4. 加载库

`system/controlDict`:
```cpp
libs ("libatmBoundaryLayerMapped.so");
```

## 关键API

### 基类新增函数

```cpp
// 从实际U值计算摩擦速度
tmp<scalarField> UstarFromU
(
    const vectorField& Uvalues,
    const vectorField& pCf
) const;

// 从给定u*计算k
tmp<scalarField> kFromUstar
(
    const scalarField& uStar,
    const vectorField& pCf
) const;

// 从给定u*计算epsilon
tmp<scalarField> epsilonFromUstar
(
    const scalarField& uStar,
    const vectorField& pCf
) const;

// 从给定u*计算omega
tmp<scalarField> omegaFromUstar
(
    const scalarField& uStar,
    const vectorField& pCf
) const;

// 返回映射速度（带ABL缩放）
tmp<vectorField> Umapped(const vectorField& pCf) const;
```

## 与传统ABL BC的区别

| 特性 | atmBoundaryLayerInlet* | atmBoundaryLayerMapped* |
|------|------------------------|-------------------------|
| U数据源 | 理论公式 (Uref) | boundaryData (时变) |
| k/ε/ω计算 | 基于Uref | 基于实际mapped U |
| 数据一致性 | U与k/ε/ω可能不一致 | U与k/ε/ω保持一致 |
| 适用场景 | 理想ABL | 实测数据驱动 |

## 注意事项

1. **数据一致性**: K/Epsilon/Omega会自动检测U的BC类型，无需手动配置
2. **z0和d**: 所有BC使用相同的z0和d参数，确保一致性
3. **时变数据**: 目前支持空间分布映射，时变功能可通过`timeVaryingMappedFixedValue`方式扩展
4. **兼容性**: 与OpenFOAM 2112及更高版本兼容

## 测试案例

见 `3_SIM/atmBCTable_forAgent/` 目录下的测试案例。

## 更新日志

### v1.0
- 初始版本：fork OpenFOAM atmBoundaryLayer家族
- 重命名所有类：atmBoundaryLayer* → atmBoundaryLayerMapped*

### v1.1
- 添加mapping基础功能：Umapped(), PatchFunction1支持
- 添加useMapping标志

### v1.2
- K/Epsilon/Omega自动检测mapped U
- 添加UstarFromU(), *FromUstar()辅助函数
- 实现k/ε/ω与mapped U的一致性计算

## 许可证

GNU GPL v3 (与OpenFOAM保持一致)

## 作者

基于OpenFOAM 2112官方实现修改
