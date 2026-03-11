# atmBoundaryLayerMapped Suite

大气边界层映射边界条件套件 - 基于实测数据的大气边界层入口边界条件实现

## 概述

本套件基于OpenFOAM 2112，实现了一套严格配套使用的映射边界条件。设计原则是：**只有U从外部读取数据，k/ε/ω通过ABL公式从U计算得出**，确保数据一致性。

### 核心特点

- **严格配套**：k/ε/ω必须与`atmBoundaryLayerMappedVelocity`一起使用，否则报错
- **自动计算湍流**：根据实际U值计算u*，再用ABL公式计算k/ε/ω
- **ABL对数律缩放**：U的mapped数据会经过对数律缩放
- **极简配置**：k/ε/ω边界只需要指定`type`和`phi`，其他参数从U复制

## 数据流

```
boundaryData/
├── East/
│   ├── points          # 采样点坐标
│   └── 0/
│       └── U           # ← 唯一的外部数据源
│
计算流程：
U (mapped + ABL缩放) → u*计算 → k/ε/ω (ABL公式)
                     ↓
              UstarFromU(): u* = κ·U/ln((z-d+z₀)/z₀)
                     ↓
        k = u*²/√Cμ · √(C₁·ln(...) + C₂)
        ε = u*³/(κ·z) · √(C₁·ln(...) + C₂)  
        ω = u*/(κ·√Cμ·z)
```

### 计算公式

#### 1. 摩擦速度 (Friction Velocity)
```
u* = κ · U / ln((z - d + z₀)/z₀)
```
其中：
- κ: von Kármán常数 (默认0.41)
- z: 高度
- d: 位移高度
- z₀: 粗糙度长度
- U: 流向速度分量

#### 2. 湍动能 (k)
```
k = (u*)²/√Cμ · √(C₁·ln((z-d+z₀)/z₀) + C₂)
```

#### 3. 湍流耗散率 (ε)
```
ε = (u*)³/(κ·z) · √(C₁·ln((z-d+z₀)/z₀) + C₂)
```

#### 4. 比耗散率 (ω)
```
ω = u*/(κ·√Cμ·z)
```

默认常数：
- Cμ = 0.09
- C₁ = 0.0
- C₂ = 1.0

## 用户设置

### 1. 准备边界数据

创建`constant/boundaryData/<patchName>/`目录：

```bash
constant/boundaryData/East/
├── points              # 采样点坐标 (必须)
└── 0/
    └── U               # 速度数据 (必须)
```

**points文件格式：**
```cpp
// Points
3
(
(10 10 10)
(-10 10 10)
(10 -10 10)
)
```

**U文件格式：**
```cpp
// Data on points
3
(
(10 0 0)
(12 0 0)
(15 0 0)
)
```

### 2. 配置边界条件

#### 0/U（必须配置）

```cpp
boundaryField
{
    "(North|West|South|East)"
    {
        type            atmBoundaryLayerMappedVelocity;
        
        // ABL参数
        zDir            (0 0 1);        // 垂直方向
        Zref            10;             // 参考高度
        z0              uniform 1;      // 粗糙度长度
        d               uniform 0;      // 位移高度
        
        // 可选：ABL常数
        kappa           0.41;           // von Kármán常数
        Cmu             0.09;
        C1              0.0;
        C2              1.0;
    }
}
```

#### 0/k, 0/epsilon, 0/omega（极简配置）

```cpp
boundaryField
{
    "(North|West|South|East)"
    {
        type            atmBoundaryLayerMappedK;  // 或MappedEpsilon/MappedOmega
        phi             phi;                      // 可选，默认为"phi"
    }
}
```

**注意**：k/ε/ω边界会自动从同patch的U边界复制zDir、Zref、z0、d等参数。

### 3. 加载库

在`system/controlDict`中添加：

```cpp
libs ("libatmBoundaryLayerMapped.so");
```

### 4. 参数说明

| 参数 | U BC | k/ε/ω BC | 说明 |
|------|------|----------|------|
| `type` | 必需 | 必需 | 边界类型 |
| `zDir` | 必需 | 自动复制 | 垂直方向向量 |
| `Zref` | 必需 | 自动复制 | 参考高度(m) |
| `z0` | 必需 | 自动复制 | 粗糙度长度(m) |
| `d` | 必需 | 自动复制 | 位移高度(m) |
| `kappa` | 可选(0.41) | 自动复制 | von Kármán常数 |
| `Cmu` | 可选(0.09) | 自动复制 | 模型常数 |
| `C1` | 可选(0.0) | 自动复制 | ABL剖面常数1 |
| `C2` | 可选(1.0) | 自动复制 | ABL剖面常数2 |

## 程序架构

### 类继承关系

```
Base Classes:
├── atmBoundaryLayerMapped (基类，包含ABL计算和映射逻辑)
│   ├── PatchFunction1<vector> UMapper_       // U的mapped数据
│   ├── PatchFunction1<scalar> scalarMapper_  // 标量的mapped数据
│   ├── Function1<vector> zDir_               // 垂直方向
│   ├── Function1<scalar> Zref_               // 参考高度
│   ├── PatchFunction1<scalar> z0_            // 粗糙度
│   └── PatchFunction1<scalar> d_             // 位移高度
│
Derived Classes:
├── atmBoundaryLayerMappedVelocityFvPatchVectorField
│   └── 继承: inletOutletFvPatchVectorField + atmBoundaryLayerMapped
│   └── updateCoeffs(): 使用Umapped()获取带ABL缩放的速度
│
├── atmBoundaryLayerMappedKFvPatchScalarField
│   └── 继承: inletOutletFvPatchScalarField + atmBoundaryLayerMapped
│   └── updateCoeffs(): 从U计算u*，再用kFromUstar()计算k
│
├── atmBoundaryLayerMappedEpsilonFvPatchScalarField
│   └── 继承: inletOutletFvPatchScalarField + atmBoundaryLayerMapped
│   └── updateCoeffs(): 从U计算u*，再用epsilonFromUstar()计算ε
│
└── atmBoundaryLayerMappedOmegaFvPatchScalarField
    └── 继承: inletOutletFvPatchScalarField + atmBoundaryLayerMapped
    └── updateCoeffs(): 从U计算u*，再用omegaFromUstar()计算ω
```

### 关键函数

#### 基类 (atmBoundaryLayerMapped)

```cpp
// 从mapped数据获取速度（带ABL对数律缩放）
tmp<vectorField> Umapped(const vectorField& pCf) const;

// 从实际U值计算摩擦速度
tmp<scalarField> UstarFromU(const vectorField& Uvalues, const vectorField& pCf) const;

// 从u*计算湍流参数
tmp<scalarField> kFromUstar(const scalarField& uStar, const vectorField& pCf) const;
tmp<scalarField> epsilonFromUstar(const scalarField& uStar, const vectorField& pCf) const;
tmp<scalarField> omegaFromUstar(const scalarField& uStar, const vectorField& pCf) const;
```

#### 派生类updateCoeffs流程

**Velocity BC:**
```cpp
void updateCoeffs()
{
    // 1. 从boundaryData读取原始U
    // 2. 应用ABL对数律缩放
    refValue() = Umapped(patch().Cf());
    inletOutletFvPatchVectorField::updateCoeffs();
}
```

**K/Epsilon/Omega BC:**
```cpp
void updateCoeffs()
{
    // 1. 严格检查：U必须使用atmBoundaryLayerMappedVelocity
    checkUCompatiblity();
    
    // 2. 获取U的patch内部场（值拷贝）
    const vectorField Uvalues(Upatch.patchInternalField());
    
    // 3. 计算摩擦速度u*
    tmp<scalarField> uStar = UstarFromU(Uvalues, patch().Cf());
    
    // 4. 使用ABL公式计算湍流参数
    refValue() = kFromUstar(uStar(), patch().Cf());
    
    inletOutletFvPatchScalarField::updateCoeffs();
}
```

### 安全机制

1. **类型检查**：k/ε/ω边界严格检查U边界类型，不匹配则报错
2. **空指针检查**：基类函数检查所有Function1/PatchFunction1是否初始化
3. **值拷贝**：使用`const vectorField Uvalues(...)`而非引用，避免内存问题

## 注意事项

1. **配套使用**：k/ε/ω必须与`atmBoundaryLayerMappedVelocity`同patch使用
2. **数据一致性**：只有U读取外部数据，湍流参数通过公式计算保证与U一致
3. **流向计算**：流向自动从mapped U的平均方向计算，无需手动指定
4. **边界数据**：只需提供U的boundaryData，k/ε/ω不需要数据文件

## 文件清单

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
├── Make/
│   ├── files
│   └── options
└── README.md
```

## 许可证

GNU GPL v3 (与OpenFOAM保持一致)

## 更新日志

### v1.0 (Current)
- 实现严格的配套检查机制
- U从boundaryData读取并应用ABL缩放
- k/ε/ω通过ABL公式从U自动计算
- 流向自动计算，无需flowDir参数
