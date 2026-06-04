# TOE超导体3D建模 — 参数详细介绍SOP

> 版本：v1.0 | 日期：2026-06-04  
> 适用文件：`porosity_calculator.py` + `superconductor_gyroid.py`  
> 理论依据：TOE_v8.8框架

---

## 目录

1. [工作流总览](#1-工作流总览)
2. [porosity_calculator.py 参数详解](#2-porosity_calculatorpy-参数详解)
3. [superconductor_gyroid.py 参数详解](#3-superconductor_gyroidpy-参数详解)
4. [材料数据库速查表](#4-材料数据库速查表)
5. [典型使用场景SOP](#5-典型使用场景sop)
6. [常见问题与排错](#6-常见问题与排错)

---

## 1. 工作流总览

```
┌─────────────────────────┐
│  Step 1: 选材料/算孔隙率  │  porosity_calculator.py
│  → 输出 target_porosity  │  --material 或 --density
│  → 输出 target_df        │  --d_f 2.973
└──────────┬──────────────┘
           │ 把输出的参数填入 ↓
┌──────────▼──────────────┐
│  Step 2: 生成3D模型       │  superconductor_gyroid.py
│  → Gyroid骨架 + 分形空隙  │  --porosity --df
│  → 输出STL文件            │  --void-mode --material
└─────────────────────────┘
```

**核心逻辑**：先算物理参数，再用物理参数驱动几何生成。

---

## 2. porosity_calculator.py 参数详解

### 2.1 核心公式链（理解参数的前提）

| 序号 | 公式 | 物理含义 |
|------|------|----------|
| 1 | η_Λ = 3 - d_f = 0.027 | 分形维度亏缺（非拟合参数，来自宇宙大尺度结构观测） |
| 2 | d_f,sup = 2.973 | 最优分形维度（由D=13大统一闭包导出） |
| 3 | d_f,sup = 2 + d_max/100 | 100倍标度关系（d_max=95→d_f=2.95中心值） |
| 4 | φ_总 = φ_本征 + φ_挖孔 = 69.7% = Ω_Λ | 超导总孔隙率 = 宇宙孔隙率 |
| 5 | φ_本征 = 1 - f_m = 1 - ρ_bulk/ρ_theoretical | 材料自带孔隙（§1.3公式） |
| 6 | φ_挖孔 = 69.7% - φ_本征 | 人工补足至宇宙孔隙率 |
| 7 | λ_eff = λ₀ + C_λ × P × √(1-P) | 分形空位耦合增强 |
| 8 | μ\*_eff = μ\*₀ / (1 + C_μ × P) | 分形结构库仑屏蔽 |
| 9 | T_c = McMillan(Θ_D, λ_eff, μ\*_eff) | 超导临界温度 |

### 2.2 命令行参数

#### 2.2.1 材料选择类

| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `--material` | str | 无 | 从材料数据库选择，如 `graphene_superlattice`、`hydrogenated_copper_oxide` |
| `--preset` | str | None | 已验证预设名（如 `hydrogenated_copper_oxide_vacuum`），自动填充全部校验参数 |
| `--density` | float | 无 | 自定义材料本征密度 (g/cm³)，不与 `--material` 同时使用 |
| `--density-theoretical` | float | 无 | X射线理论密度 (g/cm³)，用于计算 φ_0。不提供则用 `--density` 值 |
| `--list-materials` | flag | — | 列出所有可用材料及其关键参数 |

**使用原则**：
- 有预设材料 → 优先用 `--material`
- 有预设验证方案 → 用 `--preset`（自动覆盖其他参数）
- 新材料无数据库 → 用 `--density` + 手动补 `--theta_D` 等

#### 2.2.2 超导物理参数

| 参数 | 类型 | 默认值 | 量纲 | 物理含义 | 调参建议 |
|------|------|--------|------|----------|----------|
| `--theta_D` | float | 材料库/自动 | K | 德拜温度。原子振动频率的特征温度 | 越高Tc越高；轻元素（C、B、Li）>重元素 |
| `--lambda_0` | float | 材料库/自动 | 无量纲 | 电声耦合常数（本征值） | 石墨烯~1.5，CuO~1.2，金属0.1-0.5 |
| `--mu_star_0` | float | 材料库/自动 | 无量纲 | 库仑赝势（本征值） | 大多数材料0.10-0.15；越低越有利于超导 |
| `--C_lambda` | float | 材料库/自动 | 无量纲 | 分形耦合增强系数。控制 Δλ = C_λ × P × √(1-P) | 越大→空位对耦合增强越强；Θ_D越高C_λ越大 |
| `--C_mu` | float | 材料库/自动 | 无量纲 | 分形屏蔽系数。控制 μ\*_eff = μ\*₀/(1+C_μ×P) | 越大→库仑屏蔽越强→μ\*越低→越有利于超导 |
| `--d_f` | float | 2.973 | 无量纲 | 目标分形维度 | TOE框架固定值2.973，范围2.92-2.98，**不建议修改** |

**参数关联图**：

```
Θ_D 高 ──→ C_λ 大 ──→ Δλ 大 ──→ λ_eff 高 ──→ Tc 高
                                    ↑
              P 适当 ──→ P×√(1-P) 最大 ──┘

Θ_D 高 ──→ C_μ 大 ──→ μ*_eff 低 ──→ Tc 高

P → 1: Θ_D_eff → 0 (材料失去连续性, 超导消失)
P = 0: 无空位增强, 回归本征值
```

#### 2.2.3 三重修正参数（TOE框架§22）

| 参数 | 类型 | 默认值 | 物理含义 | 何时修改 |
|------|------|--------|----------|----------|
| `--phi-air` | float | 0.015 | 空气填充修正率。常压下空气占据挖孔体积的等效比例 | **真空工况设为0**；常压保持1.5% |
| `--phi-0bubble` | float | 0.003 | 0维世界泡填充修正率。极限真空下仍存在的本源真空元 | **不可消除**，保持0.3% |
| `--delta-phi-fractal` | float | 0.008 | 分形凹凸挤占修正率。孔壁分形粗糙面挤占有效空腔 | 表面越粗糙值越大；平滑表面可适当降低 |

**三重修正的逻辑**：

```
φ_挖表 (工程开孔尺寸) = φ_挖实 (有效物理空隙) + 三重修正
                                  ├── φ_空气     (1.5% 常压, 0% 真空)
                                  ├── φ_0泡      (0.3% 不可消除)
                                  └── Δφ_分形    (0.8% 粗糙度)

关键: 忽略修正 → φ_总<69.7% → 失超!
```

#### 2.2.4 输出控制参数

| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `--target-tc` | float | 200.0 | 目标超导临界温度 (K)。计算最优P后对比此目标 |
| `--json` | flag | — | 以JSON格式输出，便于脚本解析 |
| `--quiet` | flag | — | 安静模式，只输出 `target_porosity`、`target_df`、`Tc_max` |
| `--gyroid-params` | flag | — | 输出Gyroid生成器可直接复制使用的命令行 |

### 2.3 材料参数字典内部结构

每个材料条目包含以下字段：

| 字段 | 类型 | 量纲 | 说明 |
|------|------|------|------|
| `name` | str | — | 材料全名 |
| `category` | str | — | 分类标签 |
| `density` | float | g/cm³ | 实际密度（含烧结/挤出微孔） |
| `density_theoretical` | float | g/cm³ | X射线理论密度（完全致密晶格） |
| `Theta_D` | float | K | 德拜温度 |
| `lambda_0` | float | 无量纲 | 本征电声耦合常数 |
| `mu_star_0` | float | 无量纲 | 本征库仑赝势 |
| `C_lambda` | float | 无量纲 | 分形耦合增强系数 |
| `C_mu` | float | 无量纲 | 分形屏蔽系数 |
| `note` | str | — | 备注 |

**density vs density_theoretical 的关键区别**：

- `density` = 实际测量密度（烧结体/挤出体，含有天然微孔）
- `density_theoretical` = X射线衍射算出的理论密度（完全致密晶格）
- 两者之差决定了 **φ_本征**：`φ_0 = 1 - density/density_theoretical`
- 如果两者相等 → φ_本征=0 → 需要全部靠人工挖孔达到69.7%

### 2.4 已验证预设（VERIFIED_PRESETS）

| 预设名 | 材料 | 工况 | 关键结论 |
|--------|------|------|----------|
| `hydrogenated_copper_oxide_vacuum` | 氢化铜氧化物 | 真空(φ_air=0) | Tc_max=83.4K, 缺口116.6K至200K, 需提温优化 |

---

## 3. superconductor_gyroid.py 参数详解

### 3.1 生成算法流程

```
Phase 1: Gyroid骨架生成    ── 阈值控制基础孔隙率
    ↓
Phase 2: 幂律球形空隙生成  ── d_f=2.973 保证
    ↓  (或 分形层级空隙 Menger-Sponge变体)
Phase 3: 空隙雕刻          ── 不破坏Gyroid连通性
    ↓  (迭代式, 多轮直到孔隙率达标)
Phase 4: 验证              ── d_f盒计数 + 连通性检查
    ↓
Phase 5: STL导出           ── marching_cubes + 水密性
```

### 3.2 命令行参数

#### 3.2.1 几何尺寸参数

| 参数 | 类型 | 默认值 | 量纲 | 物理含义 | 调参建议 |
|------|------|--------|------|----------|----------|
| `--length` | float | 15.0 | cm | 线体长度（沿主轴方向） | 最大20cm；越长越费内存 |
| `--width` | float | 1.2 | cm | 线体宽度 | 应≥2×Gyroid波长 |
| `--height` | float | 1.2 | cm | 线体高度 | 应≥2×Gyroid波长 |

**尺寸与内存的关系**：

```
总体素数 = (L/voxel) × (W/voxel) × (H/voxel)
例: 15cm / 0.005cm = 3000 体素/轴
    3000 × 240 × 240 ≈ 173M 体素 → ~0.16GB (bool)
```

#### 3.2.2 物理参数（核心！）

| 参数 | 类型 | 默认值 | 量纲 | 物理含义 | 调参建议 |
|------|------|--------|------|----------|----------|
| `--porosity` | float | 0.6970 | 无量纲 | 目标总孔隙率 | **TOE框架定值0.6970=Ω_Λ**，一般不修改；结构验证可用0.40 |
| `--df` | float | 2.973 | 无量纲 | 目标分形维度 | **TOE框架固定值2.973**，范围2.92-2.98 |

**孔隙率与分形维度的关系**：

- `--porosity` 控制空隙的**总量**（体积占比）
- `--df` 控制空隙的**分布律**（大小分布的幂律指数）
- 两者独立控制，物理上必须同时满足：
  - φ_总 = 69.7%（五大约束锁定）
  - d_f = 2.973（D=13大统一闭包导出）

#### 3.2.3 Gyroid骨架参数

| 参数 | 类型 | 默认值 | 量纲 | 物理含义 | 调参建议 |
|------|------|--------|------|----------|----------|
| `--gyroid-porosity` | float | 0.50 | 无量纲 | Gyroid基础孔隙率 | 值越大壁越薄→空隙负担越小→但壁可能太薄不连通 |
| `--wavelength` | float | 自动 | cm | Gyroid周期波长 | 默认=min(宽,高)/2.5；手动设更大→壁更厚更稳 |

**Gyroid阈值与孔隙率的对应**：

```
threshold < 0 → 孔隙率 < 50%（更实心，壁更厚）
threshold = 0 → 孔隙率 ≈ 50%（标准Gyroid）
threshold > 0 → 孔隙率 > 50%（更多孔，壁更薄）
```

**自适应逻辑**：当 `--porosity > 0.50` 时，程序自动提高 `gyroid_porosity` 至 `target_porosity × 0.65`（上限0.50），以减轻空隙叠加负担。

#### 3.2.4 空隙模式参数

| 参数 | 类型 | 默认值 | 可选值 | 物理含义 |
|------|------|--------|--------|----------|
| `--void-mode` | str | powerlaw | `powerlaw` / `fractal` | 空隙分布算法模式 |

**两种模式对比**：

| 特性 | powerlaw（幂律随机） | fractal（分形层级） |
|------|---------------------|-------------------|
| 空间分布 | 完全随机 | 多尺度自相似网格 |
| 尺寸分布 | p(r) ∝ r^{-(d_f+1)} | 层级绑定（几何级数） |
| d_f保证 | 尺寸分布构造 | 空间递归构造 |
| 空间均匀性 | 可能聚集 | 天然均匀 |
| 物理优势 | 简单直接 | 各向同性+均匀钉扎 |
| 计算速度 | 快 | 较慢（递归） |

| 参数 | 类型 | 默认值 | 量纲 | 适用模式 | 物理含义 |
|------|------|--------|------|----------|----------|
| `--subdivision` | int | 3 | 无量纲 | fractal | 每轴细分数。3→3³=27子块, 4→4³=64子块 |
| `--min-void-mm` | float | 自动 | mm | 两种 | 最小空隙直径。默认=min_feature×1.5 |
| `--max-void-mm` | float | 自动 | mm | 两种 | 最大空隙直径。默认=min(40%波长, 30%截面) |

**subdivision参数详解**：

```
subdivision=3 (默认):
  每层分 3³=27 个子块
  移除概率 p₀ = 1 - 3^(2.973-3) ≈ 1.80%
  → 每个子块仅1.8%概率被移除 → 极稀疏 → d_f≈3 ✓

subdivision=4:
  每层分 4³=64 个子块
  移除概率 p₀ = 1 - 4^(2.973-3) ≈ 1.90%
  → 更细粒度的空隙分布

⚠️ subdivision越大 → 递归越深 → 内存越大 → 计算越慢
```

#### 3.2.5 3D打印控制参数

| 参数 | 类型 | 默认值 | 量纲 | 物理含义 | 调参建议 |
|------|------|--------|------|----------|----------|
| `--min-feature-mm` | float | 0.10 | mm | 最小可打印特征尺寸 | SLA/光固化=0.10；FDM=0.40；SLS=0.20 |
| `--voxel-mm` | float | 0.05 | mm | 体素尺寸（网格分辨率） | 越小越精细但越慢；0.05=50μm适合SLA；0.025=25μm超高精度 |
| `--shell-mm` | float | 0.5 | mm | 外壳厚度（封闭边界） | 保证STL水密和打印边界清晰 |

**体素尺寸 vs 计算时间**：

| voxel_mm | 分辨率 | 适用工艺 | 15×1.2×1.2cm 内存估算 |
|----------|--------|----------|----------------------|
| 0.10 | 100μm | FDM快速验证 | ~0.02GB |
| 0.05 | 50μm | SLA标准 | ~0.16GB |
| 0.025 | 25μm | SLA超高精度 | ~1.3GB |
| 0.01 | 10μm | 微纳尺度 | ~20GB ⚠️ |

#### 3.2.6 验证与修正参数

| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `--no-verify` | flag | — | 跳过d_f盒计数验证，大幅加速（快速调试用） |
| `--phi-air` | float | 0.015 | 空气填充修正率（同porosity_calculator，真空设0） |
| `--phi-0bubble` | float | 0.003 | 0维世界泡修正率 |
| `--delta-phi-fractal` | float | 0.008 | 分形凹凸挤占修正率 |

#### 3.2.7 其他参数

| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `--material` | str | None | 材料名称（从porosity_calculator查孔隙率） |
| `--output` | str | 自动 | 输出STL文件名（默认TOE_Wire_15cm_P70_df2.973.stl） |
| `--seed` | int | 42 | 随机种子（保证可复现性） |

### 3.3 主生成函数参数（generate_TOE_wire_gyroid）

如需在Python脚本中直接调用，以下是完整参数列表：

```python
result = generate_TOE_wire_gyroid(
    # ──── 几何 ────
    length_cm=15.0,          # 线体长度 cm
    width_cm=1.2,            # 线体宽度 cm
    height_cm=1.2,           # 线体高度 cm
    # ──── 物理 ────
    target_porosity=0.697,   # 目标孔隙率 (TOE定值0.697)
    target_df=2.973,          # 目标分形维度 (TOE定值2.973)
    # ──── Gyroid骨架 ────
    gyroid_wavelength_cm=None,# Gyroid波长 (None=自动=min_dim/2.5)
    gyroid_porosity=0.50,    # Gyroid基础孔隙率 (0.50=标准)
    # ──── 空隙 ────
    void_mode="powerlaw",    # "powerlaw" 或 "fractal"
    subdivision=3,            # 分形层级细分数 (仅fractal)
    min_void_mm=None,         # 最小空隙直径mm (None=自动)
    max_void_mm=None,         # 最大空隙直径mm (None=自动)
    # ──── 打印 ────
    min_feature_mm=0.10,     # 最小特征mm (SLA=0.10, FDM=0.40)
    voxel_size_mm=0.05,      # 体素尺寸mm
    shell_thickness_mm=0.5,  # 外壳厚度mm
    # ──── 验证 ────
    verify_df=True,           # 盒计数验证d_f
    verify_connectivity=True, # 连通性验证
    # ──── 其他 ────
    seed=42,                  # 随机种子
    verbose=True              # 详细输出
)
```

**返回值**：

| 字段 | 类型 | 含义 |
|------|------|------|
| `mesh` | trimesh.Trimesh | 生成的3D网格（可直接export为STL） |
| `porosity` | float | 实际孔隙率 |
| `size_cm` | [float×3] | 输出尺寸 [长, 宽, 高] |
| `fractal_dimension` | float | 盒计数实测d_f（或设计值） |
| `n_voids` | int | 空隙总数 |
| `voids` | list | 空隙列表 [(cx,cy,cz,r), ...] |
| `n_components_solid` | int | Solid连通分量数 |
| `connectivity_ok` | bool | 连通性是否通过 |
| `ok` | bool | **全部验证是否通过** |

---

## 4. 材料数据库速查表

### 4.1 TOE超导候选材料（核心）

| 材料键名 | 全名 | ρ (g/cm³) | ρ_th (g/cm³) | Θ_D (K) | λ₀ | μ\*₀ | C_λ | C_μ | φ_0 | Tc预估 |
|----------|------|-----------|-------------|---------|------|-------|------|------|------|--------|
| `graphene_superlattice` | 分形空位修饰石墨烯超晶格 | 2.20 | 2.267 | 1200 | 1.50 | 0.10 | 5.49 | 3.75 | 2.95% | 200-250K |
| `hydrogenated_copper_oxide` | 氢化铜氧化物 | 5.50 | 6.31 | 650 | 1.20 | 0.14 | 2.90 | 0.42 | 12.84% | ~83K |
| `lithium_graphane` | 锂插层石墨烷(Li-Gr) | 1.80 | 2.20 | 950 | 2.00 | 0.08 | 4.40 | 2.5 | 18.18% | ~350K |
| `boride_Y2B3C2` | 硼化物三角格子(Y₂B₃C₂) | 4.20 | 5.20 | 800 | 1.60 | 0.09 | 3.66 | 1.8 | 19.23% | ~150K |

### 4.2 已知超导参考材料

| 材料键名 | 全名 | 已知Tc | 用途 |
|----------|------|--------|------|
| `YBCO` | YBa₂Cu₃O₇ | 93K | 参考校准 |
| `MgB2` | MgB₂ | 39K | 双带超导参考 |
| `Nb3Sn` | Nb₃Sn | 18K | A15结构参考 |

### 4.3 3D打印结构验证材料

| 材料键名 | 全名 | ρ (g/cm³) | 用途 |
|----------|------|-----------|------|
| `PLA` | 聚乳酸 | 1.24 | FDM打印验证结构 |
| `PETG` | PETG | 1.27 | FDM打印验证结构 |
| `resin` | 光敏树脂 | 1.15 | SLA打印验证结构 |
| `nylon_PA12` | 尼龙PA12 | 1.01 | SLS打印验证结构 |
| `stainless_316L` | 316L不锈钢 | 7.99 | SLM金属打印 |
| `titanium_TI64` | Ti-6Al-4V | 4.43 | EBM/SLM金属打印 |
| `copper_printed` | 3D打印纯铜 | 8.96 | SLM高电导率 |
| `aluminum_printed` | 铝合金AlSi10Mg | 2.67 | SLM轻量化 |

---

## 5. 典型使用场景SOP

### 场景A：TOE超导材料全流程（石墨烯超晶格）

```bash
# Step 1: 计算孔隙率
python porosity_calculator.py \
  --material graphene_superlattice \
  --target-tc 250 \
  --phi-air 0.015

# 输出关键结果:
#   target_porosity = 0.6970
#   target_df = 2.973
#   Tc_max ≈ 200-250K

# Step 2: 生成3D模型
python superconductor_gyroid.py \
  --porosity 0.6970 \
  --df 2.973 \
  --material graphene_superlattice \
  --void-mode powerlaw \
  --voxel-mm 0.05 \
  --min-feature-mm 0.10
```

### 场景B：氢化铜氧化物真空工况（预设方案）

```bash
# Step 1: 使用预设（自动填充全部参数）
python porosity_calculator.py \
  --preset hydrogenated_copper_oxide_vacuum

# 输出关键结果:
#   φ_0 = 12.84%, φ_挖实 = 56.86%, φ_挖表 = 57.96%
#   λ_eff = 2.312, μ*_eff = 0.1082
#   Tc_max = 83.4K

# Step 2: 生成3D模型（真空工况 φ_air=0）
python superconductor_gyroid.py \
  --porosity 0.6970 \
  --df 2.973 \
  --phi-air 0 \
  --voxel-mm 0.05
```

### 场景C：FDM快速结构验证（用PLA测试打印可行性）

```bash
# Step 1: 查看PLA参数
python porosity_calculator.py --material PLA

# Step 2: 生成粗分辨率FDM模型
python superconductor_gyroid.py \
  --porosity 0.6970 \
  --df 2.973 \
  --min-feature-mm 0.40 \
  --voxel-mm 0.20 \
  --no-verify \
  --output PLA_test.stl
```

### 场景D：自定义新材料快速评估

```bash
# Step 1: 从密度估算（粗略）
python porosity_calculator.py \
  --density 3.5 \
  --theta_D 900

# Step 2: 精确参数覆盖
python porosity_calculator.py \
  --density 3.5 \
  --theta_D 900 \
  --lambda_0 1.80 \
  --mu_star_0 0.09 \
  --C_lambda 4.0 \
  --C_mu 2.0

# Step 3: 生成模型
python superconductor_gyroid.py \
  --porosity 0.6970 \
  --df 2.973
```

### 场景E：分形层级空隙模式

```bash
# 使用Menger-Sponge变体（空间更均匀的钉扎分布）
python superconductor_gyroid.py \
  --porosity 0.6970 \
  --df 2.973 \
  --void-mode fractal \
  --subdivision 3
```

---

## 6. 常见问题与排错

### Q1: 孔隙率达不到69.7%？

**可能原因**：
1. Gyroid骨架占位过多 → 降低 `--gyroid-porosity`
2. 空隙雕刻效率低 → 程序已内置多轮迭代，检查输出中的效率值
3. 浮岛清理移除了过多固体 → `--shell-mm` 设大一点

### Q2: d_f盒计数验证不通过？

**原因**：盒计数法测量的d_f与幂律构造的d_f可能有系统差异
- 幂律构造 d_f=2.973 由 N(>r)∝r^{-2.973} **数学保证**
- 盒计数在有限网格上测量有统计涨落
- **设计值2.973是可信的**，盒计数偏差<0.05属正常

### Q3: 连通性验证失败？

**可能原因**：
1. Gyroid壁太薄被空隙切断 → 降低 `--gyroid-porosity` 或增大 `--wavelength`
2. 最大空隙太大 → 降低 `--max-void-mm`
3. 外壳太薄 → 增大 `--shell-mm`

### Q4: 内存不足？

**解决方案**：
1. 增大 `--voxel-mm`（0.05→0.10）
2. 减小 `--length`
3. 减少 `--width` 和 `--height`

### Q5: Tc计算结果为0？

**可能原因**：
1. 材料本征耦合太弱（如PLA等聚合物λ₀<0.1）
2. 库仑赝势太高 → 检查 `--mu_star_0` 是否合理
3. 分形增强不够 → 检查 `--C_lambda` 是否太小

### Q6: φ_挖表 和 φ_挖实 的区别？

```
φ_挖实 = 有效物理空隙（真正贡献超导的净空）
φ_挖表 = 工程设计开孔尺寸（实际需要在3D模型上开的孔）

φ_挖表 = φ_挖实 + φ_空气 + φ_0泡 + Δφ_分形
         ↑           ↑──────────────────────────↑
    超导需要的      三重修正（会挤占有效空间）

⚠️ 如果按 φ_挖表 = 69.7% - φ_0 开孔（忽略修正）
   → φ_挖实 = φ_挖表 - 三重修正 < 69.7% - φ_0
   → φ_总 = φ_0 + φ_挖实 < 69.7%
   → 失超!
```

### Q7: powerlaw 和 fractal 模式怎么选？

| 需求 | 推荐 |
|------|------|
| 快速验证/调试 | powerlaw |
| 最终超导样品 | fractal（更均匀的钉扎分布） |
| 大尺寸线体 | powerlaw（fractal递归可能内存不足） |
| 物理各向同性 | fractal |

---

## 附录：参数与TOE理论章节对应

| 参数 | 理论章节 | 核心公式 |
|------|----------|----------|
| φ_总=69.7% | §1.3, §22 | φ_总=φ_本征+φ_挖孔=(C_d-1)/C_d=Ω_Λ |
| d_f=2.973 | §1.10.4 | d_f,sup=2+d_max/100, D=13→d_f=2.973 |
| η_Λ=0.027 | §2.2 | η_Λ=3-d_f |
| λ_eff | §22.5 | λ_eff=λ₀+C_λ×P×√(1-P) |
| μ\*_eff | §22.5 | μ\*_eff=μ\*₀/(1+C_μ×P) |
| φ_空气 | §22 | 常压1.5%, 真空0% |
| φ_0泡 | §22 | 0.3% (0维世界泡本底) |
| Δφ_分形 | §22 | 0.8% (孔壁粗糙度) |
| V_固%=30.3% | §2.2 | Ω_m=1/C_d=30.3% |
| McMillan公式 | §22.5 | T_c=(Θ_D/1.45)×exp[-1.04(1+λ)/(λ-μ\*(1+0.62λ))] |
