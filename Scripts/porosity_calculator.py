"""
TOE超导体孔隙率计算器 v1.0
=============================

基于TOE_v8.8框架，根据材料本征密度计算所需的挖孔孔隙率
输出可直接填入 superconductor_cube.py / superconductor_TOE_exact.py / superconductor_gyroid.py

核心公式链（TOE框架）:
  1. η_Λ = 3 - d_f = 0.027  (分形维度亏缺，非拟合参数，来自宇宙大尺度结构观测)
  2. d_f,sup = 2.973         (最优分形维度，由D=13大统一闭包导出)
  3. d_f,sup = 2 + d_max/100 (d_max=92→2.92下限, d_max=98→2.98上限, d_max=95→2.95中心)
  4. φ_总 = φ_本征 + φ_挖孔 = 69.7% = Ω_Λ = (C_d-1)/C_d  (超导总孔隙率=宇宙孔隙率)
  5. φ_本征 = 1 - f_m = 1 - ρ_bulk/ρ_theoretical  (物质体积占比→本征孔隙率, §1.3公式)
  6. φ_挖孔 = 69.7% - φ_本征  (人工补足至宇宙孔隙率)
  7. λ_eff = λ_0 + C_λ × P × √(1-P)        (分形空位耦合增强)
  8. μ*_eff = μ*_0 / (1 + C_μ × P)          (分形结构库仑屏蔽)
  9. T_c = McMillan(Θ_D, λ_eff, μ*_eff)      (超导临界温度)

物理依据：
  - D=13闭包 → 耦合稀释 → d_f∈(2.92,2.98) → 分形空位结构复制真空耦合效应
  - 超导总孔隙率必须匹配宇宙孔隙率69.7%: φ_总=φ_本征+φ_挖孔=69.7%
  - 固体基体V_固%=30.3%=Ω_m，是超导导电主体
  - 五大约束锁定69.7%: 冷媒流通、磁通钉扎、库柏对限域、热膨胀缓冲、力学形变
  - 分形空位修饰增强电声耦合：石墨烯本征λ=1.5 → 修饰后λ=3.2
  - 分形结构屏蔽库仑相互作用：常规μ*=0.10-0.15 → 屏蔽后μ*=0.04

使用方法:
  python porosity_calculator.py --material graphene_superlattice
  python porosity_calculator.py --density 2.2 --theta_D 1200
  python porosity_calculator.py --list-materials
  python porosity_calculator.py --material graphene_superlattice --target-tc 250
  python porosity_calculator.py --material graphene_superlattice --json
"""

import numpy as np
import sys
import json

# =====================================================================
# 材料数据库
# =====================================================================

MATERIALS = {
    # ──────────── 3D打印常用材料（结构验证用） ────────────
    "PLA": {
        "name": "PLA (聚乳酸)",
        "category": "3D打印塑料",
        "density": 1.24,           # g/cm³ (实际密度, 含烧结/挤出微孔)
        "density_theoretical": 1.25,  # g/cm³ (X射线理论密度, 完全致密)
        "Theta_D": 180,            # K (聚合物德拜温度较低)
        "lambda_0": 0.05,          # 聚合物几乎无电声耦合
        "mu_star_0": 0.20,         # 屏蔽差
        "C_lambda": 0.5,           # 低耦合增强
        "C_mu": 0.2,               # 低屏蔽
        "note": "FDM打印验证结构用，本身不超导"
    },
    "PETG": {
        "name": "PETG",
        "category": "3D打印塑料",
        "density": 1.27,
        "density_theoretical": 1.27,
        "Theta_D": 170,
        "lambda_0": 0.04,
        "mu_star_0": 0.22,
        "C_lambda": 0.4,
        "C_mu": 0.18,
        "note": "FDM打印验证结构用"
    },
    "resin": {
        "name": "光敏树脂",
        "category": "3D打印树脂",
        "density": 1.15,
        "density_theoretical": 1.20,
        "Theta_D": 150,
        "lambda_0": 0.03,
        "mu_star_0": 0.25,
        "C_lambda": 0.3,
        "C_mu": 0.15,
        "note": "SLA打印验证结构用"
    },
    "nylon_PA12": {
        "name": "尼龙PA12",
        "category": "3D打印塑料",
        "density": 1.01,
        "density_theoretical": 1.03,
        "Theta_D": 200,
        "lambda_0": 0.04,
        "mu_star_0": 0.22,
        "C_lambda": 0.4,
        "C_mu": 0.15,
        "note": "SLS打印验证结构用"
    },

    # ──────────── 金属3D打印材料 ────────────
    "stainless_316L": {
        "name": "316L不锈钢",
        "category": "金属3D打印",
        "density": 7.99,              # SLM打印实际密度
        "density_theoretical": 8.03,  # X射线理论密度(Fe-Cr-Ni-Mo固溶体)
        "Theta_D": 290,
        "lambda_0": 0.30,
        "mu_star_0": 0.13,
        "C_lambda": 1.8,
        "C_mu": 0.5,
        "note": "SLM打印，结构验证"
    },
    "titanium_TI64": {
        "name": "Ti-6Al-4V钛合金",
        "category": "金属3D打印",
        "density": 4.43,
        "density_theoretical": 4.43,   # α+β双相, 理论≈实测
        "Theta_D": 380,
        "lambda_0": 0.54,
        "mu_star_0": 0.12,
        "C_lambda": 2.4,
        "C_mu": 0.7,
        "note": "EBM/SLM打印"
    },
    "copper_printed": {
        "name": "3D打印纯铜",
        "category": "金属3D打印",
        "density": 8.96,              # 纯铜实际密度(致密)
        "density_theoretical": 8.96,  # Cu FCC晶格理论密度
        "Theta_D": 343,
        "lambda_0": 0.14,
        "mu_star_0": 0.10,
        "C_lambda": 2.1,
        "C_mu": 0.4,
        "note": "SLM打印纯铜，高电导率"
    },
    "aluminum_printed": {
        "name": "3D打印铝合金(AlSi10Mg)",
        "category": "金属3D打印",
        "density": 2.67,
        "density_theoretical": 2.67,  # Al-Si合金理论密度
        "Theta_D": 428,
        "lambda_0": 0.43,
        "mu_star_0": 0.10,
        "C_lambda": 2.7,
        "C_mu": 0.8,
        "note": "SLM打印铝合金"
    },

    # ──────────── TOE论文超导候选材料 ────────────
    "graphene_superlattice": {
        "name": "分形空位修饰石墨烯超晶格",
        "category": "TOE超导候选",
        "density": 2.20,              # 石墨烯超晶格实际密度(含层间空隙)
        "density_theoretical": 2.267, # 石墨晶体X射线理论密度
        "Theta_D": 1200,
        "lambda_0": 1.50,           # 石墨烯本征耦合
        "mu_star_0": 0.10,         # 石墨烯本征库仑赝势
        "C_lambda": 5.49,          # 校准: Δλ=1.7 at P=0.40, 公式C_λ×P×√(1-P)
        "C_mu": 3.75,              # 校准: μ*=0.04 at P=0.40
        "note": "TOE论文主推材料，Tc可达200-250K，Θ_D=1200K(碳极轻)"
    },
    "hydrogenated_copper_oxide": {
        "name": "氢化铜氧化物",
        "category": "TOE超导候选",
        "density": 5.50,              # CuO烧结体实际密度
        "density_theoretical": 6.31,  # CuO单斜晶系X射线理论密度
        "Theta_D": 650,
        "lambda_0": 1.20,          # CuO本征耦合
        "mu_star_0": 0.14,
        "C_lambda": 2.90,          # P=0.697: Δλ=1.113, λ_eff=2.312
        "C_mu": 0.42,              # P=0.697: μ*_eff=0.1082
        "note": "TOE验证材料，Tc_max=83.4K(P=69.7%)，真空工况φ_air=0，ΔT_c=116.6K(缺口)"
    },
    "lithium_graphane": {
        "name": "锂插层石墨烷(Li-Gr)",
        "category": "TOE超导候选",
        "density": 1.80,              # Li插层石墨烷实际密度
        "density_theoretical": 2.20,  # 理论密度(石墨+Li插层)
        "Theta_D": 950,
        "lambda_0": 2.00,
        "mu_star_0": 0.08,
        "C_lambda": 4.40,
        "C_mu": 2.5,
        "note": "TOE论文预言Tc≈350K(室温)，CVD+插层工艺"
    },
    "boride_Y2B3C2": {
        "name": "硼化物三角格子(Y₂B₃C₂)",
        "category": "TOE超导候选",
        "density": 4.20,              # 烧结体实际密度
        "density_theoretical": 5.20,  # X射线理论密度(Y-B-C三角格子)
        "Theta_D": 800,
        "lambda_0": 1.60,
        "mu_star_0": 0.09,
        "C_lambda": 3.66,
        "C_mu": 1.8,
        "note": "TOE论文备选材料，多晶块体稳定，a≈0.495nm"
    },

    # ──────────── 已知超导参考材料 ────────────
    "YBCO": {
        "name": "YBa₂Cu₃O₇ (YBCO)",
        "category": "已知高温超导",
        "density": 6.30,              # 烧结体实际密度
        "density_theoretical": 6.38,  # 正交晶系X射线理论密度
        "Theta_D": 400,
        "lambda_0": 1.50,
        "mu_star_0": 0.13,
        "C_lambda": 2.5,
        "C_mu": 0.5,
        "note": "已知Tc=93K，用作参考校准"
    },
    "MgB2": {
        "name": "MgB₂",
        "category": "已知超导",
        "density": 2.63,
        "density_theoretical": 2.63,  # AlB2型结构, 理论≈实测
        "Theta_D": 750,
        "lambda_0": 1.00,
        "mu_star_0": 0.10,
        "C_lambda": 4.7,
        "C_mu": 1.0,
        "note": "已知Tc=39K，双带超导"
    },
    "Nb3Sn": {
        "name": "Nb₃Sn",
        "category": "已知超导",
        "density": 8.87,
        "density_theoretical": 8.87,  # A15结构, 理论≈实测
        "Theta_D": 280,
        "lambda_0": 1.80,
        "mu_star_0": 0.12,
        "C_lambda": 1.75,
        "C_mu": 0.3,
        "note": "已知Tc=18K，A15结构"
    },
}

# =====================================================================
# 已验证预设 (TOE框架全套参数校验通过)
# =====================================================================

VERIFIED_PRESETS = {
    "hydrogenated_copper_oxide_vacuum": {
        "name": "氢化铜氧化物（真空工况）— TOE校验通过 2026-06-04",
        "material": "hydrogenated_copper_oxide",
        # ──── 基础密度&本征孔隙 (§1.3) ────
        "rho_bulk": 5.50,              # g/cm³ (CuO烧结体实际密度)
        "rho_theoretical": 6.31,       # g/cm³ (CuO单斜晶系X射线理论密度)
        "f_m": 0.8716,                 # = ρ_bulk/ρ_theo = 5.50/6.310
        "phi0": 0.1284,                # = 1-f_m = 12.84% (本征闭孔, 定值, 不受修正影响)
        # ──── 挖孔三重修正 (§22) ────
        "phi_air": 0.000,              # 空气填充=0% (真空工况; 常压1.5%)
        "phi_0bubble": 0.003,          # 0维世界泡=0.30% (真空本底固有, 无法消除)
        "delta_phi_fractal": 0.008,    # 分形凹凸挤占=0.80%
        "dphi_sum": 0.011,             # 三重修正总和=1.10%
        # ──── 有效空隙 vs 表观开孔 ────
        "phi_drill_effective": 0.5686, # φ_挖实 = 69.7%-φ_0 = 56.86%
        "phi_drill_apparent": 0.5796,  # φ_挖表 = φ_挖实+Δφ_Σ = 57.96% (唯一适配69.7%净空的设计开孔值)
        # ──── TOE分形维度 ────
        "target_porosity": 0.6970,     # φ_总 = Ω_Λ = (C_d-1)/C_d
        "target_df": 2.973,            # d_f,sup (2.92<d_f<2.98, 合规)
        "eta_Lambda": 0.0270,          # = 3-d_f
        # ──── 分形修饰物性 (P=0.697) ────
        "Theta_D_eff": 650,            # K (P<渗流阈值, 无修正)
        "lambda_eff": 2.312,           # = λ_0+C_λ×P×√(1-P) = 1.20+1.113
        "mu_star_eff": 0.1082,         # = μ*_0/(1+C_μ×P) = 0.14/1.293
        # ──── McMillan超导临界温度 ────
        "Tc_max": 83.4,               # K (McMillan, P=69.7%)
        "Tc_target": 200.0,            # K (目标)
        "delta_Tc": 116.6,            # K (缺口)
        # ──── 固体基体 ────
        "V_solid": 0.303,              # Ω_m = 30.3%
        # ──── 容错结论 ────
        "fault_tolerance_note": "忽略三重修正→φ_挖实=55.76%→φ_总=68.60%<69.7%→失超!",
        # ──── 提温优化方向 ────
        "optimization": [
            "替换更轻质元素抬升Θ_D",
            "增加晶格空位/人工修饰提升C_λ",
            "优化分形屏蔽结构压低μ*_eff",
        ],
    },
}


# =====================================================================
# 核心公式
# =====================================================================

def calculate_intrinsic_porosity(mat_params, phi_air=0.015, phi_0bubble=0.003, delta_phi_fractal=0.008):
    """
    从本征密度计算材料本征孔隙率 + 挖孔三重修正 (TOE框架§1.3 + §22修正)
    
    核心原理: 69.7%是净空有效物理空隙，而非几何账面空隙
    挖孔表观空隙 ≠ 有效物理空隙，存在三重挤占修正:
      1. 空气填充修正 φ_空气: 常压下空气占据挖孔体积
      2. 0维世界泡修正 φ_0泡: 极限真空仍存在的本源真空元填充
      3. 分形凹凸挤占 Δφ_分形: 孔壁分形粗糙面挤占有效空腔
    
    公式链:
      φ_0 (本征) = 1 - f_m = 1 - ρ_bulk/ρ_theoretical  (定值, 不受修正影响)
      φ_挖实 = φ_挖表 - φ_空气 - φ_0泡(D_f) - Δφ_分形  (有效物理空隙)
      φ_总 = φ_0 + φ_挖实 = 69.7% = Ω_Λ               (有效物理空隙定值)
      φ_挖表 = 69.7% - φ_0 + Δφ_分形 + φ_空气 + φ_0泡    (几何设计尺寸)
    
    参数:
        mat_params: 材料参数字典，需含 density 和 density_theoretical
        phi_air: 空气填充率 (默认1.5%, 常压环境; 真空工况=0)
        phi_0bubble: 0维世界泡填充率 (默认0.3%, 真空本底固有, 无法消除)
        delta_phi_fractal: 分形凹凸挤占率 (默认0.8%, 孔壁粗糙度)
    
    返回:
        dict: 包含 f_m, φ_0, φ_挖表, φ_挖实, 三重修正量, V_固%
    """
    rho_bulk = mat_params["density"]
    rho_theoretical = mat_params.get("density_theoretical", rho_bulk)
    
    # 物质体积占比 = 相对密度 (TOE §1.3: f_m)
    f_m = rho_bulk / rho_theoretical
    
    # 本征孔隙率 φ_0 = 1 - f_m (定值, 不受修正影响)
    # 闭孔隔绝外界气体与真空世界泡侵入, 数值严格不变
    phi_0 = 1.0 - f_m
    
    # 超导总孔隙率目标 (宇宙孔隙率, 有效物理空隙定值)
    PHI_TOTAL = 0.697  # = Ω_Λ = (C_d-1)/C_d
    
    # 三重修正量
    # φ_空气: 常压下空气占据挖孔体积 (真空工况=0)
    # φ_0泡: 0维世界泡在挖孔空间的等效体积占比 (真空本底固有, 无法消除)
    # Δφ_分形: 孔壁分形凹凸挤占的有效空腔体积
    total_correction = phi_air + phi_0bubble + delta_phi_fractal
    
    # 有效挖孔空隙 (扣除三重修正后的真实物理空隙)
    phi_drilling_effective = PHI_TOTAL - phi_0  # = φ_挖实
    
    # 几何设计挖孔空隙 (工程实际需要开孔的账面尺寸)
    # φ_挖表 = 69.7% - φ_0 + Δφ_分形 + φ_空气 + φ_0泡
    phi_drilling_geometric = phi_drilling_effective + total_correction
    
    # 固体基体体积占比
    V_solid = 1.0 - PHI_TOTAL  # = 30.3% = Ω_m
    
    return {
        "f_m": round(f_m, 6),                          # 物质体积占比 (§1.3)
        "phi_0": round(phi_0, 6),                       # φ_0 本征空隙 (定值)
        "phi_drilling_effective": round(max(0, phi_drilling_effective), 6),  # φ_挖实
        "phi_drilling_geometric": round(max(0, phi_drilling_geometric), 6),  # φ_挖表 (设计尺寸)
        "phi_air": phi_air,                             # 空气填充修正
        "phi_0bubble": phi_0bubble,                      # 0维世界泡修正
        "delta_phi_fractal": delta_phi_fractal,          # 分形凹凸挤占修正
        "total_correction": round(total_correction, 6),  # 三重修正总量
        "V_solid": round(V_solid, 4),                    # V_固% = 30.3% = Ω_m
        "rho_bulk": rho_bulk,
        "rho_theoretical": rho_theoretical,
    }

def mcmillan_tc(Theta_D, lam, mu_star):
    """
    McMillan公式计算超导临界温度

    T_c = (Θ_D/1.45) × exp[-1.04(1+λ) / (λ - μ*(1+0.62λ))]

    参数:
        Theta_D: 德拜温度 (K)
        lam: 电声耦合常数
        mu_star: 库仑赝势

    返回:
        T_c (K), 不满足超导条件返回0
    """
    denominator = lam - mu_star * (1 + 0.62 * lam)
    if denominator <= 0:
        return 0.0
    exponent = -1.04 * (1 + lam) / denominator
    Tc = (Theta_D / 1.45) * np.exp(exponent)
    return max(0.0, Tc)


def calculate_fractal_parameters(mat_params, porosity, d_f=2.973):
    """
    计算分形空位修饰后的有效参数

    TOE框架核心公式:
      η_Λ = 3 - d_f = 0.027    (分形维度亏缺，非拟合参数)

    修正公式（物理自洽，高P不发散）:
      Θ_D_eff = Θ_D_0 × (1-P)^{1/3}     (密度标度: 多孔材料声速下降→Θ_D下降)
      λ_eff = λ_0 + C_λ × P × (1-P)^{1/2} (分形空位耦合增强, P=0和P→1时均为0)
      μ*_eff = μ*_0 / (1 + C_μ × P)      (分形结构库仑屏蔽)

    物理约束:
      1. P→0: 无空位→无增强, λ=λ_0, μ*=μ*_0
      2. P→1: 无材料→无耦合, Δλ→0, 但Θ_D→0 (物理截止)
      3. 最大Δλ在 P = 2/3 处 (理论极值, 实际最优P通常在0.30-0.55)

    校准数据 (TOE论文§22.5):
      石墨烯超晶格 P=0.40: λ=1.5→3.2, μ*=0.10→0.04
      氢化铜氧化物 P=0.40: λ=1.2→2.1, μ*=0.14→0.12
    """
    P = porosity
    eta_Lambda = 3 - d_f  # = 0.027 (非拟合参数，来自宇宙大尺度结构观测)

    if P >= 1.0 or P < 0:
        return None

    # 1. 德拜温度 — 原子级不变，但超过渗流阈值后有效值急剧下降
    #    物理原因：Θ_D是晶格振动频率，空位不改变剩余原子的振动频率
    #    但当孔隙率超过渗流阈值(~0.70)，材料失去连续性，有效Θ_D→0
    #    TOE论文§22.5直接使用本征Θ_D（不做密度修正），此处相同
    P_percolation = 0.75  # 渗流阈值(高孔隙多孔超导专属区间, TOE框架φ_总=69.7%)
    if P < P_percolation:
        Theta_D_eff = mat_params["Theta_D"]
    else:
        # 超过渗流阈值，Θ_D急剧下降（材料不连续）
        Theta_D_eff = mat_params["Theta_D"] * max(0, (1 - P) / (1 - P_percolation))

    # 2. 分形空位耦合增强 (修正: P→1时不发散)
    #    物理原因：分形空位修饰改变电子态密度和声子谱，增强电声耦合
    #    但高孔隙率时材料减少，可用电子也减少，增强效应自然饱和并下降
    #    公式: Δλ = C_λ × P × (1-P)^{1/2}
    #    - P=0: Δλ=0 (无空位)
    #    - P=2/3: 最大值
    #    - P→1: Δλ→0 (无材料)
    delta_lambda = mat_params["C_lambda"] * P * (1 - P) ** 0.5
    lambda_eff = mat_params["lambda_0"] + delta_lambda

    # 3. 分形结构库仑屏蔽
    #    物理原因：分形空位结构缩短屏蔽长度，降低有效库仑赝势
    #    公式: μ*_eff = μ*_0 / (1 + C_μ × P)
    #    - P=0: μ*=μ*_0 (无屏蔽)
    #    - P→1: μ*→0 (完全屏蔽，但此时也无材料)
    mu_star_eff = mat_params["mu_star_0"] / (1 + mat_params["C_mu"] * P)

    # 4. 有效密度
    rho_eff = mat_params["density"] * (1 - P)

    return {
        "Theta_D_eff": Theta_D_eff,
        "lambda_eff": lambda_eff,
        "mu_star_eff": mu_star_eff,
        "rho_eff": rho_eff,
        "eta_Lambda": eta_Lambda,
        "delta_lambda": delta_lambda,
    }


# =====================================================================
# 最优孔隙率搜索
# =====================================================================

def find_optimal_porosity(mat_params, d_f=2.973, target_tc=200.0, verbose=True,
                          phi_air=0.015, phi_0bubble=0.003, delta_phi_fractal=0.008):
    """
    寻找最优孔隙率

    遍历孔隙率P，找到使T_c最大的P。
    如果T_c能达到target_tc，则返回该P；否则返回T_c最大时的P并给出建议。

    参数:
        mat_params: 材料参数字典
        d_f: 目标分形维度 (默认2.973，TOE框架固定值)
        target_tc: 目标临界温度 (K)
        verbose: 是否打印详细信息
        phi_air: 空气填充修正率 (默认1.5%, 常压环境; 真空工况=0)
        phi_0bubble: 0维世界泡修正率 (默认0.3%, 真空本底固有)
        delta_phi_fractal: 分形凹凸挤占修正率 (默认0.8%, 孔壁粗糙度)

    返回:
        dict: 包含最优孔隙率和相关参数
    """
    best_P = 0
    best_Tc = 0
    results = []

    # 粗搜索 (P上限0.75: 覆盖TOE框架超导孔隙率69.7%)
    P_max_search = 0.75
    for P in np.arange(0.02, P_max_search, 0.005):
        params = calculate_fractal_parameters(mat_params, P, d_f)
        if params is None:
            continue
        Tc = mcmillan_tc(params["Theta_D_eff"], params["lambda_eff"], params["mu_star_eff"])
        results.append({"P": P, "Tc": Tc, **params})
        if Tc > best_Tc:
            best_Tc = Tc
            best_P = P

    # 精搜索
    if best_P > 0:
        for P in np.arange(max(0.02, best_P - 0.02), min(P_max_search - 0.01, best_P + 0.02), 0.001):
            params = calculate_fractal_parameters(mat_params, P, d_f)
            if params is None:
                continue
            Tc = mcmillan_tc(params["Theta_D_eff"], params["lambda_eff"], params["mu_star_eff"])
            if Tc > best_Tc:
                best_Tc = Tc
                best_P = P

    # 最优点的参数
    opt_params = calculate_fractal_parameters(mat_params, best_P, d_f)

    # 约束检查
    d_f_range_ok = 2.92 <= d_f <= 2.98
    Tc_achievable = best_Tc >= target_tc

    # 非超导材料安全检查: Tc=0时不应返回搜索上限P
    if best_Tc < 1.0:
        best_P = 0.40  # 回退到默认合理值
        best_Tc = 0.0
        opt_params = calculate_fractal_parameters(mat_params, best_P, d_f)
        Tc_achievable = False

    # 孔隙率上限: TOE框架超导孔隙率为69.7%(宇宙孔隙率Ω_Λ)
    # φ_总=φ_本征+φ_挖孔=69.7%, 由冷媒流通/磁通钉扎/库柏对限域/热膨胀/力学五大约束锁定
    best_P = min(best_P, 0.697)

    result = {
        "material": mat_params["name"],
        "category": mat_params.get("category", "自定义"),
        "density": mat_params["density"],
        "target_porosity": round(best_P, 4),
        "target_df": d_f,
        "Tc_max": round(best_Tc, 1),
        "Tc_target": target_tc,
        "Tc_achievable": Tc_achievable,
        "d_f_range_ok": d_f_range_ok,
    }
    result.update(opt_params)

    # TOE框架: φ_0 + φ_挖实 = 69.7% = Ω_Λ (有效物理空隙定值)
    # 由密度直接换算: φ_0 = 1 - f_m = 1 - ρ_bulk/ρ_theoretical (§1.3)
    # 三重修正: φ_挖实 = φ_挖表 - φ_空气 - φ_0泡 - Δφ_分形
    intrinsic_info = calculate_intrinsic_porosity(mat_params, phi_air, phi_0bubble, delta_phi_fractal)
    result["framework_total_porosity"] = 0.697
    result["f_m"] = intrinsic_info["f_m"]  # 物质体积占比
    result["phi_0"] = intrinsic_info["phi_0"]  # φ_0 本征空隙 (定值)
    result["phi_drilling_effective"] = intrinsic_info["phi_drilling_effective"]  # φ_挖实
    result["phi_drilling_geometric"] = intrinsic_info["phi_drilling_geometric"]  # φ_挖表 (设计尺寸)
    result["phi_air"] = intrinsic_info["phi_air"]  # 空气填充修正
    result["phi_0bubble"] = intrinsic_info["phi_0bubble"]  # 0维世界泡修正
    result["delta_phi_fractal"] = intrinsic_info["delta_phi_fractal"]  # 分形凹凸挤占修正
    result["total_correction"] = intrinsic_info["total_correction"]  # 三重修正总量
    result["rho_theoretical"] = intrinsic_info["rho_theoretical"]  # ρ_theoretical
    result["solid_volume_fraction"] = intrinsic_info["V_solid"]  # V_固% = 30.3% = Ω_m

    if verbose:
        print_report(result, results)

    return result


# =====================================================================
# 报告输出
# =====================================================================

def print_report(result, all_results=None):
    """打印详细计算报告"""
    print("\n" + "=" * 72)
    print("          【TOE超导体孔隙率计算报告】")
    print("=" * 72)

    # 材料信息
    print(f"\n📋 材料信息:")
    print(f"   名称:     {result['material']}")
    print(f"   类别:     {result['category']}")
    print(f"   本征密度: ρ_m = {result['density']:.2f} g/cm³")
    print(f"   有效密度: ρ_eff = {result['rho_eff']:.2f} g/cm³")

    # TOE框架参数
    print(f"\n📐 TOE框架参数 (d_f,sup = 2 + d_max/100):")
    print(f"   分形维度:       d_f = {result['target_df']}")
    print(f"   d_f ∈ (2.92, 2.98):  {'✅ 满足' if result['d_f_range_ok'] else '❌ 不满足'}")
    print(f"   维度亏缺:       η_Λ = 3 - d_f = {result['eta_Lambda']:.4f}")

    # 分形修饰后参数
    print(f"\n🔬 分形修饰后有效参数 (η_Λ={result['eta_Lambda']:.4f}):")
    print(f"   Θ_D_eff = {result['Theta_D_eff']:.0f} K     (德拜温度，原子振动不变)")
    print(f"   λ_eff   = {result['lambda_eff']:.3f}       (本征→+{result['delta_lambda']:.3f} 分形增强)")
    print(f"   μ*_eff  = {result['mu_star_eff']:.4f}      (分形结构库仑屏蔽)")

    # 超导临界温度
    print(f"\n🎯 超导临界温度 (McMillan公式):")
    print(f"   T_c(max)  = {result['Tc_max']:.1f} K")
    print(f"   目标 T_c  = {result['Tc_target']:.0f} K")
    if result['Tc_achievable']:
        print(f"   ✅ 可达目标! (余量 {result['Tc_max'] - result['Tc_target']:.1f} K)")
    else:
        gap = result['Tc_target'] - result['Tc_max']
        print(f"   ⚠️  未达目标 (差 {gap:.1f} K)")
        print(f"   💡 建议: 提高Θ_D(更轻元素)、增大C_λ(更多空位修饰)、降低μ*(更好屏蔽)")

    # 推荐参数
    print(f"\n" + "─" * 72)
    print(f"💡 推荐参数 — 直接填入生成程序:")
    print(f"─" * 72)
    print(f"   target_porosity = {result['target_porosity']:.4f}")
    print(f"   target_df        = {result['target_df']}")

    # TOE框架孔隙率原理 (§1.3密度换算 + §22三重修正)
    print(f"\n📐 TOE框架孔隙率原理 (§1.3密度换算 + §22三重修正):")
    print(f"   ── 第一步: 本征空隙 (定值, 不受修正影响) ──")
    print(f"   φ_0 = 1 - f_m = 1 - ρ_bulk/ρ_theoretical")
    print(f"   ρ_bulk = {result['density']:.3f} g/cm³ (实际密度)")
    print(f"   ρ_theoretical = {result.get('rho_theoretical', result['density']):.3f} g/cm³ (X射线理论密度)")
    print(f"   f_m = ρ_bulk/ρ_theoretical = {result.get('f_m', 1.0):.4f} (物质体积占比)")
    print(f"   φ_0 = {result.get('phi_0', 0)*100:.2f}% (材料自带本征孔隙, 定值)")
    print(f"")
    print(f"   ── 第二步: 挖孔三重修正 ──")
    print(f"   修正1: φ_空气 = {result.get('phi_air', 0)*100:.2f}% (常压空气填充)")
    print(f"   修正2: φ_0泡  = {result.get('phi_0bubble', 0)*100:.2f}% (0维世界泡填充, 真空本底固有)")
    print(f"   修正3: Δφ_分形= {result.get('delta_phi_fractal', 0)*100:.2f}% (孔壁分形凹凸挤占)")
    print(f"   三重修正总量 = {result.get('total_correction', 0)*100:.2f}%")
    print(f"")
    print(f"   ── 第三步: 有效物理空隙定值69.7% ──")
    print(f"   φ_总 = φ_0 + φ_挖实 = 69.7% = Ω_Λ (净空有效空隙)")
    print(f"   φ_挖实 = {result.get('phi_drilling_effective', 0)*100:.2f}% (扣除修正后的真实有效空隙)")
    print(f"   φ_挖表 = {result.get('phi_drilling_geometric', 0)*100:.2f}% (几何设计开孔尺寸)")
    print(f"   φ_挖表 = 69.7% - φ_0 + Δφ_分形 + φ_空气 + φ_0泡")
    print(f"         = 69.7% - {result.get('phi_0', 0)*100:.2f}% + {result.get('total_correction', 0)*100:.2f}%")
    print(f"         = {result.get('phi_drilling_geometric', 0)*100:.2f}%")
    print(f"")
    print(f"   V_固% = {result['solid_volume_fraction']*100:.1f}% = Ω_m (固体基体=超导导电主体)")
    print(f"   五大约束: 冷媒流通 | 磁通钉扎 | 库柏对限域 | 热膨胀缓冲 | 力学形变")
    print(f"")
    print(f"   ── 容错校验 ──")
    phi_0_val = result.get('phi_0', 0)
    phi_drill_eff_val = result.get('phi_drilling_effective', 0)
    phi_drill_geo_val = result.get('phi_drilling_geometric', 0)
    total_corr_val = result.get('total_correction', 0)
    # 若忽略修正直接按φ_挖表=69.7%-φ_0开孔
    naive_drill_geo = 0.697 - phi_0_val
    naive_drill_eff = naive_drill_geo - total_corr_val
    naive_total = phi_0_val + naive_drill_eff
    print(f"   ❌ 忽略修正: φ_挖表=69.7%-φ_0={naive_drill_geo*100:.2f}%")
    print(f"      → φ_挖实={naive_drill_geo*100:.2f}%-{total_corr_val*100:.2f}%={naive_drill_eff*100:.2f}%")
    print(f"      → φ_总=φ_0+φ_挖实={phi_0_val*100:.2f}%+{naive_drill_eff*100:.2f}%={naive_total*100:.2f}%<69.7% → 失超!")
    print(f"   ✅ 正确设计: φ_挖表={phi_drill_geo_val*100:.2f}% (已补偿三重修正)")
    print(f"      → φ_挖实={phi_drill_geo_val*100:.2f}%-{total_corr_val*100:.2f}%={phi_drill_eff_val*100:.2f}%")
    print(f"      → φ_总=φ_0+φ_挖实={phi_0_val*100:.2f}%+{phi_drill_eff_val*100:.2f}%={(phi_0_val+phi_drill_eff_val)*100:.2f}%=69.7% ✓")
    print(f"")
    print(f"   示例调用:")
    print(f"   generate_TOE_cube_fractal(target_porosity={result['target_porosity']:.4f}, target_df={result['target_df']})")
    print(f"   generate_TOE_simple(target_porosity={result['target_porosity']:.4f}, target_df={result['target_df']})")
    print(f"   generate_fractal_gyroid(target_porosity={result['target_porosity']:.4f}, target_df={result['target_df']})")
    print(f"─" * 72)

    # McMillan验算
    Tc_check = mcmillan_tc(result['Theta_D_eff'], result['lambda_eff'], result['mu_star_eff'])
    print(f"\n✅ McMillan公式验算: T_c = {Tc_check:.1f} K")

    # T_c vs P 关系曲线
    if all_results and len(all_results) > 10:
        print(f"\n📊 T_c vs 孔隙率P 关系:")
        tc_values = [r["Tc"] for r in all_results]
        p_values = [r["P"] for r in all_results]
        tc_max = max(tc_values) if tc_values else 1
        if tc_max > 0:
            n_cols = 50
            n_rows = 20
            step = max(1, len(p_values) // n_rows)
            print(f"   {'P':>6s} │ T_c(K)")
            print(f"   {'─'*6}─┼{'─'*n_cols}")
            for i in range(0, len(p_values), step):
                p = p_values[i]
                tc = tc_values[i]
                bar_len = int(tc / max(tc_max, 1) * n_cols)
                bar_len = min(bar_len, n_cols)
                marker = " ◀最优" if abs(p - result['target_porosity']) < 0.006 else ""
                print(f"   {p:6.3f} │{'█' * bar_len}{marker}")

    print(f"\n{'=' * 72}\n")


# =====================================================================
# 从密度估算材料参数
# =====================================================================

def estimate_material_from_density(density, Theta_D=None, material_class="auto"):
    """
    从密度估算材料参数（粗略估算，建议提供更精确参数）

    经验关系:
    - Θ_D ∝ √(弹性模量/原子质量) — 轻元素材料Θ_D高
    - λ_0 ∝ N(0) × V — 与费米面态密度有关
    - μ*_0 ≈ 0.10-0.15 — 库仑赝势，大多数金属接近
    - C_λ ∝ Θ_D/160 — 德拜温度越高，分形增强越强
    - C_μ ∝ Θ_D/320 — 德拜温度越高，屏蔽越好

    参数:
        density: 本征密度 (g/cm³)
        Theta_D: 德拜温度 (K), None则自动估算
        material_class: 材料类别 ("metal"/"ceramic"/"polymer"/"auto")
    """
    # 自动判断材料类别
    if material_class == "auto":
        if density < 2.0:
            material_class = "polymer"
        elif density < 4.0:
            material_class = "ceramic"  # 可能是轻元素化合物
        else:
            material_class = "metal"

    # 估算德拜温度
    if Theta_D is None:
        if material_class == "polymer":
            Theta_D = 150 + 50 * (1.0 / max(density, 0.5)) ** 0.3
        elif material_class == "ceramic":
            # 轻元素化合物（碳、硼、锂等）Θ_D高
            Theta_D = 600 + 200 * (3.0 / max(density, 1.0)) ** 0.5
        else:  # metal
            Theta_D = 300 + 100 * (8.0 / max(density, 2.0)) ** 0.3
        Theta_D = round(Theta_D)

    # 估算耦合常数
    if material_class == "polymer":
        lambda_0 = 0.05
    elif material_class == "ceramic":
        lambda_0 = 1.0 + 0.5 * (3.0 / max(density, 1.0)) ** 0.3
    else:
        lambda_0 = 0.5 * (8.0 / max(density, 2.0)) ** 0.3
    lambda_0 = max(0.05, min(2.5, lambda_0))

    # 库仑赝势
    if material_class == "polymer":
        mu_star_0 = 0.22
    elif material_class == "ceramic":
        mu_star_0 = 0.10
    else:
        mu_star_0 = 0.12

    # 分形增强系数 (校准: Δλ=C_λ×P×√(1-P), 与Θ_D成正比)
    C_lambda = Theta_D / 220
    C_mu = Theta_D / 320

    class_names = {"polymer": "聚合物/塑料", "ceramic": "轻元素化合物/陶瓷", "metal": "金属"}

    return {
        "name": f"自定义(ρ={density:.2f}g/cm³, {class_names.get(material_class, '未知')})",
        "category": f"自定义/{class_names.get(material_class, '未知')}",
        "density": density,
        "density_theoretical": density * 1.05 if material_class == "ceramic" else density,  # 陶瓷烧结体通常有~5%致密度亏损
        "Theta_D": Theta_D,
        "lambda_0": round(lambda_0, 3),
        "mu_star_0": mu_star_0,
        "C_lambda": round(C_lambda, 2),
        "C_mu": round(C_mu, 2),
        "note": f"从密度={density}g/cm³估算，建议提供精确的Θ_D/λ_0/μ*_0"
    }


# =====================================================================
# 列表材料
# =====================================================================

def list_materials():
    """列出所有可用材料"""
    print("\n" + "=" * 72)
    print("          【可用材料列表】")
    print("=" * 72)

    categories = {}
    for key, mat in MATERIALS.items():
        cat = mat["category"]
        if cat not in categories:
            categories[cat] = []
        categories[cat].append((key, mat))

    for cat, materials in categories.items():
        print(f"\n📂 {cat}:")
        for key, mat in materials:
            tc_quick = mcmillan_tc(
                mat["Theta_D"],
                mat["lambda_0"] + mat["C_lambda"] * 0.40 * (1 - 0.40)**0.5,
                mat["mu_star_0"] / (1 + mat["C_mu"] * 0.40)
            )
            # 计算本征孔隙率 + 三重修正
            rho_theo = mat.get("density_theoretical", mat["density"])
            phi_0 = 1 - mat["density"] / rho_theo
            phi_drilling_eff = 0.697 - phi_0
            # 三重修正默认值
            phi_air = 0.015
            phi_0b = 0.003
            dphi_f = 0.008
            phi_drilling_geo = phi_drilling_eff + phi_air + phi_0b + dphi_f
            print(f"   {key:30s} │ ρ={mat['density']:.2f}g/cm³ "
                  f"│ Θ_D={mat['Theta_D']:4.0f}K "
                  f"│ φ_0={phi_0*100:5.1f}% "
                  f"│ φ_挖实={phi_drilling_eff*100:5.1f}% "
                  f"│ φ_挖表={phi_drilling_geo*100:5.1f}% "
                  f"│ T_c(P=40%)≈{tc_quick:6.1f}K "
                  f"│ {mat['note']}")

    print(f"\n💡 使用方法:")
    print(f"   python porosity_calculator.py --material <材料名>")
    print(f"   python porosity_calculator.py --density <密度> [--theta_D <Θ_D>]")
    print(f"   python porosity_calculator.py --list-materials")
    print()


# =====================================================================
# 命令行接口
# =====================================================================

def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="TOE超导体孔隙率计算器 v1.0 — 基于TOE_v8.8框架",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:
  %(prog)s --material graphene_superlattice
  %(prog)s --density 2.2 --theta_D 1200
  %(prog)s --list-materials
  %(prog)s --material graphene_superlattice --target-tc 250
  %(prog)s --material graphene_superlattice --json
  %(prog)s --density 8.96 --theta_D 343 --lambda_0 0.14 --C_lambda 2.1
        """
    )

    parser.add_argument("--material", type=str, help="材料名称(从数据库选择)")
    parser.add_argument("--preset", type=str, default=None,
                        help="已验证预设名(如hydrogenated_copper_oxide_vacuum)，自动填充全部校验参数")
    parser.add_argument("--density", type=float, help="自定义材料本征密度(g/cm³)")
    parser.add_argument("--density-theoretical", type=float, help="X射线理论密度(g/cm³)，用于计算φ_0")
    parser.add_argument("--phi-air", type=float, default=0.015, help="空气填充修正率(默认0.015=1.5%%, 真空工况设0)")
    parser.add_argument("--phi-0bubble", type=float, default=0.003, help="0维世界泡填充修正率(默认0.003=0.3%%)")
    parser.add_argument("--delta-phi-fractal", type=float, default=0.008, help="分形凹凸挤占修正率(默认0.008=0.8%%)")
    parser.add_argument("--theta_D", type=float, help="德拜温度(K)")
    parser.add_argument("--lambda_0", type=float, help="电声耦合常数")
    parser.add_argument("--mu_star_0", type=float, help="库仑赝势")
    parser.add_argument("--C_lambda", type=float, help="分形耦合增强系数")
    parser.add_argument("--C_mu", type=float, help="分形屏蔽系数")
    parser.add_argument("--d_f", type=float, default=2.973,
                        help="分形维度(默认2.973, 范围2.92-2.98)")
    parser.add_argument("--target-tc", type=float, default=200.0,
                        help="目标临界温度K(默认200)")
    parser.add_argument("--list-materials", action="store_true",
                        help="列出所有可用材料")
    parser.add_argument("--json", action="store_true",
                        help="JSON格式输出")
    parser.add_argument("--quiet", action="store_true",
                        help="安静模式，只输出关键数值")
    parser.add_argument("--gyroid-params", action="store_true",
                        help="输出Gyroid生成器可直接使用的参数")

    args = parser.parse_args()

    if args.list_materials:
        list_materials()
        return

    # ──── 预设处理 ────
    if args.preset:
        if args.preset not in VERIFIED_PRESETS:
            print(f"❌ 未知预设: {args.preset}")
            print(f"   可用预设: {', '.join(VERIFIED_PRESETS.keys())}")
            return
        preset = VERIFIED_PRESETS[args.preset]
        mat_key = preset["material"]
        mat_params = MATERIALS[mat_key].copy()
        # 用预设中的修正参数覆盖CLI默认值
        args.phi_air = preset["phi_air"]
        args.phi_0bubble = preset["phi_0bubble"]
        args.delta_phi_fractal = preset["delta_phi_fractal"]
        args.d_f = preset["target_df"]
        args.target_tc = preset["Tc_target"]
        if not args.material:
            args.material = mat_key
        print(f"\n📋 已加载预设: {preset['name']}")
        print(f"   φ_0={preset['phi0']*100:.2f}% | φ_挖实={preset['phi_drill_effective']*100:.2f}% | φ_挖表={preset['phi_drill_apparent']*100:.2f}%")
        print(f"   λ_eff={preset['lambda_eff']} | μ*_eff={preset['mu_star_eff']} | Tc_max={preset['Tc_max']}K")
        print(f"   三重修正: φ_空气={preset['phi_air']*100:.2f}% + φ_0泡={preset['phi_0bubble']*100:.2f}% + Δφ_分形={preset['delta_phi_fractal']*100:.2f}% = {preset['dphi_sum']*100:.2f}%")

    # 确定材料参数
    if args.material:
        if args.material not in MATERIALS:
            print(f"❌ 未知材料: {args.material}")
            print(f"   使用 --list-materials 查看可用材料")
            return
        mat_params = MATERIALS[args.material].copy()
    elif args.density:
        mat_params = estimate_material_from_density(args.density, args.theta_D)
    else:
        print("❌ 请指定材料(--material)或密度(--density)")
        print("   使用 --help 查看帮助")
        return

    # 覆盖自定义参数
    if args.theta_D:
        mat_params["Theta_D"] = args.theta_D
    if args.lambda_0:
        mat_params["lambda_0"] = args.lambda_0
    if args.mu_star_0:
        mat_params["mu_star_0"] = args.mu_star_0
    if args.C_lambda:
        mat_params["C_lambda"] = args.C_lambda
    if args.C_mu:
        mat_params["C_mu"] = args.C_mu
    if args.density_theoretical:
        mat_params["density_theoretical"] = args.density_theoretical

    # 计算最优孔隙率
    result = find_optimal_porosity(
        mat_params,
        d_f=args.d_f,
        target_tc=args.target_tc,
        verbose=not args.quiet,
        phi_air=args.phi_air,
        phi_0bubble=args.phi_0bubble,
        delta_phi_fractal=args.delta_phi_fractal,
    )

    # JSON输出
    if args.json:
        json_result = {}
        for k, v in result.items():
            if isinstance(v, (np.integer, np.floating)):
                json_result[k] = float(v)
            elif isinstance(v, np.bool_):
                json_result[k] = bool(v)
            else:
                json_result[k] = v
        print(json.dumps(json_result, indent=2, ensure_ascii=False))

    # 安静模式
    if args.quiet:
        print(f"target_porosity = {result['target_porosity']:.4f}")
        print(f"target_df = {result['target_df']}")
        print(f"Tc_max = {result['Tc_max']:.1f} K")

    # Gyroid生成器参数输出
    if args.gyroid_params:
        print(f"\n{'─'*60}")
        print(f"🔧 Gyroid生成器参数 (直接复制):")
        print(f"{'─'*60}")
        print(f"python superconductor_gyroid.py \\")
        print(f"  --porosity {result['target_porosity']:.4f} \\")
        print(f"  --df {result['target_df']} \\")
        print(f"  --material {args.material if args.material else 'custom'}")
        print(f"{'─'*60}")


if __name__ == "__main__":
    main()
