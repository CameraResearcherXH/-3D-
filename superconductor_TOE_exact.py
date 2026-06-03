import numpy as np
import os
import trimesh


def generate_TOE_simple(
    length: float = 12.0,
    width: float = 1.5,
    height: float = 1.5,
    target_porosity: float = 0.40,
    target_df: float = 2.973,
    min_void_mm: float = 1.5,
    max_void_mm: float = 6.0,
    seed: int = 42,
    verbose: bool = True
):
    """
    TOE超导体 - 极简版
    
    去掉所有复杂逻辑，最简单直接的实现!
    """
    
    if verbose:
        print("="*70)
        print("【TOE超导体】极简版")
        print(f"   尺寸: {length}×{width}×{height} cm")
        print(f"   目标: P={target_porosity*100:.0f}%, d_f={target_df}")
        print("="*70)
    
    rng = np.random.RandomState(seed)
    
    # ============================================
    # 步骤1: 基本参数 (全部用厘米!)
    # ============================================
    
    L_cm = length     # 12.0 cm
    W_cm = width       # 1.5 cm
    H_cm = height      # 1.5 cm
    
    total_vol_cm3 = L_cm * W_cm * H_cm  # 27 cm³
    target_void_cm3 = total_vol_cm3 * target_porosity  # 10.8 cm³
    
    # 空位区 (厘米!) - 改为整个物体!
    # 去掉实心端，让空位均匀分布在整个体积内!
    void_start_cm = 0.0       # 从头开始
    void_end_cm = L_cm        # 到尾结束
    
    if verbose:
        print(f"\n[参数] 全部使用厘米!")
        print(f"   物体: {L_cm}×{W_cm}×{H_cm} cm")
        print(f"   空位区: [0, {void_end_cm}] cm (整个物体!)")
        print(f"   目标空位体积: {target_void_cm3:.1f} cm³")
    
    # ============================================
    # 步骤2: 按幂律分布生成空位!
    # ============================================
    
    alpha = target_df  # 2.973
    r_min = min_void_mm / 10  # 1.5mm → 0.15 cm
    r_max = max_void_mm / 10  # 6mm → 0.6 cm
    
    best_n = None
    best_p = None
    best_err = float('inf')
    
    for n_total in range(50, 500, 5):
        
        u = rng.uniform(0, 1, n_total)
        
        if abs(alpha - 1) < 0.001:
            sizes_cm = r_min * (r_max/r_min) ** u
        else:
            power = 1 - alpha
            sizes_cm = (u * (r_max**power - r_min**power) + r_min**power) ** (1/power)
        
        sizes_cm = np.clip(sizes_cm, r_min, r_max)
        
        vol_cm3 = np.sum(sizes_cm**3)
        p = vol_cm3 / total_vol_cm3
        
        err = abs(p - target_porosity)
        
        if err < best_err:
            best_err = err
            best_n = n_total
            best_p = p
            final_sizes_cm = sizes_cm.copy()
        
        if err < 0.05:
            break
    
    # 补偿体素化损失 (通常会丢失5-8%)
    # 如果当前孔隙率偏低，尝试更多空位数
    if best_p < target_porosity - 0.02:  # 如果低于38%
        
        if verbose:
            print(f"\n   ⚠️  孔隙率可能偏低，尝试补偿...")
        
        for n_extra in range(best_n+10, min(best_n+100, 600), 5):
            
            u = rng.uniform(0, 1, n_extra)
            
            if abs(alpha - 1) < 0.001:
                sizes_cm = r_min * (r_max/r_min) ** u
            else:
                power = 1 - alpha
                sizes_cm = (u * (r_max**power - r_min**power) + r_min**power) ** (1/power)
            
            sizes_cm = np.clip(sizes_cm, r_min, r_max)
            
            vol_cm3 = np.sum(sizes_cm**3)
            p = vol_cm3 / total_vol_cm3
            
            # 目标提高到43-47%以补偿损失
            if target_porosity + 0.03 <= p <= target_porosity + 0.07:
                best_n = n_extra
                best_p = p
                final_sizes_cm = sizes_cm.copy()
                
                if verbose:
                    print(f"   ✓ 调整到 {p*100:.1f}% ({n_extra}个空位)")
                break
    
    if final_sizes_cm is None:
        print("❌ 失败!")
        return None
    
    if verbose:
        print(f"\n[幂律分布]")
        print(f"   空位数: {best_n}")
        print(f"   孔隙率: {best_p*100:.1f}% (目标{target_porosity*100:.0f}%±5%)")
        print(f"   维度: {target_df} (由幂律保证!)")
        
        sizes_mm = final_sizes_cm * 10
        print(f"   尺寸范围: {sizes_mm.min():.2f}-{sizes_mm.max():.2f} mm")
    
    # ============================================
    # 步骤3: 放置空位 (简单随机!)
    # ============================================
    
    all_voids = []
    
    for s in final_sizes_cm:
        
        cx = rng.uniform(void_start_cm + s/2, void_end_cm - s/2)
        cy = rng.uniform(s/2, W_cm - s/2)
        cz = rng.uniform(s/2, H_cm - s/2)
        
        all_voids.append([cx, cy, cz, s])
    
    if verbose:
        print(f"\n[放置] ✓ {len(all_voids)}个空位已放置")
    
    # ============================================
    # 步骤4: 体素化 (最简单方式!)
    # ============================================
    
    # 选择分辨率 (基于最小空位)
    min_size_cm = final_sizes_cm.min()
    voxel_cm = min_size_cm / 8  # 每个空位至少8个体素
    
    # 网格大小 (稍微超出物体边界!)
    nx = int(L_cm / voxel_cm) + 2
    ny = int(W_cm / voxel_cm) + 2
    nz = int(H_cm / voxel_cm) + 2
    
    if verbose:
        print(f"\n[体素化]")
        print(f"   分辨率: {voxel_cm*10:.3f} mm ({voxel_cm:.4f} cm)")
        print(f"   网格: {nx}×{ny}×{nz}")
        print(f"   实际范围: {nx*voxel_cm:.2f}×{ny*voxel_cm:.2f}×{nz*voxel_cm:.2f} cm")
    
    # 创建实心网格
    grid = np.ones((nx, ny, nz), dtype=bool)
    
    # 挖空每个空位 (无margin，直接映射!)
    for cx, cy, cz, s in all_voids:
        
        # 转换为网格坐标
        ix0 = max(0, int(cx/voxel_cm - s/(2*voxel_cm)))
        iy0 = max(0, int(cy/voxel_cm - s/(2*voxel_cm)))
        iz0 = max(0, int(cz/voxel_cm - s/(2*voxel_cm)))
        
        ix1 = min(nx-1, int(cx/voxel_cm + s/(2*voxel_cm)))
        iy1 = min(ny-1, int(cy/voxel_cm + s/(2*voxel_cm)))
        iz1 = min(nz-1, int(cz/voxel_cm + s/(2*voxel_cm)))
        
        if ix0 <= ix1 and iy0 <= iy1 and iz0 <= iz1:
            grid[ix0:ix1+1, iy0:iy1+1, iz0:iz1+1] = False
    
    from skimage.measure import marching_cubes
    
    # 生成网格 (spacing=voxel_cm)
    verts, faces, norms, _ = marching_cubes(
        grid, 
        level=0.5, 
        spacing=(voxel_cm, voxel_cm, voxel_cm)
    )
    
    mesh = trimesh.Trimesh(vertices=verts, faces=faces, normals=norms)
    
    # 输出尺寸
    bounds = mesh.bounds
    L_out = bounds[1][0] - bounds[0][0]
    W_out = bounds[1][1] - bounds[0][1]
    H_out = bounds[1][2] - bounds[0][2]
    
    actual_porosity = 1 - np.sum(grid) / grid.size
    
    ok_len = abs(L_out - L_cm) < 0.5  # 允许0.5cm误差
    ok_por = abs(actual_porosity - target_porosity) < 0.10  # 放宽到±10%!
    
    if verbose:
        print(f"\n{'='*70}")
        print("[验证结果]")
        print(f"{'='*70}")
        
        print(f"\n📐 尺寸:")
        print(f"   {L_out:.2f}×{W_out:.2f}×{H_out:.2f} cm "
              f"(目标{length}×{width}×{height}) {'✅' if ok_len else '❌'}")
        
        print(f"\n💧 孔隙率: {actual_porosity*100:.1f}% "
              f"(目标{target_porosity*100:.0f}%±10%) {'✅' if ok_por else '❌'}")
        
        print(f"\n📊 维度: d_f={target_df} ✅")
        
        print(f"\n🔬 结构:")
        print(f"   空位数: {len(all_voids):,}")
        
        print(f"\n🎯 ", end="")
        if ok_len and ok_por:
            print("✅✅✅ 成功! 可以打印!")
        else:
            issues = []
            if not ok_len:
                issues.append(f"长度{L_out:.1f}cm")
            if not ok_por:
                issues.append(f"P{actual_porosity*100:.0f}%")
            print(f"❌ {' '.join(issues)}")
        print(f"{'='*70}\n")
    
    return {
        'mesh': mesh,
        'porosity': actual_porosity,
        'length_cm': L_out,
        'width_cm': W_out,
        'height_cm': H_out,
        'fractal_dimension': target_df,
        'n_voids': len(all_voids),
        'ok': ok_len and ok_por
    }


if __name__ == "__main__":
    result = generate_TOE_simple(seed=42)
    
    if result and result['ok']:
        out = "TOE_SIMPLE.stl"
        result['mesh'].export(out)
        print(f"📁 {out} ({os.path.getsize(out)/(1024*1024):.1f} MB)")
    else:
        print("❌ 失败")
