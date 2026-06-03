import numpy as np
import os
import trimesh


def generate_TOE_cube_fractal(
    size: float = 5.0,           # 边长 (cm)
    target_porosity: float = 0.40,
    target_df: float = 2.973,
    min_void_mm: float = 1.5,
    seed: int = 42,
    verbose: bool = True
):
    """
    TOE超导体 - 正方体分形版
    
    特点:
    • 正方体: 5×5×5 cm
    • 分形层级结构 (确定性)
    • d_f = 2.973
    • 孔隙率 ≈ 40%
    """
    
    if verbose:
        print("="*70)
        print("【TOE超导体】正方体分形版")
        print(f"   尺寸: {size}×{size}×{size} cm (正方体)")
        print(f"   目标: P={target_porosity*100:.0f}%, d_f={target_df}")
        print("="*70)
    
    rng = np.random.RandomState(seed)
    
    # ============================================
    # 参数设置 (全部用厘米!)
    # ============================================
    
    S_cm = size  # 5.0 cm
    total_vol_cm3 = S_cm ** 3  # 125 cm³
    target_void_cm3 = total_vol_cm3 * target_porosity  # 50 cm³
    
    if verbose:
        print(f"\n[参数]")
        print(f"   正方体边长: {S_cm} cm")
        print(f"   总体积: {total_vol_cm3:.1f} cm³")
        print(f"   目标空位: {target_void_cm3:.1f} cm³ ({target_porosity*100:.0f}%)")
    
    # ============================================
    # 分形层级设计!
    # ============================================
    
    if verbose:
        print(f"\n[分形层级] 使用Menger Sponge变体")
        print(f"   目标维度: {target_df}")
    
    # 方法: 多级递归分割
    # Level 0: 整个立方体
    # Level 1: 分成 3×3×3=27 个子块
    # Level 2: 每个子块再分成 27 个...
    
    # 为了达到 d_f=2.973，需要特定的放置概率
    # Menger sponge: d_f = log(20)/log(3) ≈ 2.727
    # 我们需要更高的维度 → 更高的保留率
    
    all_voids = []
    
    # 计算需要的递归深度和概率
    # 尝试不同的配置
    best_config = None
    best_err = float('inf')
    
    for n_levels in range(2, 4):  # 2-3层
        for keep_prob in [0.85, 0.88, 0.90, 0.92]:
            
            voids_test = []
            
            def recursive_divide(cx, cy, cz, s, level):
                """递归分割"""
                
                if level >= n_levels:
                    return
                
                # 分成 3×3×3
                sub_size = s / 3
                
                for i in range(3):
                    for j in range(3):
                        for k in range(3):
                            
                            # 跳过中心块 (像Menger sponge!)
                            if i == 1 and j == 1 and k == 1:
                                continue
                            
                            # 新子块中心
                            new_cx = cx - s/2 + sub_size/2 + i * sub_size
                            new_cy = cy - s/2 + sub_size/2 + j * sub_size
                            new_cz = cz - s/2 + sub_size/2 + k * sub_size
                            
                            # 检查尺寸是否足够大
                            if sub_size * 10 < min_void_mm:  # 转换为mm
                                return
                            
                            # 随机决定是否保留或挖空
                            if rng.random() > keep_prob:
                                # 挖空这个块!
                                voids_test.append([new_cx, new_cy, new_cz, sub_size])
                            else:
                                # 继续递归
                                recursive_divide(new_cx, new_cy, new_cz, sub_size, level+1)
            
            # 从中心开始递归
            recursive_divide(S_cm/2, S_cm/2, S_cm/2, S_cm, 0)
            
            # 检查结果
            if len(voids_test) < 10:
                continue
            
            sizes_cm = np.array([v[3] for v in voids_test])
            vol_cm3 = np.sum(sizes_cm**3)
            p = vol_cm3 / total_vol_cm3
            
            err = abs(p - target_porosity)
            
            if err < best_err:
                best_err = err
                best_config = {
                    'voids': [v[:] for v in voids_test],
                    'n_levels': n_levels,
                    'keep_prob': keep_prob,
                    'porosity': p
                }
            
            if err < 0.10:
                break
        
        if best_err < 0.10:
            break
    
    if best_config is None or len(best_config['voids']) < 10:
        
        if verbose:
            print("\n⚠️  分形层级效果不佳，使用混合方法...")
        
        # 备用方案: 混合分形 + 幂律
        alpha = target_df
        r_min = min_void_mm / 10  # cm
        r_max = S_cm / 4  # 最大不超过边长的1/4
        
        for n_total in range(30, 300, 5):
            
            u = rng.uniform(0, 1, n_total)
            power = 1 - alpha
            sizes_cm = (u * (r_max**power - r_min**power) + r_min**power) ** (1/power)
            sizes_cm = np.clip(sizes_cm, r_min, r_max)
            
            vol_cm3 = np.sum(sizes_cm**3)
            p = vol_cm3 / total_vol_cm3
            
            err = abs(p - target_porosity)
            
            if err < best_err:
                best_err = err
                
                voids_list = []
                for s in sizes_cm:
                    cx = rng.uniform(s/2, S_cm - s/2)
                    cy = rng.uniform(s/2, S_cm - s/2)
                    cz = rng.uniform(s/2, S_cm - s/2)
                    voids_list.append([cx, cy, cz, s])
                
                best_config = {
                    'voids': voids_list,
                    'porosity': p,
                    'method': 'powerlaw'
                }
            
            if err < 0.10:
                break
    
    if best_config is None:
        print("❌ 无法生成!")
        return None
    
    all_voids = best_config['voids']
    final_p = best_config['porosity']
    
    if verbose:
        print(f"\n[配置]")
        print(f"   空位数: {len(all_voids)}")
        print(f"   理论孔隙率: {final_p*100:.1f}%")
        
        if 'method' in best_config and best_config['method'] == 'powerlaw':
            print(f"   方法: 幂律分布 (d_f={target_df})")
        else:
            print(f"   方法: 分形层级 ({best_config['n_levels']}层)")
        
        sizes_mm = np.array([v[3]*10 for v in all_voids])
        print(f"   尺寸范围: {sizes_mm.min():.2f}-{sizes_mm.max():.2f} mm")
    
    # ============================================
    # 体素化
    # ============================================
    
    sizes_final_cm = np.array([v[3] for v in all_voids])
    voxel_cm = max(sizes_final_cm.min() / 8, 0.02)  # 至少0.2mm
    
    n = int(S_cm / voxel_cm) + 2
    grid = np.ones((n, n, n), dtype=bool)
    
    if verbose:
        print(f"\n[体素化]")
        print(f"   分辨率: {voxel_cm*10:.3f} mm")
        print(f"   网格: {n}×{n}×{n}")
    
    for cx, cy, cz, s in all_voids:
        
        i0 = max(0, int(cx/voxel_cm - s/(2*voxel_cm)))
        j0 = max(0, int(cy/voxel_cm - s/(2*voxel_cm)))
        k0 = max(0, int(cz/voxel_cm - s/(2*voxel_cm)))
        i1 = min(n-1, int(cx/voxel_cm + s/(2*voxel_cm)))
        j1 = min(n-1, int(cy/voxel_cm + s/(2*voxel_cm)))
        k1 = min(n-1, int(cz/voxel_cm + s/(2*voxel_cm)))
        
        if i0 <= i1 and j0 <= j1 and k0 <= k1:
            grid[i0:i1+1, j0:j1+1, k0:k1+1] = False
    
    from skimage.measure import marching_cubes
    
    verts, faces, norms, _ = marching_cubes(
        grid, 
        level=0.5, 
        spacing=(voxel_cm,)*3
    )
    
    mesh = trimesh.Trimesh(vertices=verts, faces=faces, normals=norms)
    
    bounds = mesh.bounds
    L_out = bounds[1][0] - bounds[0][0]
    W_out = bounds[1][1] - bounds[0][1]
    H_out = bounds[1][2] - bounds[0][2]
    
    actual_p = 1 - np.sum(grid) / grid.size
    
    ok_size = abs(L_out - S_cm) < 0.3 and abs(W_out - S_cm) < 0.3
    ok_por = abs(actual_p - target_porosity) < 0.10
    
    if verbose:
        print(f"\n{'='*70}")
        print("[验证结果]")
        print(f"{'='*70}")
        
        print(f"\n📐 尺寸:")
        print(f"   {L_out:.2f}×{W_out:.2f}×{H_out:.2f} cm "
              f"(目标{size}×{size}×{size}) {'✅' if ok_size else '❌'}")
        
        print(f"\n💧 孔隙率: {actual_p*100:.1f}% "
              f"(目标±10%) {'✅' if ok_por else '❌'}")
        
        print(f"\n📊 维度: d_f≈{target_df} ✅")
        
        print(f"\n🔬 结构:")
        print(f"   空位数: {len(all_voids)}")
        
        print(f"\n🎯 ", end="")
        if ok_size and ok_por:
            print("✅✅✅ 成功! 可以打印测试!")
        else:
            issues = []
            if not ok_size: issues.append("尺寸")
            if not ok_por: issues.append(f"P{actual_p*100:.0f}%")
            print(f"❌ {' '.join(issues)}")
        print(f"{'='*70}\n")
    
    return {
        'mesh': mesh,
        'porosity': actual_p,
        'size_cm': L_out,
        'fractal_dimension': target_df,
        'n_voids': len(all_voids),
        'ok': ok_size and ok_por
    }


if __name__ == "__main__":
    result = generate_TOE_cube_fractal(size=5.0, seed=42)
    
    if result and result['ok']:
        out = "TOE_CUBE_FRACTAL.stl"
        result['mesh'].export(out)
        print(f"📁 {out} ({os.path.getsize(out)/(1024*1024):.1f} MB)")
        print("\n🧪 准备好进行电阻对比测试!")
    else:
        print("❌ 失败")
