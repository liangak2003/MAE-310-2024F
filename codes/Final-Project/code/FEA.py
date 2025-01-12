import numpy as np
import meshio
import scipy.sparse as sp
import scipy.sparse.linalg
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
import os
import subprocess
import shutil

# 配置参数部分
# -----------------

# Gmsh 配置
geo_file = "geometry.geo"    # .geo 文件路径
mesh_file = "geometry.msh"   # 生成的 .msh 文件路径

# 模型类型：'plane_stress' 或 'plane_strain'
model_type = 'plane_stress'

# 材料属性（均匀各向同性弹性体）
E = 1e9      # 杨氏模量 (Pa)
nu = 0.3     # 泊松比

# 圆孔半径
R = 0.5      # 根据网格定义调整

# 远场拉伸应力
T = 1e4      # 应力 (Pa)

# 变形放大倍数
deformation_scale = 5400  # 调整此值以放大变形

# 是否绘图
plot_config = {
    'deformed_mesh': True,       # Deformed mesh plot
    'von_mises_stress': True,    # Von Mises stress distribution plot
    'von_mises_strain': True,    # Von Mises strain distribution plot
    'exact_von_mises_stress': False  # Exact Von Mises stress distribution plot
}

# 是否计算误差
compute_errors_flag = True

# 边界条件定义
boundary_conditions = {
    'SymX': {
        'type': 'Dirichlet',
        'values': {'uy': lambda x, y: 0.0}
    },
    'SymY': {
        'type': 'Dirichlet',
        'values': {'ux': lambda x, y: 0.0}
    },
    'OuterX': {
        'type': 'Neumann',
        'values': {
            'tx': lambda x, y: T ,
            'ty': lambda x, y: 0
        }
    },
    'OuterY': {
        'type': 'Neumann',
        'values': {
            'tx': lambda x, y: 0,
            'ty': lambda x, y: 0
        }
    },
 'Hole': {
        'type': 'Neumann',
        'values': {
            'tx': lambda x, y: 0,
            'ty': lambda x, y: 0
        }
    }
}


# 精确解定义（用于误差计算）
def exact_solution(x, y):
    r = np.sqrt(x**2 + y**2)
    theta = np.arctan2(y, x)
    r = np.where(r == 0, 1e-10, r)
    u_r = (T * R**2 / (2 * E)) * ((1 + nu) / r + (1 - nu) * (R**2) / r**3) * np.cos(theta)
    u_theta = (T * R**2 / (2 * E)) * ((1 + nu) / r - (1 + 3 * nu) * (R**2) / r**3) * np.sin(theta)
    u_x = u_r * np.cos(theta) - u_theta * np.sin(theta)
    u_y = u_r * np.sin(theta) + u_theta * np.cos(theta)
    return np.array([u_x, u_y])

# -----------------
# 配置结束

def read_mesh(mesh_file):
    mesh = meshio.read(mesh_file)
    # 节点坐标
    nodes = mesh.points[:, :2]
    # 单元连接关系
    elements = None
    for cell in mesh.cells:
        if cell.type == "triangle":
            elements = cell.data
            break
    if elements is None:
        raise ValueError("Mesh does not contain triangular elements.")
    # 读取物理组信息
    boundaries = {}
    if "line" in mesh.cell_data_dict.get("gmsh:physical", {}):
        boundary_data = mesh.cell_data_dict["gmsh:physical"]["line"]
        lines = [cell.data for cell in mesh.cells if cell.type == "line"]
        if not lines:
            raise ValueError("No line elements found in the mesh.")
        lines = lines[0]
        for i, phys in enumerate(boundary_data):
            if phys not in boundaries:
                boundaries[phys] = []
            boundaries[phys].append(lines[i])
    # 读取物理组标签名称
    phys_tags = {}
    for name, val in mesh.field_data.items():
        phys_tags[val[0]] = name
    return nodes, elements, boundaries, phys_tags

def plane_stress_matrix(E, nu):
    D = (E / (1 - nu**2)) * np.array([
        [1, nu, 0],
        [nu, 1, 0],
        [0, 0, (1 - nu) / 2]
    ])
    return sp.csr_matrix(D)

def plane_strain_matrix(E, nu):
    D = (E / ((1 + nu) * (1 - 2 * nu))) * np.array([
        [1 - nu, nu, 0],
        [nu, 1 - nu, 0],
        [0, 0, (1 - 2 * nu) / 2]
    ])
    return sp.csr_matrix(D)

def assemble_system(nodes, elements, D):
    num_nodes = nodes.shape[0]
    K = sp.lil_matrix((2*num_nodes, 2*num_nodes))
    F = np.zeros(2*num_nodes)
    for elem in elements:
        node_indices = elem
        x = nodes[node_indices, 0]
        y = nodes[node_indices, 1]
        # 计算三角形面积
        area_matrix = np.array([
            [1, x[0], y[0]],
            [1, x[1], y[1]],
            [1, x[2], y[2]]
        ])
        area = 0.5 * np.abs(np.linalg.det(area_matrix))
        if area < 1e-12:
            raise ValueError("Element with zero or near-zero area detected.")
        
        beta = np.array([y[1] - y[2], y[2] - y[0], y[0] - y[1]])
        gamma = np.array([x[2] - x[1], x[0] - x[2], x[1] - x[0]])
        B = np.zeros((3, 6))
        B[0, 0::2] = beta
        B[1, 1::2] = gamma
        B[2, 0::2] = gamma
        B[2, 1::2] = beta
        B /= (2 * area)
        
        # 计算 ke = B^T * D * B * area
        BT = B.T
        BD = BT @ D.toarray()
        ke = BD @ B * area
        
        # 组装到全局刚度矩阵
        dof = []
        for ni in node_indices:
            dof.extend([2*ni, 2*ni+1])
        for i in range(6):
            for j in range(6):
                K[dof[i], dof[j]] += ke[i, j]
    return K.tocsr(), F

def apply_boundary_conditions(K, F, nodes, boundaries, phys_tags, boundary_conditions):
    num_nodes = nodes.shape[0]
    fixed_dofs = []
    U = np.zeros(2*num_nodes)  # 初始化位移向量
    # 处理每个物理组的边界条件
    for phys_name, bc in boundary_conditions.items():
        # 查找物理组编号
        phys_id = None
        for key, name in phys_tags.items():
            if name == phys_name:
                phys_id = key
                break
        if phys_id is None:
            raise ValueError(f"Physical group '{phys_name}' not found in mesh.")
        # 获取边界边
        if phys_id not in boundaries:
            continue
        edges = boundaries[phys_id]
        # 获取节点索引
        node_indices = set()
        for edge in edges:
            node_indices.update(edge)
        node_indices = list(node_indices)
        if bc['type'] == 'Dirichlet':
            if bc['values'].get('ux') is not None:
                ux_func = bc['values']['ux']
                for ni in node_indices:
                    U[2*ni] = ux_func(nodes[ni, 0], nodes[ni, 1])
                fixed_dofs.extend([2*ni for ni in node_indices])
            if bc['values'].get('uy') is not None:
                uy_func = bc['values']['uy']
                for ni in node_indices:
                    U[2*ni+1] = uy_func(nodes[ni, 0], nodes[ni, 1])
                fixed_dofs.extend([2*ni+1 for ni in node_indices])
        elif bc['type'] == 'Neumann':
            tx_func = bc['values'].get('tx', lambda x, y: 0.0)
            ty_func = bc['values'].get('ty', lambda x, y: 0.0)
            # 定义2点高斯积分
            xi = [-1.0 / np.sqrt(3.0), 1.0 / np.sqrt(3.0)]
            weights = [1.0, 1.0]
            for edge in edges:
                n1, n2 = edge
                x1, y1 = nodes[n1]
                x2, y2 = nodes[n2]
                l = np.sqrt((x2 - x1)**2 + (y2 - y1)**2)
                for point, weight in zip(xi, weights):
                    # 计算全局坐标
                    x_g = (1 - point)/2 * x1 + (1 + point)/2 * x2
                    y_g = (1 - point)/2 * y1 + (1 + point)/2 * y2
                    # 计算面力
                    tx = tx_func(x_g, y_g)
                    ty = ty_func(x_g, y_g)
                    # 计算形状函数
                    N1 = (1 - point)/2
                    N2 = (1 + point)/2
                    # 计算节点力贡献
                    F[2*n1] += N1 * tx * weight * l / 2
                    F[2*n2] += N2 * tx * weight * l / 2
                    F[2*n1 + 1] += N1 * ty * weight * l / 2
                    F[2*n2 + 1] += N2 * ty * weight * l / 2
    # 去重并排序
    fixed_dofs = sorted(list(set(fixed_dofs)))
    free_dofs = np.setdiff1d(np.arange(K.shape[0]), fixed_dofs)
    # 应用 Dirichlet 边界条件
    F[free_dofs] -= K[free_dofs, :][:, fixed_dofs] @ U[fixed_dofs]
    # 组装缩减系统
    K_ff = K[free_dofs, :][:, free_dofs]
    F_f = F[free_dofs]
    return K_ff, F_f, free_dofs, fixed_dofs, U

def solve_system(K_ff, F_f):
    # 使用稀疏线性求解器
    U_f = sp.linalg.spsolve(K_ff, F_f)
    return U_f

def compute_strains_stresses(nodes, elements, U, D):
    strains = []
    stresses = []
    von_mises_strain = []
    for elem in elements:
        node_indices = elem
        x = nodes[node_indices, 0]
        y = nodes[node_indices, 1]
        # 计算三角形面积
        area_matrix = np.array([
            [1, x[0], y[0]],
            [1, x[1], y[1]],
            [1, x[2], y[2]]
        ])
        area = 0.5 * np.abs(np.linalg.det(area_matrix))
        if area < 1e-12:
            raise ValueError("Element with zero or near-zero area detected during stress computation.")
        
        beta = np.array([y[1] - y[2], y[2] - y[0], y[0] - y[1]])
        gamma = np.array([x[2] - x[1], x[0] - x[2], x[1] - x[0]])
        B = np.zeros((3, 6))
        B[0, 0::2] = beta
        B[1, 1::2] = gamma
        B[2, 0::2] = gamma
        B[2, 1::2] = beta
        B /= (2 * area)
        
        # 单元位移向量
        Ue = np.zeros(6)
        for i, ni in enumerate(node_indices):
            Ue[2*i] = U[2*ni]
            Ue[2*i+1] = U[2*ni+1]
        
        # 计算应变和应力
        strain = B @ Ue
        stress = D @ strain
        strains.append(strain)
        stresses.append(stress)
        
        # 计算等效应变
        eps_vm = np.sqrt((strain[0] - strain[1])**2 + 3*(strain[2]/2)**2)
        von_mises_strain.append(eps_vm)
    return np.array(strains), np.array(stresses), np.array(von_mises_strain)

def visualize_results(nodes, elements, U, strains, stresses, von_mises_strain, boundaries, phys_tags, boundary_conditions, deformation_scale):
    # 应用变形放大倍数
    U_scaled = U * deformation_scale
    
    # 创建变形后的节点坐标
    deformed = nodes + U_scaled.reshape(-1, 2)
    
    # 创建 Triangulation 对象
    triang_original = Triangulation(nodes[:,0], nodes[:,1], elements)
    triang_deformed = Triangulation(deformed[:,0], deformed[:,1], elements)
    
    # 初始化一个标志，判断是否有figure被创建
    plots_generated = False
    
    # Plot deformed mesh with original mesh
    if plot_config['deformed_mesh']:
        fig, ax = plt.subplots(figsize=(10, 10))
        # 绘制原始网格
        ax.triplot(triang_original, color='lightgrey', label='Original Mesh')
        # 绘制变形后的网格
        ax.triplot(triang_deformed, color='blue', alpha=0.7, label='Deformed Mesh')
        
        # 使用 plt.cm.tab20 colormap 生成颜色循环
        cmap = plt.cm.get_cmap('tab20')
        colors = [cmap(i) for i in range(len(boundaries))]
        
        # 记录图例项
        legend_handles = []
        legend_labels = []
        
        # 绘制边界边，并根据物理组名称设置颜色
        for idx, (phys_id, edges) in enumerate(boundaries.items()):
            if phys_id not in phys_tags:
                continue
            phys_name = phys_tags[phys_id]
            if phys_name not in boundary_conditions:
                continue
            bc = boundary_conditions[phys_name]
            bc_type = bc['type']
            color = colors[idx % len(colors)]
            for edge in edges:
                n1, n2 = edge
                x = [nodes[n1,0], nodes[n2,0]]
                y = [nodes[n1,1], nodes[n2,1]]
                line, = ax.plot(x, y, color=color, linewidth=2)
            
            # 添加图例项
            label = f"{phys_name} ({bc_type})"
            legend_handles.append(line)
            legend_labels.append(label)
        
        # 绘制 Neumann 边界条件的力箭头在原始网格上
        for phys_id, edges in boundaries.items():
            if phys_id not in phys_tags:
                continue
            phys_name = phys_tags[phys_id]
            if phys_name not in boundary_conditions:
                continue
            bc = boundary_conditions[phys_name]
            if bc['type'] == 'Neumann':
                tx_func = bc['values'].get('tx', lambda x, y: 0.0)
                ty_func = bc['values'].get('ty', lambda x, y: 0.0)
                for edge in edges:
                    n1, n2 = edge
                    x1, y1 = nodes[n1]
                    x2, y2 = nodes[n2]
                    x_centroid = (x1 + x2) / 2
                    y_centroid = (y1 + y2) / 2
                    tx = tx_func(x_centroid, y_centroid)
                    ty = ty_func(x_centroid, y_centroid)
                    ax.quiver(x_centroid, y_centroid, tx, ty, color='red', scale=10*T, width=0.001)
        
        # 添加图例
        ax.legend(handles=legend_handles, labels=legend_labels, loc='upper right')
        
        # 设置图形属性
        ax.set_title("Original and Deformed Mesh with Boundary Conditions")
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.axis('equal')
        
        plots_generated = True
    
    # Plot Von Mises stress distribution
    if plot_config['von_mises_stress']:
        plt.figure(figsize=(7, 7))
        von_mises_elem = np.sqrt(stresses[:,0]**2 - stresses[:,0]*stresses[:,1] + stresses[:,1]**2 + 3*stresses[:,2]**2)
        von_mises_node = np.zeros(nodes.shape[0])
        count = np.zeros(nodes.shape[0])
        for i, elem in enumerate(elements):
            for node in elem:
                von_mises_node[node] += von_mises_elem[i]
                count[node] += 1
        count[count == 0] = 1
        von_mises_node /= count
        contour = plt.tricontourf(triang_deformed, von_mises_node, levels=20, cmap='jet')
        plt.colorbar(label='Von Mises Stress (Pa)')
        plt.title("Von Mises Stress Distribution")
        plt.xlabel("X")
        plt.ylabel("Y")
        plt.axis('equal')
        
        # 找到最大值和最小值及其位置
        max_vm_stress = np.max(von_mises_node)
        min_vm_stress = np.min(von_mises_node)
        max_node_stress = np.argmax(von_mises_node)
        min_node_stress = np.argmin(von_mises_node)
        
        # 在图上标出最大值和最小值
        plt.scatter(deformed[max_node_stress,0], deformed[max_node_stress,1], color='black', marker='*', s=100)
        plt.scatter(deformed[min_node_stress,0], deformed[min_node_stress,1], color='black', marker='*', s=100)
        
        # 添加数值标签
        plt.text(deformed[max_node_stress,0]+0.01, deformed[max_node_stress,1]+0.01, f'{max_vm_stress:.2e}', color='black', fontsize=10)
        plt.text(deformed[min_node_stress,0]+0.01, deformed[min_node_stress,1]+0.01, f'{min_vm_stress:.2e}', color='black', fontsize=10)
        
        plots_generated = True
    
    # Plot Von Mises strain distribution
    if plot_config['von_mises_strain']:
        plt.figure(figsize=(7, 7))
        von_mises_strain_node = np.zeros(nodes.shape[0])
        count = np.zeros(nodes.shape[0])
        for i, elem in enumerate(elements):
            for node in elem:
                von_mises_strain_node[node] += von_mises_strain[i]
                count[node] += 1
        count[count == 0] = 1
        von_mises_strain_node /= count
        contour = plt.tricontourf(triang_deformed, von_mises_strain_node, levels=20, cmap='jet')
        plt.colorbar(label='Von Mises Strain')
        plt.title("Von Mises Strain Distribution")
        plt.xlabel("X")
        plt.ylabel("Y")
        plt.axis('equal')
        
        # 找到最大值和最小值及其位置
        max_vm_strain = np.max(von_mises_strain_node)
        min_vm_strain = np.min(von_mises_strain_node)
        max_node_strain = np.argmax(von_mises_strain_node)
        min_node_strain = np.argmin(von_mises_strain_node)
        
        # 在图上标出最大值和最小值
        plt.scatter(deformed[max_node_strain,0], deformed[max_node_strain,1], color='black', marker='*', s=100)
        plt.scatter(deformed[min_node_strain,0], deformed[min_node_strain,1], color='black', marker='*', s=100)
        
        # 添加数值标签
        plt.text(deformed[max_node_strain,0]+0.01, deformed[max_node_strain,1]+0.01, f'{max_vm_strain:.2e}', color='black', fontsize=10)
        plt.text(deformed[min_node_strain,0]+0.01, deformed[min_node_strain,1]+0.01, f'{min_vm_strain:.2e}', color='black', fontsize=10)
        
        plots_generated = True
    
    # Plot exact Von Mises stress distribution
    if plot_config['exact_von_mises_stress']:
        plt.figure(figsize=(7, 7))
        exact_stresses = np.zeros((nodes.shape[0], 3))
        for i in range(nodes.shape[0]):
            x, y = nodes[i]
            r = np.sqrt(x**2 + y**2)
            theta = np.arctan2(y, x)
            if r == 0:
                r = 1e-200  # 避免除以零
            sigma_rr = (T / 2) * (1 - (R**2 / r**2)) + (T / 2) * (1 - 4 * (R**2 / r**2) + 3 * (R**4 / r**4)) * np.cos(2 * theta)
            sigma_tt = (T / 2) * (1 + (R**2 / r**2)) - (T / 2) * (1 + 3 * (R**4 / r**4)) * np.cos(2 * theta)
            sigma_rt = -(T / 2) * (1 + 2 * (R**2 / r**2) - 3 * (R**4 / r**4)) * np.sin(2 * theta)
            exact_stresses[i, 0] = sigma_rr * np.cos(theta)**2 + sigma_tt * np.sin(theta)**2 - 2 * sigma_rt * np.sin(theta) * np.cos(theta)
            exact_stresses[i, 1] = sigma_rr * np.sin(theta)**2 + sigma_tt * np.cos(theta)**2 + 2 * sigma_rt * np.sin(theta) * np.cos(theta)
            exact_stresses[i, 2] = (sigma_rr - sigma_tt) * np.sin(theta) * np.cos(theta) + sigma_rt * (np.cos(theta)**2 - np.sin(theta)**2)
        
        von_mises_exact = np.sqrt(exact_stresses[:,0]**2 - exact_stresses[:,0]*exact_stresses[:,1] + exact_stresses[:,1]**2 + 3*exact_stresses[:,2]**2)
        contour = plt.tricontourf(triang_deformed, von_mises_exact, levels=20, cmap='jet')
        plt.colorbar(label='Von Mises Stress (Pa)')
        plt.title("Exact Von Mises Stress Distribution")
        plt.xlabel("X")
        plt.ylabel("Y")
        plt.axis('equal')
        
        # 找到最大值和最小值及其位置
        max_vm_exact = np.max(von_mises_exact)
        min_vm_exact = np.min(von_mises_exact)
        max_node_exact = np.argmax(von_mises_exact)
        min_node_exact = np.argmin(von_mises_exact)
        
        # 在图上标出最大值和最小值
        plt.scatter(deformed[max_node_exact,0], deformed[max_node_exact,1], color='black', marker='*', s=100)
        plt.scatter(deformed[min_node_exact,0], deformed[min_node_exact,1], color='black', marker='*', s=100)
        
        # 添加数值标签
        plt.text(deformed[max_node_exact,0]+0.01, deformed[max_node_exact,1]+0.01, f'{max_vm_exact:.2e}', color='black', fontsize=10)
        plt.text(deformed[min_node_exact,0]+0.01, deformed[min_node_exact,1]+0.01, f'{min_vm_exact:.2e}', color='black', fontsize=10)
        
        plots_generated = True
    
    # 如果有任何一个图被绘制，则显示所有figure
    if plots_generated:
        plt.show()

def compute_errors(nodes, elements, U, exact_solution):
    num_nodes = nodes.shape[0]
    U_exact = np.zeros(2 * num_nodes)
    grad_U_exact = np.zeros((num_nodes, 2, 2))  # 存储精确解的梯度

    # 高斯积分点和权重
    gauss_points = np.array([[1/6, 1/6], [2/3, 1/6], [1/6, 2/3]])
    gauss_weights = np.array([1/3, 1/3, 1/3])

    # 初始化误差积分
    L2_error_numerator = 0.0
    L2_error_denominator = 0.0
    H1_error_numerator = 0.0
    H1_error_denominator = 0.0

    for elem in elements:
        node_indices = elem
        x = nodes[node_indices, 0]
        y = nodes[node_indices, 1]

        # 计算三角形面积
        area_matrix = np.array([
            [1, x[0], y[0]],
            [1, x[1], y[1]],
            [1, x[2], y[2]]
        ])
        area = 0.5 * np.abs(np.linalg.det(area_matrix))

        # 单元位移向量
        Ue = np.zeros(6)
        for i, ni in enumerate(node_indices):
            Ue[2 * i] = U[2 * ni]
            Ue[2 * i + 1] = U[2 * ni + 1]

        # 计算形状函数梯度
        beta = np.array([y[1] - y[2], y[2] - y[0], y[0] - y[1]])
        gamma = np.array([x[2] - x[1], x[0] - x[2], x[1] - x[0]])
        B = np.zeros((3, 6))
        B[0, 0::2] = beta
        B[1, 1::2] = gamma
        B[2, 0::2] = gamma
        B[2, 1::2] = beta
        B /= (2 * area)

        # 计算单元应变
        grad_Ue = B @ Ue  # 单元应变 (3,)

        # 高斯积分
        for point, weight in zip(gauss_points, gauss_weights):
            # 计算全局坐标
            xi, eta = point
            N1 = 1 - xi - eta
            N2 = xi
            N3 = eta
            x_g = N1 * x[0] + N2 * x[1] + N3 * x[2]
            y_g = N1 * y[0] + N2 * y[1] + N3 * y[2]

            # 计算精确解
            u_exact, v_exact = exact_solution(x_g, y_g)

            # 插值有限元解
            u_h = N1 * Ue[0] + N2 * Ue[2] + N3 * Ue[4]
            v_h = N1 * Ue[1] + N2 * Ue[3] + N3 * Ue[5]

            # 计算 L2 误差积分
            L2_error_numerator += weight * area * ((u_h - u_exact)**2 + (v_h - v_exact)**2)
            L2_error_denominator += weight * area * (u_exact**2 + v_exact**2)

            # 计算精确解的梯度
            r = np.sqrt(x_g**2 + y_g**2)
            theta = np.arctan2(y_g, x_g)
            if r == 0:
                r = 1e-20  # 避免除以零

            # 计算极坐标下的位移梯度
            du_r_dr = (T * R / (2 * E)) * (-(1 + nu) / r**2 - 3 * (1 - nu) * (R**2) / r**4) * np.cos(theta)
            du_theta_dr = (T * R / (2 * E)) * (-(1 + nu) / r**2 + 3 * (1 + 3 * nu) * (R**2) / r**4) * np.sin(theta)
            du_r_dtheta = (T * R / (2 * E)) * (-(1 + nu) / r - (1 - nu) * (R**2) / r**3) * np.sin(theta)
            du_theta_dtheta = (T * R / (2 * E)) * ((1 + nu) / r - (1 + 3 * nu) * (R**2) / r**3) * np.cos(theta)

            # 转换为笛卡尔坐标下的梯度
            grad_u_exact = np.zeros((2, 2))
            grad_u_exact[0, 0] = du_r_dr * np.cos(theta)**2 - (du_r_dtheta / r - du_theta_dr) * np.sin(theta) * np.cos(theta) + du_theta_dtheta * np.sin(theta)**2 / r
            grad_u_exact[0, 1] = du_r_dr * np.sin(theta) * np.cos(theta) + (du_r_dtheta / r - du_theta_dr) * np.cos(theta)**2 - du_theta_dtheta * np.sin(theta) * np.cos(theta) / r
            grad_u_exact[1, 0] = du_r_dr * np.sin(theta) * np.cos(theta) + (du_r_dtheta / r - du_theta_dr) * np.sin(theta)**2 - du_theta_dtheta * np.sin(theta) * np.cos(theta) / r
            grad_u_exact[1, 1] = du_r_dr * np.sin(theta)**2 + (du_r_dtheta / r - du_theta_dr) * np.sin(theta) * np.cos(theta) + du_theta_dtheta * np.cos(theta)**2 / r

            # 计算 H1 误差积分
            H1_error_numerator += weight * area * (
                (u_h - u_exact)**2 + (v_h - v_exact)**2 +
                (grad_Ue[0] - grad_u_exact[0, 0])**2 +
                (grad_Ue[1] - grad_u_exact[1, 1])**2 +
                (grad_Ue[2] - (grad_u_exact[0, 1] + grad_u_exact[1, 0]) / 2)**2
            )
            H1_error_denominator += weight * area * (
                u_exact**2 + v_exact**2 +
                grad_u_exact[0, 0]**2 +
                grad_u_exact[1, 1]**2 +
                ((grad_u_exact[0, 1] + grad_u_exact[1, 0]) / 2)**2
            )

    # 计算相对误差
    L2_error = np.sqrt(L2_error_numerator / L2_error_denominator)
    H1_error = np.sqrt(H1_error_numerator / H1_error_denominator)

    return L2_error, H1_error

def generate_mesh(geo_file, mesh_file):
    """
    使用 Gmsh 从 .geo 文件生成指定版本的 .msh 文件。
    如果目标文件已存在，则覆盖它。
    """
    if not os.path.exists(geo_file):
        raise FileNotFoundError(f".geo 文件 '{geo_file}' 未找到。")
    
    # 尝试在 PATH 中查找 gmsh
    gmsh_path = shutil.which("gmsh")
    if gmsh_path is None:
        # 如果没有找到，可以指定默认路径或抛出错误
        # 根据实际安装路径修改以下路径
        gmsh_path = r"D:\Desktop\下载\gmsh-4.0.4-Windows64\gmsh-4.0.4-Windows64\gmsh.exe"
        if not os.path.exists(gmsh_path):
            raise FileNotFoundError("Gmsh 未安装或未添加到系统 PATH。请安装 Gmsh 并确保其可执行文件在 PATH 中。")
    
    # 构建 Gmsh 命令
    cmd = [
        gmsh_path,
        geo_file,
        "-2",               # 生成二维网格
        "-o", mesh_file,    # 输出文件
        "-format", "msh2"   # 指定 MshFileVersion 为 2.2
    ]
    
    try:
        print(f"正在生成网格文件 '{mesh_file}' 使用 '{geo_file}'...")
        # 运行命令并捕获输出，明确指定编码
        result = subprocess.run(
            cmd,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            encoding='utf-8',    # 指定编码为 UTF-8
            errors='replace'     # 替换无法解码的字节
        )
        print("网格文件生成成功。")
    except subprocess.CalledProcessError as e:
        print("Gmsh 生成网格时出错:")
        print(e.stderr)
        raise

def main():
    # 直接生成网格文件，覆盖原有文件
    generate_mesh(geo_file, mesh_file)
    
    # 读取网格数据
    try:
        nodes, elements, boundaries, phys_tags = read_mesh(mesh_file)
    except Exception as e:
        print(f"读取网格时出错: {e}")
        return
    num_nodes = nodes.shape[0]
    # 定义材料刚度矩阵
    if model_type == 'plane_stress':
        D = plane_stress_matrix(E, nu)
    elif model_type == 'plane_strain':
        D = plane_strain_matrix(E, nu)
    else:
        raise ValueError("model_type必须是'plane_stress'或'plane_strain'")
    # 组装系统矩阵
    try:
        K, F = assemble_system(nodes, elements, D)
    except Exception as e:
        print(f"组装系统时出错: {e}")
        return
    # 应用边界条件
    try:
        K_ff, F_f, free_dofs, fixed_dofs, U = apply_boundary_conditions(K, F, nodes, boundaries, phys_tags, boundary_conditions)
    except Exception as e:
        print(f"应用边界条件时出错: {e}")
        return
    # 求解位移
    try:
        U_f = solve_system(K_ff, F_f)
    except Exception as e:
        print(f"求解系统时出错: {e}")
        return
    # 组装完整的位移向量
    U[free_dofs] = U_f
    # 计算应变和应力
    try:
        strains, stresses, von_mises_strain = compute_strains_stresses(nodes, elements, U, D)
    except Exception as e:
        print(f"计算应变和应力时出错: {e}")
        return
    
    if compute_errors_flag:
        try:
            L2_err, H1_err = compute_errors(nodes, elements, U, exact_solution)
            print(f"Relative L2 Error: {L2_err:.6e}")
            print(f"Relative H1 Error: {H1_err:.6e}")
        except Exception as e:
            print(f"计算误差时出错: {e}")
    else:
        print("跳过误差计算。")
        
    # 可视化结果
    try:
        visualize_results(nodes, elements, U, strains, stresses, von_mises_strain, boundaries, phys_tags, boundary_conditions, deformation_scale)
    except Exception as e:
        print(f"可视化结果时出错: {e}")
        return
    # 根据开关决定是否计算误差
    

if __name__ == "__main__":
    main()