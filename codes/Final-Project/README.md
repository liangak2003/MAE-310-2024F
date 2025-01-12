README
概述：

这是一个求解在 x-y 平面上定义的平面弹性问题的代码

功能介绍：

1. 网格生成：
使用 Gmsh 从 .geo 文件生成二维三角形网格。
支持定义物理组（如边界条件）。

2. 有限元求解：
组装全局刚度矩阵和载荷向量。
支持平面应力（plane_stress）和平面应变（plane_strain）模型。
应用 Dirichlet 和 Neumann 边界条件。
使用稀疏矩阵求解器计算节点位移。

3. 后处理与可视化：
计算单元应变.应力和 Von Mises 等效应力/应变。
支持绘制变形后的网格、Von Mises 应力分布、Von Mises 应变分布。
可选绘制精确解的 Von Mises 应力分布。

4. 误差分析：
计算数值解与精确解之间的相对 L2 误差和 H1 误差。

前期准备：
1. 安装 Python 环境：
确保你已经安装了 Python 3.6 或更高版本。你可以从 Python 官方网站 下载并安装 Python。

2. 安装所需的 Python 库：
代码中使用了多个 Python 库，包括 numpy, meshio, scipy, matplotlib 等。你可以通过以下命令安装这些库：
pip install numpy meshio scipy matplotlib

3. 装 Gmsh：
Gmsh 是一个用于生成网格的开源工具。你需要安装 Gmsh 并确保它可以在命令行中调用。
你可以从 Gmsh 官方网站 下载并安装 Gmsh。
安装完成后，确保 Gmsh 的可执行文件路径已添加到系统的 PATH 环境变量中。
如果 Gmsh 没有自动添加到 PATH，你可以手动添加路径。例如，在 Windows 上，你可以将 Gmsh 的安装路径（如 C:\Program Files\Gmsh）添加到 PATH 中。

4. 准备 .geo 文件：
代码中使用了 .geo 文件来定义几何模型。你需要准备一个 .geo 文件，定义你要分析的几何形状和边界条件。
确保它与 Python 脚本在同一目录下。
这里有提供一个geometry.geo文件。

配置代码参数：
在代码的 配置参数部分，根据你的问题调整以下参数：

Gmsh 配置
geo_file：.geo 文件的路径（默认："geometry.geo"）。
mesh_file：生成的 .msh 文件的路径（默认："geometry.msh"）。

模型类型
model_type：选择模型类型，支持 'plane_stress'（平面应力）或 'plane_strain'（平面应变）。

材料属性
E：杨氏模量（默认：1e9 Pa）。
nu：泊松比（默认：0.3）。

几何参数
R：圆孔半径 0.5（默认）。

载荷条件
T：远场拉伸应力1e4 Pa（默认）。

变形放大倍数
deformation_scale：变形放大倍数，用于可视化（默认：5400）。

绘图配置
plot_config：控制绘图内容的字典，包括：
deformed_mesh：是否绘制变形后的网格（默认：True）。
von_mises_stress：是否绘制 Von Mises 应力分布（默认：True）。
von_mises_strain：是否绘制 Von Mises 应变分布（默认：True）。
exact_von_mises_stress：是否绘制精确解的 Von Mises 应力分布（默认：False）。

误差计算
compute_errors_flag：是否计算数值解与精确解之间的误差（默认：True）。

边界条件
boundary_conditions：定义边界条件的字典，包括：
默认：
SymX：X 对称边界（固定 Y 方向位移）。

SymY：Y 对称边界（固定 X 方向位移）。

OuterX：X 方向远场边界（施加远场应力）。

OuterY：Y 方向远场边界（无载荷）。

Hole：圆孔边界（无载荷）。

字典的格式：
boundary_conditions = {
    "物理组名称": {
        "type": "边界条件类型",  # Dirichlet 或 Neumann
        "values": {
            "自由度": 函数或值,  # 例如 ux, uy, tx, ty
            ...
        }
    },
    ...
}

说明：
1. 物理组名称
对应 .geo 文件中定义的物理组名称（如 SymX, OuterX, Hole 等），必须与 .geo 文件中的物理组名称一致。

2. 边界条件类型
Dirichlet：固定位移边界条件。
例如：固定某个方向的位移。
Neumann：施加力或应力边界条件。
例如：施加远场应力或分布力。

3. 自由度
Dirichlet 条件：
ux：固定 X 方向位移。
uy：固定 Y 方向位移。
Neumann 条件：
tx：施加 X 方向的面力（应力）。
ty：施加 Y 方向的面力（应力）。

4. 值
可以是一个常量值，也可以是一个函数（如 lambda 函数）。
如果是函数，函数的输入是节点的坐标 (x, y)，输出是对应的值。

运行代码
将这些参数配置好以后，便可以运行代码了。

附：
关于带圆孔平板，解析解的边界条件设置
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
            'tx': lambda x, y: T * (1 - (R**2 * (4 * x**2 - 2 * y**2)) / (2 * (x**2 + y**2)**2) + (3 * R**4 * (x**2 - y**2)) / (2 * (x**2 + y**2)**3)),
            'ty': lambda x, y: -T * ((3 * R**2 * x * y) / (x**2 + y**2)**2 - (3 * R**4 * x * y) / (x**2 + y**2)**3)
        }
    },
    'OuterY': {
        'type': 'Neumann',
        'values': {
            'tx': lambda x, y: -T * ((3 * R**2 * x * y) / (x**2 + y**2)**2 - (3 * R**4 * x * y) / (x**2 + y**2)**3),
            'ty': lambda x, y: T * ((R**2 * (2 * x**2 - 4 * y**2)) / (2 * (x**2 + y**2)**2) - (3 * R**4 * (x**2 - y**2)) / (2 * (x**2 + y**2)**3))
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





