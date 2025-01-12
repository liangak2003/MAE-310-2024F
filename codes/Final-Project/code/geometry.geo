// 参数定义
R = 0.5;
L = 2;

// 定义几何点（原点移到圆心）
Point(1) = {L, 0, 0};                     // (2, 0, 0)
Point(2) = {L, L, 0};                   // (2, 2, 0)
Point(3) = {0, L, 0};                     // (0, 2, 0)
Point(4) = {0, 0, 0};                       // 圆心 (0, 0, 0)
Point(5) = {R, 0, 0};                       // (0.5, 0, 0)
Point(6) = {0, R, 0};                       // (0, 0.5, 0)
Point(7) = {Cos(Pi/4)*R, Sin(Pi/4)*R, 0};   // (≈0.3536, ≈0.3536, 0)

// 定义圆弧
Circle(1) = {5, 4, 7};
Circle(2) = {7, 4, 6};

// 定义直线
Line(3) = {6, 3};
Line(4) = {3, 2};
Line(5) = {2, 1};
Line(6) = {1, 5};

// 定义曲线环
Curve Loop(1) = {4, 5, 6, 1, 2, 3};
Plane Surface(1) = {1};

// 定义物理组
Physical Surface("Body") = {1};
Physical Curve("SymX") = {6};   // 对称边界 x=0
Physical Curve("SymY") = {3};   // 对称边界 y=0
Physical Curve("OuterX") = {5}; // 右边界
Physical Curve("OuterY") = {4}; // 上边界
Physical Curve("Hole") = {1, 2};// 圆孔边界

// 网格划分
Transfinite Line{1, 2, 3, 4, 5, 6} = 50;
Mesh.ElementOrder = 1;
Mesh.Algorithm = 8;

// EOF