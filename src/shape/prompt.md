# shape

```rust
/// 泛型凸多面体：由顶点和一组 (DIM-1) 维面构成。
/// 保证为凸多面体。
/// 每个面是 DIM 个顶点确定的单纯形，`normal` 为外向单位法线，
/// `offset` 满足 `normal · x = offset` 对该面上所有点成立。
#[derive(Debug, Clone)]
pub struct ConvexPolyhedron<Num: RealField + Copy, const DIM: usize> {
    vertexes: Vec<SVector<Num, DIM>>,
    surfaces: Vec<Surface<Num, DIM>>,

    volume:Num
}

#[derive(Debug, Clone)]
pub struct Surface<Num: RealField + Copy, const DIM: usize> {
    /// 该面的 DIM 个顶点索引。
    indices: [usize; DIM],
    /// 外向单位法线。
    normal_vector: SVector<Num, DIM>,
    /// 平面偏移：normal · x = offset。
    offset: Num,
    /// 该面当前关联的外部点集合（仅 Quickhull 内部使用）。
    outside: Vec<usize>,

    /// 面向量，方向为法向量，长度为面积
    edge_vector:SVector<Num,DIM>,
    area:Num,
    center:SVector<Num,DIM>
}
```

utils:

```rust
/// 计算一个面的“未归一化”法向量。
/// DIM=3 时就是三角形叉积；
/// DIM=2 时就是线段的垂直向量；
/// 一般 DIM 下使用余子式/广义叉积。
pub fn generalized_cross<Num: RealField, const DIM: usize>(
    vertexes: impl IntoIterator<Item = SVector<Num,DIM>>,
    // edge: &[usize; DIM],
) -> SVector<Num, DIM> {
    assert!(DIM >= 2);
    let mut vertexes=vertexes.into_iter();
    let p0 = vertexes.next().unwrap();

    // DIM 维空间中，(DIM-1)-单纯形有 DIM 个顶点，
    // 因此有 DIM-1 条边向量。
    let m = DIM - 1;
    let mut edges = Vec::with_capacity(m);

    for i in 0..m {
        edges.push(vertexes.next().unwrap() - p0.clone());
    }

    // n[k] = (-1)^k * det(去掉第 k 列的 (DIM-1)x(DIM) 边矩阵)
    let mut n = SVector::<Num, DIM>::zeros();

    for k in 0..DIM {
        let mut mat = nalgebra::DMatrix::<Num>::zeros(m, m);

        for i in 0..m {
            let mut col = 0;
            for j in 0..DIM {
                if j == k {
                    continue;
                }
                mat[(i, col)] = edges[i][j].clone();
                col += 1;
            }
        }

        let mut d = mat.determinant();
        if k % 2 == 1 {
            d = -d;
        }
        n[k] = d;
    }

    n
}

pub fn surface_center<Num: RealField, const DIM: usize>(
    vertexes: &[SVector<Num, DIM>],
    edge: &[usize; DIM],
) -> SVector<Num, DIM> 
{
    let a:SMatrix<Num,DIM,DIM>=SMatrix::from( array::from_fn::<_,DIM,_>(|i|vertexes[edge[i]].clone().into())  );
    a.column_sum().scale(Num::one()/Num::from_isize(DIM as isize).unwrap())
}

pub fn volume<'a,Num:RealField,const DIM:usize>(vertexes:&'a [SVector<Num,DIM>],edges:impl IntoIterator<Item = &'a [usize;DIM]>) -> Num 
where Const<DIM>:DimSquare
{
    let mut sum = Num::zero();

    for edges in edges {
        let mat = SMatrix::<Num, DIM, DIM>::from_fn(|r, c| {
            vertexes[edges[r]][c].clone()
        });

        sum += mat.determinant();
    }

    let fact = Num::from_usize(factorial(DIM).unwrap()).unwrap();

    (sum / fact).abs()
}
```

## 任务

为他们实现

判断一个点是否在图形内部

计算一条直线与图形相交的部分

计算两个图形的碰撞

## 要求

必须为DIM>=2的所有情况实现
