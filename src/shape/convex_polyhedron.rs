

use nalgebra::{DMatrix, RealField, SVector};
use std::collections::HashMap;

use crate::shape::generalized_cross;

// ────────────────────────────────────────────────────────────
//  面结构
// ────────────────────────────────────────────────────────────

pub struct Face<Num: RealField+Copy, const DIM: usize> {
    /// DIM 个顶点索引（逆时针/一致定向），法线朝外。
    pub vertices: [usize; DIM],
    /// 单位外法线
    pub normal: SVector<Num, DIM>,
    /// 平面偏移：normal · x = offset
    pub offset: Num,
    
}

impl<Num: RealField+Copy, const DIM: usize> Face<Num, DIM> {
    #[inline]
    pub fn signed_distance(&self, p: &SVector<Num, DIM>) -> Num {
        self.normal.dot(p) - self.offset
    }
}

// ────────────────────────────────────────────────────────────
//  凸多面体
// ────────────────────────────────────────────────────────────

pub struct ConvexPolyhedron<Num: RealField+Copy, const DIM: usize> {
    /// 凸包上的顶点（去重后保留原始坐标）
    pub points: Vec<SVector<Num, DIM>>,
    /// 所有面，每面有 DIM 个顶点
    pub faces: Vec<Face<Num, DIM>>,
}

impl<Num: RealField+Copy, const DIM: usize> ConvexPolyhedron<Num, DIM> {
    pub fn from_convex_hull(input: &[SVector<Num, DIM>]) -> Option<Self> {
        assert!(DIM >= 2, "DIM must be >= 2");
        let eps = Num::default_epsilon().sqrt();

        // ① 去重
        let points = dedup(input, eps);
        if points.len() < DIM + 1 {
            return None;
        }

        // ② 找初始单纯形（DIM+1 个仿射独立的点）
        let simplex = find_initial_simplex(&points, eps)?;

        // ③ 内部参考点 —— 单纯形重心，严格位于单纯形内部，因此严格位于凸包内部。
        let mut centroid = SVector::<Num, DIM>::zeros();
        for &i in &simplex {
            centroid += points[i];
        }
        centroid /= Num::from_usize(DIM + 1).unwrap();

        // ④ 初始单纯形的 DIM+1 个面（每个面 = 去掉一个单纯形顶点）
        let mut faces: Vec<Face<Num, DIM>> = Vec::new();
        for skip in 0..=DIM {
            let mut verts = [0usize; DIM];
            let mut idx = 0;
            for k in 0..=DIM {
                if k != skip {
                    verts[idx] = simplex[k];
                    idx += 1;
                }
            }
            if let Some(face) = build_face(&points, &verts, &centroid, eps) {
                faces.push(face);
            }
        }

        // ⑤ 增量添加剩余点
        let mut in_hull = vec![false; points.len()];
        for &i in &simplex {
            in_hull[i] = true;
        }

        for i in 0..points.len() {
            if in_hull[i] {
                continue;
            }
            let p = &points[i];

            // 找出所有"可见面"：p 在其外法线正侧
            let visible: Vec<usize> = faces
                .iter()
                .enumerate()
                .filter(|(_, f)| f.signed_distance(p) > eps)
                .map(|(idx, _)| idx)
                .collect();

            if visible.is_empty() {
                // p 在凸包内部
                continue;
            }

            // 收集 horizon ridges：ridge 是每面去掉一个顶点得到的 DIM-1 元组
            // 出现两次 → 内部 ridge（两相邻可见面共享）；出现一次 → 边界 ridge
            let mut ridge_count: HashMap<Vec<usize>, (usize, Vec<usize>)> = HashMap::new();
            for &fi in &visible {
                let f = &faces[fi];
                for skip in 0..DIM {
                    let ordered: Vec<usize> = (0..DIM)
                        .filter(|&k| k != skip)
                        .map(|k| f.vertices[k])
                        .collect();
                    let mut sorted = ordered.clone();
                    sorted.sort();
                    ridge_count
                        .entry(sorted)
                        .and_modify(|e| e.0 += 1)
                        .or_insert((1, ordered));
                }
            }

            // 删除可见面（从大到小，swap_remove 安全）
            let mut vis_sorted = visible.clone();
            vis_sorted.sort();
            for &fi in vis_sorted.iter().rev() {
                faces.swap_remove(fi);
            }

            // 每条 horizon ridge 加 p → 一个新面
            for (_, (count, ordered)) in ridge_count {
                if count % 2 != 1 {
                    continue; // 内部 ridge
                }
                let mut new_verts = [0usize; DIM];
                for (k, &v) in ordered.iter().enumerate() {
                    new_verts[k] = v;
                }
                new_verts[DIM - 1] = i; // 追加新点
                if let Some(face) = build_face(&points, &new_verts, &centroid, eps) {
                    faces.push(face);
                }
            }

            in_hull[i] = true;
        }

        Some(ConvexPolyhedron { points, faces })
    }
}

// ────────────────────────────────────────────────────────────
//  辅助：去重
// ────────────────────────────────────────────────────────────

fn dedup<Num: RealField+Copy, const DIM: usize>(
    pts: &[SVector<Num, DIM>],
    eps: Num,
) -> Vec<SVector<Num, DIM>> {
    let eps2 = eps * eps;
    let mut out: Vec<SVector<Num, DIM>> = Vec::with_capacity(pts.len());
    for &p in pts {
        if !out.iter().any(|q| (*q - p).norm_squared() < eps2) {
            out.push(p);
        }
    }
    out
}

// ────────────────────────────────────────────────────────────
//  辅助：找 DIM+1 个仿射独立的初始单纯形
//        通过 Gram–Schmidt 迭代正交化，每一轮找"距离当前子空间最远"的点
// ────────────────────────────────────────────────────────────

fn find_initial_simplex<Num: RealField+Copy, const DIM: usize>(
    points: &[SVector<Num, DIM>],
    eps: Num,
) -> Option<Vec<usize>> {
    let n = points.len();
    let origin = points[0];
    let mut simplex: Vec<usize> = vec![0];
    let mut ortho_basis: Vec<SVector<Num, DIM>> = Vec::new();

    for _k in 0..DIM {
        let mut best_dist = -Num::one();
        let mut best_idx: Option<usize> = None;
        let mut best_ortho: Option<SVector<Num, DIM>> = None;

        for i in 0..n {
            if simplex.contains(&i) {
                continue;
            }
            let v = points[i] - origin;
            let mut w = v;
            for b in &ortho_basis {
                let bb = b.dot(b);
                if bb > eps * eps {
                    w -= *b * (w.dot(b) / bb);
                }
            }
            let dist = w.norm();
            if dist > best_dist {
                best_dist = dist;
                best_idx = Some(i);
                best_ortho = Some(w);
            }
        }

        if best_dist < eps {
            return None; // 退化（所有点落在低维仿射子空间上）
        }
        simplex.push(best_idx.unwrap());
        ortho_basis.push(best_ortho.unwrap());
    }
    Some(simplex)
}

// // ────────────────────────────────────────────────────────────
// //  辅助：广义叉积 —— 在 DIM 维中，给定 DIM-1 个向量，返回垂直于它们的向量
// //  用余子式展开：n_j = (-1)^j · det(M 去掉第 j 列)
// // ────────────────────────────────────────────────────────────

// fn generalized_cross<Num: RealField, const DIM: usize>(
//     vectors: &[SVector<Num, DIM>],
// ) -> SVector<Num, DIM> {
//     debug_assert_eq!(vectors.len(), DIM - 1);
//     let n = DIM - 1;
//     let mut result = SVector::<Num, DIM>::zeros();
//     for j in 0..DIM {
//         let mut mat = DMatrix::<Num>::zeros(n, n);
//         for i in 0..n {
//             let mut col = 0;
//             for c in 0..DIM {
//                 if c != j {
//                     mat[(i, col)] = vectors[i][c];
//                     col += 1;
//                 }
//             }
//         }
//         let det = mat.determinant();
//         result[j] = if j % 2 == 0 { det } else { -det };
//     }
//     result
// }

// ────────────────────────────────────────────────────────────
//  辅助：构建一个面（法线单位化 + 重心校正朝外 + 顶点顺序一致）
// ────────────────────────────────────────────────────────────

fn build_face<Num: RealField+Copy, const DIM: usize>(
    points: &[SVector<Num, DIM>],
    verts_in: &[usize; DIM],
    centroid: &SVector<Num, DIM>,
    eps: Num,
) -> Option<Face<Num, DIM>> {
    let mut verts = *verts_in;
    let origin = points[verts[0]];

    // DIM-1 个边向量（都从 verts[0] 出发）
    let edge_vectors: Vec<SVector<Num, DIM>> =
        (1..DIM).map(|k| points[verts[k]] - origin).collect();

    let mut normal = generalized_cross(edge_vectors);
    let norm = normal.norm();
    if norm < eps {
        return None; // 退化面
    }
    normal /= norm;
    let mut offset = normal.dot(&origin);

    // 重心必须在负侧（内部侧）
    if normal.dot(centroid) - offset > Num::zero() {
        normal = -normal;
        offset = -offset;
        verts.swap(0, 1); // 翻转前两个顶点，保持顶点顺序与法线一致
    }

    Some(Face {
        vertices: verts,
        normal,
        offset,
    })
}