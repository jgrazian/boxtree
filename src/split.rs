use crate::traits::{Bounded, BoundingBox};

pub enum SplitMethod {
    SurfaceAreaHeuristic,
}

impl SplitMethod {
    #[inline]
    pub fn split<'a, B: BoundingBox, T: Bounded<B>>(
        &self,
        bounds: &B,
        objects: &[T],
        centroids: &'a [B::Vector],
        indexes: &'a mut [usize],
    ) -> (&'a mut [usize], &'a mut [usize]) {
        match self {
            Self::SurfaceAreaHeuristic => sah_split(bounds, objects, centroids, indexes),
        }
    }
}

#[inline]
fn sah_split<'a, B: BoundingBox, T: Bounded<B>>(
    bounds: &B,
    objects: &[T],
    centroids: &'a [B::Vector],
    indexes: &'a mut [usize],
) -> (&'a mut [usize], &'a mut [usize]) {
    const MAX_DIM: usize = 3;
    const NUM_BUCKETS: usize = 16;
    let mut cuts = [[0.0; NUM_BUCKETS]; MAX_DIM];
    let mut split_idx = [[0 as usize; NUM_BUCKETS]; MAX_DIM];
    let mut scores = [[f32::INFINITY; NUM_BUCKETS]; MAX_DIM];
    let bounds_sa = bounds.surface_area();

    for axis in 0..B::DIM {
        // Sort objects by axis
        indexes.sort_unstable_by(|a, b| {
            let centroid1 = centroids[*a][axis];
            let centroid2 = centroids[*b][axis];
            centroid1.partial_cmp(&centroid2).unwrap()
        });

        for bucket in 0..NUM_BUCKETS {
            cuts[axis][bucket] = bounds.min()[axis]
                + bounds.axis_length(axis) * ((bucket + 1) as f32 / (NUM_BUCKETS + 1) as f32);

            split_idx[axis][bucket] =
                indexes.partition_point(|o| centroids[*o][axis] < cuts[axis][bucket]);

            let (left, right) = indexes.split_at(split_idx[axis][bucket]);

            // Bad cut location one of the sides has no shapes.
            if left.len() == 0 || right.len() == 0 {
                continue;
            }

            let mut left_bounds = bounds.clone();
            left_bounds.max_mut()[axis] = objects[left[left.len() - 1]].bounds().max()[axis];
            let mut right_bounds = bounds.clone();
            right_bounds.min_mut()[axis] = objects[right[0]].bounds().min()[axis];

            let left_sa = left_bounds.surface_area();
            let right_sa = right_bounds.surface_area();

            let left_ratio = left_sa / bounds_sa;
            let right_ratio = right_sa / bounds_sa;

            scores[axis][bucket] =
                1.0 + left_ratio * left.len() as f32 + right_ratio * right.len() as f32;
        }
    }

    // Select minimum score
    let mut min_score = f32::INFINITY;
    let mut min_axis = 0;
    let mut min_bucket = 0;
    for axis in 0..B::DIM {
        for bucket in 0..NUM_BUCKETS {
            if scores[axis][bucket] < min_score {
                min_score = scores[axis][bucket];
                min_axis = axis;
                min_bucket = bucket;
            }
        }
    }

    // Decide where to split
    let split_index = if min_score == f32::INFINITY {
        // Centroids are all very close just split in half
        indexes.len() / 2
    } else {
        indexes.sort_unstable_by(|a, b| {
            let centroid1 = objects[*a].bounds().centroid()[min_axis];
            let centroid2 = objects[*b].bounds().centroid()[min_axis];
            centroid1.partial_cmp(&centroid2).unwrap()
        });

        split_idx[min_axis][min_bucket]
    };

    let (left, right) = indexes.split_at_mut(split_index);
    (left, right)
}
