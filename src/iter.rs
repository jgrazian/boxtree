use crate::bounds::{Bounds2, Bounds3, Bounds3A};
use crate::bvh::{Bvh2, Bvh3, Bvh3A, BvhNode2, BvhNode3, BvhNode3A, BvhObjKey};
use crate::traits::*;

pub struct BvhLeafIterator2<'a, T, F, const N: usize>
where
    T: Bounded2,
    F: Fn(&Bounds2) -> bool,
{
    pub(crate) bvh: &'a Bvh2<T, N>,
    pub(crate) predicate: F,
    pub(crate) stack: Vec<&'a BvhNode2<N>>,
}
pub struct BvhLeafIterator3<'a, T, F, const N: usize>
where
    T: Bounded3,
    F: Fn(&Bounds3) -> bool,
{
    pub(crate) bvh: &'a Bvh3<T, N>,
    pub(crate) predicate: F,
    pub(crate) stack: Vec<&'a BvhNode3<N>>,
}
pub struct BvhLeafIterator3A<'a, T, F, const N: usize>
where
    T: Bounded3A,
    F: Fn(&Bounds3A) -> bool,
{
    pub(crate) bvh: &'a Bvh3A<T, N>,
    pub(crate) predicate: F,
    pub(crate) stack: Vec<&'a BvhNode3A<N>>,
}

macro_rules! impl_leaf_iter {
    ($iter:ident, $bounds:ty, $bounded:ident, $node:ident, $bvh:ident) => {
        impl<'a, T, F, const N: usize> Iterator for $iter<'a, T, F, N>
        where
            T: $bounded,
            F: Fn(&$bounds) -> bool,
        {
            type Item = (BvhObjKey, &'a T);

            fn next(&mut self) -> Option<Self::Item> {
                let node = match self.stack.pop() {
                    Some(node) => node,
                    None => return None,
                };

                match node {
                    $node::Leaf(data) => Some((*data, &self.bvh[data])),
                    $node::Node { bounds, children } => {
                        if (self.predicate)(bounds) {
                            for child_key in children.iter().rev() {
                                match child_key {
                                    Some(key) => self.stack.push(&self.bvh.nodes[key.0]),
                                    None => (),
                                }
                            }
                        }
                        self.next()
                    }
                }
            }
        }
    };
}

impl_leaf_iter!(BvhLeafIterator2, Bounds2, Bounded2, BvhNode2, Bvh2);
impl_leaf_iter!(BvhLeafIterator3, Bounds3, Bounded3, BvhNode3, Bvh3);
impl_leaf_iter!(BvhLeafIterator3A, Bounds3A, Bounded3A, BvhNode3A, Bvh3A);
