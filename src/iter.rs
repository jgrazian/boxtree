use crate::bvh::{Bvh, BvhNode, BvhObjKey};
use crate::traits::*;

pub struct BvhLeafIterator<'a, B, T, F, const N: usize>
where
    B: BoundingBox,
    T: Bounded<B>,
    F: Fn(&B) -> bool,
{
    pub(crate) bvh: &'a Bvh<B, T, N>,
    pub(crate) predicate: F,
    pub(crate) stack: Vec<&'a BvhNode<B, N>>,
}

impl<'a, B, T, F, const N: usize> Iterator for BvhLeafIterator<'a, B, T, F, N>
where
    B: BoundingBox,
    T: Bounded<B>,
    F: Fn(&B) -> bool,
{
    type Item = (BvhObjKey, &'a T);

    fn next(&mut self) -> Option<Self::Item> {
        let node = match self.stack.pop() {
            Some(node) => node,
            None => return None,
        };

        match node {
            BvhNode::Leaf(data) => Some((*data, &self.bvh[data])),
            BvhNode::Node { bounds, children } => {
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
