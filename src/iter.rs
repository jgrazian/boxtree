use crate::bounds::Bounds;
use crate::bvh::{Bvh, BvhNode, BvhObjKey};
use crate::traits::*;

pub struct BvhLeafIterator<'a, T, F, const D: usize, const N: usize>
where
    T: Bounded<D>,
    F: Fn(&Bounds<D>) -> bool,
{
    pub(crate) bvh: &'a Bvh<T, D, N>,
    pub(crate) predicate: F,
    pub(crate) stack: Vec<&'a BvhNode<D, N>>,
}

impl<'a, T, F, const D: usize, const N: usize> Iterator for BvhLeafIterator<'a, T, F, D, N>
where
    T: Bounded<D>,
    F: Fn(&Bounds<D>) -> bool,
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
