use crate::bounds::Bounds;
use crate::bvh::{Bvh, BvhIndex, BvhNode};
use crate::traits::*;

pub struct BvhLeafIterator<'a, T, F, const D: usize>
where
    T: Bounded<D>,
    F: Fn(&Bounds<D>) -> bool,
{
    pub(crate) bvh: &'a Bvh<T, D>,
    pub(crate) predicate: F,
    pub(crate) stack: Vec<&'a BvhNode<D>>,
}

impl<'a, T, F, const D: usize> Iterator for BvhLeafIterator<'a, T, F, D>
where
    T: Bounded<D>,
    F: Fn(&Bounds<D>) -> bool,
{
    type Item = (BvhIndex, &'a T);

    fn next(&mut self) -> Option<Self::Item> {
        let node = match self.stack.pop() {
            Some(node) => node,
            None => return None,
        };

        match node {
            BvhNode::Leaf { data } => Some((*data, &self.bvh[data])),
            BvhNode::Node { bounds, children } => {
                if (self.predicate)(bounds) {
                    self.stack.push(&children[1]);
                    self.stack.push(&children[0]);
                }
                self.next()
            }
        }
    }
}
