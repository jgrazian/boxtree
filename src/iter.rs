use crate::bvh::{self, Bvh, BvhNode, BvhObjKey};
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
    obj_key_stack: Vec<BvhObjKey>,
}

impl<'a, B, T, F, const N: usize> BvhLeafIterator<'a, B, T, F, N>
where
    B: BoundingBox,
    T: Bounded<B>,
    F: Fn(&B) -> bool,
{
    pub fn new(bvh: &'a Bvh<B, T, N>, predicate: F) -> Self {
        let mut stack = Vec::with_capacity(16);
        stack.push(&bvh.nodes[0]);
        let obj_key_stack = Vec::with_capacity(N);
        Self {
            bvh,
            predicate,
            stack,
            obj_key_stack,
        }
    }
}

impl<'a, B, T, F, const N: usize> Iterator for BvhLeafIterator<'a, B, T, F, N>
where
    B: BoundingBox,
    T: Bounded<B>,
    F: Fn(&B) -> bool,
{
    type Item = (BvhObjKey, &'a T);

    fn next(&mut self) -> Option<Self::Item> {
        // We have seen a leaf so there are objects on the stack
        match self.obj_key_stack.pop() {
            Some(obj_key) => {
                let obj = &self.bvh[obj_key];
                return Some((obj_key, obj));
            }
            None => (),
        }

        // No objects so we iterate on node stack
        let node = match self.stack.pop() {
            Some(node) => node,
            None => return None,
        };

        match node {
            BvhNode::Leaf { bounds, objects } => {
                if (self.predicate)(bounds) {
                    self.obj_key_stack
                        .extend(objects.iter().rev().filter_map(|o| *o))
                }
            }
            BvhNode::Node { bounds, children } => {
                if (self.predicate)(bounds) {
                    for child_key in children.iter().rev() {
                        match child_key {
                            Some(key) => self.stack.push(&self.bvh.nodes[key.0]),
                            None => (),
                        }
                    }
                }
            }
        };
        self.next()
    }
}
