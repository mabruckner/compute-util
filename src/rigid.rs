use std::vec::IntoIter;
use std::iter::Chain;
use std::ops::{Add,Mul};
use std::marker::PhantomData;

pub type Real = f64;

pub trait Space : Add<Self,Output=Self> + Mul<Real,Output=Self> + Clone {
    type BIter : Iterator<Item=Self>;
    type CIter : Iterator<Item=Real>;
    fn basis() -> Self::BIter;
    fn coords(&self) -> Self::CIter;
}

pub trait Function {
    type Domain : Space;
    type Range : Space;
    fn eval(&self, Self::Domain) -> Self::Range;
}

pub trait Differentiable: Function {
    fn partial(&self, <Self as Function>::Domain, <Self as Function>::Domain) -> Self::Range;
}

#[derive(Clone,Copy,Debug)]
pub struct C2<A:Space,B:Space> {
    pub a: A,
    pub b: B,
}

impl <A:Space,B:Space> Add<C2<A,B>> for C2<A,B> {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        C2 {
            a: self.a + rhs.a,
            b: self.b + rhs.b
        }
    }
}

impl <A:Space,B:Space> Mul<Real> for C2<A,B> {
    type Output = Self;
    fn mul(self, rhs: Real) -> Self {
        C2 {
            a: self.a * rhs,
            b: self.b * rhs
        }
    }
}

impl <A:Space,B:Space> Space for C2<A,B> {
    type BIter = IntoIter<Self>;
    type CIter = IntoIter<Real>;
    fn basis() -> Self::BIter {
        let mut output = Vec::new();
        let a_basis:Vec<A> = A::basis().collect();
        let b_basis:Vec<B> = B::basis().collect();
        let a_null = a_basis[0].clone() * 0.0;
        let b_null = b_basis[0].clone() * 0.0;
        for i in 0..a_basis.len() {
            output.push(C2{
                a:a_basis[i].clone(),
                b:b_null.clone()
            })
        }
        for i in 0..b_basis.len() {
            output.push(C2{
                a:a_null.clone(),
                b:b_basis[i].clone()
            })
        }
        output.into_iter()
    }
    fn coords(&self) -> Self::CIter {
        let a_basis:Vec<A> = A::basis().collect();
        let b_basis:Vec<B> = B::basis().collect();
        let a_null = a_basis[0].clone() * 0.0;
        let b_null = b_basis[0].clone() * 0.0;
        let out:Vec<Real> = self.a.coords().chain(self.b.coords()).collect();
        out.into_iter()
    }
}

/*pub trait WorldState : Space {
    type Input : Space;
    fn direction(&self, &Self::Input) -> Self;
}*/

impl Space for f64 {
    type BIter = IntoIter<Self>;
    type CIter = IntoIter<Real>;
    fn basis() -> IntoIter<Self> {
        vec![1.0].into_iter()
    }
    fn coords(&self) -> IntoIter<Real> {
        vec![self.clone()].into_iter()
    }
}

#[derive(Clone,Debug)]
pub struct LinearMap<D:Space,R:Space> {
    scales: Vec<Vec<Real>>,
    offsets: Vec<Real>,
    domain: PhantomData<D>,
    range: PhantomData<R>
}

impl <D:Space,R:Space> Function for LinearMap<D,R> {
    type Domain=D;
    type Range=R;
    fn eval(&self, value:D) -> R {
        let d_basis: Vec<D> = D::basis().collect();
        let r_basis: Vec<R> = R::basis().collect();
        let coords: Vec<Real> = value.coords().collect();
        let mut output = r_basis[0].clone() * 0.0;
        for i in 0..r_basis.len() {
            let mut val = self.offsets[i];
            for j in 0..d_basis.len() {
                val = val + coords[j] * self.scales[i][j];
            }
            output = output + r_basis[i].clone() * val
        }
        output
    }
}

impl <D:Space,R:Space> Mul<Real> for LinearMap<D,R> {
    type Output = LinearMap<D,R>;
    fn mul(self, other:Real) -> LinearMap<D,R> {
        LinearMap {
            scales: self.scales.into_iter().map(|row|{
                row.into_iter().map(|val| val * other).collect()
            }).collect(),
            offsets: self.offsets.into_iter().map(|val| val*other).collect(),
            domain: PhantomData,
            range: PhantomData
        }
    }
}

impl <D:Space,R:Space> Add<LinearMap<D,R>> for LinearMap<D,R> {
    type Output = LinearMap<D,R>;
    fn add(self, other:LinearMap<D,R>) -> LinearMap<D,R> {
        LinearMap {
            scales: self.scales.into_iter().zip(other.scales.into_iter()).map(|(row_a,row_b)|{
                row_a.into_iter().zip(row_b.into_iter()).map(|(a,b)| a+b).collect()
            }).collect(),
            offsets: self.offsets.into_iter().zip(other.offsets.into_iter()).map(|(a,b)| a+b).collect(),
            domain: PhantomData,
            range: PhantomData
        }
    }
}

impl <D:Space,R:Space> Space for LinearMap<D,R> {
    type BIter=IntoIter<LinearMap<D,R>>;
    type CIter=IntoIter<Real>;
    fn basis() -> IntoIter<LinearMap<D,R>> {
        let d_size = D::basis().count();
        let r_size = R::basis().count();
        let mut output = Vec::new();
        let mut template:LinearMap<D,R> = LinearMap{
            scales: Vec::new(),
            offsets: Vec::new(),
            domain: PhantomData,
            range:PhantomData
        };
        for i in 0..r_size {
            template.offsets.push(0.0);
            template.scales.push(Vec::new());
            for j in 0..d_size {
                template.scales[i].push(0.0);
            }
        }
        for i in 0..r_size {
            for j in 0..d_size {
                let mut thing = template.clone();
                thing.scales[i][j] = 1.0;
                output.push(thing);
            }
            let mut thing = template.clone();
            thing.offsets[i] = 1.0;
            output.push(thing);
        }
        output.into_iter()
    }
    fn coords(&self) -> IntoIter<Real> {
        let mut output = Vec::new();
        for thing in &self.scales {
            output.extend(thing.into_iter());
        }
        output.extend(self.offsets.clone().into_iter());
        output.into_iter()
    }
}
