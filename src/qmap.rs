use std::marker::PhantomData;
use std::vec::IntoIter;
use rigid::{
    Space,
    Real,
    Function
};
use std::ops::{
    Add,
    Mul
};

fn fact(n:usize) -> usize {
    if n == 0 {
        1
    } else {
        fact(n-1)*n
    }
}

#[derive(Clone,Debug)]
pub struct PolyBase<D:Space,R:Space> {
    values: Vec<Vec<Vec<Real>>>,
    domain: PhantomData<D>,
    range: PhantomData<R>
}

impl <D:Space,R:Space> PolyBase<D,R> {
    pub fn constant_facts(rank: usize,index: usize,base:usize) -> Vec<usize> {
        if rank == 0 {
            vec![]
        }else if rank == 1 {
            vec![index]
        } else {
            let mut rem = index as i32;
            let mut res = Vec::new();
            for i in 0..base {
                let sub = Self::rank_len(rank-1,i+1) as i32;
                if rem - sub < 0 {
                    res.push(i);
                    res.extend(Self::constant_facts(rank-1,rem as usize,base).iter().cloned());
                    break;
                }
                rem = rem - sub;
            }
            res
        }
    }
    pub fn rank_len(rank:usize,base:usize) -> usize {
        if rank == 0 {
            1
        }else if rank == 1 {
            base
        } else {
            let mut sum = 0;
            for i in 0..base {
                sum = sum + Self::rank_len(rank-1,i+1);
            }
            sum
/*            let mut p = 1;
            for i in 0..rank {
                p = p*base
            }
            let v = fact(rank);
            (p + Self::rank_len(rank-1,base)*(v-1))/v*/
        }
    }
    pub fn basis(rank:usize) -> IntoIter<Self> {
        let mut template = Vec::new();
        let d_basis = D::basis().count();
        println!("{}",d_basis);
        let r_basis = R::basis().count();
        for i in 0..rank {
            let mut file = Vec::new();
            let mut thing = Vec::new();
            for j in 0..Self::rank_len(i,d_basis) {
                thing.push(0.0);
            }
            for j in 0..r_basis {
                file.push(thing.clone());
            }
            template.push(file);
        }
        let mut out = Vec::new();
        for i in 0..template.len() {
            for j in 0..template[i].len() {
                for k in 0..template[i][j].len() {
                    let mut thing = template.clone();
                    thing[i][j][k] = 1.0;
                    out.push(PolyBase {
                        values: thing,
                        range: PhantomData,
                        domain: PhantomData
                    })
                }
            }
        }
        out.into_iter()
    }
    pub fn coords(&self) -> IntoIter<Real> {
        let mut out = Vec::new();
        for rank in &self.values {
            for row in rank {
                for &x in row {
                    out.push(x);
                }
            }
        }
        out.into_iter()
    }
}

impl <D:Space,R:Space> Function for PolyBase<D,R> {
    type Domain=D;
    type Range=R;
    fn eval(&self, value:D) -> R {
        let d_basis: Vec<D> = D::basis().collect();
        let r_basis: Vec<R> = R::basis().collect();
        let coords: Vec<Real> = value.coords().collect();
        let mut output = r_basis[0].clone()*0.0;
        for i in 0..self.values.len() {
            for j in 0..self.values[i].len() {
                for k in 0..self.values[i][j].len() {
                    let mut mul = 1.0;
                    for x in Self::constant_facts(i,k,d_basis.len()) {
                        mul = mul * coords[x];
                    }
                    output = output + r_basis[j].clone()*(mul*self.values[i][j][k]);
                }
            }
        }
        output 
    }
}

impl <D:Space,R:Space> Mul<Real> for PolyBase<D,R> {
    type Output = Self;
    fn mul(self, other:Real) -> Self {
        PolyBase {
            values: self.values.into_iter().map(|rank|{
                rank.into_iter().map(|row|{
                    row.into_iter().map(|val| val * other).collect()
                }).collect()
            }).collect(),
            domain: PhantomData,
            range: PhantomData
        }
    }
}

impl <D:Space,R:Space> Add<PolyBase<D,R>> for PolyBase<D,R> {
    type Output = Self;
    fn add(self, other:Self) -> Self {
        PolyBase {
            values: self.values.into_iter().zip(other.values.into_iter()).map(|(rank_a,rank_b)|{
                rank_a.into_iter().zip(rank_b.into_iter()).map(|(row_a,row_b)|{
                    row_a.into_iter().zip(row_b.into_iter()).map(|(a,b)| a+b).collect()
                }).collect()
            }).collect(),
            domain: PhantomData,
            range: PhantomData
        }
    }
}

pub type QuinticMap<D:Space,R:Space> = PolyBase<D,R>;

impl <D:Space,R:Space> QuinticMap<D,R> {
    fn thing()->Self{
        let b:Vec<PolyBase<D,R>> = PolyBase::basis(5).collect();
        b[0].clone() as Self
    }
}

impl <D:Space,R:Space> Space for QuinticMap<D,R> {
    type BIter=IntoIter<Self>;
    type CIter=IntoIter<Real>;
    fn basis() -> IntoIter<Self> {
        PolyBase::<D,R>::basis(5) as IntoIter<Self>
    }
    fn coords(&self) -> IntoIter<Real> {
        (self as &PolyBase<D,R>).coords()
    }
}
