extern crate rustc_serialize;

use self::rustc_serialize::{Decoder,Decodable,Encoder,Encodable,json};
use rigid::{
    Space,
    Real,
    Function,
    C2
};
use std::fmt::Debug;

pub fn read_space<D:Decoder,T:Space>(d:&mut D) -> Result<T, D::Error> {
    let basis: Vec<T> = T::basis().collect();
    let l: Result<Vec<Real>, D::Error> = Decodable::decode(d);
    match l {
        Ok(v) => if v.len() == basis.len() {
            let mut out = basis[0].clone()*0.0;
            for i in 0..basis.len() {
                out = out + basis[i].clone()*v[i];
            }
            Ok(out)
        } else {
            Err(d.error("length mismatch"))
        },
        Err(e) => Err(e)
    }
}

pub fn read_pair<D:Decoder,A:Space,B:Space>(d:&mut D) -> Result<(A,B), D::Error> {
    d.read_seq(|d,n| {
        if n != 2 {
            Err(d.error("not a pair"))
        } else {
            let a = d.read_seq_elt(0,read_space);
            let b = d.read_seq_elt(1,read_space);
            match (a,b) {
                (Ok(x),Ok(y)) => Ok((x,y)),
                (Ok(x),Err(y)) => Err(y),
                (Err(x),_) => Err(x)
            }
        }
    })
}

pub fn read_list<D:Decoder,A:Space,B:Space>(d:&mut D) -> Result<Vec<(A,B)>, D::Error> where D::Error:Debug{
    d.read_seq(|d,size|{
        let mut out = Vec::new();
        for i in 0..size {
            let res = d.read_seq_elt(i,|d|{
                read_pair(d)
            });
            match(res) {
                Ok(val) => out.push(val),
                Err(err) => println!("{:?}",err)
            }
        }
        Ok(out)
    })
}

