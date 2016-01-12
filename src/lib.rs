mod rigid;
mod spaceread;
mod qmap;
use rigid::{
    Space,
    Real,
    Function,
    Differentiable,
    LinearMap,
    C2
};
use qmap::{PolyBase,QuinticMap};
extern crate rustc_serialize;

use self::rustc_serialize::{Decoder,Decodable,Encoder,Encodable,json};
use std::f64;
use std::fmt::Debug;
use std::marker::PhantomData;
use spaceread::{read_list,read_space};
use std::io;
use std::vec::IntoIter;
use std::collections::vec_deque::VecDeque;
//Everything over (0,1)

#[derive(Debug)]
struct Fourier {
    values: Vec<Vec<Real>>
}
fn what_solve(mat : &Vec<Vec<Real>>, targ: &Vec<Real>) -> Vec<Real> {
    let (w,h) = (mat[0].len(),mat.len());
    let mut wat = Vec::new();
    for i in 0..w {
        let mut sum = 0.0;
        for j in 0..h {
            sum = sum + mat[j][i]*targ[j];
        }
        wat.push(sum);
    }
    let mut hermit = Vec::new();
    for y in 0..w {
        let mut row = Vec::new();
        for x in 0..w {
            let mut sum = 0.0;
            for i in 0..h {
                sum += mat[i][x] * mat[i][y];
            }
            row.push(sum);
        }
        hermit.push(row);
    }
    solve_mat(&hermit,&wat)
}
fn wide_solve(imat : &Vec<Vec<Real>>, targ: &Vec<Real>) -> Vec<Real> {
    let (w,h) = (imat[0].len(),imat.len());
    let mut mat = imat.clone();
    let mut out = targ.clone();
    let mut pairs = Vec::new();
    let mut y = 0;
    for i in 0..w {
        let (mut m,mut mi) = (0.0,0);
        for j in y..h {
            if mat[j][i] > m {
                mi = j;
                m = mat[j][i];
            }
        }
        if m == 0.0 {
            continue
        }
        mat.swap(y,mi);
        out.swap(y,mi);
        for j in i..w {
            mat[y][j] = mat[y][j]/m;
        }
        mat[y][i] = 1.0;
        out[y] = out[y]/m;
        for j in 0..h {
            if j == y {
                continue
            }
            let mul = mat[j][i];
            for k in 0..w {
                mat[j][k] = mat[j][k] - mat[y][k]*mul;
            }
            out[j] = out[j] - out[y]*mul;
            mat[j][i] = 0.0;
        }
        pairs.push((y,i));
        y = y+1;
        if y == h {
            break
        }
    }
    /*for i in 0..h {
        println!("{:?} [{}]",mat[i],out[i]);
    }
    println!("{:?}",pairs);*/
    let mut empty =Vec::new();
    for i in 0..w {
        empty.push(0.0);
    }
    for i in h..w {
        mat.push(empty.clone());
        out.push(0.0);
    }
    for &(y,x) in pairs.iter().rev() {
        mat.swap(y,x);
        out.swap(y,x);
    }
    /*for i in 0..w {
        println!("{:?} [{}]",mat[i],out[i]);
    }*/
    for i in 0..w {
        if mat[i][i] != 0.0 {
            continue
        }
        mat[i][i] = 1.0;
        for &(_,r) in &pairs {
            let m = mat[r][i];
            for j in 0..w {
                if j != r {
                    mat[i][j] = mat[i][j] + mat[r][j]*m;
                }
            }
            out[i] = out[i] + out[r]*m;
        }
    }
    /*println!("--");
    for i in 0..w {
        println!("{:?} [{}]",mat[i],out[i]);
    }
    println!("--");*/
    solve_mat(&mat,&out)
}
fn approx_solve_mat(mat : &Vec<Vec<Real>>, targ: &Vec<Real>) -> Vec<Real> {
    let (w,h) = (mat[0].len(),mat.len());
    let mut tab = mat.clone();
    let mut out = targ.clone();
    vec![]
}

//TODO!!!: take max of row for selection
fn solve_mat(mat : &Vec<Vec<Real>>, targ : &Vec<Real>) -> Vec<Real> {
    let (w,h) = (mat[0].len(),mat.len());
    let mut tab = mat.clone();
    let mut out = targ.clone();
    for r in 0..h {
        let mut m = tab[r][r];
        let mut mi = r;
        for i in r+1..h {
            if tab[i][r]*tab[i][r] > m*m {
                m = tab[i][r];
                mi = i;
            }
        }
        if mi != r {
            tab.swap(mi,r);
            out.swap(mi,r);
        }
        if m == 0.0 {
            continue
        }
        for sr in (r+1)..h {
            let val = - tab[sr][r] / tab[r][r];
            tab[sr][r] = 0.0;
            for i in (r+1)..w {
                tab[sr][i] += val * tab[r][i];
            }
            out[sr] += val * out[r];
        }
    }
    /*for i in 0..h {
        println!("{:?} [{}]",tab[i],out[i]);
    }*/
    for r in 0..h {
        if tab[r][r]*tab[r][r] < 1e-10 {
            for i in (r+1)..w {
                tab[r][i] = 0.0;
            }
            out[r] = 0.0;
        } else {
            for i in (r+1)..w {
                tab[r][i] /= tab[r][r];
            }
            out[r] /= tab[r][r];
            tab[r][r] = 1.0;
        }
    }
/*    for x in &tab {
        println!("{:?}",x);
    }*/
    for r in (0..h).rev() {
        for sr in (0..r) {
            out[sr] -= out[r] * tab[sr][r];
        }
    }
    out
}

struct MappingCost<Domain:Space+Debug,Range:Space+Debug,Mapping:Function<Domain=Domain,Range=Range>> {
    points: Vec<(Domain,Range)>,
    map: PhantomData<Mapping>
}

impl <D:Space + Debug,R:Space+Debug,M:Space + Function<Domain=D,Range=R>> Function for MappingCost<D,R,M> {
    type Domain=M;
    type Range= Real;
    fn eval(&self, func: M) -> Real {
        let basis: Vec<R> = R::basis().collect();
        let mut cost = 0.0;
        for x in &self.points {
            let (aim,target) = x.clone();
            let res = func.eval(aim.clone());
            let rcoords: Vec<Real> = res.coords().collect();
            let tcoords: Vec<Real> = target.coords().collect();
//            println!("{:?} | {:?}",x,rcoords);
            for j in 0..basis.len() {
                let diff = rcoords[j] - tcoords[j] ;
                cost = cost + diff*diff;
//                println!("{}",cost);
            }
        }
        cost
    }
}

impl <D:Space+Debug,R:Space+Debug,M:Space + Function<Domain=D,Range=R>> Differentiable for MappingCost<D,R,M> {
    fn partial(&self, func: M, dir: M) -> Real {
        let basis: Vec<R> = R::basis().collect();
        let m_basis: Vec<M> = M::basis().collect();
        let dir_coords: Vec<Real> = dir.coords().collect();
        let func_coords: Vec<Real> = func.coords().collect();
        let mut value = 0.0;
        for x in &self.points {
            let (aim,target) = x.clone();
            let res = func.eval(aim.clone());
            let dres = dir.eval(aim.clone());
            let rcoords: Vec<Real> = res.coords().collect();
            let dcoords: Vec<Real> = dres.coords().collect();
            let tcoords: Vec<Real> = target.coords().collect();
            for j in 0..basis.len() {
                value = value - 0.5 * (tcoords[j] - rcoords[j]) * dcoords[j];
            }
        }
        value
    }
}

fn minimize_simple<F:Differentiable + Function<Range=Real>>(func:F, start: F::Domain, stepsize: f64, stepnum: u32) -> F::Domain where F::Domain:Debug{
    let basis: Vec<F::Domain> = F::Domain::basis().collect();
    let mut current = start;
    println!("STARTING AT : {:?}",&current);
    for i in 0..stepnum {
        let mut things = Vec::new();
        for j in 0..basis.len() {
            let component = func.partial(current.clone(),basis[j].clone());
            things.push(-component * stepsize);
        }
        for j in 0..basis.len() {
            current = current + basis[j].clone() * things[j];
        }
        println!("{} : {} : {:?}",i,func.eval(current.clone()),&current);
    }
    current
}

fn newtonian_zero<F:Differentiable + Function<Range=Real>>(func:F,start:F::Domain, stepnum: u32) -> F::Domain where F::Domain:Debug {
    let basis: Vec<F::Domain> = F::Domain::basis().collect();
    let mut current = start;
    for i in 0..stepnum {
        let value = func.eval(current.clone());
        let mut things = Vec::new();
        let mut sum = 0.0;
        for j in 0..basis.len() {
            let component = func.partial(current.clone(),basis[j].clone());
            things.push(-component * value);
            sum = sum + component*component;
        }
        for j in 0..basis.len() {
            current = current + basis[j].clone() * things[j];
        }
        println!("{:?}",current)
    }
    current
}



impl Fourier {
    fn eval<O:Space>(&self, time:Real) -> O {
        let basis: Vec<O> = O::basis().collect();
        let mut out =  basis[0].clone() * 0.0;
        for i in 0..self.values.len()  {
            let mut sum = 0.0;
            for j in 0..self.values[i].len() {
                sum += self.values[i][j] * f64::sin(time * ((j+1) as f64) * f64::consts::PI)
            }
            out = out + basis[i].clone() * sum;
        }
        out
    }
}

fn regression<F:Function+Space>(data:&Vec<(F::Domain,F::Range)>) -> F {
    let basis:Vec<F> = F::basis().collect();
    let r_basis:Vec<F::Range> = F::Range::basis().collect();
    let mut columns = Vec::new();
    let mut targ = Vec::new();
    for x in 0..basis.len() {
        columns.push(Vec::new());
    }
    for x in data {
        let (inp,out) = x.clone();
        for i in 0..basis.len() {
            let res = basis[i].eval(inp.clone());
            columns[i].extend(res.coords());
        }
        targ.extend(out.coords());
    }
    let mut mat = Vec::new();
    for i in 0..columns[0].len() {
        mat.push(Vec::new());
        for j in 0..columns.len() {
            mat[i].push(columns[j][i]);
        }
        println!("{:?} [{}]",&mat[i],targ[i]);
    }
    println!("--");
    let endval = what_solve(&mat,&targ);
    let mut output = basis[0].clone() * 0.0;
    for i in 0..basis.len() {
        output = output + basis[i].clone() * endval[i];
    }
    output
}

fn smooth<I,S:Space>(choppy: I, size: usize)-> IntoIter<S> where I:Iterator<Item=S> { 
    let mut ring = VecDeque::new();
    let mut out = Vec::new();
    for item in choppy {
        ring.push_back(item);
        if ring.len() == size {
            let mut thing = ring.pop_front().unwrap();
            for x in &ring {
                thing = thing + x.clone();
            }
            thing = thing * (1.0/(size as f64));
            out.push(thing);
        }
    }
    out.into_iter()
}
fn pair_to_cell<A:Space,B:Space>((a,b):(A,B))->C2<A,B> {
    C2{
        a:a,
        b:b
    }
}
fn cell_to_pair<A:Space,B:Space>(x:C2<A,B>) -> (A,B) {
    (x.a,x.b)
}
