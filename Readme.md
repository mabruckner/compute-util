# My personal mathematical utilities

First, an apology. This library as is currently stands is the most comprehensive, powerful, and cohesive set of mathematical utilities I have ever written. It is also horrible code.
If you find yourself in need of some of this code, if you have no alternative, I am deeply sorry.

Now, a basic overview.
This library defines the notion of a real vector space in rust. All operations are 64-bit floating point (though if you redefine Real it will probably work).
Matrix operations assume row-major order (at least, I think thats what it would be called).
I **DID NOT** write this code to be understandable, fast, or safe. I wrote this code to work. Maybe someday I will find the time or be forced to make this code better in those ways.

## Why?

I threw together the base of this framework in order to prototype a novel machine learning algorithm.
(Ever since a frusturating encounter with BLAS in high school I have always written my own matrix solving functions.)
A couple of weeks later I had a completely seperate idea that would benefit from the matrix operations and vector space abstractions that I had just finished implementing.
These abstractions have problems and holes, and eventually I may have to refactor.
However, rusts typing system makes the current iteration extremely powerful, and is more than sufficient for my purposes.

