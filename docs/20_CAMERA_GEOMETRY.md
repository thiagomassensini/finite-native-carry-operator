# Camera geometry

## C2 aligned chart

At cutoff \(M\):

- seed: \(\psi_t(1)\);
- centers: \(c=4,8,\ldots,4M\);
- one radius: \(r=1\);
- exactly \(M\) brackets.

Each center can be classified by its dyadic depth, but this classification does
not insert an extra numerical weight into the operator.

## Natural saturated cameras

For \(b\ge3\), write \(h=\lfloor b/2\rfloor\).  The camera uses seeds
\(1,\ldots,h\), centers \(b,2b,\ldots,Mb\), and all radii \(1,\ldots,h\).
It therefore emits \(Mh\) brackets.

## Even antipodal channel

For even \(b\), \(r=b/2\) is antipodal modulo \(b\), but the integers
\(c-r\) and \(c+r\) are distinct legs.  In C4, `r=1` reads odd neighbors and
`r=2` reads even integers with dyadic valuation one.  The full C4 camera is not
the C2 chart; C2 matches only its radius-one sector at centers \(4m\).
