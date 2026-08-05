# Precision and cutoff record

## Fixed-cutoff arithmetic

At camera 3 and \(M=16\,364\), requested stages 40, 60, 100, and 160 decimal
digits converged to the same finite stationary minimum with 118 significant
digits of reported arithmetic agreement.

This establishes that binary64 was not the limiting arithmetic representation
for that finite problem.  It does not remove finite-cutoff drift.

## Cutoff ladder

The preserved ledgers contain the following late C3 minima:

| \(M\) | finite minimum |
|---:|---:|
| 65,536 | 92.4918991964453292025855476777… |
| 131,072 | 92.4918992692865555029556892732… |
| 262,144 | 92.4918992795639263431084313891… |
| 524,288 | 92.4918992725523190067519554385… |
| 1,048,576 | 92.4918992698389339609338560010… |
| 2,097,152 | 92.4918992701626961792035134570… |
| 4,194,304 | 92.4918992705678195408206862095… |

The approach is oscillatory rather than monotone.  Direct agreement of the
last two cutoffs is weaker than the internal arithmetic agreement, which is
why a separate tail model was investigated.
