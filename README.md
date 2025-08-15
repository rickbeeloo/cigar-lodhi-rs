# Lodhi on cigars

## Intro
Most scoring schemes, such as unit cost edit distnace, ignore the alignment 
structure. Take the two alignments below:

```python
A            B
 ATAGCA        ATAGCA
 --||||        |-||-|
 CCAGCA        AGAGTA
```
Many would say that `A` looks more "natural" because of having a single 
"mutation" instead of scattered mutations as in `B`. However, both have 
the same edit distance of 2. 

## Lodhi
A 2002 paper, [Lodhi et. al. 2002](https://www.jmlr.org/papers/volume2/lodhi02a/lodhi02a.pdf) introduced  <i>"Text Classification using String Kernels"</i>. The idea is that instead of treating strings as flat sequences of symbols, we can represent them as weighted sets of all subsequences of length k, where longer consecutive matches contribute more.

todo: finish