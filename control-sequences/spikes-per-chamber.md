

Expected number of spike molecules per capture chamber
======================================================

Summary
-------

Each capture chamber contains synthetic RNA molecules, to _spike_ the
libraries.  These RNA spikes are from the Ambion ArrayControl set (cat. num.
[AM1780M](http://www.lifetechnologies.com/order/catalog/product/AM1780M)).
Their base composition is as follows.  It differs from the manufacturer's
information, that miss some base at the 5′ end and the poly-A tail.


```r
s <- data.frame(
 row.names = c( "A", "C", "G", "U")
 ,      s1 = c( 191, 177, 215, 201)
 ,      s2 = c( 181, 226, 199, 178)
 ,      s3 = c( 194, 232, 347, 261)
 ,      s4 = c( 298, 236, 277, 219)
 ,      s5 = c( 272, 288, 274, 233)
 ,      s6 = c( 306, 347, 355, 276)
 ,      s7 = c( 361, 412, 405, 304)
 ,      s8 = c( 501, 545, 559, 429)
)
```


```r
base.molecular.weight <- c(A=347.2, C=323.2, G=363.2, U=324.2)
base.molecular.weight <- base.molecular.weight - 18
```

The molecular weight of a spike is the sum of the molecular weight of its
monophosphate bases, minus one molecule of water (18 _g/mol_) per bond:
329.2, 305.2, 345.2, 306.2 _g/mol_ for A, C, G, and U.


```r
molecular.weight <- colSums ( s * base.molecular.weight )
```

The molecular weight for each spike (from 1 to 8) is 2.527 &times; 10<sup>5</sup>, 2.518 &times; 10<sup>5</sup>, 3.344 &times; 10<sup>5</sup>, 3.328 &times; 10<sup>5</sup>, 3.434 &times; 10<sup>5</sup>, 4.137 &times; 10<sup>5</sup>, 4.775 &times; 10<sup>5</sup>, 6.556 &times; 10<sup>5</sup> _g/mol_.

The original C1 RNA-seq protocol (PN 100-5950 B1) recommends a serial dilution
of 3 spikes.  For a better calibration, we added 2 more diluted spikes, as
follows.

```
||=Tube                                =||=   A  =||=   B  =||=   C  =||=   D  =||=   E   =||
||THE RNA Storage Solution              || 13.5 μL|| 12.0 μL|| 12.0 μL|| 12.0 μL|| 148.5 μL||
||RNA Spikes  (1.5 μL, 100 ng/μL)       || spike 6|| spike 3|| spike 7|| spike 4||  spike 1||
||Add 1.5 μL from tube:                 ||    —   ||    A   ||    B   ||    C   ||    D    ||
```



```r
dilution <- c(s6=10000, s3=1000, s7=100, s4=10, s1=1) * 100
```

The dilution factors are therefore: 10<sup>6</sup>, 10<sup>5</sup>, 10<sup>4</sup>, 1000, 100.

The stock solution for all spikes is 100 _ng/μL_, that is, 0.1 _g/L_.


```r
stock <- 0.1
mass.concentration <- stock / dilution
```

The mass concentration is therefore 10<sup>-7</sup>, 10<sup>-6</sup>, 10<sup>-5</sup>, 10<sup>-4</sup>, 0.001 _g/L_.


```r
molar.concentration <- mass.concentration / molecular.weight[c('s6', 's3', 's7', 's4', 's1')]
```

The molar concentration (molarity) is the mass concentration divided by the
molecular weight: 2.417 &times; 10<sup>-13</sup>, 2.991 &times; 10<sup>-12</sup>, 2.094 &times; 10<sup>-11</sup>, 3.005 &times; 10<sup>-10</sup>, 3.958 &times; 10<sup>-9</sup> _mol/L_.


```r
N <- 6.02e23
number.concentration <- molar.concentration * N
```

The number of molecules in one liter is the molar concentration multiplied by
the Avogadro number: 1.455 &times; 10<sup>11</sup>, 1.8 &times; 10<sup>12</sup>, 1.261 &times; 10<sup>13</sup>, 1.809 &times; 10<sup>14</sup>, 2.383 &times; 10<sup>15</sup>.


```r
chamber.volume <- 9e-9
spike.mix.dilution.factor <- 22.5 / 1.5 * 100
number.per.chamber <- number.concentration * chamber.volume / spike.mix.dilution.factor
```

Each chamber contains 9 _nL_ of lysis mix, where the spike mix was diluted 100
times in C1 loading reagent, which is then diluted 15.013 times.

**The number of spike molecules per chamber is therefore 0.873, 10.802, 75.648, 1085.314, 1.43 &times; 10<sup>4</sup>
for s6, s3, s7, s4, s1 respectively.**

Note that the spike counts are not exactly following the 10-fold dilution
because of difference of length and base composition.


```r
plot(number.per.chamber, dilution, log='xy', type='b')
```

![plot of chunk spike-dilution-plot](figure/spike-dilution-plot.png) 
