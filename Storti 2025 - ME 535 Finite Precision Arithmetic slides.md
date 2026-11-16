---

theme: default #gaia
paginate: true
# class:
#  - lead
size: 16:9

marp: true
footer: "ME 535 Winter 2025, Prof. Storti, UW-Seattle"

---

# Finite Precision Arithmetic

Based on: __Lessons in Scientific Computing: Numerical Mathematics, Computer Technology, and Scientific Discovery__ ByNorbert Schorghofer
Ch. 3 - _Roundoff & Number Representation_

E-book available from UW Libraries. Ch. 3 on Canvas (Week 2).

For _much_ more detail on IEEE754 Standard  Floating-Point Arithmetic, see "What Every Computer Scientist Should Know About Floating-Point Arithmetic," by David Goldberg, ACM Computing Surveys. Available in the Canvas "readings" folder.

---

## A bit of context

_What is this class all about?_   __Practical__ numerical methods for obtaining good approximate solutions to engineering problems (for which analytic solutions may not be readily available)

What are __numerical methods__ about?
Using digital computers to obtain those good approximate solutions...

---

Computers are both __fast/capable__ and __dumb/literal__.

The computer is __not__ smarter than you, so you have to tell it __very explicitly__ what you want it to do. (_Exception?_)

The specification of what to do is an __algorithm__: 
Unambiguous finite rule that, after a finite number of steps, provides a solution to a class of problems.

We will look at __algorithms__ and __implementations__ as codes/programs.

Algorithms have inputs (often in $ℤ^n$ or $ℝ^n$) that are converted into outputs by a sequence of __elementary operations__ including `+,-,*,/`.

---

## Classes of problems

__Ill-conditioned problem__:
Small change in input can cause LARGE change in output

__Well-conditioned problem__:
Small change in input $\iff$ small change in output

Identify well-conditioned problems we can solve with __numerically stable algorithms__:
Actual computational error comparable to __unavoidable error__.

---

Why are there unavoidable errors?
- Exact representation of $ℝ$ requires infinitely many digits.
- Infinite data storage requires infinite computing resources.
- Memory chip prices have come down,  
but infinite memory is impossible (physically and economically).

Typical compromise: use finite precision numeric representations. 

---

#### $\S$ 3.1 - Number Representation

Real numbers in an algorithm get approximated for digital computing purposes as __floating-point__ numbers of the form:
$$(-1)^s \; (d_0.d_1 d_2 \ldots d_{p-1}) \; \beta^e$$
- $s\in \{0,1\}$ indicates the sign
- $p$ is the number of digits or __precision__
- each digit $d_i \in \{0, 1, \ldots, \beta-1\}$
- $\beta$ is the base and $e$ is the exponent.

Humans may prefer base $\beta = 10$  
Computers almost universally use a __binary representation with $\beta = 2$__.

---

A floating-point representation corresponds to a unique real number:
$\bar{x} = (-1)^s \big( d_0 + \frac{d_1} {\beta} + \frac{d_2}{\beta^2} + \ldots + \frac{d_{p-1}}{\beta^{p-1}} \big) \beta^e$

Some real numbers have an exact (but not unique) floating-point representation; e.g. $x=0.5$ has the following representations:
$\bar{x} = (-1)^0 \; (1.0) \; 2^{-1} = (-1)^0 \; (1+\frac{0}{2}) \; 2^{-1} = \frac{1}{2}$
and
$\bar{x} = (-1)^0 \; (0.1) \; 2^0 = (-1)^0 \; (0+\frac{1}{2}) \; 2^0 = \frac{1}{2}$

For uniqueness, choose the __normalized__ (first) version with $d_0 \neq 0$.

---

A floating point number system is defined by:
- $\beta$ : base 
- $p$ : precision (number of digits) 
- $\{e_{min}, e_{max}\}$ :  range of (int) exponents 

Let's take a look at an example: $\beta=2, p=2, e \in [-2,3]$.
Positive normalized numbers:
$\bar{x} = (-1)^0 \; (1.0) \; 2^{-2} = \frac{1}{4}$, $\,\bar{x} = (-1)^0 \; (1.1) \; 2^{-2} = \frac{3}{8}$
$\bar{x} = (-1)^0 \; (1.0) \; 2^{-1} = \frac{1}{2}$, $\,\bar{x} = (-1)^0 \; (1.1) \; 2^{-1} = \frac{3}{4}$
$\bar{x} = (-1)^0 \; (1.0) \; 2^0 = 1$, $\quad \bar{x} = (-1)^0 \; (1.1) \; 2^0 = \frac{3}{2}$
$\bar{x} = (-1)^0 \; (1.0) \; 2^1 = 2$, $\quad \bar{x} = (-1)^0 \; (1.1) \; 2^1 = 3$
$\bar{x} = (-1)^0 \; (1.0) \; 2^2 = 4$, $\quad \bar{x} = (-1)^0 \; (1.1) \; 2^2 = 6$
$\bar{x} = (-1)^0 \; (1.0) \; 2^3 = 8$, $\quad \bar{x} = (-1)^0 \; (1.1) \; 2^3 = 12$

---

## Gamut $\bf{(X)}$ and observations

$\bar{x} \in {\bf{X}} = ± \{1/4, 3/8, 1/2, 3/4, 1, 3/2, 2, 3, 4, 6, 8, 12 \}$

- Not all reals can be represented

- No representation for zero (although one can be designated)

- Uneven spacing: __gapsize roughly proportional to magnitude__

Typical usage: round real number $x$ to nearest element of $\bf{X}$.
- Rounding incurs an __absolute error__: $\lvert x-\bar{x} \rvert = E(x)$
- More often consider __relative error__: $\frac{E(x)}{\lvert x \rvert} = R(x)$

---

## Floating-point representation and error

Floating point aims to even out __relative error__:  larger gaps between larger numbers
$\implies E(x) \overset{\sim} \propto |x| \implies R(x) \text{ has a constant bound called \bf{unit roundoff}}$

Floating point systems have arithmetic limitations:
For $\bar{x} \in {\bf{X}} = ± \{1/4, 3/8, 1/2, 3/4, 1, 3/2, 2, 3, 4, 6, 8, 12 \}$
$12 * 2 = 24 \notin {\bf X} \implies$ overflow ($±$ `Inf`)
$\frac{1}{4} / 2 = \frac{1}{8} \notin {\bf X} \implies$ underflow ( ±`0`, `NaN`)
<br>IEEE754 __de-normalizes__ for gradual underflow: $(-1)^0 \; (0.1) \; 2^{-2} = \frac{1}{8}$

---

#### Section 3.2 - IEEE Standardization

Most common floating point systems: __bit__= Binary digIT
- __Double precision__: 52 mantissa/fraction bits, 1 sign bit, 11 exponent bits
- __Single precision__: 24 fraction bits, 1 sign bit, 7 exponent bits
- __Half precision__: 10 fraction bits, 1 sign bit, 5 exponent bits


Note that bit counts are multiples of 8 because hardware organizes the bits into __bytes__ (1 byte = 8 bits, double = 8 bytes, single = 4 bytes, half = 2 bytes)

__Double precision:__ 8 bytes, ~15-17 significant decimal digits
__Single precision:__ 4 bytes, ~6-9 significant decimal digits
__Half precision:__ 2 bytes, ~2-3 significant decimal digits (When/why is this relevant?)

See Table 3.1 for details of representation range.

---

Consider error for $x \in ℝ$ that rounds to $\bar{x}$ so they agree to $p$ digits: 
$$\bar{x} = x \; (1 + \epsilon)$$

$R(x) =  \lvert \epsilon \rvert \leq (1/2) \beta^{1-p} = u$ (__unit roundoff__: upper bound on $R(x)$ )
Spacing between normalized floats is $\approx 2 u \lvert x \rvert$
For $x=1$, the gap is $2u$, so the next f.p. number available is $1 + 2u$.

Alternate precision spec of machine epsilon, $\epsilon_{M}$:
Smallest number you can add to 1 without producing the result 1: 
$$\overline{1+x} = 1 \; \forall x \; s.t. \; |x|<\epsilon_M$$

Others define $u$ to differ by a factor of 2, but either $u$ or $\epsilon_M$ characterize f.p. precision.

---

### Section 3.3 - Roundoff Sensitivity

Consider how roundoff error in inputs propagates when performing floating point arithmetic operations:

Consider real inputs $x_1$ and $x_2$ that round to $\bar{x}_1 = x_1 (1+\epsilon_1)$, $\bar{x}_2 = x_2 (1+\epsilon_2)$

$\epsilon$ indicates signed relative error bounded by $u < 10^{-16} <<1$ (double precision)

Determine bounds on relative error for result of basic arithmetic operations.

---

### Multiplication: 

$$\bar{x}_1 * \bar{x}_2 = x_1 (1+\epsilon_1) * x_2 (1+\epsilon_2)$$
$$\bar{x}_1 * \bar{x}_2 = x_1 * x_2 \; (1 + \epsilon_1 +\epsilon_2 + \epsilon_1 * \epsilon_2)$$

Product of relative errors is VERY small $(\epsilon_1 * \epsilon_2 << 1)$, so ignore...
$$\bar{x}_1 * \bar{x}_2 = x_1 * x_2 \; (1 + \epsilon_1 +\epsilon_2 + \ldots) = (x_1*x_2)(1+\epsilon_*)$$
$\implies$ __when multiplying, relative errors add__.

---

### Division:
$$\bar{x}_1 / \bar{x}_2 = x_1 (1+\epsilon_1) / (x_2 (1+\epsilon_2))$$
$$\bar{x}_1 / \bar{x}_2 = x_1 / x_2 (1 + \epsilon_1) * \frac{1}{1+\epsilon_2}$$ 
$$\approx x_1 / x_2 * (1 + \epsilon_1) * (1 - \epsilon_2 + \epsilon^2 + \ldots)$$
$$\approx x_1 / x_2 * (1 + \epsilon_1 - \epsilon_2 + \ldots)$$
$$\approx x_1 / x_2 * (1 + \epsilon_\div)$$

Again ignore product of relative errors,  
so output error is bounded by sum of input relative errors.
$\implies$__Mutiply or divide: relative errors add.__ 
(Output relative error bounded by sum of input relative errors.)

---

### Addition: 
$$\bar{x}_1 + \bar{x}_2 = x_1 (1+\epsilon_1) + x_2 (1+\epsilon_2)$$
$$= x_1 + x_2 + (x_1 * \epsilon_1 + x_2 * \epsilon_2) = (x_1 + x_2) * (1 + \epsilon_+)$$
$$\implies \lvert \epsilon_+ \rvert = \frac{ \lvert x_1 * \epsilon_1 + x_2 * \epsilon_2 \rvert}{\lvert x_1 + x_2 \rvert} \lessapprox u$$


__Addition:__ Relative error near unit roundoff 
__if $x_1$ and $x_2$ have the same sign!__

---

### Addition/Subtraction

$$\lvert \epsilon_\pm \rvert = \frac{ \lvert x_1 * \epsilon_1 + x_2 * \epsilon_2 \rvert}{\lvert x_1 + x_2 \rvert} $$

When adding $x_1 \approx -x_2$ (or subtracting $x_1-x_2$ with $x_1 \approx x_2$) 
- denominator becomes arbitrarily small
- error is not well bounded.

## ___Catastrophic cancellation___: 
__Subtracting nearly equal numbers wipes out the significant digits and leaves noise...__

---

### Dealing with limitations of fixed precision (example 1)

Compute $x^2 - y^2$ in example number system with $x=2$ and $y=1$
$\bar{x} \in {\bf{X}} = ± \{1/2, 3/4, 1, 2, 3\}$
Blindly computing $x*x = 2*2$ causes overflow. 
Instead, rewrite the expression:  $x^2 - y^2 = (x-y) * (x+y)$
Computation becomes $(2-1) * (2+1) = 1*3 = 3$ 
Avoids overflow and produces a result (accidentally exact)
More usual case: $x=2, y=\frac{3}{4}$ 
Exact (and rounded) result: 
$(2-\frac{3}{4})*(2+\frac{3}{4}) \rightarrow 1*3 = 3$ (differs from exact result $\frac{7}{4}$)
__More general approach to overflow__: 
Normalize with units that make your numerical values close to 1

---

### Dealing with finite precision (Example 2: Quadratic formula)


Evaluate roots of $a x^2 + b x + c = 0 = a(x-x_+)(x-x_-)$ using quadratic formula: 
$$x_\pm = \frac{-b \pm \sqrt{b^2- 4 a c}}{2 a}$$
What about this could be problematic?

---

### Quadratic formula (cont'd): $\quad a x^2 + b x + c = a(x-x_+)(x-x_-) = 0 \implies x_\pm = \frac{-b \pm \sqrt{b^2- 4 a c}}{2 a}$
When $4 a c \ll b^2$ the square root $r = \sqrt{b^2-4ac}$ is nearly equal to $b$.  
If $b>0$ the root $x_+$ is subject to catastrophic cancellation (when computing $(-b+r)$). 
However the root $x_-$ should be OK, and we know that $a* x_+ * x_- = c$.
$\implies$ Effective plan for $b>0$ : 
1) Compute $q = b + \sqrt{b^2 -4 a c}$
2) Compute $x_1 = -\frac{q}{2a}$, $x_2 = \frac{c}{a x_1} = -2 \frac{c}{q}$
Possible improvement: 
Compute $q = (-1/2) (b + \mathrm{sgn(b)} \sqrt{b^2 -4 a c}) \implies$ $x_1 = \frac{q}{a}$, $x_2 = \frac{c}{a x_1} = \frac{c}{q}$
Note that basic python does not include a "sign" function.  
Define your own or import from python's symbolics package `sympy`.

---

### Other "tricks"

Sometimes you can rearrange a compound computation to avoid cancellation.
Consider a sum where cancellation would be of concern (Why?):
$$ S = 1 -\frac{1}{2}+\frac{1}{3}- \frac{1}{4}+\ldots$$
Collect pairs of terms over common denominator:
$$ S = \frac{1}{1 \cdot 2}+\frac{1}{3 \cdot 4} +\ldots$$
Eliminates subtraction and avoids cancellation trouble.

---

#### Interval Arithmetic (IA): Alternative to dealing with roundoff errors (aside)

- Instead of trying to represent a particular real number (which we cannot do with finite precision), keep track of the lower and upper bounds of an interval that is guaranteed to contain the exact number

- Operation output must guarantee inclusion of exact answer

- Useful property for root-finding/isolation: finite convergence

--- 

### IA example due to Rump
E. Loh, G. W. Walster, “Rump’s example revisited”, Reliable Computing, vol. 8 (2002), n. 2, pp. 245–248.

$$f(x,y) = (\frac{1335}{4} - x^2) y^6 +x^2 (11 x^2 y^2 - 121 y^4 -2) + \frac{11}{2} y^8 + \frac{x}{2 y}$$

With rational inputs, this function with rational coefficients can be evaluated exactly (e.g. with Mathematica): 
$$f(77617, 33096) = -\frac{54767}{66192} \approx −0.827396\ldots$$

Now let's try floating point evaluation:

---

## Floating point evaluation of Rump's example
```python
def f(x,y):
    return ((333.75 - x**2)* y**6 + x**2 * (11* x**2 * y**2 - 121 * y**4 - 2) + 5.5 * y**8 + x/(2*y))

f(77617.0, 33096.0)
```
$\rightarrow 1.1726039400531787 \neq -0.827396\ldots$
Can't even evaluate, so would there be any chance of locating roots?
Answer turns out to be YES (but stick to  simpler case with 1 variable)

---

## IA rootfinding/exclusion
Basic plan: 
- Define __interval extensions__ of functions that reliably contain correct results
- Evaluate interval extension function over input interval to obtain an output interval
- If output interval includes 0, then input interval can contain a root
- Subdivide and evaluate on subinterval until the output interval excludes 0  
$\implies$ NO ROOTS in the subinterval
- In finite steps, obtain set of narrow intervals that are candidate root locations.

---

![title](intervalRoots.png)


