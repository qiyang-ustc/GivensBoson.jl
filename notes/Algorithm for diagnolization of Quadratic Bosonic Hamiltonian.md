## Algorithm for Diagnolizing Quadratic Bosonic Hamiltonian

What we want to diagnolize is:
$$
H_1 = |J|S\sum_{\langle ij\rangle}[b^{\dagger}_ib_i+b^{\dagger}_jb_j-b_i^{\dagger}b_j^{\dagger}-b_ib_j]
$$
In matrix form, $H_1$ is real, sparse and symmetric . These properties allow us using (generalized) Givens Rotation to diagnolize it.

##### Why previous work are based on eigen problem?

For a quadratic form $H$, we want to find a matrix $G$,  which satisfy $G^{\dagger}\eta G=\eta$  and $G^{\dagger}HG =\Omega$ where $\Omega$ is diagnolized.

They do not know how to find $G$ directly, so they tried to convert this problem to a eigen problem:

As $G^{\dagger}\eta G = \eta$ ,$\eta G^{\dagger}\eta G = 1$, $G^{-1} = \eta G^{\dagger} \eta$

so, $G^{\dagger}HG = \eta (G^{-1}\eta HG)$. They convert this problem to an eigen problem. When they try to solve eigen problem, they should suffer many problem such as degenerate space and zero modes.

**We can deal with this "quadratic form reduction problem" directly by using Givens Rotation**

#### Givens Rotation in normal metric situation

A Givens rotation can be regarded as a rotation which mix two eigen-vectors.
$$
G(i, j, \theta)=\left[\begin{array}{ccccccc}
1 & \cdots & 0 & \cdots & 0 & \cdots & 0 \\
\vdots & \ddots & \vdots & & \vdots & & \vdots \\
0 & \cdots & c & \cdots & -s & \cdots & 0 \\
\vdots & & \vdots & \ddots & \vdots & & \vdots \\
0 & \cdots & s & \cdots & c & \cdots & 0 \\
\vdots & & \vdots & & \vdots & \ddots & \vdots \\
0 & \cdots & 0 & \cdots & 0 & \cdots & 1
\end{array}\right]
$$
If $c^2+s^2=1$, in another word: $c = \cos(\theta),s=\sin(\theta)$. $G(i,j,\theta)$ is unitary. And $G^{-1}(i,j,\theta)=G(i,j,-\theta)$

In $S' = G(i,j,\theta)SG(i,j,-\theta)$, $i<j$
$$
\begin{align}
&S'_{ii}=G_{ik}S_{kl}G^{-1}_{li}= c^2S_{ii}-2csS_{ij}+s^2S_{jj}\\
&S'_{jj}=G_{jk}S_{kl}G^{-1}_{lj}= c^2S_{jj}+2csS_{ij}+s^2S_{ii}\\
&S'_{ij}=S'_{ji}=G_{ik}S_{kl}G^{-1}_{lj}= (c^2-s^2)S_{ij}+cs(S_{ii}-S_{jj})\\
&S'_{ik}=S'_{ki}=G_{il}S_{lm}G^{-1}_{mk}=cS_{ik}-sS_{jk}   &k\neq i,j\\
&S'_{jk}=S'_{kj}=G_{jl}S_{lm}G^{-1}_{mk}=cS_{jk}+sS_{ik}   &k\neq i,j\\
&S'_{kl} = S_{kl} &k,l\neq i,j
\end{align}
$$
We can choose $\theta$ to let $S'_{ij}=0$, where we need to read:
$$
\cos(2\theta)S_{ij}+\frac{1}{2}\sin(2\theta)(S_{ii}-S_{jj})=0\longrightarrow \tan(2\theta) = \frac{2S_{ij}}{S_{jj}-S_{ii}}
$$
If $S_{ii}=S_{jj}$, $\theta = \frac{\pi}{4}$.

#### Givens Rotation in bosonic metric situation

In bosonic problem, the canonical relation for $G$ become $G\eta G^{-1}=\eta$ , where $\eta = \text{diag}(I_{n\times n},-I_{n\times n})$.

We need to modify the origin Givens Rotation, for $2n\times 2n$ real symmetric matrix $S$:
$$
G_a(i, j, \theta)=\left[\begin{array}{ccccccc}
1 & \cdots & 0 & \cdots & 0 & \cdots & 0 \\
\vdots & \ddots & \vdots & & \vdots & & \vdots \\
0 & \cdots & c & \cdots & s & \cdots & 0 \\
\vdots & & \vdots & \ddots & \vdots & & \vdots \\
0 & \cdots & s & \cdots & c & \cdots & 0 \\
\vdots & & \vdots & & \vdots & \ddots & \vdots \\
0 & \cdots & 0 & \cdots & 0 & \cdots & 1
\end{array}\right]  ,\text{if}: i<n,j>n
$$

$$
G_n(i, j, \theta)=\left[\begin{array}{ccccccc}
1 & \cdots & 0 & \cdots & 0 & \cdots & 0 \\
\vdots & \ddots & \vdots & & \vdots & & \vdots \\
0 & \cdots & c & \cdots & -s & \cdots & 0 \\
\vdots & & \vdots & \ddots & \vdots & & \vdots \\
0 & \cdots & s & \cdots & c & \cdots & 0 \\
\vdots & & \vdots & & \vdots & \ddots & \vdots \\
0 & \cdots & 0 & \cdots & 0 & \cdots & 1
\end{array}\right],\text{if:} i,j<n,\text{or}: i,j>n
$$

There are two types of Givens Rotation in this picture. We could call them **normal** and **abnormal** Givens Rotation.

The abnormal Givens Rotation is (10), its inverse is
$$
G_a^{-1}(i, j, \theta)=\left[\begin{array}{ccccccc}
1 & \cdots & 0 & \cdots & 0 & \cdots & 0 \\
\vdots & \ddots & \vdots & & \vdots & & \vdots \\
0 & \cdots & c & \cdots & -s & \cdots & 0 \\
\vdots & & \vdots & \ddots & \vdots & & \vdots \\
0 & \cdots & -s & \cdots & c & \cdots & 0 \\
\vdots & & \vdots & & \vdots & \ddots & \vdots \\
0 & \cdots & 0 & \cdots & 0 & \cdots & 1
\end{array}\right]  ,\text{if}: i<n,j>n
$$
It can be verified all types of Givens Rotation satisfy: (need $c^2-s^2=1$)
$$
\begin{align}&G^{\dagger}\eta G = \eta\\
&G^{-1}G=I_{2n\times 2n}

\end{align}
$$
We also need to modify the real symmetric hamiltonian to a special form which owns a set of new symmetry:
$$
H\longrightarrow \eta H=S=\begin{pmatrix}A&B\\-B &-A\end{pmatrix}
$$
This $H$ has such a set of symmetry:
$$
\begin{align}
&S_{ij} = S_{ji} , \text{if}\ i,j<n, \text{or}, i,j>n\\
&S_{ij} = -S_{ji} , \text{if}\ i>n,j<n, \text{or}, i<n,j>n\\
&S_{ij}=-S_{i+n,j+n},\text{if}\ i,j<n
\end{align}
$$
Suppose $S$ has such symmetry, we can verified the effect of abnormal Givens Rotation is (only one of $i,j$ > n): 
$$
\begin{align}
&S'_{ii}=G_{ik}S_{kl}G^{-1}_{li}= c^2S_{ii}-2csS_{ij}-s^2S_{jj}\\
&S'_{jj}=G_{jk}S_{kl}G^{-1}_{lj}= c^2S_{jj}-2csS_{ji}-s^2S_{ii}\\
&S'_{ij}=-S'_{ji}=G_{ik}S_{kl}G^{-1}_{lj}= cs(S_{jj}-S_{ii})+(c^2+s^2)S_{ij}\\
&S'_{ik}=G_{il}S_{lm}G^{-1}_{mk}=cS_{ik}+sS_{jk}   &k\neq i,j\\
&S'_{jk}=G_{jl}S_{lm}G^{-1}_{mk}=cS_{jk}+sS_{ik}   &k\neq i,j\\
&S'_{kl} = S_{kl} &k,l\neq i,j
\end{align}
$$
It can be verified $S'_{ki}$ and $S'_{kj}$ follows the top two symmetry of $S$, but will break the last symmetry.

We can let $S'_{ij}=0$, while $c= \cosh(\theta),s=\sinh(\theta)$
$$
\frac{1}{2}\sinh(2\theta)(S_{jj}-S_{ii})=\cosh(2\theta)S_{ij}\rightarrow \tanh(2\theta) = \frac{2S_{ij}}{S_{jj}-S_{ii}}
$$
Obviously, abnormal Givens Rotation broke down when $S_{jj}=S_{ii}$, but we have a trick to avoid this problem:

Suppose $i<n<j$, if $i+n\neq j$, after an abnormal Givens Rotation $G_a(i,j,\theta)$, We could immediately do another abnormal Givens Rotation:$G_{a}(j-n,i+n,\theta)$,

In another way, let $a=i,b=j-n$, $\tilde{b} =b+n,\tilde{a}=a+n $.  The composite Givens Rotation is $G_{ca}(a,b,\theta)=G_a(a+n,b,\theta)G_a(a,b+n,\theta)$, whose effect is
$$
\begin{align}
&S'_{aa}=G_{ak}S_{kl}G^{-1}_{la}= c^2S_{aa}-2csS_{a,\tilde{b}}-s^2S_{\tilde{b}\tilde{b}}\\
&S'_{\tilde{b}\tilde{b}}=G_{\tilde{b}k}S_{kl}G^{-1}_{l\tilde{b}}= c^2S_{\tilde{b}\tilde{b}}-2csS_{\tilde{b}a}-s^2S_{aa}\\
&S'_{\tilde{a}\tilde{a}}=G_{\tilde{a}k}S_{kl}G^{-1}_{l\tilde{a}}= c^2S_{\tilde{a}\tilde{a}}-2csS_{\tilde{a},b}-s^2S_{bb}\\
&S'_{bb}=G_{bk}S_{kl}G^{-1}_{lb}= c^2S_{bb}-2csS_{b\tilde{a}}-s^2S_{\tilde{a}\tilde{a}}\\
&S_{\tilde{a}b}=S_{a\tilde{b}}=S_{\tilde{b}a}=S_{\tilde{a}b}=0\\
&S'_{ak}=G_{al}S_{lm}G^{-1}_{mk}=cS_{ak}+sS_{\tilde{b}k}   &k\neq a,b,\tilde{a},\tilde{b}\\
&S'_{\tilde{b}k}=G_{\tilde{b}l}S_{lm}G^{-1}_{mk}=cS_{\tilde{b}k}+sS_{ak}   &k\neq a,b,\tilde{a},\tilde{b}\\
&S'_{\tilde{a}k}=G_{\tilde{a}l}S_{lm}G^{-1}_{mk}=cS_{\tilde{a}k}+sS_{bk}   &k\neq a,b,\tilde{a},\tilde{b}\\
&S'_{bk}=G_{bl}S_{lm}G^{-1}_{mk}=cS_{bk}+sS_{\tilde{a}k}   &k\neq a,b,\tilde{a},\tilde{b}\\
&S'_{kl} = S_{kl} &k,l\neq a,b,\tilde{a},\tilde{b}
\end{align}
$$
Similarly, $G_{cn}=G_n(a,b,\theta)G_n(a+n,b+n,\theta)$ for $a,b<n$.

$\{G_{cn},G_{ca}\}$ will keep all symmetry of $S$, which can prevent $S_{ii}=S_{jj}$ in abnormal Givens Rotation.

**It is also very straightforward that Bogoliubov transformation is an abnormal Givens Rotation! **



#### Directly reduct quadratic form by using Givens Rotation

Now $H$ is Hermitian, we will not try to solve eigen problem of $\eta H$.

We modify every steps in abnormal givens rotation:$G^{-1}HG\rightarrow G^{\dagger}HG$, that's all. $G^{\dagger}HG$'s effect is:

For 
$$
G_a(i, j, \theta)=\left[\begin{array}{ccccccc}
1 & \cdots & 0 & \cdots & 0 & \cdots & 0 \\
\vdots & \ddots & \vdots & & \vdots & & \vdots \\
0 & \cdots & c & \cdots & s & \cdots & 0 \\
\vdots & & \vdots & \ddots & \vdots & & \vdots \\
0 & \cdots & s & \cdots & c & \cdots & 0 \\
\vdots & & \vdots & & \vdots & \ddots & \vdots \\
0 & \cdots & 0 & \cdots & 0 & \cdots & 1
\end{array}\right]  ,\text{if}: i<n,j>n
$$

$$
\begin{align}
&S'_{ii}=G_{ik}S_{kl}G^{\dagger}_{li}= c^2S_{ii}+2csS_{ij}+s^2S_{jj}\\
&S'_{jj}=G_{jk}S_{kl}G^{\dagger}_{lj}= c^2S_{jj}+2csS_{ji}+s^2S_{ii}\\
&S'_{ij}=S'_{ji}=G_{ik}S_{kl}G^{\dagger}_{lj}= cs(S_{jj}+S_{ii})+(c^2+s^2)S_{ij}\\
&S'_{ik}=G_{il}S_{lm}G^{\dagger}_{mk}=cS_{ik}+sS_{jk}   &k\neq i,j\\
&S'_{jk}=G_{jl}S_{lm}G^{\dagger}_{mk}=cS_{jk}+sS_{ik}   &k\neq i,j\\
&S'_{kl} = S_{kl} &k,l\neq i,j
\end{align}
$$

so 
$$
\frac{1}{2}\sinh(2\theta)(S_{jj}+S_{ii})=\cosh(2\theta)S_{ij}\rightarrow \tanh(2\theta) = \frac{2S_{ij}}{S_{jj}+S_{ii}}
$$
For  
$$
G_n(i, j, \theta)=\left[\begin{array}{ccccccc}
1 & \cdots & 0 & \cdots & 0 & \cdots & 0 \\
\vdots & \ddots & \vdots & & \vdots & & \vdots \\
0 & \cdots & c & \cdots & s & \cdots & 0 \\
\vdots & & \vdots & \ddots & \vdots & & \vdots \\
0 & \cdots & -s & \cdots & c & \cdots & 0 \\
\vdots & & \vdots & & \vdots & \ddots & \vdots \\
0 & \cdots & 0 & \cdots & 0 & \cdots & 1
\end{array}\right],\text{if:} i,j<n,\text{or}: i,j>n
$$

$$
\begin{align}
&S'_{ii}=G_{ik}^{\dagger}S_{kl}G_{li}= c^2S_{ii}-2csS_{ij}+s^2S_{jj}\\
&S'_{jj}=G_{jk}^{\dagger}S_{kl}G_{lj}= c^2S_{jj}+2csS_{ij}+s^2S_{ii}\\
&S'_{ij}=S'_{ji}=G_{ik}^{\dagger}S_{kl}G_{lj}= (c^2-s^2)S_{ij}+cs(S_{ii}-S_{jj})\\
&S'_{ik}=S'_{ki}=G_{il}^{\dagger}S_{lm}G_{mk}=cS_{ik}-sS_{jk}   &k\neq i,j\\
&S'_{jk}=S'_{kj}=G_{jl}^{\dagger}S_{lm}G_{mk}=cS_{jk}+sS_{ik}   &k\neq i,j\\
&S'_{kl} = S_{kl} &k,l\neq i,j
\end{align}
$$

#### Zero mode compatibility

##### zero mode is out of SP(2,R)

The kernel to solve zero mode divergence problem is finding the answer of diagonlization of such a matrix
$$
\begin{pmatrix}1&1\\1&1\end{pmatrix}
$$
The result of how to diagonalize such matrix by $\text{SP(2,R)}$  is impossible. The eigen vector of this matrix is 

$\begin{pmatrix}-1 & 1\\ 1 & 1\end{pmatrix}$ . So that under the metric $\eta$, the norm of these eigen vectors are zero. This means that our givens rotation method could not be compatible with these vectors because these vectors correspond to such operators $[d_i,d_i^{\dagger}]=0$ which is not canonical. We need to patch the algorithm in zero mode cases.

##### Patch in Algorithm

In the patched algorithm, we used a strategy to patch such problem. 

In previous algorithm, we always try to eliminate the largest element in iteration. Now, we always try to eliminate the largest element which will not **raise zero mode problem**. Obviously, in iteration, the rotation parameter $\theta$ should converge to zero.  So if a rotation step propose a quite large $\theta$ (in our program, $\tanh(2\theta)>0.9$ means $\theta$ is too large) , the program will realize this step will cause a **zero mode problem**. Then the program will try to put off such rotation and try to use the second largest off-diagonal element.

In this way, the program can happily eliminate all elements which will not raise zero mode problem. And at last, the matrix we want to diagonlize will become such form.
$$
\left[\begin{array}{ccccccc}
\lambda_1 & \cdots & 0 & \cdots & 0 & \cdots & 0 \\
\vdots & \ddots & \vdots & & \vdots & & \vdots \\
0 & \cdots & c & \cdots & c & \cdots & 0 \\
\vdots & & \vdots & \ddots & \vdots & & \vdots \\
0 & \cdots & c & \cdots & c & \cdots & 0 \\
\vdots & & \vdots & & \vdots & \ddots & \vdots \\
0 & \cdots & 0 & \cdots & 0 & \cdots & \lambda_N
\end{array}\right] or \left[\begin{array}{ccccccc}
\lambda_1 & \cdots & 0 & \cdots & 0 & \cdots & 0 \\
\vdots & \ddots & \vdots & & \vdots & & \vdots \\
0 & \cdots & c & \cdots & -c & \cdots & 0 \\
\vdots & & \vdots & \ddots & \vdots & & \vdots \\
0 & \cdots & -c & \cdots & c & \cdots & 0 \\
\vdots & & \vdots & & \vdots & \ddots & \vdots \\
0 & \cdots & 0 & \cdots & 0 & \cdots & \lambda_N
\end{array}\right]
$$
Now we can simply use $\begin{pmatrix}-1 & 1\\ 1 & 1\end{pmatrix}$ to diagonalize such matrix.

This party of function is correlated to GivensBoson.deal_zeromodes(type=:normal)

##### Another Patch.

This rotation is mathematically strict. Though this is not what we want.

The matrix we used in the previous patch is  $\begin{pmatrix}-1 & 1\\ 1 & 1\end{pmatrix}$.  There is a important problem for this matrix. The final matrix $T$ does not have the form  $\begin{pmatrix}U & V\\ V & U\end{pmatrix}$. This means the system can not be described by bosonic generator and annihilator. If we want the final modes are "not canonical but indeed generator and annilator"(sorry, I know this is weird...). We need to try another patch.

Another patch came from a mathematically ill-defined rotation. Which takes $\theta\rightarrow \pm\infty$.

while $c= \cosh(\theta),s=\sinh(\theta)$
$$
G_a(i, j, \theta)=\left[\begin{array}{ccccccc}
1 & \cdots & 0 & \cdots & 0 & \cdots & 0 \\
\vdots & \ddots & \vdots & & \vdots & & \vdots \\
0 & \cdots & c & \cdots & s & \cdots & 0 \\
\vdots & & \vdots & \ddots & \vdots & & \vdots \\
0 & \cdots & s & \cdots & c & \cdots & 0 \\
\vdots & & \vdots & & \vdots & \ddots & \vdots \\
0 & \cdots & 0 & \cdots & 0 & \cdots & 1
\end{array}\right]  ,\text{if}: i<n,j>n
$$
Then it seems $c$ and $s$ tend to infinity. We need to claim. let $c$ and $s$ equals one.

This party of function is correlated to GivensBoson.deal_zeromodes(type=:abnormal)